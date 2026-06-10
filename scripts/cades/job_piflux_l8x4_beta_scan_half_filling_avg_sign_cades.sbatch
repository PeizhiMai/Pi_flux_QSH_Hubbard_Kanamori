#!/bin/bash
#SBATCH -A ccsd
#SBATCH -p burst
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --ntasks-per-node=32
#SBATCH -c 1
#SBATCH --mem=100G
#SBATCH -t 02:30:00
#SBATCH --chdir=/home/9pm/Pi_flux_QSH_Hubbard_Kanamori
#SBATCH --job-name=pfQSH_sign
#SBATCH --output=/home/9pm/Pi_flux_QSH_Hubbard_Kanamori/logs/%x_%A_%a_%j.out
#SBATCH --error=/home/9pm/Pi_flux_QSH_Hubbard_Kanamori/logs/%x_%A_%a_%j.err
#SBATCH --array=0-15%8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=peizhimai@gmail.com

# CADES MPI/checkpoint average-phase scan for the clean pi-flux QSH
# Hubbard-Kanamori model. Compact grid only: four interaction tiers at their
# physical half-filled chemical potentials for beta=7,10,12,16. No mu scan.

set -u -o pipefail

PROJECT=${PROJECT:-/home/9pm/Pi_flux_QSH_Hubbard_Kanamori}
LOG_DIR=${PROJECT}/logs
mkdir -p "${LOG_DIR}"
cd "${PROJECT}"

LX=${LX:-8}
LY=${LY:-4}
T=${T:-1.0}
LAMBDA=${LAMBDA:-0.2}
U=${U:-1.0}
JH=${JH:-0.25}
DTAU=${DTAU:-0.1}
NTHERM=${NTHERM:-2000}
NMEAS=${NMEAS:-10000}
NUPDATES=${NUPDATES:-2}
N_STAB=${N_STAB:-10}
DGMAX=${DGMAX:-1e-6}
OPEN_X=${OPEN_X:-false}
BC=${BC:-periodic}

BETAS=(7.0 10.0 12.0 16.0)
TIERS=(hubbard density spinflip full)
NBETA=${#BETAS[@]}
NTIER=${#TIERS[@]}
TOTAL=$((NBETA * NTIER))
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
if (( TASK_ID < 0 || TASK_ID >= TOTAL )); then
  echo "Invalid SLURM_ARRAY_TASK_ID=${TASK_ID}; expected 0..$((TOTAL-1))" >&2
  exit 2
fi

tier_idx=$((TASK_ID % NTIER))
beta_idx=$((TASK_ID / NTIER))
TIER=${TIERS[$tier_idx]}
BETA=${BETAS[$beta_idx]}

case "${TIER}" in
  hubbard)
    DENSITY_KANAMORI=false
    SPIN_FLIP_HUND=false
    PAIR_HOPPING=false
    MU=$(python3 - <<PY
U=float('${U}')
print(f'{U/2:.16g}')
PY
)
    ;;
  density)
    DENSITY_KANAMORI=true
    SPIN_FLIP_HUND=false
    PAIR_HOPPING=false
    MU=$(python3 - <<PY
U=float('${U}'); JH=float('${JH}')
print(f'{(3*U-5*JH)/2:.16g}')
PY
)
    ;;
  spinflip)
    DENSITY_KANAMORI=true
    SPIN_FLIP_HUND=true
    PAIR_HOPPING=false
    MU=$(python3 - <<PY
U=float('${U}'); JH=float('${JH}')
print(f'{(3*U-5*JH)/2:.16g}')
PY
)
    ;;
  full)
    DENSITY_KANAMORI=true
    SPIN_FLIP_HUND=true
    PAIR_HOPPING=true
    MU=$(python3 - <<PY
U=float('${U}'); JH=float('${JH}')
print(f'{(3*U-5*JH)/2:.16g}')
PY
)
    ;;
  *)
    echo "Unknown TIER=${TIER}" >&2
    exit 2
    ;;
esac

RANKS=${SLURM_NTASKS:-32}
CHECKPOINT_FREQ_HOURS=${CHECKPOINT_FREQ_HOURS:-2.0}
RUNTIME_LIMIT_HOURS=${RUNTIME_LIMIT_HOURS:-2.0}
MAX_RESUBMITS=${MAX_RESUBMITS:-12}
RESUBMIT_COUNT=${RESUBMIT_COUNT:-0}
AUTO_RESUBMIT=${AUTO_RESUBMIT:-true}
SEED_BASE=${SEED_BASE:-$((2026061000 + 100000 * TASK_ID))}
SID_BASE=${SID_BASE:-$((920000 + 1000 * TASK_ID))}

mu_tag=$(printf '%+.3f' "${MU}" | sed 's/+//; s/-/m/; s/\./p/g')
beta_tag=$(printf '%.1f' "${BETA}" | sed 's/\./p/g')
run_name="L${LX}x${LY}_${BC}_${TIER}_mu${mu_tag}_U${U}_JH${JH}_b${BETA}_dt${DTAU}_w${NTHERM}_m${NMEAS}perrank"
run_name=$(echo "${run_name}" | sed 's/\./p/g')
OUT_ROOT=${OUT_ROOT:-${PROJECT}/runs/piflux_L8x4_beta7_10_12_16_half_filling_avg_sign/${run_name}}

success_outputs_present() {
  [[ -d "${OUT_ROOT}/ranks" ]] || return 1
  local n
  n=$(find "${OUT_ROOT}/ranks" -mindepth 3 -maxdepth 3 -name summary.txt | wc -l | tr -d ' ')
  [[ "${n}" -ge "${RANKS}" ]]
}

if success_outputs_present; then
  echo "[$(date -Is)] run already complete: ${OUT_ROOT}"
  exit 0
fi

mkdir -p "${OUT_ROOT}/ranks"

# Avoid CADES user-env/system-OpenMPI contamination; let mpiexecjl use Julia's OpenMPI_jll.
unset LD_LIBRARY_PATH MPI_PATH MPI_ROOT MPICC MPICXX MPIF77 MPIF90 MPIFC
unset OMPI_MCA_pml OMPI_MCA_btl
export PATH="$HOME/.julia/bin:$HOME/.juliaup/bin:$PATH"
export JULIA_DEPOT_PATH=${JULIA_DEPOT_PATH:-${PROJECT}/.julia_depot:$HOME/.julia}
export JULIA_NUM_THREADS=${JULIA_NUM_THREADS:-1}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export MKL_DYNAMIC=${MKL_DYNAMIC:-FALSE}

export PROJECT JULIA_PROJECT=${JULIA_PROJECT:-${PROJECT}/vendor/SmoQyDQMC.jl}
export DQMC_DRIVER=${DQMC_DRIVER:-${PROJECT}/DQMC/SmoqyDQMC/scripts/run_piflux_qsh_smoqy.jl}
export OUT_ROOT LX LY T LAMBDA U JH MU OPEN_X
export DENSITY_KANAMORI SPIN_FLIP_HUND PAIR_HOPPING BETA DTAU NTHERM NMEAS NUPDATES N_STAB DGMAX
export CHECKPOINT_FREQ_HOURS RUNTIME_LIMIT_HOURS SEED_BASE SID_BASE

START_EPOCH=$(date +%s)
echo "[$(date -Is)] starting pi-flux QSH avg-sign task=${TASK_ID}/${TOTAL} job=${SLURM_JOB_ID:-manual} array=${SLURM_ARRAY_JOB_ID:-manual}"
echo "host=$(hostname) ranks=${RANKS} resubmit=${RESUBMIT_COUNT}/${MAX_RESUBMITS} out=${OUT_ROOT}"
echo "params: L=${LX}x${LY} BC=${BC} open_x=${OPEN_X} tier=${TIER} density=${DENSITY_KANAMORI} spin_flip=${SPIN_FLIP_HUND} pair=${PAIR_HOPPING}"
echo "one-body: t=${T} lambda=${LAMBDA}; interactions: U=${U} JH=${JH}; mu=${MU}; beta=${BETA} dtau=${DTAU}"
echo "statistics per rank: Ntherm=${NTHERM} Nmeas=${NMEAS} Nupdates=${NUPDATES}; total chains=${RANKS}"
echo "checkpoint: freq_hours=${CHECKPOINT_FREQ_HOURS} runtime_limit_hours=${RUNTIME_LIMIT_HOURS} walltime=02:30:00 keep=true"
echo "JULIA=$(command -v julia || true)"
echo "MPIEXECJL=$(command -v mpiexecjl || true)"

rc=0
mpiexecjl -n "${RANKS}" bash "${PROJECT}/scripts/cades/run_piflux_rank_checkpoint.sh" || rc=$?

END_EPOCH=$(date +%s)
echo "[$(date -Is)] launcher exit code: ${rc} runtime_seconds=$((END_EPOCH - START_EPOCH))"

if success_outputs_present; then
  echo "[$(date -Is)] complete: ${OUT_ROOT}"
  find "${OUT_ROOT}/ranks" -mindepth 3 -maxdepth 3 -name summary.txt | sort | sed 's#^#summary: #' | head -20
  exit 0
fi

checkpoint_count=0
if [[ -d "${OUT_ROOT}/ranks" ]]; then
  checkpoint_count=$(find "${OUT_ROOT}/ranks" -mindepth 2 -maxdepth 2 -name checkpoint.jls | wc -l | tr -d ' ')
fi
echo "[$(date -Is)] checkpoint_count=${checkpoint_count}/${RANKS}"

checkpoint_stop=false
if [[ "${checkpoint_count}" -ge "${RANKS}" ]]; then
  if [[ "${rc}" == "0" || "${rc}" == "1" || "${rc}" == "13" ]]; then
    checkpoint_stop=true
  fi
fi

if [[ "${AUTO_RESUBMIT}" == "true" && "${checkpoint_stop}" == "true" ]]; then
  if (( RESUBMIT_COUNT < MAX_RESUBMITS )); then
    next_count=$((RESUBMIT_COUNT + 1))
    script_path=$(readlink -f "${BASH_SOURCE[0]}")
    resubmit_account=${SLURM_JOB_ACCOUNT:-ccsd}
    resubmit_job_name=${SLURM_JOB_NAME:-pfQSH_resume}
    echo "[$(date -Is)] incomplete with valid checkpoints; resubmitting task ${TASK_ID} as account=${resubmit_account} script=${script_path} attempt ${next_count}/${MAX_RESUBMITS}"
    RESUBMIT_COUNT="${next_count}" sbatch -A "${resubmit_account}" --array="${TASK_ID}" --job-name="${resubmit_job_name}" "${script_path}"
    exit 0
  else
    echo "[$(date -Is)] incomplete and MAX_RESUBMITS reached: ${OUT_ROOT}" >&2
    exit 13
  fi
fi

echo "[$(date -Is)] failed/incomplete without valid checkpoint stop. rc=${rc} out=${OUT_ROOT}" >&2
exit "${rc}"
