#!/usr/bin/env bash
# Rank-local wrapper launched by mpiexecjl. Each MPI rank runs one independent
# pi-flux QSH Hubbard-Kanamori DQMC Markov chain with its own seed, output
# directory, and serialized checkpoint file.
set -u -o pipefail

rank=${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-${SLURM_PROCID:-0}}}
world=${OMPI_COMM_WORLD_SIZE:-${PMI_SIZE:-${SLURM_NTASKS:-1}}}
rank_pad=$(printf "%05d" "${rank}")

: "${PROJECT:?PROJECT must point to the Pi_flux_QSH_Hubbard_Kanamori checkout}"
: "${OUT_ROOT:?OUT_ROOT must be set by the sbatch script}"

JULIA_PROJECT=${JULIA_PROJECT:-${PROJECT}/vendor/SmoQyDQMC.jl}
DQMC_DRIVER=${DQMC_DRIVER:-${PROJECT}/DQMC/SmoqyDQMC/scripts/run_piflux_qsh_smoqy.jl}

rank_dir="${OUT_ROOT}/ranks/rank_${rank_pad}"
mkdir -p "${rank_dir}"
checkpoint_file="${rank_dir}/checkpoint.jls"

sid=$(( ${SID_BASE:-920000} + rank ))
seed=$(( ${SEED_BASE:-20260610} + rank ))

printf '[%s] rank=%s/%s host=%s sid=%s seed=%s out=%s\n' \
  "$(date -Is)" "${rank}" "${world}" "$(hostname)" "${sid}" "${seed}" "${rank_dir}"

set +e
julia --project="${JULIA_PROJECT}" "${DQMC_DRIVER}" \
  --Lx="${LX}" \
  --Ly="${LY}" \
  --t="${T}" \
  --lambda="${LAMBDA}" \
  --U="${U}" \
  --JH="${JH}" \
  --mu="${MU}" \
  --open_x="${OPEN_X}" \
  --density_kanamori="${DENSITY_KANAMORI}" \
  --spin_flip_hund="${SPIN_FLIP_HUND}" \
  --pair_hopping="${PAIR_HOPPING}" \
  --beta="${BETA}" \
  --dtau="${DTAU}" \
  --Ntherm="${NTHERM}" \
  --Nmeas="${NMEAS}" \
  --Nupdates="${NUPDATES}" \
  --n_stab="${N_STAB}" \
  --dGmax="${DGMAX}" \
  --symmetric=false \
  --checkerboard=false \
  --seed="${seed}" \
  --sID="${sid}" \
  --outdir="${rank_dir}" \
  --checkpoint_enable=true \
  --checkpoint_file="${checkpoint_file}" \
  --checkpoint_freq_hours="${CHECKPOINT_FREQ_HOURS}" \
  --runtime_limit_hours="${RUNTIME_LIMIT_HOURS}" \
  --checkpoint_keep=true \
  --checkpoint_resume=true
rc=$?
set -e

# The Julia driver exits 13 after a clean soft-stop checkpoint. Convert that
# rank-local condition to success so mpiexecjl waits for other independent ranks.
if [[ "${rc}" == "13" ]]; then
  if [[ -f "${checkpoint_file}" ]]; then
    printf '[%s] rank=%s checkpoint stop converted to wrapper success: %s\n' \
      "$(date -Is)" "${rank}" "${checkpoint_file}"
    exit 0
  fi
fi

exit "${rc}"
