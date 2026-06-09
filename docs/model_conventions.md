# Model conventions: square π-flux QSH Hubbard-Kanamori

## One-body parent

Each unit cell has two orbitals/sublattices `(A,B)`.  The two-orbital unit cell is
the doubled magnetic unit cell of a square π-flux lattice.  Spin is conserved.
The spin-up block is a Chern insulator and the spin-down block is its time-reversal
partner, giving a Kane-Mele/BHZ-style QSH parent.

With target hopping coefficient `h c†_{r,a,σ} c_{r+δ,b,σ} + h.c.`, the implemented
nonzero hoppings are

| term | `(a,b,δ)` | coefficient for spin σ |
|---|---|---|
| π-flux x, intra-cell | `(A,B,(0,0))` | `t` |
| π-flux x, inter-cell | `(B,A,(1,0))` | `t` |
| π-flux y on A | `(A,A,(0,1))` | `t` |
| π-flux y on B | `(B,B,(0,1))` | `-t` |
| QSH diagonal | `(A,B,(0,1))` | `-i λ sσ` |
| QSH diagonal | `(A,B,(0,-1))` | `-i λ sσ` |
| QSH diagonal | `(B,A,(1,1))` | `-i λ sσ` |
| QSH diagonal | `(B,A,(1,-1))` | `+i λ sσ` |

where `s_up=+1` and `s_dn=-1`.  For `λ != 0`, the π-flux Dirac cones are gapped.
The H0 check computes the spin Chern numbers and expects `C_up=+1`, `C_dn=-1`
for positive `λ` with the default convention.

## Interaction convention: physical, not PH-shifted

ED and DQMC use the same physical interaction convention as
`BHZ_Hubbard_Kanamori`; the interaction operators are **not** written as
`(n-1/2)(n-1/2)` in either code path.

The Hubbard-only tier is

```math
H_U = U \sum_{i,a=A,B} n_{ia↑} n_{ia↓}.
```

The density-density Kanamori part is

```math
H_{dd} = \sum_i \Big[
  (U-3J_H)(n_{iA↑} n_{iB↑} + n_{iA↓} n_{iB↓})
 + (U-2J_H)(n_{iA↑} n_{iB↓} + n_{iA↓} n_{iB↑})
\Big].
```

The transverse terms are

```math
H_{sf} = -J_H \sum_i (S_{iA}^+S_{iB}^- + S_{iA}^-S_{iB}^+),
```

and, for the full tier,

```math
H_{pair} = +J_H \sum_i (P_{iA}^†P_{iB} + P_{iB}^†P_{iA}).
```

Because the interactions are physical, the physical half-filling chemical
potential is not generally `mu=0`:

```text
Hubbard-only:      mu_half = U/2
Kanamori tiers:    mu_half = (3U - 5JH)/2
```

For the default `U=1.0, JH=0.25`, this is `mu=0.5` for Hubbard-only and
`mu=0.875` for the density, spin-flip, and full Kanamori tiers.

## Sign-free Hubbard parent

The one-body π-flux/QSH Hamiltonian has the relevant bipartite/time-reversal
structure.  In the original square-lattice parity convention,
`η_{x,y,A}=(-1)^y`, `η_{x,y,B}=-(-1)^y`.  Real π-flux hoppings connect opposite
η sectors; imaginary QSH hoppings connect equal η sectors and are complex
conjugate between spins.

With the physical Hubbard interaction above, DQMC uses
`HubbardModel(ph_sym_form=false)`.  SmoQyDQMC internally adds the one-body shift
needed to decouple the physical `U n↑ n↓` term.  Therefore the Hubbard-only
sign-free half-filled point is reached by passing the physical `mu=U/2`.

## Interaction tiers

The tiers use the same local two-orbital labels `(A,B)`:

- `hubbard`: onsite intra-orbital Hubbard only.
- `density`: `hubbard` plus density-density/Ising Kanamori.
- `spinflip`: `density` plus `-JH(S_A^+ S_B^- + S_A^- S_B^+)`.
- `full`: `spinflip` plus `+JH(P_A†P_B + P_B†P_A)`, with `Jpair=JH`.

Only the Hubbard-only parent at its physical half-filling chemical potential is
expected to be exactly sign-problem-free; the Kanamori tiers are included to
measure controlled sign/phase degradation.

Implementation note: the DQMC spin-flip-only HST uses the same compensated
density coupling as the BHZ project internally, but the represented Hamiltonian
is still the physical spin-flip Kanamori tier above.
