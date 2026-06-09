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

## Particle-hole/sign convention

At half filling (`mu=0`) the Hubbard-only model is simulated in particle-hole
symmetric form

```math
U \sum_{i,a} (n_{ia↑}-1/2)(n_{ia↓}-1/2).
```

The bipartite parity in the original square lattice is
`η_{x,y,A}=(-1)^y`, `η_{x,y,B}=-(-1)^y`.  Real π-flux hoppings connect opposite
η sectors; imaginary QSH hoppings connect equal η sectors and are complex conjugate
between spins.  This gives the DQMC Hubbard-only sign-free parent at half filling.

## Interaction tiers

The tiers use the same local two-orbital labels `(A,B)`:

- `hubbard`: onsite intra-orbital Hubbard only.
- `density`: `hubbard` plus density-density/Ising Kanamori with
  `V_opposite=U-2JH`, `V_same=U-3JH`.
- `spinflip`: `density` plus `-JH(S_A^+ S_B^- + S_A^- S_B^+)`.
- `full`: `spinflip` plus `+JH(P_A†P_B + P_B†P_A)`, with `Jpair=JH`.

Only the Hubbard-only parent is expected to be exactly sign-problem-free; the
Kanamori tiers are included to measure controlled sign/phase degradation.
