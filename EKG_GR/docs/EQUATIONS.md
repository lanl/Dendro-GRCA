# Formulation and generated sources

Conventions: signature (-,+,+,+), Pi=-n^a partial_a phi, K_ij=-L_n gamma_ij/2,
gamma_ij=gt_ij/chi, G=c=1. The conformal metric is denoted gt; its evolved
contracted connection Gt is a DIFFERENT variable. The inverse is computed as an
actual matrix inverse rather than replacing its determinant by one by fiat.
The native algebraic-constraint treatment still enforces the chosen BSSN
normalization during evolution.

## Matter

$$
\partial_t\phi=\beta^i\partial_i\phi-\alpha\Pi,
$$
$$
\partial_t\Pi=\beta^i\partial_i\Pi+\alpha K\Pi
-\alpha D^iD_i\phi-D^i\alpha D_i\phi+\alpha V'(\phi).
$$

$$
\rho=\tfrac12(\Pi^2+D_i\phi D^i\phi)+V,\qquad S_i=\Pi\partial_i\phi,
$$
$$
S_{ij}=\partial_i\phi\partial_j\phi+
\gamma_{ij}[\tfrac12\Pi^2-\tfrac12 D_k\phi D^k\phi-V],
\quad S=\tfrac32\Pi^2-\tfrac12 D_k\phi D^k\phi-3V.
$$

## Einstein and gauge coupling

$$
\delta\dot K=4\pi\alpha(\rho+S)=8\pi\alpha(\Pi^2-V),
$$
$$
\delta\dot{\tilde A}_{ij}=-8\pi\alpha\chi
[\partial_i\phi\partial_j\phi]^{\mathrm{TF}},\qquad
\delta\dot{\tilde\Gamma}^{i}=-16\pi\alpha\tilde\gamma^{ij}S_j.
$$

The standard second-order Gamma driver uses the FULL connection RHS in its B
RHS. Therefore its source correction is the same connection correction for the
implemented puncture gauge. Synchronous gauge fixes lapse and shift and does
not use that driver. Gauge constants and advection switches are visible in the
Python source, not hidden in generated C++.

The conformal Ricci term follows the symbolic BSSN structure in the uploaded
`bssn_eqns_config.py`. Scalar Laplacians use geometrical Christoffels. Constraint
Ricci is calculated from metric derivatives, without replacing the geometric
connection by an independently evolved field.

$$
\mathcal H=R+\tfrac23 K^2-\tilde A_{ij}\tilde A^{ij}-16\pi\rho,
\qquad
\mathcal M_i=D_j A^j{}_i-\tfrac23\partial_iK-8\pi S_i.
$$

Matter terms in these diagnostics follow the compile-time backreaction choice.
In passive-field mode the diagnostic is the vacuum Einstein constraint on the
background, not the residual of a coupled Einstein system that is not evolved.

## ID103: time-symmetric radial Hamiltonian solve

Take gt=identity, K=At=Pi=0, chi=psi^(-4). The compact C-infinity scalar profile is

$$
\phi(r)=A\exp[-r^2/(2\sigma^2)-r^2/(R_c^2-r^2)]\quad(0\le r<R_c),
\qquad \phi=0\quad(r\ge R_c).
$$

The Hamiltonian constraint reduces to

$$
\psi''+\frac{2}{r}\psi'=-\pi\psi(\phi')^2-2\pi\psi^5 V(\phi).
$$

Use regularity at the origin and, for R>Rc, the exact vacuum exterior condition
psi(R)+R psi'(R)=1. The exterior mass is M=-2R^2 psi'(R). The implementation uses
radial RK4 shooting on the weak positive branch and Hermite interpolation.
The scalar potential comes from the same generated potential expression used
in evolution. Radial and octree discretization errors must be tested separately.
This does not construct a black hole or solve binary puncture constraints.

## Scope

The algebraic tests check local expressions. Their success does not establish
finite-difference convergence, boundary accuracy, remeshing conservation or
native MPI correctness. Those require the run sequence in README.md.
