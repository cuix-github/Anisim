# Experimental PF-FLIP phase contract

2026-10-23. This is a uniform-grid experiment for the unresolved NumPy/Bifrost
2D stability gate, not a new production default. The old graph (`phase_model=0`)
keeps its weighted phase-tag path, and the old `Sim` keeps its original path.

The two explicit experimental modes share the following mathematical definition
between an independent `PhaseContractSim` and the actual fused C++ graph:

1. Deposit **raw sums** of primary-particle mass with PF-FLIP Eq.6,
   `w=max(1-distance²/r²,0)³`, with `r=dx`. Auxiliary spray has zero mass and
   contributes to no phase or velocity channel. Momentum uses the same kernel.
2. Use Eq.7, `phi=min(sqrt(max(raw-log(rho_l/rho_g)*rho_g*rho0,0) /
   (alpha_phi*rho0*rho_l)),1)`. `rho0` is the calibrated kernel sum of a uniform
   reference seeding, not a per-face normalization by the instantaneous particle
   count. `alpha_phi=1` is the paper's default.
3. Mode 1 stops here. Mode 2 also deposits mass at cell centres with the same
   kernel and evaluates Eq.7 there. Cells with phi >= 0.5 are liquid. Interior
   faces between two liquid cells are set to 1; between two air cells to 0;
   differing neighbours retain their face-sampled Eq.7 value. Closed boundary
   faces retain their raw phase value but have zero pressure coupling/velocity.

The final pure-face operation is required by PF-FLIP §3.3. **The separate
cell-centred splat is our explicit sampling choice**: the paper states the cell
tagging rule but does not specify how to derive cell samples from face samples.
We do not claim this choice is the authors' unpublished implementation. The
tested 2D seeding is the same regular 2x2 subcell lattice used for face calibration.

The momentum and mean-mass float32 ports retain their existing semantics.
Experimental phase ports instead carry the compressed/snap-adjusted phase,
computed from double raw accumulators before a float32 port conversion. The
independent Python implementation rounds particle momentum and the same field
boundaries to float32 explicitly when `port_rounding=True`. Disabling that switch
is a separate precision experiment, never an unrecorded change to a comparison.

ST weights, adaptive particle volumes and kernel support, interior solid support,
moving boundaries, forcing and adaptive pressure grids require separate contracts
and gates. This first implementation **refuses** those combinations rather than
applying the uniform calibration silently. Their intended raw deposit is a sum of
`kernel * rho_phase * represented_volume * temporal_weight`, but the calibration,
kernel support near solid/refinement boundaries and pure-cell tagging must be
verified before these modes are enabled. Existing legacy features remain available
with phase_model=0.

Acceptance starts with Eq.7 analytic points, the pure-face truth table, actual
P2G field comparison, full graph single steps, and unsupported-mode refusal. Then
compare mode 1 and mode 2 on identical still-water and plunge seeds, without
emission first, measuring water/air speeds, particle energy and distribution.
Tolerance and timestep changes are separate arms. Agreement is necessary but is
not sufficient: both implementations must also satisfy the physical controls.

The run ledger and retained control outputs live under
`out/runs/phase_contract_20261023/`. Existing 40-second simulation caches and
render presets are not inputs to this job and are not changed.
