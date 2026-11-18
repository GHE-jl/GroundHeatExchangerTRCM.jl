# Design notes

This package is the packaged version of the research code in the `TRCM` folder (`TRCM.jl`,
`MTRCM.jl`, `Quasi3D.jl`, and `gTRCM.m`, the MATLAB quasi-3D model of Pasquier & Marcotte, 2014).
The physics and the verification findings of that work were kept; this page lists what changed.
The model itself is described in the Modeling theory pages.

## Corrections to the research code

1. **Convective film resistance.** The research code imposed `Rf = 0` (fluid at the pipe wall
   temperature), also in the quasi-3D response model where the fluid flows. The films are now
   computed from the flow rate (`BoreholeResistance.resistance_fluid`, Gnielinski, and the annulus
   correlation of Lamarche, 2021) and enter both the multipole `Rb`, `Ra` and the network.
   `films = false` reproduces the fixed pipe-wall temperature of Pasquier & Marcotte (2012).
2. **Grout core resistance at the wrong location.** The MTRCM core chain (`Rgg`) was attached to the
   first grout node, at `x = 1/ng` of the grout resistance, but `Rgg` was computed from Eq. 10 at
   `x = 0`. The network then reproduced `Rb` exactly but not `Ra` (the thermal short-circuit); only
   `Rb` was checked. `Rgg` is now derived at the actual attachment point (`_core_conductances`), and both
   `Rb` and `Ra` are tested to machine precision.
3. **Capacities of open chains.** In `add_chain!` (open chains: pipe walls, casing), the first node
   received half of shell 1 already credited to the previous node, and half of the last shell was
   lost. The capacity of each layer now equals the capacity of its materials exactly (tested).
4. **Energy balance of the advection.** The advective term `v/dz·(T_up - T)` was added to `dT/dt`
   directly, although the fluid node also holds half of the first pipe shell: the advected heat was
   overestimated by the ratio of the node capacity to the fluid capacity.
   The advection is now assembled in energy units before dividing by the capacities; the steady-state
   heat balance is exact.
5. **Inlet condition.** The research port enforced `Tin = Tout + α` with a symmetric near-zero
   resistance `Rloop`, which instantly split the inlet and the outlet to `±α/2`. In the 2014 paper
   (Fig. 5) the outlet stays at zero until the fluid residence time and the inlet starts at `α`. The
   inlet temperature is now the upstream temperature of the inlet channels, `Tin = Tout + α`, which
   reproduces that behaviour and removes the stiff penalty resistance.
6. **Jacobian.** The research code let the solver build a dense Jacobian by automatic
   differentiation; the exact sparse Jacobian is now provided.

## Other changes

- Single U-tube, double U-tube (diagonal or adjacent inlets), coaxial and standing column well
  cross-sections share one network builder; the U-tube grout network is derived for any number of
  pipes from the Δ-circuit.
- Negative Δ-circuit resistances: the non-physical `rgg_fix = :abs` workaround of `gTRCM.m` was
  dropped (see [Negative resistances](@ref)).
- The far-field radius `rs` and the number of soil shells `ns` are settings
  ([`TRCMSettings`](@ref)).
- The short soil column below the toe of `gTRCM.m` (Pasquier & Marcotte, 2014, Fig. 3) is not
  included without a buried depth; with a buried depth, [`trcm`](@ref) has the complete vertical
  conduction of the ground instead (see [Vertical conduction in the ground](@ref)).
- The U-turn at the toe is an upwind link from the down-going to the up-going channel, instead of a
  merged node.
- With a partial casing, the uncased layers carry a ring of ground of the same thickness, so the
  soil discretization is the same in every layer.

## Verification findings

- The steady state of the network is best verified with a direct linear solve, not by running the
  ODE for a long time (the far field equilibrates slowly).
- A steel casing lowers the outlet response, mostly at short times.
