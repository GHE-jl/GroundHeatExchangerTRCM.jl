# Short and long term

Without a buried depth, the TRCM of [`trcm`](@ref) is built for the short term. Two of its
assumptions limit it in the long term:

- the ground is held at its initial temperature at the far-field radius ``r_s``: once the thermal
  front reaches ``r_s``, the response levels off towards an artificial steady state;
- the depth layers are coupled only through the fluid: the TRCM behaves like a stack of infinite
  cylindrical sources and misses the axial effects (ground surface, borehole toe, neighbouring
  boreholes) of the finite line source.

[`trcm`](@ref) with a buried depth `D` lifts both limitations at a higher cost (see
[Vertical conduction in the ground](@ref)). The faster alternative joins the TRCM to a long-term
ground model: the method `outlet_transfer_function(t, ex, ks, Cs, H, D, V̇, fluid; ...)` of
[GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl), defined in its package
extension and available after `using GroundHeatExchanger, ThermalResistanceCapacitanceModel`. With
`long = :fls` (default), it uses each model where it is valid. Up to a contact time ``t_c`` it
returns the TRCM; after ``t_c`` it adds the increase of the borehole wall response of a ground model
``g_m`` (keyword `m`), the finite line source `FLSModel(H, D, ks, Cs)` by default:

```math
g(t) = g_\text{TRCM}(t_c) + \frac{\dot V C_f}{H}\left[g_m(t) - g_m(t_c)\right], \qquad t > t_c.
```

The response model injects heat at the constant rate ``\dot V C_f`` from ``t = 0``, so once the
borehole capacities are saturated the outlet temperature follows the wall temperature of a
constant-flux source, shifted by a constant (the effective borehole resistance). That constant
cancels out of the junction, which is continuous by construction: the same vertical shift that
`outlet_transfer_function` applies to the short-term ANNs.

With `long = :trcm`, `outlet_transfer_function` returns [`trcm`](@ref) with the buried depth `D`
for a single borehole, and the BC-IV borefield built from [`trcm_pairs`](@ref) for several
boreholes (see [Borehole-to-borehole responses](@ref)).

## Contact time and far-field radius

The reference is [`trcm`](@ref) with the buried depth `D`: the same TRCM columns with vertical
conduction in the ground, the buried depth, the ground surface and the ground below the toe, valid
from seconds to decades. Single U-tube, ``H = 150`` m, ``D = 4`` m, ``k_s = 2.5`` W/mK,
``C_s = 2.2`` MJ/m³K, ``\dot V = 0.5`` L/s (`script/script_contact_time.jl`).

Maximum relative error of the hybrid between 1 hour and 50 years:

| ``t_c`` | 1 h | 6 h | 1 d | 3 d | 7 d | 21 d | 60 d | 180 d |
|---|---|---|---|---|---|---|---|---|
| ``r_s = 5`` m | 12.1 % | 1.8 % | 0.62 % | 0.46 % | 0.46 % | 0.46 % | 0.46 % | — |
| ``r_s = 20`` m | 12.2 % | 1.8 % | 0.65 % | 0.52 % | 0.52 % | 0.52 % | 0.52 % | 0.52 % |
| ``r_s = 50`` m | 12.2 % | 1.8 % | 0.65 % | 0.52 % | 0.52 % | 0.52 % | 0.52 % | 0.52 % |

- Before about one day the borehole capacities are not saturated, and the finite line source
  (no capacity) underestimates the increase that follows.
- From 3 days on the error is flat: the remaining half percent is the difference between the two
  TRCMs at short times (vertical and radial discretization), not the junction. Much later contact
  times eventually accumulate the axial effects that the TRCM lacks.
- The far-field radius does not matter as long as the far field is not felt at ``t_c`` (with
  ``r_s = 5`` m it is felt before 180 days); it costs almost nothing (log-spaced shells).

The defaults are ``t_c = 7`` days (keyword `tc` of `outlet_transfer_function`), the horizon of the
published short-term ANN of GroundHeatExchangerANN.jl, and ``r_s = 20`` m ([`TRCMSettings`](@ref)).

## With or without vertical conduction

Same single borehole, 160 times from 10 s to 50 years: relative difference with the reference
(`trcm` with `D`, vertical conduction in the ground) and computing time of one call (warm,
indicative: it varies by about a factor 2 between runs).

| Model | 1 d | 1 yr | 10 yr | 50 yr | Time |
|---|---|---|---|---|---|
| [`trcm`](@ref) with `D` (reference): vertical conduction | — | — | — | — | 1.5–3 s |
| [`trcm`](@ref) without `D`: no vertical conduction (``r_s = 500`` m) | −0.36 % | +0.40 % | +2.07 % | +5.28 % | 0.1–0.2 s |
| Hybrid TRCM + FLS, `long = :fls` (``t_c = 7`` d, ``r_s = 20`` m) | −0.41 % | −0.29 % | −0.13 % | +0.07 % | 0.08 s |

- Without vertical conduction, the TRCM drifts upwards in the long term: it misses the heat lost
  through the ground surface and below the toe.
- The hybrid stays within about half a percent of the reference at the cost of the TRCM up to
  7 days, 20 to 40 times faster.
- `outlet_transfer_function(...; long = :trcm)` with a single borehole is the reference itself.

## Borefields

With several boreholes (`xy` with more than one row) and `long = :fls`, the long-term part is the
field response of `ground_response` of GroundHeatExchanger.jl under the boundary condition `bc`:
BC-II by default (the same heat rate and the same mean wall temperature in every borehole), or
`:I`, `:III` (`nseg` segments per borehole). `bc = :IV` (single U-tube, finite line source) uses
`bc_IV` of GroundHeatExchanger.jl with the ``R_b`` and ``R_a`` of the Δ-circuit of the TRCM. The
contact time must precede the thermal interactions between the boreholes.

With `long = :trcm`, the field has the same inlet temperature in every borehole (BC-IV) and is built
from the borehole-to-borehole responses of [`trcm_pairs`](@ref) (Rose et al., 2026). For a 5×5
field, the BC-II hybrid departs from it by −0.40 %, −0.65 %, −0.35 % and +0.16 % at 1, 5, 10 and
20 years, in 0.04 s instead of about 13 s (see [Borehole-to-borehole responses](@ref)).

## Groundwater flow

Any ground model re-exported by GroundHeatExchanger.jl can be passed as `m`, e.g. the moving
finite line source `MFLSModel(H, rb, D, ks, Cs, Cw, vD)` for a Darcy flux `vD`. The TRCM part
before ``t_c`` has no groundwater flow.
