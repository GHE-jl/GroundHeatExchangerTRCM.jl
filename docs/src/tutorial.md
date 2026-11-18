# Tutorial

This tutorial computes the outlet transfer function of a single U-tube, extends it to decades,
compares exchanger types and computes borehole-to-borehole responses.

## 1. Exchanger, ground and fluid

An exchanger stores its own geometry and materials; the ground, the length and the flow rate are
passed to the functions.

```julia
using ThermalResistanceCapacitanceModel

ex = SingleUTube(rb = 0.075, ro = 0.0167, ri = 0.0137, s = 0.06,   # radii and shank spacing [m]
                 kg = 1.5, Cg = 2.0e6,                              # grout [W/mK], [J/m³K]
                 kp = 0.4, Cp = 1.8e6)                              # pipe  [W/mK], [J/m³K]
ks, Cs = 2.5, 2.2e6            # ground conductivity [W/mK] and volumetric capacity [J/m³K]
H, D, V̇ = 150.0, 4.0, 0.5e-3   # borehole length [m], buried depth [m] and flow rate [m³/s]
fluid = FluidProperties(10.0, :water)      # water at 10 °C (CoolProp, via BoreholeResistance.jl)
```

`FluidProperties(kf, cf, ρf, μf)` gives the properties directly (conductivity, mass-specific heat
capacity, density, viscosity).

## 2. Short term: the TRCM

```julia
t = exp10.(range(0, 7; length = 150))   # 1 s to ~4 months
g = trcm(t, ex, ks, Cs, H, V̇, fluid)
```

`g` is the outlet transfer function: the outlet temperature for `Tin − Tout = 1`, i.e. for a
constant heat injection rate `V̇·ρf·cf` [W] from `t = 0`. It stays at zero until the fluid has
travelled the loop. For a heat injection rate `Q(t)`, the outlet temperature follows by convolution
with the increments of `Q/(ṁ·cf)` (Pasquier & Marcotte, 2014, Eqs. 19-22), e.g. with
`GroundHeatExchanger.convolution`.

## 3. From seconds to decades

Without a buried depth, the TRCM ground has no vertical conduction and is bounded at 20 m. With the
buried depth `D`, [`trcm`](@ref) adds the vertical conduction in the ground, the ground surface and
the ground below the toe (see [Vertical conduction in the ground](@ref)):

```julia
yr = 365.25 * 86400
t = exp10.(range(0, log10(25yr); length = 200))
g = trcm(t, ex, ks, Cs, H, V̇, fluid; D)
```

The faster alternative is the method `outlet_transfer_function` of GroundHeatExchanger.jl, which
keeps the TRCM up to 7 days and continues with the increase of the finite line source:

```julia
using GroundHeatExchanger

g = outlet_transfer_function(t, ex, ks, Cs, H, D, V̇, fluid)
```

The contact time is the keyword `tc` (7 days by default). For a borefield, pass the coordinates
(the long-term part is then the BC-II field response; `bc = :I`, `:III` or `:IV` for the other
boundary conditions), and any other ground model, e.g. with groundwater flow:

```julia
xy = borefield(:rectangle, 3, 3, 6.0)
g_field = outlet_transfer_function(t, ex, ks, Cs, H, D, V̇, fluid; xy)
g_gw = outlet_transfer_function(t, ex, ks, Cs, H, D, V̇, fluid;
                                m = MFLSModel(H, ex.rb, D, ks, Cs, 4.18e6, 1e-7))  # Darcy flux 1e-7 m/s
```

The TRCM part (before the contact time) has no groundwater flow (see [Short and long term](@ref)).

## 4. Other exchangers, casing, discretization

```julia
ex2 = DoubleUTube(rb = 0.075, ro = 0.0167, ri = 0.0137, s = 0.09, kg = 1.5, Cg = 2.0e6,
                  kp = 0.4, Cp = 1.8e6, network = :adjacent)
cx = Coaxial(rb = 0.075, rii = 0.0177, rio = 0.0201, roi = 0.0545, roo = 0.0572,
             kg = 1.5, Cg = 2.0e6, kpi = 0.4, Cpi = 1.8e6, kpo = 0.4, Cpo = 1.8e6)
scw = StandingColumn(rb = 0.076, ro = 0.05, ri = 0.045, kp = 0.4, Cp = 1.8e6)

g = trcm(t, ex2, ks, Cs, H, V̇, fluid;
         casing = Casing(10.0, 0.01),                  # 10 m of 1 cm steel casing
         settings = TRCMSettings(nz = 60, ns = 30))    # finer discretization
g_scw = trcm(t, scw, ks, Cs, H, V̇, fluid; D)          # every exchanger, with vertical conduction
```

The casing is only available without a buried depth.

## 5. Borehole-to-borehole responses and BC-IV borefields

[`trcm_pairs`](@ref) models an emitter (`Tin = Tout + 1`) and receivers in recirculation
(`Tin = Tout`) at the distances `d`, with vertical conduction in the ground (Rose et al., 2026):

```julia
g_self, g_pair = trcm_pairs(t, ex, ks, Cs, H, D, V̇, fluid, [3.0, 6.0, 12.0])
```

`g_self` is `trcm(t, ex, ks, Cs, H, V̇, fluid; D)` and `g_pair[:, j]` the outlet transfer function
of the receiver at `d[j]`. `outlet_transfer_function` with `long = :trcm` combines them into the
outlet transfer function of a borefield connected in parallel, with the same inlet temperature in
every borehole (BC-IV):

```julia
g_field = outlet_transfer_function(t, ex, ks, Cs, H, D, V̇, fluid; xy, long = :trcm)
```

## Scripts

The `script/` folder (own environment, `julia --project=script script/<name>.jl`) reproduces the
figures and numbers of this documentation:

| Script | Content |
|---|---|
| `script_exchanger_types.jl` | Outlet transfer functions of every exchanger type, TRCM and hybrid TRCM + FLS |
| `script_pasquier2012.jl` | Transient borehole resistance, scenarios a and e of Pasquier & Marcotte (2012) |
| `script_contact_time.jl` | Contact time and far-field radius of the hybrid TRCM + FLS, against the TRCM with vertical conduction |
| `script_rose2026.jl` | Borehole-to-borehole response functions and BC-IV borefield g-function (Rose et al., 2026) |
