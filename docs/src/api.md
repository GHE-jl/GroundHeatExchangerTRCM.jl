# API reference

```@index
Modules = [ThermalResistanceCapacitanceModel]
```

## Exchangers, fluid and discretization

```@docs
AbstractExchanger
SingleUTube
DoubleUTube
Coaxial
StandingColumn
Casing
FluidProperties
TRCMSettings
```

## Outlet transfer functions

```@docs
trcm
trcm_pairs
```

`outlet_transfer_function(t, ex, ks, Cs, H, D, V̇, fluid; ...)` is documented with
[GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl).

## Verification tools

Not exported (call them as `ThermalResistanceCapacitanceModel.name`).

```@docs
ThermalResistanceCapacitanceModel.delta_circuit
ThermalResistanceCapacitanceModel.transient_resistance
```
