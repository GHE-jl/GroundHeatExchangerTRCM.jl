# ThermalResistanceCapacitanceModel.jl

[![CI](https://github.com/GHE-jl/ThermalResistanceCapacitanceModel.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/GHE-jl/ThermalResistanceCapacitanceModel.jl/actions/workflows/CI.yml)
[![Docs: dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GHE-jl.github.io/ThermalResistanceCapacitanceModel.jl/dev)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPLv3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)

A Julia package for the **outlet transfer function of ground heat exchangers** with quasi-3D
thermal resistance and capacity models (TRCM; Pasquier & Marcotte, 2012, 2014): single and double
U-tubes, coaxial exchangers and standing column wells, from seconds to months, or to decades with
vertical conduction in the ground. It also computes the borehole-to-borehole responses of Rose et
al. (2026).

The package is the borehole layer of the GHE.jl ecosystem: it uses
[BoreholeResistance.jl](https://github.com/GHE-jl/BoreholeResistance.jl) for the steady-state
resistances, and [GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl) joins
it to the long-term ground models (`outlet_transfer_function`).

## Installation

Not registered yet. From the `GHE.jl` folder:

```julia
using Pkg
Pkg.develop(path = "ThermalResistanceCapacitanceModel.jl")
```

## Example

```julia
using ThermalResistanceCapacitanceModel

ex = SingleUTube(0.075, 0.0167, 0.0137, 0.06, 1.5, 2.0e6, 0.4, 1.8e6)   # rb, ro, ri, s, kg, Cg, kp, Cp
fluid = FluidProperties(10.0, :water)
t = exp10.(range(0, 7; length = 150))                                     # 1 s to ~4 months
g = trcm(t, ex, 2.5, 2.2e6, 150.0, 0.5e-3, fluid)                         # ks, Cs, H, V̇
```

## Documentation

The [documentation](https://GHE-jl.github.io/ThermalResistanceCapacitanceModel.jl/dev) has the
choice of function, a tutorial, the modeling theory and the API reference. The `script/` folder
(own environment, `julia --project=script script/<name>.jl`) reproduces its figures and numbers.
