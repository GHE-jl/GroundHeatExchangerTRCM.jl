module ThermalResistanceCapacitanceModel

using LinearAlgebra
using SparseArrays
using SciMLBase: ODEFunction, ODEProblem, solve, successful_retcode
using OrdinaryDiffEqRosenbrock: Rodas5P
using OrdinaryDiffEqBDF: FBDF
using LinearSolve: KLUFactorization
using BoreholeResistance: resistance_pipe, resistance_fluid, resistance_ULoop_borehole,
    resistance_ULoop_total_internal, Nusselt, Nusselt_annulus, convection_coefficient,
    fluid_property

# Ground heat exchangers, casing, fluid and discretization settings
include("exchangers.jl")
# Steady-state Δ-circuit of each exchanger and transient borehole resistance
include("borehole_resistance.jl")
# One depth layer of the TRCM: resistor/capacitor network of each exchanger
include("cross_section.jl")
# Network builder, linear ODE solve and small helpers
include("utils.jl")

# Quasi-3D TRCM response model (Pasquier & Marcotte, 2014), from seconds to months, or to decades
# with vertical conduction in the ground
include("trcm.jl")
# TRCM columns with vertical conduction, borehole-to-borehole responses (Rose et al., 2026)
include("trcm_pairs.jl")

# Ground heat exchangers, casing, fluid and discretization
export AbstractExchanger, SingleUTube, DoubleUTube, Coaxial, StandingColumn
export Casing, FluidProperties, TRCMSettings

# Outlet transfer functions
export trcm, trcm_pairs

end
