# Ground heat exchangers (geometry and materials), casing, heat carrier fluid and discretization of
# the networks. Every constructor checks its inputs.

"""
    AbstractExchanger

Abstract type for the ground heat exchangers of the TRCM. A concrete exchanger only stores its
geometry and the properties of its own materials (pipes and grout); the ground, the borehole length,
the flow rate and the fluid are arguments of the response functions.
"""
abstract type AbstractExchanger end

"""
    SingleUTube(rb, ro, ri, s, kg, Cg, kp, Cp, ϵ = 5e-6)
    SingleUTube(; rb, ro, ri, s, kg, Cg, kp, Cp, ϵ = 5e-6)

Grouted single U-tube borehole: two pipes placed symmetrically at a distance `s/2` from the
borehole center. The down-going leg is pipe 1, the up-going leg is pipe 2.
# Arguments
    - `rb`: Borehole radius [m]
    - `ro`: Pipe outer radius [m]
    - `ri`: Pipe inner radius [m]
    - `s`: Shank spacing (center-to-center distance between the two legs) [m]
    - `kg`: Grout thermal conductivity [W/mK]
    - `Cg`: Grout volumetric heat capacity [J/m³K]
    - `kp`: Pipe thermal conductivity [W/mK]
    - `Cp`: Pipe volumetric heat capacity [J/m³K]
    - `ϵ`: Pipe roughness [m] (default 5e-6 for HDPE pipes)
"""
struct SingleUTube <: AbstractExchanger
    rb::Float64
    ro::Float64
    ri::Float64
    s::Float64
    kg::Float64
    Cg::Float64
    kp::Float64
    Cp::Float64
    ϵ::Float64
    function SingleUTube(rb, ro, ri, s, kg, Cg, kp, Cp, ϵ = 5e-6)
        all(>(0), (rb, ro, ri, s, kg, Cg, kp, Cp)) || throw(ArgumentError(
            "SingleUTube: all dimensions and properties must be strictly positive"))
        ri < ro || throw(ArgumentError(
            "SingleUTube: the pipe inner radius ri must be smaller than ro"))
        s > 2ro || throw(ArgumentError("SingleUTube: the two pipes overlap (s ≤ 2ro)"))
        s / 2 + ro < rb || throw(ArgumentError(
            "SingleUTube: the pipes cross the borehole wall (s/2 + ro ≥ rb)"))
        return new(rb, ro, ri, s, kg, Cg, kp, Cp, ϵ)
    end
end
SingleUTube(; rb, ro, ri, s, kg, Cg, kp, Cp, ϵ = 5e-6) =
    SingleUTube(rb, ro, ri, s, kg, Cg, kp, Cp, ϵ)

"""
    DoubleUTube(rb, ro, ri, s, kg, Cg, kp, Cp, network = :diagonal, ϵ = 5e-6)
    DoubleUTube(; rb, ro, ri, s, kg, Cg, kp, Cp, network = :diagonal, ϵ = 5e-6)

Grouted double U-tube borehole: four pipes placed symmetrically on a circle of radius `s/2` (pipes 1
to 4 at 0°, 90°, 180° and 270°), the two U-tubes being connected in parallel (each carries half of
the borehole flow rate).
- `network = :diagonal`: the inlet pipes are diagonally opposite (1 and 3 down, 2 and 4 up);
- `network = :adjacent`: the inlet pipes are adjacent (1 and 2 down, 3 and 4 up).
# Arguments
    - `rb`: Borehole radius [m]
    - `ro`: Pipe outer radius [m]
    - `ri`: Pipe inner radius [m]
    - `s`: Shank spacing (center-to-center distance between two opposite pipes) [m]
    - `kg`: Grout thermal conductivity [W/mK]
    - `Cg`: Grout volumetric heat capacity [J/m³K]
    - `kp`: Pipe thermal conductivity [W/mK]
    - `Cp`: Pipe volumetric heat capacity [J/m³K]
    - `network`: Flow configuration, `:diagonal` or `:adjacent` (default `:diagonal`)
    - `ϵ`: Pipe roughness [m] (default 5e-6 for HDPE pipes)
# Reference
    - Claesson, J., & Javed, S. (2019). Explicit multipole formulas and thermal network models for
        calculating thermal resistances of double U-pipe borehole heat exchangers. Science and
        Technology for the Built Environment, 25(8), 980–992.
        https://doi.org/10.1080/23744731.2019.1620565
"""
struct DoubleUTube <: AbstractExchanger
    rb::Float64
    ro::Float64
    ri::Float64
    s::Float64
    kg::Float64
    Cg::Float64
    kp::Float64
    Cp::Float64
    network::Symbol
    ϵ::Float64
    function DoubleUTube(rb, ro, ri, s, kg, Cg, kp, Cp, network = :diagonal, ϵ = 5e-6)
        all(>(0), (rb, ro, ri, s, kg, Cg, kp, Cp)) || throw(ArgumentError(
            "DoubleUTube: all dimensions and properties must be strictly positive"))
        ri < ro || throw(ArgumentError(
            "DoubleUTube: the pipe inner radius ri must be smaller than ro"))
        s / sqrt(2) > 2ro || throw(ArgumentError(
            "DoubleUTube: adjacent pipes overlap (s/√2 ≤ 2ro)"))
        s / 2 + ro < rb || throw(ArgumentError(
            "DoubleUTube: the pipes cross the borehole wall (s/2 + ro ≥ rb)"))
        network in (:diagonal, :adjacent) || throw(ArgumentError(
            "DoubleUTube: network must be :diagonal or :adjacent, got :$network"))
        return new(rb, ro, ri, s, kg, Cg, kp, Cp, network, ϵ)
    end
end
DoubleUTube(; rb, ro, ri, s, kg, Cg, kp, Cp, network = :diagonal, ϵ = 5e-6) =
    DoubleUTube(rb, ro, ri, s, kg, Cg, kp, Cp, network, ϵ)

"""
    Coaxial(rb, rii, rio, roi, roo, kg, Cg, kpi, Cpi, kpo, Cpo, inlet = :annulus, ϵ = 5e-6)
    Coaxial(; rb, rii, rio, roi, roo, kg, Cg, kpi, Cpi, kpo, Cpo, inlet = :annulus, ϵ = 5e-6)

Coaxial (concentric-tube) borehole: a center pipe inside an outer pipe, the outer pipe being grouted
in the borehole (`roo < rb`) or in direct contact with the ground (`roo == rb`).
- `inlet = :annulus`: the fluid goes down the annulus and comes back up the center pipe;
- `inlet = :center`: the fluid goes down the center pipe and comes back up the annulus.
# Arguments
    - `rb`: Borehole radius [m]
    - `rii`: Center pipe inner radius [m]
    - `rio`: Center pipe outer radius [m]
    - `roi`: Outer pipe inner radius [m]
    - `roo`: Outer pipe outer radius [m]
    - `kg`: Grout thermal conductivity [W/mK] (unused if `roo == rb`)
    - `Cg`: Grout volumetric heat capacity [J/m³K] (unused if `roo == rb`)
    - `kpi`: Center pipe thermal conductivity [W/mK]
    - `Cpi`: Center pipe volumetric heat capacity [J/m³K]
    - `kpo`: Outer pipe thermal conductivity [W/mK]
    - `Cpo`: Outer pipe volumetric heat capacity [J/m³K]
    - `inlet`: Flow direction, `:annulus` or `:center` (default `:annulus`)
    - `ϵ`: Pipe roughness [m] (default 5e-6)
# Reference
    - Bauer, D., Heidemann, W., Müller-Steinhagen, H., & Diersch, H.-J. G. (2011). Thermal
        resistance and capacity models for borehole heat exchangers. International Journal of
        Energy Research, 35(4), 312–320. https://doi.org/10.1002/er.1689
    - Lamarche, L. (2021). Analytic models and effective resistances for coaxial ground heat
        exchangers. Geothermics, 97, 102224. https://doi.org/10.1016/j.geothermics.2021.102224
"""
struct Coaxial <: AbstractExchanger
    rb::Float64
    rii::Float64
    rio::Float64
    roi::Float64
    roo::Float64
    kg::Float64
    Cg::Float64
    kpi::Float64
    Cpi::Float64
    kpo::Float64
    Cpo::Float64
    inlet::Symbol
    ϵ::Float64
    function Coaxial(rb, rii, rio, roi, roo, kg, Cg, kpi, Cpi, kpo, Cpo, inlet = :annulus, ϵ = 5e-6)
        all(>(0), (rb, rii, rio, roi, roo, kpi, Cpi, kpo, Cpo)) || throw(ArgumentError(
            "Coaxial: all dimensions and pipe properties must be strictly positive"))
        rii < rio < roi < roo <= rb || throw(ArgumentError(
            "Coaxial: radii must satisfy rii < rio < roi < roo ≤ rb"))
        roo == rb || (kg > 0 && Cg > 0) || throw(ArgumentError(
            "Coaxial: grout properties kg and Cg must be strictly positive when roo < rb"))
        inlet in (:annulus, :center) || throw(ArgumentError(
            "Coaxial: inlet must be :annulus or :center"))
        return new(rb, rii, rio, roi, roo, kg, Cg, kpi, Cpi, kpo, Cpo, inlet, ϵ)
    end
end
Coaxial(; rb, rii, rio, roi, roo, kg, Cg, kpi, Cpi, kpo, Cpo, inlet = :annulus, ϵ = 5e-6) =
    Coaxial(rb, rii, rio, roi, roo, kg, Cg, kpi, Cpi, kpo, Cpo, inlet, ϵ)

"""
    StandingColumn(rb, ro, ri, kp, Cp, inlet = :annulus, ϵ = 5e-6, ϵb = 5e-6)
    StandingColumn(; rb, ro, ri, kp, Cp, inlet = :annulus, ϵ = 5e-6, ϵb = 5e-6)

Standing column well (SCW) without bleed: an open, water-filled borehole with a central dip tube.
The groundwater in the annulus is in direct contact with the rock (no grout). As in Nguyen et al.
(2015), the water of each channel is perfectly mixed radially.
- `inlet = :annulus` (usual operation): the water is reinjected at the top of the annulus, flows
  down along the borehole wall and is pumped back up through the dip tube;
- `inlet = :center`: the reverse flow direction.
The groundwater flow induced by a bleed is not modeled (pure conduction in the rock).
# Arguments
    - `rb`: Borehole radius [m]
    - `ro`: Dip tube outer radius [m]
    - `ri`: Dip tube inner radius [m]
    - `kp`: Dip tube thermal conductivity [W/mK]
    - `Cp`: Dip tube volumetric heat capacity [J/m³K]
    - `inlet`: Flow direction, `:annulus` or `:center` (default `:annulus`)
    - `ϵ`: Dip tube roughness [m] (default 5e-6)
    - `ϵb`: Borehole wall roughness [m] (default 5e-6; a rough rock wall can be much higher)
# Reference
    - Nguyen, A., Pasquier, P., & Marcotte, D. (2015). Thermal resistance and capacity model for
        standing column wells operating under a bleed control. Renewable Energy, 76, 743–756.
        https://doi.org/10.1016/j.renene.2014.11.080
"""
struct StandingColumn <: AbstractExchanger
    rb::Float64
    ro::Float64
    ri::Float64
    kp::Float64
    Cp::Float64
    inlet::Symbol
    ϵ::Float64
    ϵb::Float64
    function StandingColumn(rb, ro, ri, kp, Cp, inlet = :annulus, ϵ = 5e-6, ϵb = 5e-6)
        all(>(0), (rb, ro, ri, kp, Cp)) || throw(ArgumentError(
            "StandingColumn: all dimensions and properties must be strictly positive"))
        ri < ro < rb || throw(ArgumentError("StandingColumn: radii must satisfy ri < ro < rb"))
        inlet in (:annulus, :center) || throw(ArgumentError(
            "StandingColumn: inlet must be :annulus or :center"))
        return new(rb, ro, ri, kp, Cp, inlet, ϵ, ϵb)
    end
end
StandingColumn(; rb, ro, ri, kp, Cp, inlet = :annulus, ϵ = 5e-6, ϵb = 5e-6) =
    StandingColumn(rb, ro, ri, kp, Cp, inlet, ϵ, ϵb)

"""
    Casing(L, e, kc = 30.0, Cc = 3.9e6)
    Casing(; L, e, kc = 30.0, Cc = 3.9e6)

Casing ring between the borehole wall (`rb`) and the ground (`rb + e`) over the top `L` meters of
the borehole, as in the quasi-3D model of Pasquier & Marcotte (2014). Below `L`, the same ring is
filled with ground, so that the soil discretization is identical in every layer.
# Arguments
    - `L`: Casing length from the ground surface [m]
    - `e`: Casing thickness [m]
    - `kc`: Casing thermal conductivity [W/mK] (default 30, steel)
    - `Cc`: Casing volumetric heat capacity [J/m³K] (default 3.9e6, steel)
"""
struct Casing
    L::Float64
    e::Float64
    kc::Float64
    Cc::Float64
    function Casing(L, e, kc = 30.0, Cc = 3.9e6)
        L >= 0 && e > 0 && kc > 0 && Cc > 0 || throw(ArgumentError(
            "Casing: L must be non-negative and e, kc, Cc strictly positive"))
        return new(L, e, kc, Cc)
    end
end
Casing(; L, e, kc = 30.0, Cc = 3.9e6) = Casing(L, e, kc, Cc)

"""
    FluidProperties(kf, cf, ρf, μf)
    FluidProperties(T, fluid = :water; percentage = 100)

Thermophysical properties of the heat carrier fluid. The second form evaluates them at the fluid
temperature `T` [°C] with `BoreholeResistance.fluid_property` (CoolProp), for pure water (`:water`)
or a water–antifreeze mixture (`:MPG`, `:MEG`, `:MMA`, `:MEA`, `:MKA`, `:MKF`).
# Arguments
    - `kf`: Fluid thermal conductivity [W/mK]
    - `cf`: Fluid specific heat capacity [J/kgK]
    - `ρf`: Fluid density [kg/m³]
    - `μf`: Fluid dynamic viscosity [kg/m⋅s]
    - `T`: Fluid temperature [°C]
    - `fluid`: Fluid symbol (default `:water`)
    - `percentage`: Mass fraction of the additive for mixtures [%m] (default 100)
"""
struct FluidProperties
    kf::Float64
    cf::Float64
    ρf::Float64
    μf::Float64
end
FluidProperties(T::Real, fluid::Symbol = :water; percentage::Real = 100) =
    FluidProperties(fluid_property(T, fluid; percentage)...)

"""
    TRCMSettings(; nz = nothing, np = 6, ng = 6, ngg = 6, nc = 6, ns = nothing, rs = nothing,
        nbelow = 20)

Discretization of the TRCM networks. The defaults give an accuracy of a few tenths of a percent on
the outlet transfer function: there is usually no need to change them. Only the fields to change are
given, e.g. `TRCMSettings(nz = 25)`; a field left to `nothing` takes the default of the model.
# Keywords
    - `nz`: Number of depth layers along the borehole (default `max(round(Int, H/4), 50)`)
    - `np`: Number of sub-resistances of each pipe wall (default 6)
    - `ng`: Number of sub-resistances of the grout between a pipe and the borehole wall (default 6)
    - `ngg`: Number of sub-resistances of each grout chain between pipes (default 6)
    - `nc`: Number of sub-resistances of the casing ring (default 6, only with a casing)
    - `ns`: Number of log-spaced soil shells (default 24, 100 with vertical conduction:
        [`trcm`](@ref) with a buried depth `D`, and [`trcm_pairs`](@ref))
    - `rs`: Far-field radius, where the ground is held at its initial temperature [m] (default 20,
        300 with vertical conduction)
    - `nbelow`: Number of soil layers below the borehole toe with vertical conduction (default 20)
"""
Base.@kwdef struct TRCMSettings
    nz::Union{Nothing,Int} = nothing
    np::Int = 6
    ng::Int = 6
    ngg::Int = 6
    nc::Int = 6
    ns::Union{Nothing,Int} = nothing
    rs::Union{Nothing,Float64} = nothing
    nbelow::Int = 20
    function TRCMSettings(nz, np, ng, ngg, nc, ns, rs, nbelow)
        something(nz, 1) >= 1 && min(np, ng, ngg, nc, nbelow) >= 1 && something(ns, 2) >= 2 &&
            something(rs, 1.0) > 0 || throw(ArgumentError(
            "TRCMSettings: nz, np, ng, ngg, nc, nbelow ≥ 1, ns ≥ 2 and rs > 0 are required"))
        return new(nz, np, ng, ngg, nc, ns, rs, nbelow)
    end
end

# The same settings with the defaults of a model filled in (`nz` from the borehole length `H`)
_resolve(s::TRCMSettings, H; ns, rs) = TRCMSettings(something(s.nz, max(round(Int, H / 4), 50)),
    s.np, s.ng, s.ngg, s.nc, something(s.ns, ns), something(s.rs, rs), s.nbelow)
