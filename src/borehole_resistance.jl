# Steady-state and transient borehole resistances of one cross-section. A TRCM must reduce exactly
# to the steady-state resistances of the borehole (Bauer et al., 2011, Section 2: "models that
# differ from the Delta-circuit still need to satisfy the Ra and Rb values"). They are computed with
# BoreholeResistance.jl, and the networks of `cross_section.jl` are sized to reproduce them.

"""
    delta_circuit(ex, ks, V̇, fluid; films = true)

Steady-state resistances [mK/W] of the cross-section of `ex`, i.e. the Δ-circuit that every layer of
the TRCM reduces to at steady state.
- `SingleUTube`: `Rb`, `Ra` from the first-order multipole method (Javed & Spitler, 2017), and the
  Δ-circuit `R1b = 2Rb` (pipe to borehole wall) and `R12` (pipe to pipe, `1/R12 = 1/Ra - 1/(4Rb)`);
- `DoubleUTube`: `Rb`, `Rda`, `Raa` (Claesson & Javed, 2019, Eqs. 14, 19, 23) and the Δ-circuit
  `R1b = 4Rb`, `R12` (adjacent pipes) and `R13` (opposite pipes) of their Eq. 26;
- `Coaxial`: `R12` (center fluid to annulus fluid) and `R1 = Rb` (annulus fluid to borehole wall),
  as in Lamarche (2021, Eqs. 1-2) and Bauer et al. (2011, Eqs. 4-14);
- `StandingColumn`: `R12` (dip tube water to annulus water) and `R1 = Rb` (annulus water to the
  rock, a convective film only, Nguyen et al., 2015).
The parts (`Rf…` convective films, `Rp…` pipe conduction, `Rg` grout) are also returned. A Δ-circuit
resistance can be negative (Bauer et al., 2011, Section 2.6) or infinite. Not exported.
# Arguments
    - `ex`: Ground heat exchanger
    - `ks`: Ground thermal conductivity [W/mK]
    - `V̇`: Total volumetric flow rate in the borehole [m³/s]
    - `fluid`: Heat carrier fluid ([`FluidProperties`](@ref))
# Keywords
    - `films`: Include the convective film resistances (default `true`); `false` reproduces the
        fixed pipe-wall temperature of Pasquier & Marcotte (2012)
# Output
    - A `NamedTuple` of resistances [mK/W]
# Reference
    - Bauer, D., Heidemann, W., Müller-Steinhagen, H., & Diersch, H.-J. G. (2011). Thermal
        resistance and capacity models for borehole heat exchangers. International Journal of
        Energy Research, 35(4), 312–320. https://doi.org/10.1002/er.1689
    - Claesson, J., & Javed, S. (2019). Explicit multipole formulas and thermal network models for
        calculating thermal resistances of double U-pipe borehole heat exchangers. Science and
        Technology for the Built Environment, 25(8), 980–992.
    - Lamarche, L. (2021). Analytic models and effective resistances for coaxial ground heat
        exchangers. Geothermics, 97, 102224.
"""
function delta_circuit(ex::SingleUTube, ks, V̇, fluid::FluidProperties; films::Bool = true)
    v = V̇ / (π * ex.ri^2)                                  # mean fluid speed in a leg [m/s]
    Rf = films ? resistance_fluid(v, ex.ri, fluid.kf, fluid.cf, fluid.ρf, fluid.μf, ex.ϵ) : 0.0
    Rp = resistance_pipe(ex.ro, ex.ri, ex.kp)
    Rb = resistance_ULoop_borehole(ex.s, ex.rb, ex.ro, ks, ex.kg, Rp, Rf; nLoop = 1)
    Ra = resistance_ULoop_total_internal(ex.s, ex.rb, ex.ro, ks, ex.kg, Rp, Rf; nLoop = 1)
    # Δ-circuit: Ra = R12 ∥ (R1b + R2b) with R1b = R2b = 2Rb
    R12 = 1 / (1 / Ra - 1 / (4Rb))
    return (; Rf, Rp, Rb, Ra, R1b = 2Rb, R12)
end

function delta_circuit(ex::DoubleUTube, ks, V̇, fluid::FluidProperties; films::Bool = true)
    v = V̇ / 2 / (π * ex.ri^2)                              # each U-tube carries half the flow
    Rf = films ? resistance_fluid(v, ex.ri, fluid.kf, fluid.cf, fluid.ρf, fluid.μf, ex.ϵ) : 0.0
    Rp = resistance_pipe(ex.ro, ex.ri, ex.kp)
    Rb = resistance_ULoop_borehole(ex.s, ex.rb, ex.ro, ks, ex.kg, Rp, Rf; nLoop = 2)
    Rda = resistance_ULoop_total_internal(ex.s, ex.rb, ex.ro, ks, ex.kg, Rp, Rf; nLoop = 2,
        network = "diagonal")
    Raa = resistance_ULoop_total_internal(ex.s, ex.rb, ex.ro, ks, ex.kg, Rp, Rf; nLoop = 2,
        network = "adjacent")
    # Δ-circuit of Claesson & Javed (2019, Eq. 26): R1b = 4Rb, adjacent (R12) and opposite (R13)
    # pipes
    R12 = 1 / (1 / (2Rda) - 1 / (16Rb))
    R13 = 1 / (1 / Raa - 1 / (2Rda) - 1 / (16Rb))
    return (; Rf, Rp, Rb, Rda, Raa, R1b = 4Rb, R12, R13)
end

function delta_circuit(ex::Coaxial, ks, V̇, fluid::FluidProperties; films::Bool = true)
    Rpi = resistance_pipe(ex.rio, ex.rii, ex.kpi)
    Rpo = resistance_pipe(ex.roo, ex.roi, ex.kpo)
    Rg = ex.roo < ex.rb ? resistance_pipe(ex.rb, ex.roo, ex.kg) : 0.0
    Rf_in, Rf_ann_in, Rf_ann_out = 0.0, 0.0, 0.0
    if films
        # Center pipe (hydraulic diameter 2rii) and annulus (hydraulic diameter 2(roi - rio)), as in
        # BoreholeResistance.resistance_coaxial
        (; kf, cf, ρf, μf) = fluid
        v_in = V̇ / (π * ex.rii^2)
        h_in = convection_coefficient(Nusselt(v_in, ex.rii, kf, cf, ρf, μf, ex.ϵ), ex.rii, kf)
        v_ann = V̇ / (π * (ex.roi^2 - ex.rio^2))
        Nu_ann = Nusselt_annulus(v_ann, ex.roi, ex.rio, kf, cf, ρf, μf, ex.ϵ, ex.ϵ)
        h_ann = convection_coefficient(Nu_ann, ex.roi - ex.rio, kf)
        Rf_in = 1 / (2π * ex.rii * h_in)
        Rf_ann_in = 1 / (2π * ex.rio * h_ann)
        Rf_ann_out = 1 / (2π * ex.roi * h_ann)
    end
    R12 = Rf_in + Rpi + Rf_ann_in
    R1 = Rf_ann_out + Rpo + Rg
    return (; Rf_in, Rf_ann_in, Rf_ann_out, Rpi, Rpo, Rg, R12, R1, Rb = R1)
end

function delta_circuit(ex::StandingColumn, ks, V̇, fluid::FluidProperties; films::Bool = true)
    Rp = resistance_pipe(ex.ro, ex.ri, ex.kp)
    # Without films the annulus water would sit exactly at the rock temperature: a vanishing (but
    # non-zero) wall resistance keeps the two nodes distinct
    Rf_in, Rf_ann_in, Rf_wall = 0.0, 0.0, 1e-6
    if films
        (; kf, cf, ρf, μf) = fluid
        v_in = V̇ / (π * ex.ri^2)
        h_in = convection_coefficient(Nusselt(v_in, ex.ri, kf, cf, ρf, μf, ex.ϵ), ex.ri, kf)
        v_ann = V̇ / (π * (ex.rb^2 - ex.ro^2))
        Nu_ann = Nusselt_annulus(v_ann, ex.rb, ex.ro, kf, cf, ρf, μf, ex.ϵb, ex.ϵ)
        h_ann = convection_coefficient(Nu_ann, ex.rb - ex.ro, kf)
        Rf_in = 1 / (2π * ex.ri * h_in)
        Rf_ann_in = 1 / (2π * ex.ro * h_ann)
        Rf_wall = 1 / (2π * ex.rb * h_ann)
    end
    R12 = Rf_in + Rp + Rf_ann_in
    return (; Rf_in, Rf_ann_in, Rf_wall, Rp, R12, R1 = Rf_wall, Rb = Rf_wall)
end

"""
    transient_resistance(t, ex, ks, Cs, V̇, fluid; films = true, settings = TRCMSettings())

Transient borehole thermal resistance `Rb(t)` of one cross-section of the TRCM (no depth), the
verification setting of Bauer et al. (2011) and Pasquier & Marcotte (2012, Eq. 21): the network
starts at 0, the fluid temperatures are held at 1 from `t = 0` and the far field (`rs`) at 0. With
`q(t)` the heat rate leaving the fluids and `Tb(t)` the borehole wall temperature,

    Rb(t) = (1 - Tb(t)) / q(t),

which tends to the steady-state `Rb` of
[`delta_circuit`](@ref ThermalResistanceCapacitanceModel.delta_circuit) once the capacities of the
borehole are saturated. The departure from that limit at short times measures the thermal capacity
of the borehole. Not exported.
# Arguments
    - `t`: Time vector, strictly positive and increasing [s]
    - `ex`: Ground heat exchanger
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric heat capacity [J/m³K]
    - `V̇`: Total volumetric flow rate in the borehole [m³/s] (only used for the films)
    - `fluid`: Heat carrier fluid ([`FluidProperties`](@ref))
# Keywords
    - `films`: Include the convective film resistances (default `true`). Pasquier & Marcotte (2012)
        impose the temperature of the pipe inner wall, i.e. `films = false`.
    - `settings`: Discretization of the cross-section ([`TRCMSettings`](@ref))
# Output
    - `Rb`: Transient borehole thermal resistance at `t` [mK/W]
# Reference
    - Pasquier, P., & Marcotte, D. (2012). Short-term simulation of ground heat exchanger with an
        improved TRCM. Renewable Energy, 46, 92–99. https://doi.org/10.1016/j.renene.2012.03.014
"""
function transient_resistance(t, ex::AbstractExchanger, ks, Cs, V̇, fluid::FluidProperties;
    films::Bool = true, settings::TRCMSettings = TRCMSettings())

    tv = _check_times(t)
    s = _resolve(settings, 1.0; ns = 24, rs = 20.0)
    ch = _channels(ex, V̇)
    lay = _layer(ex, delta_circuit(ex, ks, V̇, fluid; films), ch, fluid.ρf * fluid.cf, ks, Cs, s)

    # The fluid nodes (1 to nf) are held at 1: the other nodes form the system
    A = _laplacian(lay.net)
    nf = length(ch.area)
    free = (nf + 1):size(A, 1)
    C = lay.net.C[free]
    M = Diagonal(1 ./ C) * A[free, free]
    b = A[free, 1:nf] * ones(nf) ./ C
    sol = _solve(M, b, tv[end], _TRCM_SOLVER; saveat = tv, reltol = 1e-8, abstol = 1e-10)

    # Heat rate leaving the fluids, q = -Σ_fluid rows (A T), and wall node Tb (local node nf + 1)
    q(u) = -sum(A[1:nf, 1:nf] * ones(nf) + A[1:nf, free] * u)
    return [(1 - u[lay.Tb - nf]) / q(u) for u in sol.u]
end
