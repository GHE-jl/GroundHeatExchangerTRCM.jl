using Test
using LinearAlgebra
using ThermalResistanceCapacitanceModel
using BoreholeResistance: resistance_coaxial

const TRCM = ThermalResistanceCapacitanceModel
const delta_circuit = TRCM.delta_circuit

# Shared physical parameters
const ks = 2.5          # ground thermal conductivity [W/mK]
const Cs = 2.2e6        # ground volumetric heat capacity [J/m³K]
const H = 100.0         # borehole length [m]
const V̇ = 0.4e-3        # flow rate [m³/s]
const fluid = FluidProperties(0.6, 4180.0, 1000.0, 1.0e-3)   # water-like, no CoolProp call
const Cf = fluid.ρf * fluid.cf

const single = SingleUTube(rb = 0.075, ro = 0.0167, ri = 0.0137, s = 0.06, kg = 1.5, Cg = 2.0e6,
    kp = 0.4, Cp = 1.8e6)
const double_d = DoubleUTube(rb = 0.075, ro = 0.0167, ri = 0.0137, s = 0.09, kg = 1.5, Cg = 2.0e6,
    kp = 0.4, Cp = 1.8e6, network = :diagonal)
const double_a = DoubleUTube(0.075, 0.0167, 0.0137, 0.09, 1.5, 2.0e6, 0.4, 1.8e6, :adjacent)
const coax = Coaxial(rb = 0.075, rii = 0.0177, rio = 0.0201, roi = 0.0545, roo = 0.0572, kg = 1.5,
    Cg = 2.0e6, kpi = 0.4, Cpi = 1.8e6, kpo = 0.4, Cpo = 1.8e6)
const coax_bare = Coaxial(0.0572, 0.0177, 0.0201, 0.0545, 0.0572, 0.0, 0.0, 0.4, 1.8e6, 0.4,
    1.8e6, :center)
const scw = StandingColumn(rb = 0.076, ro = 0.05, ri = 0.045, kp = 0.4, Cp = 1.8e6)
const exchangers = (single, double_d, double_a, coax, coax_bare, scw)

# Cheap discretization for the transient tests
const small = TRCMSettings(nz = 10, np = 3, ng = 3, ngg = 3, ns = 12)

"One layer, and its steady conductance matrix seen from the fluid nodes and the wall node."
function boundary_conductance(ex; ng = 6, np = 6, ngg = 6, films = true)
    s = TRCMSettings(nz = 1, np = np, ng = ng, ngg = ngg, nc = 1, ns = 12, rs = 20.0)
    ch = TRCM._channels(ex, V̇)
    lay = TRCM._layer(ex, delta_circuit(ex, ks, V̇, fluid; films), ch, Cf, ks, Cs, s)
    A = Matrix(TRCM._laplacian(lay.net))
    bnd = vcat(eachindex(ch.area), lay.Tb)
    inn = setdiff(1:size(A, 1), bnd)
    S = A[bnd, bnd] - A[bnd, inn] * (A[inn, inn] \ A[inn, bnd])
    return -S, lay   # -S: heat leaving each boundary node for unit temperatures
end

"Steady state of the TRCM system and the heat rate per unit length reaching the far field."
function steady_far_field(sys)
    u = -(sys.M \ sys.b)
    net = sys.layer.net
    q = sum(net.Gg[e] * u[(k - 1) * sys.nloc + net.Ig[e]] for k in 1:sys.nz, e in eachindex(net.Gg))
    return u, q / sys.nz
end

@testset "ThermalResistanceCapacitanceModel.jl" begin

@testset "Exchanger construction and validation" begin
    @test single == SingleUTube(0.075, 0.0167, 0.0137, 0.06, 1.5, 2.0e6, 0.4, 1.8e6)
    @test double_d.network === :diagonal && double_a.network === :adjacent
    @test coax.inlet === :annulus && coax_bare.inlet === :center && scw.inlet === :annulus
    @test Casing(50.0, 0.01) == Casing(L = 50.0, e = 0.01)
    # Pipes crossing the borehole wall, overlapping pipes, wrong symbols, wrong radii
    @test_throws ArgumentError SingleUTube(0.075, 0.0167, 0.0137, 0.12, 1.5, 2e6, 0.4, 1.8e6)
    @test_throws ArgumentError SingleUTube(0.075, 0.0167, 0.0137, 0.03, 1.5, 2e6, 0.4, 1.8e6)
    @test_throws ArgumentError DoubleUTube(0.075, 0.0167, 0.0137, 0.09, 1.5, 2e6, 0.4, 1.8e6,
        :crossed)
    @test_throws ArgumentError Coaxial(0.075, 0.03, 0.02, 0.05, 0.06, 1.5, 2e6, 0.4, 1.8e6, 0.4,
        1.8e6)
    @test_throws ArgumentError StandingColumn(0.076, 0.08, 0.045, 0.4, 1.8e6)
    @test_throws ArgumentError Casing(10.0, 0.0)
    @test_throws ArgumentError TRCMSettings(ns = 1)
    @test_throws ArgumentError trcm([1.0, 0.5], single, ks, Cs, H, V̇, fluid)
    @test_throws ArgumentError trcm([0.0, 1.0], single, ks, Cs, H, V̇, fluid)
    @test_throws ArgumentError trcm([1.0], single, ks, Cs, H, V̇, fluid;
        settings = TRCMSettings(rs = 0.05))
    # Model defaults: nz from the borehole length, ns and rs from the model
    s = TRCM._resolve(TRCMSettings(rs = 30.0), 100.0; ns = 24, rs = 20.0)
    @test s.nz == 50 && s.ns == 24 && s.rs == 30.0
    @test TRCM._resolve(TRCMSettings(), 400.0; ns = 24, rs = 20.0).nz == 100
end

@testset "FluidProperties from CoolProp" begin
    w = FluidProperties(20.0)
    @test 990 < w.ρf < 1000 && 4100 < w.cf < 4250 && 0.55 < w.kf < 0.65 && 5e-4 < w.μf < 2e-3
end

@testset "Δ-circuit inputs" begin
    d1 = delta_circuit(single, ks, V̇, fluid)
    @test d1.R1b ≈ 2 * d1.Rb
    @test 1 / d1.Ra ≈ 1 / d1.R12 + 1 / (2 * d1.R1b)
    @test d1.Rf > 0
    @test delta_circuit(single, ks, V̇, fluid; films = false).Rf == 0.0

    d2 = delta_circuit(double_d, ks, V̇, fluid)
    @test d2.R1b ≈ 4 * d2.Rb
    @test 1 / d2.Rda ≈ 2 / d2.R12 + 1 / (2 * d2.R1b)                 # Claesson & Javed, Eq. 25
    @test 1 / d2.Raa ≈ 1 / d2.R12 + 1 / d2.R13 + 1 / (2 * d2.R1b)

    # Coaxial: identical to BoreholeResistance.resistance_coaxial (Lamarche, 2021)
    dc = delta_circuit(coax, ks, V̇, fluid)
    R1, R12 = resistance_coaxial(V̇, coax.rii, coax.rio, coax.roi, coax.roo, coax.rb, coax.kg,
        coax.kpi, coax.kpo, fluid.kf, fluid.cf, fluid.ρf, fluid.μf, coax.ϵ)
    @test dc.R1 ≈ R1 && dc.R12 ≈ R12 && dc.Rb == dc.R1
    @test delta_circuit(coax_bare, ks, V̇, fluid).Rg == 0.0

    ds = delta_circuit(scw, ks, V̇, fluid)
    @test ds.R1 == ds.Rf_wall > 0 && ds.R12 > ds.Rp
end

@testset "Grout core conductances (Bauer et al., 2011, Eq. 21)" begin
    d = delta_circuit(single, ks, V̇, fluid)
    Rfp, Rg = d.Rf + d.Rp, d.R1b - d.Rf - d.Rp
    Gd = [0 1/d.R12; 1/d.R12 0]
    for x in (0.0, 1 / 6, 0.5)
        Rgb, Rar = (1 - x) * Rg, d.Ra - 2Rfp
        Rgg_bauer = 2Rgb * (Rar - 2x * Rg) / (2Rgb - Rar + 2x * Rg)
        G = TRCM._core_conductances(d.R1b, Gd, Rfp + x * Rg, (1 - x) * Rg)
        @test 1 / G[1, 2] ≈ Rgg_bauer rtol = 1e-10
    end
end

@testset "Steady state reproduces the Δ-circuit exactly" begin
    for ng in (1, 2, 6)
        # Single U-tube: Rb (all fluids at 1, wall at 0) and Ra (antisymmetric fluids)
        S, _ = boundary_conductance(single; ng)
        d = delta_circuit(single, ks, V̇, fluid)
        @test 1 / sum(S[1:2, :] * [1.0, 1.0, 0.0]) ≈ d.Rb rtol = 1e-10
        @test 1 / (S[1, :] ⋅ [0.5, -0.5, 0.0]) ≈ d.Ra rtol = 1e-10

        # Double U-tube: Rb, Rda (diagonal inlets) and Raa (adjacent inlets)
        S, _ = boundary_conductance(double_d; ng)
        d = delta_circuit(double_d, ks, V̇, fluid)
        @test 1 / sum(S[1:4, :] * [1.0, 1.0, 1.0, 1.0, 0.0]) ≈ d.Rb rtol = 1e-10
        @test 2 / (S[1, :] ⋅ [1.0, -1.0, 1.0, -1.0, 0.0]) ≈ d.Rda rtol = 1e-10
        @test 2 / (S[1, :] ⋅ [1.0, 1.0, -1.0, -1.0, 0.0]) ≈ d.Raa rtol = 1e-10
    end
    for ex in (coax, coax_bare, scw)
        S, _ = boundary_conductance(ex)
        d = delta_circuit(ex, ks, V̇, fluid)
        @test 1 / sum(S[1:2, :] * [1.0, 1.0, 0.0]) ≈ d.R1 rtol = 1e-10    # both fluids at 1
        @test 1 / (S[1, :] ⋅ [1.0, 0.0, 0.0]) ≈ d.R12 rtol = 1e-10       # center to annulus
    end
end

@testset "Capacity conservation" begin
    rs, ns = 20.0, 12
    vlast = TRCM._log_shells(single.rb, rs, ns)[end]                   # far-field half shell
    soil = Cs * (π * (rs^2 - single.rb^2) - vlast / 2)
    for ex in (single, double_d)
        n = ex isa SingleUTube ? 2 : 4
        _, lay = boundary_conductance(ex)
        expected = n * Cf * π * ex.ri^2 + n * ex.Cp * π * (ex.ro^2 - ex.ri^2) +
                   ex.Cg * π * (ex.rb^2 - n * ex.ro^2) + soil
        @test sum(lay.net.C) ≈ expected rtol = 1e-12
    end
    _, lay = boundary_conductance(coax)
    expected = Cf * π * (coax.rii^2 + coax.roi^2 - coax.rio^2) +
               coax.Cpi * π * (coax.rio^2 - coax.rii^2) + coax.Cpo * π * (coax.roo^2 - coax.roi^2) +
               coax.Cg * π * (coax.rb^2 - coax.roo^2) + soil
    @test sum(lay.net.C) ≈ expected rtol = 1e-12
end

@testset "Response model: energy balance and boundary conditions" begin
    t = 10.0 .^ range(1, 7; length = 25)
    for ex in exchangers
        sys = TRCM._trcm_system(ex, ks, Cs, H, V̇, fluid, nothing,
            TRCMSettings(ns = 12, np = 3, ng = 3, ngg = 3, nc = 1, nz = 10))
        # At steady state, all the injected heat V̇·Cf·(Tin - Tout) reaches the far field
        u, q = steady_far_field(sys)
        @test q ≈ V̇ * Cf / H rtol = 1e-8

        g = trcm(t, ex, ks, Cs, H, V̇, fluid; settings = small)
        @test length(g) == length(t) && all(isfinite, g)
        @test all(diff(g) .> -1e-6)                                     # heating response
        @test all(g .>= 0)                                              # negatives set to 0
        @test all(g .<= sys.w ⋅ u[sys.outlet] + 1e-6)                   # below the steady state
    end
end

@testset "Outlet unchanged before the fluid residence time" begin
    # Pasquier & Marcotte (2014, Fig. 5): Tout starts deviating from zero after the residence time
    τ = 2H * π * single.ri^2 / V̇
    g = trcm([0.1τ, 0.3τ, 2τ], single, ks, Cs, H, V̇, fluid; settings = small)
    @test g[1] < 1e-6 && g[2] < 1e-3
    @test g[3] > 0.1
end

@testset "Casing" begin
    t = [1e4, 1e5, 1e6]
    # Same discretization, ring of ground (reference) or of steel (conductive and capacitive)
    g0 = trcm(t, single, ks, Cs, H, V̇, fluid; settings = small,
        casing = Casing(H, 0.01, ks, Cs))
    gc = trcm(t, single, ks, Cs, H, V̇, fluid; settings = small, casing = Casing(H, 0.01))
    gh = trcm(t, single, ks, Cs, H, V̇, fluid; settings = small,
        casing = Casing(H / 2, 0.01))
    @test all(gc .< g0)                         # a steel ring lowers the response
    @test all(gc .< gh .< g0)                   # half-cased lies in between
    g_none = trcm(t, single, ks, Cs, H, V̇, fluid; settings = small)
    @test g0 ≈ g_none rtol = 5e-2               # a ground ring changes only the discretization
    sys = TRCM._trcm_system(single, ks, Cs, H, V̇, fluid, Casing(H / 2, 0.01),
        TRCMSettings(ns = 12, np = 3, ng = 3, ngg = 3, nc = 2, nz = 10))
    _, q = steady_far_field(sys)
    @test q ≈ V̇ * Cf / H rtol = 1e-8
end

@testset "Transient borehole resistance" begin
    t = [10.0, 600.0, 3600.0, 1e8]
    for ex in (single, double_d, coax)
        Rb_t = TRCM.transient_resistance(t, ex, ks, Cs, V̇, fluid;
            settings = TRCMSettings(rs = 1.0, ns = 8, np = 3, ng = 3, ngg = 3))
        @test Rb_t[end] ≈ delta_circuit(ex, ks, V̇, fluid).Rb rtol = 1e-4
        @test all(Rb_t .> 0)
    end
end

@testset "Unstable cross-section is rejected" begin
    net = TRCM._Network()
    TRCM._add_node!(net, 1.0)
    TRCM._add_node!(net, 1.0)
    TRCM._connect!(net, 1, 2, -0.1)                 # strongly negative resistance
    TRCM._connect!(net, 1, nothing, 1.0)
    TRCM._connect!(net, 2, nothing, 1.0)
    @test_throws ArgumentError TRCM._check_stability(net)
end

@testset "Borehole-to-borehole response functions (Rose et al., 2026)" begin
    yr = 365.25 * 86400
    t = [3600.0, 86400.0, 30 * 86400.0, yr, 10yr]
    set = TRCMSettings(nz = 10, np = 3, ng = 3, ngg = 3, ns = 40, rs = 100.0, nbelow = 8)
    g_self, g_pair = trcm_pairs(t, single, ks, Cs, H, 4.0, V̇, fluid, [3.0, 6.0, 12.0];
        settings = set)
    @test size(g_pair) == (5, 3) && length(g_self) == 5
    @test all(g_pair .>= 0) && all(g_pair .< g_self)
    @test all(diff(g_pair[end, :]) .< 0)                      # weaker with distance
    @test all(diff(g_self) .> 0) && all(diff(g_pair[:, 1]) .>= 0)
    # Before the ground reacts at depth, the emitter equals the TRCM without vertical conduction
    g0 = trcm(t[1:2], single, ks, Cs, H, V̇, fluid;
        settings = TRCMSettings(nz = 10, np = 3, ng = 3, ngg = 3, ns = 40, rs = 100.0))
    @test g_self[1:2] ≈ g0 rtol = 1e-3
    # Scalar distance: a vector response
    r1 = trcm_pairs(t[1:2], single, ks, Cs, H, 4.0, V̇, fluid, 6.0; settings = set)
    @test r1.g_pair isa Vector && length(r1.g_pair) == 2
    @test_throws ArgumentError trcm_pairs(t, single, ks, Cs, H, 4.0, V̇, fluid, 150.0;
        settings = set)
end

@testset "Vertical conduction in the ground (buried depth D)" begin
    yr = 365.25 * 86400
    t = [3600.0, 86400.0, 30 * 86400.0, yr, 10yr]
    set = TRCMSettings(nz = 10, np = 3, ng = 3, ngg = 3, ns = 40, rs = 100.0, nbelow = 8)
    # The emitter of the borehole-to-borehole responses
    g = trcm(t, single, ks, Cs, H, V̇, fluid; D = 4.0, settings = set)
    @test g ≈ trcm_pairs(t, single, ks, Cs, H, 4.0, V̇, fluid, 6.0; settings = set).g_self rtol = 1e-3
    # Every exchanger: increasing, and the short-term model before the ground reacts at depth
    for ex in exchangers
        gv = trcm(t, ex, ks, Cs, H, V̇, fluid; D = 4.0, settings = set)
        g0 = trcm(t[1:2], ex, ks, Cs, H, V̇, fluid;
            settings = TRCMSettings(nz = 10, np = 3, ng = 3, ngg = 3, ns = 40, rs = 100.0))
        @test all(isfinite, gv) && all(diff(gv) .> 0)
        @test gv[1:2] ≈ g0 rtol = 1e-3
    end
    # Long term: below the short-term model without vertical conduction and the same rs (no ground
    # surface, no finite length there)
    @test g[end] < trcm(t, single, ks, Cs, H, V̇, fluid;
        settings = TRCMSettings(nz = 10, np = 3, ng = 3, ngg = 3, ns = 40, rs = 100.0))[end]
    @test_throws ArgumentError trcm(t, single, ks, Cs, H, V̇, fluid; D = 4.0,
        casing = Casing(L = 10.0, e = 0.01))
    @test_throws MethodError trcm_pairs(t, coax, ks, Cs, H, 4.0, V̇, fluid, 6.0)
end

end
