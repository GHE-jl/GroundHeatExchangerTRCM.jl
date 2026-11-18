# Quasi-3D response model (Pasquier & Marcotte, 2014): `nz` copies of the cross-section of
# `cross_section.jl`, stacked along the depth and coupled only by the advection of the fluid, under
# the response-model boundary conditions:
#   - T(t = 0) = 0 everywhere (Eq. 14) and far field held at 0 in every layer (Eq. 15);
#   - Tin(t) = Tout(t) + 1 at the borehole head (Eqs. 13 and 16), imposed as the upstream
#     temperature of the inlet channels.
# The system is linear with constant coefficients, du/dt = M u + b. Its exact sparse Jacobian M is
# handed to Rodas5P with the KLU sparse LU: every layer is a radial chain, so the LU has almost no
# fill-in. Rodas5P + KLU was the fastest of the solvers compared at equal accuracy.
# With a buried depth `D`, the layers are the TRCM columns of `trcm_pairs.jl` instead: vertical
# conduction in the ground, buried depth, ground surface and ground below the toe (long term).

const _TRCM_SOLVER = Rodas5P(linsolve = KLUFactorization())

"""
    trcm(t, ex, ks, Cs, H, V̇, fluid; D = nothing, casing = nothing, settings = TRCMSettings())

Outlet transfer function of a ground heat exchanger from the quasi-3D thermal resistance and
capacity model (TRCM) of Pasquier & Marcotte (2014), solved as a response model: the inlet
temperature is held at `Tin = Tout + 1` from `t = 0`, so that the outlet temperature is the response
to a constant heat injection rate `V̇·Cf`, including the fluid transit, the capacities of the fluid,
pipes and grout, and the thermal short-circuit between the pipes along the depth.

Each of the `nz` depth layers is a TRCM cross-section (Pasquier & Marcotte, 2012) whose
sub-resistances reproduce the steady-state Δ-circuit of the exchanger exactly, followed by a radial
soil chain to the far-field radius `rs`.
- `D = nothing` (default, short term): the layers are coupled by the fluid only, the ground has no
    vertical conduction and is held at its initial temperature at `rs` (20 m). The response is valid
    from seconds to a few months; it levels off once the far field is felt. Joined to a long-term
    ground model by `outlet_transfer_function` of GroundHeatExchanger.jl.
- `D` given (long term): the ground conducts heat vertically between the layers, with soil layers
    over the buried depth `D` (ground surface held at 0) and below the borehole toe, and `rs`
    defaults to 300 m with 100 soil shells. The response is valid from seconds to decades, at a
    cost 15 to 30 times higher; it is the emitter of [`trcm_pairs`](@ref). The casing is not
    available in this mode.
# Arguments
    - `t`: Time vector, strictly positive and increasing [s]
    - `ex`: Ground heat exchanger ([`SingleUTube`](@ref), [`DoubleUTube`](@ref),
        [`Coaxial`](@ref) or [`StandingColumn`](@ref))
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric heat capacity [J/m³K]
    - `H`: Borehole length [m]
    - `V̇`: Total volumetric flow rate in the borehole [m³/s]
    - `fluid`: Heat carrier fluid ([`FluidProperties`](@ref))
# Keywords
    - `D`: Buried depth of the borehole head [m], or `nothing` (default) for the short-term model
        without vertical conduction
    - `casing`: [`Casing`](@ref) of the top of the borehole, or `nothing` (default; short-term
        model only)
    - `settings`: Discretization ([`TRCMSettings`](@ref))
# Output
    - `g`: Outlet transfer function at `t` [-] (add 1 for the inlet transfer function)
# Reference
    - Pasquier, P., & Marcotte, D. (2012). Short-term simulation of ground heat exchanger with an
        improved TRCM. Renewable Energy, 46, 92–99. https://doi.org/10.1016/j.renene.2012.03.014
    - Pasquier, P., & Marcotte, D. (2014). Joint use of quasi-3D response model and spectral method
        to simulate borehole heat exchanger. Geothermics, 51, 281–299.
        https://doi.org/10.1016/j.geothermics.2014.02.001
"""
function trcm(t, ex::AbstractExchanger, ks, Cs, H, V̇, fluid::FluidProperties;
    D::Union{Nothing,Real} = nothing, casing::Union{Nothing,Casing} = nothing,
    settings::TRCMSettings = TRCMSettings())

    tv = _check_times(t)

    # Long term: the emitter column with vertical conduction in the ground
    if D !== nothing
        casing === nothing || throw(ArgumentError(
            "the casing is not available with vertical conduction (D given)"))
        sys = _column_system(ex, ks, Cs, H, D, V̇, fluid, Float64[], settings)
        sol = _solve(sys.M, sys.b, tv[end], _PAIR_SOLVER; saveat = tv, save_idxs = sys.outlets,
            reltol = 1e-4, abstol = 1e-6)
        return [sys.outlet(u, 0) for u in sol.u]
    end

    # Short term: layers coupled by the fluid only
    sys = _trcm_system(ex, ks, Cs, H, V̇, fluid, casing, settings)
    sol = _solve(sys.M, sys.b, tv[end], _TRCM_SOLVER; saveat = tv, save_idxs = sys.outlet,
        reltol = 1e-4, abstol = 1e-6)
    # A heating response cannot be negative, but a negative Δ-circuit resistance (e.g. between the
    # opposite pipes of a double U-tube with adjacent inlets) can produce a tiny undershoot of the
    # outlet before the fluid residence time: it is set to 0
    return [max(0.0, sys.w ⋅ u) for u in sol.u]
end

"""
    _trcm_system(ex, ks, Cs, H, V̇, fluid, casing, settings)

Assemble the linear system `du/dt = M u + b` of the quasi-3D TRCM. The conduction matrix of the
layer (with or without the casing material) is replicated `nz` times; the advection is added in
energy units, W/(m⋅K), *before* dividing by the node capacities, so that a fluid node carrying part
of the pipe capacity still exchanges exactly `ṁ·cf·ΔT` with its neighbour. Returns `M`, `b`, the
outlet nodes and their mixing weights `w`, and the layer of the network.
"""
function _trcm_system(ex::AbstractExchanger, ks, Cs, H, V̇, fluid::FluidProperties, casing,
    settings)
    ks > 0 && Cs > 0 && H > 0 && V̇ > 0 || throw(ArgumentError(
        "ks, Cs, H and V̇ must be strictly positive"))
    s = _resolve(settings, H; ns = 24, rs = 20.0)
    r0 = casing === nothing ? ex.rb : ex.rb + casing.e
    s.rs > r0 || throw(ArgumentError(
        "the far-field radius rs = $(s.rs) m must exceed the borehole (and casing) radius $r0 m"))

    Cf = fluid.ρf * fluid.cf
    δ = delta_circuit(ex, ks, V̇, fluid)
    ch = _channels(ex, V̇)
    nz, dz = s.nz, H / s.nz

    # Layers below the casing (ring of ground) and along it (ring of casing material)
    if casing === nothing
        lay = _layer(ex, δ, ch, Cf, ks, Cs, s)
        layer = fill(lay, nz)
    else
        lay = _layer(ex, δ, ch, Cf, ks, Cs, s; ring = (; casing.e, k = ks, C = Cs))
        cased = _layer(ex, δ, ch, Cf, ks, Cs, s; ring = (; casing.e, k = casing.kc, C = casing.Cc))
        layer = [(k - 1) * dz < casing.L ? cased : lay for k in 1:nz]
    end

    # Conduction: block diagonal, one block per layer
    nloc = length(lay.net.C)
    N = nz * nloc
    rows, cols, vals = Int[], Int[], Float64[]
    for k in 1:nz
        _append_conduction!(rows, cols, vals, layer[k].net, (k - 1) * nloc)
    end
    C = reduce(vcat, [layer[k].net.C for k in 1:nz])

    # Advection (first-order upwind) in energy units, ṁ·cf/dz = V̇c·Cf/dz per unit length. The fluid
    # of channel c is the local node c of every layer.
    node(k, c) = (k - 1) * nloc + c
    up = findall(!, ch.down)
    w = ch.flow[up] ./ sum(ch.flow[up])                     # mixing weights of the outlet channels
    bE = zeros(N)
    for k in 1:nz, c in eachindex(ch.flow)
        i = node(k, c)
        a = ch.flow[c] * Cf / dz
        push!(rows, i); push!(cols, i); push!(vals, -a)
        if ch.down[c] && k == 1                             # inlet: Tin = Tout + 1 (Eqs. 13, 16)
            append!(rows, fill(i, length(up))); append!(cols, node.(1, up)); append!(vals, a .* w)
            bE[i] += a
        else                                                # upstream node: above for a down-going
            j = ch.down[c] ? node(k - 1, c) :               # channel, below for an up-going one,
                k == nz ? node(nz, ch.pair[c]) : node(k + 1, c)   # fed at the toe by its pair
            push!(rows, i); push!(cols, j); push!(vals, a)
        end
    end

    M = Diagonal(1 ./ C) * sparse(rows, cols, vals, N, N)
    return (; M, b = bE ./ C, outlet = node.(1, up), w, nz, nloc, layer = lay)
end
