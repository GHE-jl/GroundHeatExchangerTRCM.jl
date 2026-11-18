# Single-borehole and borehole-to-borehole response functions measured at the borehole outlet
# (Rose et al., 2026), and the TRCM columns with vertical conduction in the ground that compute them.
#
# Rose et al. (2026) simulate two complete boreholes in a 3D finite element model: the emitter
# operates under the response-model boundary condition Tin = Tout + 1, the receiver in recirculation
# (Tin = Tout). The emitter outlet gives the single-borehole response function and the receiver
# outlet the borehole-to-borehole response function; the block matrix algorithm (Dusseault et al.,
# 2018) then combines them into the BC-IV g-function of a borefield (`outlet_transfer_function` of
# GroundHeatExchanger.jl, `long = :trcm`).
#
# Here the same configuration is built with TRCM columns:
#   - every column (emitter and receivers) is a stack of TRCM layers with vertical conduction in the
#     ground, soil-only layers above the borehole (buried depth D) and below its toe, the ground
#     surface and the bottom of the domain being held at 0;
#   - the receivers are driven one way by the emitter: the ground temperature of the emitter at the
#     receiver location (interpolated in ln r, per layer) is superposed on the receiver, whose own
#     soil nodes carry the disturbance it creates. The effect of a receiver on the emitter (second
#     order, a receiver exchanges no net heat at steady state) is neglected.
# All the receivers are solved in a single system. It is block lower triangular (each receiver
# depends on the emitter only), a structure that the KLU sparse LU detects: the emitter is solved
# once. FBDF reuses its factorization over many steps, which pays off on these 2D (radius × depth)
# columns. The emitter column alone is the long-term TRCM of `trcm(...; D)`, for every exchanger.

const _PAIR_SOLVER = FBDF(linsolve = KLUFactorization())

"""
    trcm_pairs(t, ex::SingleUTube, ks, Cs, H, D, V̇, fluid, d; settings = TRCMSettings())

Single-borehole and borehole-to-borehole response functions measured at the borehole outlet, as in
Rose et al. (2026), with TRCM columns instead of a 3D finite element model. The emitter borehole
operates under `Tin = Tout + 1` and a receiver at each distance `d[j]` in recirculation
(`Tin = Tout`). As in [`trcm`](@ref) with a buried depth `D`, the ground conducts heat vertically
and includes the buried depth, the ground surface (held at 0) and the ground below the borehole
toe, so that the responses are valid from seconds to decades; `g_self` is the outlet transfer
function of `trcm(...; D)`.

The heat flux of the emitter varies along its length and builds up progressively (borehole
capacities, fluid residence time), the two effects that the finite line source neglects. The
receivers see the ground temperature of the emitter at their location and redistribute it along
their length through the circulating fluid. `outlet_transfer_function` of GroundHeatExchanger.jl
(`long = :trcm`) combines these responses into the g-function of a borefield.
# Arguments
    - `t`: Time vector, strictly positive and increasing [s]
    - `ex`: Single U-tube ([`SingleUTube`](@ref)), identical for the emitter and the receivers
    - `ks`: Ground thermal conductivity [W/mK]
    - `Cs`: Ground volumetric heat capacity [J/m³K]
    - `H`: Borehole length [m]
    - `D`: Buried depth of the borehole heads [m]
    - `V̇`: Volumetric flow rate in each borehole [m³/s]
    - `fluid`: Heat carrier fluid ([`FluidProperties`](@ref))
    - `d`: Distance or vector of distances between the emitter and the receivers [m]
# Keywords
    - `settings`: Discretization ([`TRCMSettings`](@ref); here `ns` defaults to 100 shells and `rs`
        to 300 m, which must exceed the distances `d`)
# Output
    - `g_self`: Outlet transfer function of the emitter at `t` [-]
    - `g_pair`: Outlet transfer function of the receivers at `t` [-] (a vector for a scalar `d`, a
        `length(t) × length(d)` matrix otherwise)
# Reference
    - Rose, C., Pasquier, P., Nguyen, A., & Labib, R. (2026). A unified framework for the
        construction of short- and long-term g-functions for parallel-connected ground heat
        exchangers. International Journal of Heat and Mass Transfer, 261, 128565.
        https://doi.org/10.1016/j.ijheatmasstransfer.2026.128565
"""
function trcm_pairs(t, ex::SingleUTube, ks, Cs, H, D, V̇, fluid::FluidProperties,
    d::Union{Real,AbstractVector{<:Real}}; settings::TRCMSettings = TRCMSettings())

    tv = _check_times(t)
    sys = _column_system(ex, ks, Cs, H, D, V̇, fluid, d isa Real ? [d] : d, settings)
    sol = _solve(sys.M, sys.b, tv[end], _PAIR_SOLVER; saveat = tv, save_idxs = sys.outlets,
        reltol = 1e-4, abstol = 1e-6)
    g_self = [sys.outlet(u, 0) for u in sol.u]
    g_pair = [sys.outlet(u, j) for u in sol.u, j in eachindex(d)]
    return (; g_self, g_pair = d isa Real ? vec(g_pair) : g_pair)
end

"""
    _column_system(ex, ks, Cs, H, D, V̇, fluid, d, settings)

Assemble `du/dt = M u + b` for one emitter column and one receiver column per distance in `d`
(`d` empty: the emitter alone, any exchanger; receivers: single U-tube only). Every column has the
same layers (top to bottom): soil layers over the buried depth `D`, `nz` borehole layers, and
`nbelow` soil layers below the toe. Every layer has the same radial soil nodes (index 0 at `rb`:
the borehole-wall node, or the soil core of a soil layer). Returns `M`, `b`, the saved outlet nodes
and `outlet(u, c)`, the outlet temperature of column `c` (0 for the emitter) from the saved values
`u`.
"""
function _column_system(ex::AbstractExchanger, ks, Cs, H, D, V̇, fluid::FluidProperties, d,
    settings)
    ks > 0 && Cs > 0 && H > 0 && V̇ > 0 && D >= 0 || throw(ArgumentError(
        "ks, Cs, H and V̇ must be strictly positive and D non-negative"))
    isempty(d) || ex isa SingleUTube || throw(ArgumentError(
        "the receivers of the borehole-to-borehole responses are single U-tubes only"))
    s = _resolve(settings, H; ns = 100, rs = 300.0)
    ns, rs, nz = s.ns, s.rs, s.nz
    all(ex.rb .< d .< _log_radii(ex.rb, rs, ns)[ns]) || throw(ArgumentError(
        "the distances must exceed the borehole radius and lie within the soil nodes: increase rs"))

    Cf = fluid.ρf * fluid.cf
    ch = _channels(ex, V̇)
    lay = _layer(ex, delta_circuit(ex, ks, V̇, fluid), ch, Cf, ks, Cs, s)
    # Soil-only layer: a core node (r < rb, at rb) and the radial soil chain to rs
    soil = _Network()
    _add_node!(soil, Cs * π * ex.rb^2)
    soil_nodes = _add_chain!(soil, 1, nothing, log(rs / ex.rb) / (2π * ks), Cs,
        _log_shells(ex.rb, rs, ns))
    radii = vcat(ex.rb, lay.radii)                           # radius of soil index 0, 1, …, ns - 1

    # Layer thicknesses: buried depth, borehole, below the toe (geometric growth down to rs)
    dz = H / nz
    ntop = D > 0 ? max(1, round(Int, D / dz)) : 0
    q = _growth_factor(dz, rs, s.nbelow)
    dzs = vcat(fill(D / max(ntop, 1), ntop), fill(dz, nz), [dz * q^j for j in 1:s.nbelow])
    isbh = vcat(falses(ntop), trues(nz), falses(s.nbelow))
    nl = length(dzs)

    # Local node of the radial soil index 0…ns-1 in a borehole and in a soil layer, and its area for
    # the vertical conduction [m²] (the borehole interior does not conduct vertically, a soil core
    # does)
    sidx_bh, sidx_s = vcat(lay.Tb, lay.soil), vcat(1, soil_nodes)
    v = _log_shells(ex.rb, rs, ns)
    area_bh = [i == 0 ? v[1] / 2 : (v[i] + v[i+1]) / 2 for i in 0:(ns - 1)]
    area_s = vcat(area_bh[1] + π * ex.rb^2, area_bh[2:end])

    # Column layout (identical for every column)
    loff = cumsum(vcat(0, [isbh[k] ? length(lay.net.C) : length(soil.C) for k in 1:nl]))
    ncol = loff[end]
    ncols = 1 + length(d)
    N = ncols * ncol
    node(c, k, i) = (c - 1) * ncol + loff[k] + i             # local node i of layer k, column c
    snode(c, k, j) = node(c, k, isbh[k] ? sidx_bh[j+1] : sidx_s[j+1])
    area(k, j) = isbh[k] ? area_bh[j+1] : area_s[j+1]

    rows, cols, vals = Int[], Int[], Float64[]
    C, bE = zeros(N), zeros(N)
    add!(i, j, x) = (push!(rows, i); push!(cols, j); push!(vals, x))

    up = findall(!, ch.down)
    w = ch.flow[up] ./ sum(ch.flow[up])
    kb = findall(isbh)                                       # borehole layers, top to bottom
    for c in 1:ncols
        # Radial conduction and capacities, in energy units (times the layer thickness)
        for k in 1:nl
            net = isbh[k] ? lay.net : soil
            _append_conduction!(rows, cols, vals, net, node(c, k, 0), dzs[k])
            C[node(c, k, 1):node(c, k, length(net.C))] .= net.C .* dzs[k]
        end
        # Vertical conduction in the ground; ground surface and bottom held at 0
        for j in 0:(ns - 1)
            add!(snode(c, 1, j), snode(c, 1, j), -ks * area(1, j) / (dzs[1] / 2))
            add!(snode(c, nl, j), snode(c, nl, j), -ks * area(nl, j) / (dzs[nl] / 2))
            for k in 1:(nl - 1)
                # Through the smaller of the two areas, except at the head and toe of the borehole,
                # where the whole borehole section conducts
                a = j == 0 && isbh[k] != isbh[k+1] ? area_s[1] : min(area(k, j), area(k + 1, j))
                g = ks * a / ((dzs[k] + dzs[k+1]) / 2)
                i1, i2 = snode(c, k, j), snode(c, k + 1, j)
                add!(i1, i2, g); add!(i2, i1, g); add!(i1, i1, -g); add!(i2, i2, -g)
            end
        end
        # Advection (first-order upwind) in energy units, V̇c·Cf
        for (m, k) in enumerate(kb), cc in eachindex(ch.flow)
            i = node(c, k, cc)
            a = ch.flow[cc] * Cf
            add!(i, i, -a)
            if ch.down[cc] && m == 1                     # inlet: Tin = Tout (+1 for the emitter)
                for (mm, u) in enumerate(up)
                    add!(i, node(c, kb[1], u), a * w[mm])
                end
                c == 1 && (bE[i] += a)
            elseif ch.down[cc]
                add!(i, node(c, kb[m-1], cc), a)
            else
                add!(i, m == length(kb) ? node(c, k, ch.pair[cc]) : node(c, kb[m+1], cc), a)
            end
        end
    end

    # Receivers: the emitter ground temperature at distance d acts on the borehole interior of the
    # receiver through the edges that reach its wall node (the soil nodes of a receiver carry the
    # disturbance it creates, its interior nodes the actual temperature)
    net = lay.net
    inner = [(e, net.I[e] == lay.Tb ? net.J[e] : net.I[e]) for e in eachindex(net.G)
             if lay.Tb in (net.I[e], net.J[e]) && !(net.I[e] in lay.soil || net.J[e] in lay.soil)]
    Cg_half = net.C[lay.Tb] - Cs * v[1] / 2                  # grout part of the wall capacity
    Wr, Wc, Wv = Int[], Int[], Float64[]                     # d(T_amb)/dt correction of the wall
    for (jr, dj) in enumerate(d)
        i = searchsortedlast(radii, dj)
        wgt = log(dj / radii[i]) / log(radii[i+1] / radii[i])
        for k in kb
            amb = ((snode(1, k, i - 1), 1 - wgt), (snode(1, k, i), wgt))
            Tb = node(jr + 1, k, lay.Tb)
            for (e, g_in) in inner, (ne, we) in amb
                G = net.G[e] * dzs[k]
                add!(node(jr + 1, k, g_in), ne, G * we)
                add!(Tb, ne, -G * we)
            end
            for (ne, we) in amb
                push!(Wr, Tb); push!(Wc, ne); push!(Wv, Cg_half * dzs[k] * we)
            end
        end
    end

    M = Diagonal(1 ./ C) * sparse(rows, cols, vals, N, N)
    b = bE ./ C
    W = Diagonal(1 ./ C) * sparse(Wr, Wc, Wv, N, N)
    M, b = M - W * M, b - W * b

    # Saved outlet nodes, column after column, and the mixed outlet of column c (0: emitter)
    outlets = [node(c, kb[1], u) for c in 1:ncols for u in up]
    outlet(u, c) = max(0.0, sum(w[m] * u[c*length(up)+m] for m in eachindex(up)))
    return (; M, b, outlets, outlet)
end
