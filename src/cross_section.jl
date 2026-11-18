# One depth layer of the TRCM, per unit length of borehole (Pasquier & Marcotte, 2012, 2014):
#   1. one fluid node per flow channel (perfectly mixed radially): the fluid of channel `c` is node
#      `c` of the layer, and the borehole wall is the next node `Tb`;
#   2. the network of the pipes and grout, from the fluid nodes to `Tb`, sized to reproduce the
#      steady-state Δ-circuit of `delta_circuit` exactly;
#   3. an optional ring (casing) and the radial soil chain, from `Tb` to the far field at `rs`.

"""
    _channels(ex, V̇)

Flow channels of an exchanger: cross-section `area` [m²], volumetric `flow` [m³/s], flow direction
(`down[c] == true` for a down-going channel) and, for an up-going channel `c`, the down-going
channel `pair[c]` that feeds it at the bottom of the borehole.
"""
_channels(ex::SingleUTube, V̇) =
    (; area = fill(π * ex.ri^2, 2), flow = [V̇, V̇], down = [true, false], pair = [0, 1])
function _channels(ex::DoubleUTube, V̇)
    area, flow = fill(π * ex.ri^2, 4), fill(V̇ / 2, 4)
    # Diagonal: pipes 1 and 3 down, U-tubes 1→2 and 3→4; adjacent: 1 and 2 down, 1→4 and 2→3
    ex.network === :diagonal &&
        return (; area, flow, down = [true, false, true, false], pair = [0, 1, 0, 3])
    return (; area, flow, down = [true, true, false, false], pair = [0, 0, 2, 1])
end
function _channels(ex::Union{Coaxial,StandingColumn}, V̇)
    # Channel 1: center pipe (dip tube), channel 2: annulus
    area = ex isa Coaxial ? [π * ex.rii^2, π * (ex.roi^2 - ex.rio^2)] :
           [π * ex.ri^2, π * (ex.rb^2 - ex.ro^2)]
    ex.inlet === :center && return (; area, flow = [V̇, V̇], down = [true, false], pair = [0, 1])
    return (; area, flow = [V̇, V̇], down = [false, true], pair = [2, 0])
end

"""
    _layer(ex, δ, ch, Cf, ks, Cs, s; ring = nothing)

One depth layer (per unit length) of exchanger `ex`, from its Δ-circuit `δ` (`delta_circuit`), its
flow channels `ch` (`_channels`), the fluid volumetric heat capacity `Cf`, the ground `ks`, `Cs` and
the discretization `s` (`TRCMSettings` with every default filled in). `ring = (; e, k, C)` adds a
ring of thickness `e`, conductivity `k` and capacity `C` around the borehole (casing). Returns the
network `net`, the borehole-wall node `Tb`, and the soil nodes `soil` located at `radii`.
"""
function _layer(ex::AbstractExchanger, δ, ch, Cf, ks, Cs, s::TRCMSettings; ring = nothing)
    net = _Network()
    for area in ch.area
        _add_node!(net, Cf * area)
    end
    Tb = _add_node!(net, 0.0)
    _cross_section!(net, ex, δ, Tb, s)

    # Ring (casing) and soil chain to the far field, held at 0
    r0, inner = ex.rb, Tb
    if ring !== nothing
        r0, inner = ex.rb + ring.e, _add_node!(net, 0.0)
        _add_chain!(net, Tb, inner, log(r0 / ex.rb) / (2π * ring.k), ring.C,
            _log_shells(ex.rb, r0, s.nc))
    end
    soil = _add_chain!(net, inner, nothing, log(s.rs / r0) / (2π * ks), Cs,
        _log_shells(r0, s.rs, s.ns))
    _check_stability(net)
    return (; net, Tb, soil, radii = _log_radii(r0, s.rs, s.ns)[2:s.ns])
end

"""
    _cross_section!(net, ex, δ, Tb, s)

Network of the pipes and grout of a U-tube (Pasquier & Marcotte, 2012, Fig. 3, extended to `n`
pipes), between the fluid nodes `1:n` and the borehole-wall node `Tb`. Each pipe `i` is linked to
the wall by a chain: film and pipe (`np` sub-resistances), then grout (`ng` sub-resistances), with
the grout resistance `Rg = R1b - Rf - Rp` (Bauer et al., 2011, Eqs. 34/36). The grout-to-grout
(core) chains start at the first grout node, at the fraction `x = 1/ng` of `Rg` from the pipe
(`x = 0` for `ng = 1`), and their resistances follow from `_core_conductances` at that `x`: the
network reduces exactly to `Rb` and `Ra` at steady state for any discretization.

The grout area is split as in `gTRCM.m`: a core area between the pipes (`2ro·s - πro²` for two
pipes, the square between the four pipe centers minus the pipe quarters for four pipes) shared by
the core chains, and the remaining annulus between an equivalent radius `r_geq` and `rb`, shared by
the `n` outer grout chains.
"""
function _cross_section!(net, ex::Union{SingleUTube,DoubleUTube}, δ, Tb, s)
    n = ex isa SingleUTube ? 2 : 4
    Rfp = δ.Rf + δ.Rp
    Rg = δ.R1b - Rfp
    Rg > 0 || throw(ArgumentError("the multipole borehole resistance is smaller than the fluid + " *
        "pipe resistance (R1b = $(δ.R1b), Rf + Rp = $Rfp mK/W): the grout resistance is negative"))

    # Grout areas: core between the pipes, outer annulus between r_geq and rb
    rc = ex.s / 2
    A_core = n == 2 ? 4 * ex.ro * rc - π * ex.ro^2 : 2 * rc^2 - π * ex.ro^2
    A_grout = π * ex.rb^2 - n * π * ex.ro^2
    0 < A_core < A_grout || throw(ArgumentError("invalid pipe layout: the grout core area " *
        "between the pipes ($A_core m²) must be positive and smaller than the grout area " *
        "($A_grout m²)"))
    r_geq = sqrt(n * ex.ro^2 + A_core / π)

    # Fluid → (film + pipe) → pipe outer wall → (grout) → Tb
    attach = zeros(Int, n)
    for i in 1:n
        wall = _add_node!(net, 0.0)
        _add_chain!(net, i, wall, δ.Rp, ex.Cp, _log_shells(ex.ri, ex.ro, s.np); Rl = δ.Rf)
        grout = _add_chain!(net, wall, Tb, Rg, ex.Cg, _log_shells(r_geq, ex.rb, s.ng) ./ n)
        attach[i] = isempty(grout) ? wall : grout[1]
    end

    # Core chains between the grout attachment nodes of every pair of pipes coupled by the
    # Δ-circuit. A negative resistance (Bauer et al., 2011, Section 2.6) is never subdivided: an
    # interior node reached only through negative sub-resistances would be unconditionally unstable.
    x = s.ng >= 2 ? 1 / s.ng : 0.0
    Gd = n == 2 ? [0 1/δ.R12; 1/δ.R12 0] :
         [0 1/δ.R12 1/δ.R13 1/δ.R12; 1/δ.R12 0 1/δ.R12 1/δ.R13;
          1/δ.R13 1/δ.R12 0 1/δ.R12; 1/δ.R12 1/δ.R13 1/δ.R12 0]
    Ggg = _core_conductances(δ.R1b, Gd, Rfp + x * Rg, (1 - x) * Rg)
    tol = 1e-12 * maximum(abs, Ggg[i, j] for i in 1:n, j in 1:n if i != j)
    pairs = [(i, j, 1 / Ggg[i, j]) for i in 1:n for j in (i + 1):n if abs(Ggg[i, j]) > tol]
    for (i, j, R) in pairs
        m = R > 0 ? s.ngg : 1
        _add_chain!(net, attach[i], attach[j], R, ex.Cg, fill(A_core / (length(pairs) * m), m))
    end
    if isempty(pairs)                                  # no core coupling: keep the core capacity
        net.C[attach] .+= ex.Cg * A_core / n
    end
    return nothing
end

function _cross_section!(net, ex::Coaxial, δ, Tb, s)
    # Center fluid (1) → (film, center pipe, film) → annulus fluid (2)
    _add_chain!(net, 1, 2, δ.Rpi, ex.Cpi, _log_shells(ex.rii, ex.rio, s.np); Rl = δ.Rf_in,
        Rr = δ.Rf_ann_in)
    # Annulus fluid → (film, outer pipe) → (grout, if any) → Tb
    if ex.roo < ex.rb
        wall = _add_node!(net, 0.0)
        _add_chain!(net, 2, wall, δ.Rpo, ex.Cpo, _log_shells(ex.roi, ex.roo, s.np);
            Rl = δ.Rf_ann_out)
        _add_chain!(net, wall, Tb, δ.Rg, ex.Cg, _log_shells(ex.roo, ex.rb, s.ng))
    else
        _add_chain!(net, 2, Tb, δ.Rpo, ex.Cpo, _log_shells(ex.roi, ex.roo, s.np); Rl = δ.Rf_ann_out)
    end
    return nothing
end

function _cross_section!(net, ex::StandingColumn, δ, Tb, s)
    # Dip tube water (1) → (film, dip tube, film) → annulus water (2) → (wall film) → rock
    _add_chain!(net, 1, 2, δ.Rp, ex.Cp, _log_shells(ex.ri, ex.ro, s.np); Rl = δ.Rf_in,
        Rr = δ.Rf_ann_in)
    _connect!(net, 2, Tb, δ.Rf_wall)
    return nothing
end

"""
    _core_conductances(R1b, Gd, a, b)

Grout-to-grout conductances that make a U-tube TRCM reproduce a target Δ-circuit exactly at steady
state, for any number of symmetrically placed pipes.

Each pipe `i` is linked to its grout attachment node `gᵢ` by the series resistance `a` (fluid film,
pipe and the first fraction `x` of the grout resistance `Rg`), and `gᵢ` to the borehole wall by
`b = (1 - x)Rg`. The grout nodes are linked together by conductances `Gᵍᵍᵢⱼ`. Eliminating the grout
nodes (Schur complement), the conductance matrix seen from the fluids is `S = I/a - K⁻¹/a²`, with
`K = (1/a + 1/b)I + L(Gᵍᵍ)` and `L` the graph Laplacian of `Gᵍᵍ`. Imposing `S` equal to the
Δ-circuit (`Sᵢᵢ = 1/R1b + Σⱼ Gdᵢⱼ`, `Sᵢⱼ = -Gdᵢⱼ`) gives `K = (aI - a²S)⁻¹` and
`Gᵍᵍ = -(K - (1/a + 1/b)I)`.
For two pipes this is Eq. 21 of Bauer et al. (2011) (Eq. 10 of Pasquier & Marcotte, 2012) at the
grout location `x`; the same algebra covers the double U-tube Δ-circuit of Claesson & Javed (2019).
Returns the full `n × n` matrix (only its off-diagonal entries are meaningful).
"""
function _core_conductances(R1b, Gd, a, b)
    S = -Matrix{Float64}(Gd)
    for i in axes(S, 1)
        S[i, i] = 1 / R1b + sum(Gd[i, :]) - Gd[i, i]
    end
    K = inv(a * I - a^2 * S)
    return -(K - (1 / a + 1 / b) * I)
end

"""
    _check_stability(net)

A layer is stable if and only if its conduction matrix `A` is negative definite (`C⁻¹A` then has
only negative real eigenvalues, `C` being positive). This only fails when the Δ-circuit contains a
negative resistance large enough to dominate the positive paths (Bauer et al., 2011, Eqs. 44-45), a
genuine growing mode. The layers are only coupled by the advection, which is dissipative, so
checking every distinct layer is sufficient.
"""
function _check_stability(net)
    isposdef(Symmetric(-Matrix(_laplacian(net)))) || throw(ArgumentError(
        "the TRCM cross-section is unstable: a negative grout-to-grout resistance of the " *
        "Δ-circuit dominates the network (Bauer et al., 2011, Section 2.6). Bauer et al. " *
        "recommend moving the grout capacity closer to the pipes: the grout attachment point " *
        "is at x = 1/ng of the grout resistance, so increase `ng`, or use `ng = 1` (x = 0)."))
    return nothing
end
