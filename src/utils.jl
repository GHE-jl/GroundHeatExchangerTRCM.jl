# Generic tools: a resistor/capacitor network builder, the solve of the linear ODE system and small
# numerical helpers.
#
# Every network quantity is per unit length of borehole (R in mK/W, C in J/(m K)), so that the
# conduction part of the system matrix does not depend on the layer thickness. An edge stores its
# conductance g = 1/R; a grounded edge links a node to the far field, held at 0 (Eq. 15 of Pasquier
# & Marcotte, 2014).

"""
    _Network()

Growable resistor/capacitor network: `C[i]` is the capacity of node `i` [J/(m K)], `(I[k], J[k],
G[k])` an edge of conductance `G[k]` between nodes `I[k]` and `J[k]`, and `(Ig[k], Gg[k])` an edge
of conductance `Gg[k]` between node `Ig[k]` and the far field.
"""
struct _Network
    C::Vector{Float64}
    I::Vector{Int}
    J::Vector{Int}
    G::Vector{Float64}
    Ig::Vector{Int}
    Gg::Vector{Float64}
end
_Network() = _Network(Float64[], Int[], Int[], Float64[], Int[], Float64[])

"Add a node of capacity `C` [J/(m K)] and return its index."
function _add_node!(net::_Network, C)
    push!(net.C, C)
    return length(net.C)
end

"Connect node `i` to node `j`, or to the far field if `j === nothing`, by a resistance `R` [mK/W]."
function _connect!(net::_Network, i::Int, j::Union{Int,Nothing}, R)
    if j === nothing
        push!(net.Ig, i)
        push!(net.Gg, 1 / R)
    else
        push!(net.I, i)
        push!(net.J, j)
        push!(net.G, 1 / R)
    end
    return nothing
end

"""
    _add_chain!(net, left, right, R, Cv, v; Rl = 0.0, Rr = 0.0)

Connect the node `left` to the node `right` (or to the far field if `right === nothing`) by a chain
of `n = length(v)` equal sub-resistances `R/n` [mK/W], through `n - 1` new interior nodes. `v` are
the areas of the `n` shells of the chain [m²], of a material of volumetric heat capacity `Cv`
[J/m³K]: each node receives half of each of its two neighbouring shells (Pasquier & Marcotte, 2012,
Eq. 8), and the half shell next to the far field is not tracked. `Rl` and `Rr` are capacity-free
resistances (convective films) added to the first and last sub-resistances. Returns the interior
nodes.
"""
function _add_chain!(net::_Network, left::Int, right::Union{Int,Nothing}, R, Cv, v; Rl = 0.0,
    Rr = 0.0)
    n = length(v)
    Rseg = fill(R / n, n)
    Rseg[1] += Rl
    Rseg[n] += Rr
    net.C[left] += Cv * v[1] / 2
    right === nothing || (net.C[right] += Cv * v[n] / 2)
    nodes = [_add_node!(net, Cv * (v[i] + v[i+1]) / 2) for i in 1:(n - 1)]
    path = vcat(left, nodes)
    for k in 1:n
        _connect!(net, path[k], k < n ? path[k+1] : right, Rseg[k])
    end
    return nodes
end

"""
    _log_radii(r₁, r₂, n)

The `n + 1` boundaries of `n` annular shells between `r₁` and `r₂` with a constant radius ratio: for
a uniform conductivity, `n` shells of equal resistance `ln(r₂/r₁)/(2πkn)` (Pasquier & Marcotte,
2012, Eq. 6).
"""
_log_radii(r₁, r₂, n) = [r₁ * (r₂ / r₁)^(i / n) for i in 0:n]

"Areas of the `n` log-spaced shells between `r₁` and `r₂` [m²]."
function _log_shells(r₁, r₂, n)
    r = _log_radii(r₁, r₂, n)
    return [π * (r[i+1]^2 - r[i]^2) for i in 1:n]
end

"""
    _append_conduction!(rows, cols, vals, net, offset, scale = 1.0)

Append the conduction matrix of `net` (a weighted graph Laplacian, `A[i, j] = gᵢⱼ`,
`A[i, i] = -Σⱼ gᵢⱼ - g_far`), multiplied by `scale` and shifted by `offset`, to coordinate lists.
"""
function _append_conduction!(rows, cols, vals, net::_Network, offset, scale = 1.0)
    for (i, j, g) in zip(net.I, net.J, net.G)
        append!(rows, (i, j, i, j) .+ offset)
        append!(cols, (j, i, i, j) .+ offset)
        append!(vals, (g, g, -g, -g) .* scale)
    end
    append!(rows, net.Ig .+ offset)
    append!(cols, net.Ig .+ offset)
    append!(vals, -net.Gg .* scale)
    return nothing
end

"Sparse conduction matrix of a network."
function _laplacian(net::_Network)
    rows, cols, vals = Int[], Int[], Float64[]
    _append_conduction!(rows, cols, vals, net, 0)
    return sparse(rows, cols, vals, length(net.C), length(net.C))
end

"""
    _solve(M, b, tend, solver; kwargs...)

Solve the linear system `du/dt = M u + b` from `u(0) = 0` to `tend`, with the exact sparse Jacobian
`M` (and `∂f/∂t = 0`). The keywords (`saveat`, `save_idxs`, `dense`, `reltol`, `abstol`) are passed
to `solve`. Returns the ODE solution.
"""
function _solve(M, b, tend, solver; kwargs...)
    f!(du, u, _, _) = (mul!(du, M, u); du .+= b)
    jac!(J, _, _, _) = copyto!(J, M)
    tgrad!(dT, _, _, _) = fill!(dT, 0)
    f = ODEFunction(f!; jac = jac!, jac_prototype = copy(M), tgrad = tgrad!)
    sol = solve(ODEProblem(f, zeros(length(b)), (0.0, tend)), solver; kwargs...)
    successful_retcode(sol) || error("the TRCM ODE solve failed (retcode $(sol.retcode))")
    return sol
end

"Check a time vector (non-empty, strictly positive, strictly increasing) and return it as floats."
function _check_times(t)
    tv = collect(Float64, t)
    !isempty(tv) && tv[1] > 0 && all(>(0), diff(tv)) || throw(ArgumentError(
        "the time vector must be non-empty, strictly positive and strictly increasing"))
    return tv
end

"Growth factor `q` such that `n` layers of thickness `dz·q^j` (j = 1…n) span the length `L`."
function _growth_factor(dz, L, n)
    L <= n * dz && return 1.0
    f(q) = dz * q * (q^n - 1) / (q - 1) - L
    lo, hi = 1.0 + 1e-9, 2.0
    while f(hi) < 0
        hi *= 2
    end
    for _ in 1:100                                          # bisection
        mid = (lo + hi) / 2
        f(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end
