# Borehole-to-borehole response functions measured at the borehole outlet (Rose et al., 2026),
# computed with TRCM columns instead of a 3D finite element model, and their use in the spatial
# superposition of a borefield.
#   Figure 1: borehole-to-borehole response functions of the TRCM vs the finite line source (FLS)
#             at several distances, and their relative difference
#   Figure 2: single-borehole response function with vertical conduction in the ground, vs the
#             TRCM without vertical conduction and the hybrid TRCM + FLS
#   Figure 3: BC-IV g-function of a 5×5 borefield (bc_II of GroundHeatExchanger.jl) built with
#             TRCM pairs, FLS pairs, TRCM pairs up to 12 m and FLS beyond, and the current hybrid
#             approach (BC-II FLS field shifted onto the single-borehole response)
# Reference
#   Rose, C., Pasquier, P., Nguyen, A., & Labib, R. (2026). A unified framework for the construction
#   of short- and long-term g-functions for parallel-connected ground heat exchangers.
#   International Journal of Heat and Mass Transfer, 261, 128565.
#   Dusseault, B., Pasquier, P., & Marcotte, D. (2018). A block matrix formulation for efficient
#   g-function construction. Renewable Energy, 121, 249–260.

import Pkg; Pkg.activate(@__DIR__)

using CairoMakie
using GroundHeatExchanger: bc_II, borefield, fls, outlet_transfer_function
using Printf
using ThermalResistanceCapacitanceModel

# Ground, borehole and flow rate
ks = 2.5                  # Ground thermal conductivity [W/mK]
Cs = 2.2e6                # Ground volumetric heat capacity [J/m³K]
H  = 100.0                # Borehole length [m]
D  = 4.0                  # Buried depth [m]
V̇  = 0.4e-3               # Flow rate in each borehole [m³/s]
fluid = FluidProperties(0.6, 4180.0, 1000.0, 1.0e-3)
ex = SingleUTube(rb = 0.075, ro = 0.0167, ri = 0.0137, s = 0.06, kg = 1.5, Cg = 2.0e6, kp = 0.4,
    Cp = 1.8e6)
nz = 25                   # Borehole layers (4 m; nz = 50 changes the results by ~0.5 %)
yr = 365.25 * 86400       # One year [s]
q  = V̇ * fluid.ρf * fluid.cf / H      # Heat rate per unit length for Tin - Tout = 1 [W/(m⋅K)]

# Figure 1: borehole-to-borehole response functions vs FLS
t = exp10.(range(1, log10(30yr), length = 120))
d = [3.0, 6.0, 12.0, 24.0, 48.0]
el = @elapsed g_self, g_pair = trcm_pairs(t, ex, ks, Cs, H, D, V̇, fluid, d;
    settings = TRCMSettings(nz = nz))
println("Pair responses (", length(d), " receivers): ", round(el, digits = 1), " s")

# FLS driven by the mean heat rate of the emitter
g_fls_3  = q .* fls(t, 3.0, H, D, ks, Cs)
g_fls_6  = q .* fls(t, 6.0, H, D, ks, Cs)
g_fls_12 = q .* fls(t, 12.0, H, D, ks, Cs)
g_fls_24 = q .* fls(t, 24.0, H, D, ks, Cs)
g_fls_48 = q .* fls(t, 48.0, H, D, ks, Cs)

# Relative difference TRCM/FLS, shown only once the interaction is noticeable (g > 0.01)
reldiff(g, gf) = [gf[i] > 0.01 ? 100 * (g[i] / gf[i] - 1) : NaN for i in eachindex(g)]

# TRCM / FLS - 1 at 1, 10 and 30 years
k = [argmin(abs.(t .- yr)), argmin(abs.(t .- 10yr)), lastindex(t)]
println("TRCM / FLS - 1 at 1, 10 and 30 years")
@printf("   d =  3 m  %+6.2f %%  %+6.2f %%  %+6.2f %%\n", reldiff(g_pair[k, 1], g_fls_3[k])...)
@printf("   d =  6 m  %+6.2f %%  %+6.2f %%  %+6.2f %%\n", reldiff(g_pair[k, 2], g_fls_6[k])...)
@printf("   d = 12 m  %+6.2f %%  %+6.2f %%  %+6.2f %%\n", reldiff(g_pair[k, 3], g_fls_12[k])...)
@printf("   d = 24 m  %+6.2f %%  %+6.2f %%  %+6.2f %%\n", reldiff(g_pair[k, 4], g_fls_24[k])...)
@printf("   d = 48 m  %+6.2f %%  %+6.2f %%  %+6.2f %%\n", reldiff(g_pair[k, 5], g_fls_48[k])...)

fig1 = Figure(size = (1100, 450))
ax1 = Axis(fig1[1, 1], xlabel = "Time (yr)", ylabel = "g (-)", xscale = log10,
    title = "Borehole-to-borehole response functions: TRCM (solid) vs FLS (dashed)")
lines!(ax1, t ./ yr, g_pair[:, 1], color = Cycled(1), label = "d = 3 m")
lines!(ax1, t ./ yr, g_fls_3, color = Cycled(1), linestyle = :dash)
lines!(ax1, t ./ yr, g_pair[:, 2], color = Cycled(2), label = "d = 6 m")
lines!(ax1, t ./ yr, g_fls_6, color = Cycled(2), linestyle = :dash)
lines!(ax1, t ./ yr, g_pair[:, 3], color = Cycled(3), label = "d = 12 m")
lines!(ax1, t ./ yr, g_fls_12, color = Cycled(3), linestyle = :dash)
lines!(ax1, t ./ yr, g_pair[:, 4], color = Cycled(4), label = "d = 24 m")
lines!(ax1, t ./ yr, g_fls_24, color = Cycled(4), linestyle = :dash)
lines!(ax1, t ./ yr, g_pair[:, 5], color = Cycled(5), label = "d = 48 m")
lines!(ax1, t ./ yr, g_fls_48, color = Cycled(5), linestyle = :dash)
axislegend(ax1, position = :lt)

ax2 = Axis(fig1[1, 2], xlabel = "Time (yr)", ylabel = "TRCM / FLS - 1 (%)", xscale = log10,
    title = "Relative difference (g > 0.01)")
lines!(ax2, t ./ yr, reldiff(g_pair[:, 1], g_fls_3), color = Cycled(1), label = "d = 3 m")
lines!(ax2, t ./ yr, reldiff(g_pair[:, 2], g_fls_6), color = Cycled(2), label = "d = 6 m")
lines!(ax2, t ./ yr, reldiff(g_pair[:, 3], g_fls_12), color = Cycled(3), label = "d = 12 m")
lines!(ax2, t ./ yr, reldiff(g_pair[:, 4], g_fls_24), color = Cycled(4), label = "d = 24 m")
lines!(ax2, t ./ yr, reldiff(g_pair[:, 5], g_fls_48), color = Cycled(5), label = "d = 48 m")
hlines!(ax2, 0.0, color = :gray, linestyle = :dot)
axislegend(ax2, position = :lb)
display(fig1)

# Figure 2: single-borehole response function
# TRCM with vertical conduction: g_self, the emitter of the pair model
# TRCM without vertical conduction (same radial discretization)
g_trcm = trcm(t, ex, ks, Cs, H, V̇, fluid;
    settings = TRCMSettings(nz = nz, rs = 300.0, ns = 100))
# Hybrid TRCM + FLS (contact time 7 days)
g_hyb = outlet_transfer_function(t, ex, ks, Cs, H, D, V̇, fluid; settings = TRCMSettings(nz = nz))

fig2 = Figure(size = (650, 450))
ax3 = Axis(fig2[1, 1], xlabel = "Time (yr)", ylabel = "g (-)", xscale = log10,
    title = "Single-borehole outlet response function")
lines!(ax3, t ./ yr, g_self, linewidth = 3, color = :black, label = "TRCM, vertical conduction")
lines!(ax3, t ./ yr, g_trcm, linewidth = 2, color = :steelblue, linestyle = :dash,
    label = "TRCM, no vertical conduction")
lines!(ax3, t ./ yr, g_hyb, linewidth = 2, color = :tomato, linestyle = :dot,
    label = "Hybrid TRCM + FLS (tc = 7 d)")
axislegend(ax3, position = :lt)
display(fig2)

# Figure 3: BC-IV g-function of a 5×5 borefield
# The BC-II solver (bc_II of GroundHeatExchanger.jl) needs a constant time step: monthly steps over
# 20 years (outlet_transfer_function with long = :trcm does the same on its own time blocks)
tm = yr / 12 .* (1:240)
xy = borefield(:rectangle, 5, 5, 6.0)
nb = size(xy, 1)
dist = [hypot(xy[i, 1] - xy[j, 1], xy[i, 2] - xy[j, 2]) for i in 1:nb, j in 1:nb]
du = sort(unique(round.(dist[dist .> 0], digits = 6)))       # 14 unique distances

# Pair responses at every unique distance of the borefield (one solve)
el = @elapsed g_self_m, g_pair_m = trcm_pairs(tm, ex, ks, Cs, H, D, V̇, fluid, du;
    settings = TRCMSettings(nz = nz))
println("Borefield pair responses (", length(du), " receivers): ", round(el, digits = 1), " s")
trcm_pair(r) = g_pair_m[:, findfirst(x -> isapprox(x, r, atol = 1e-5), du)]
fls_pair(r) = q .* fls(tm, r, H, D, ks, Cs)

# Response matrix (nt × nb × nb): inlet response on the diagonal (Tin = Tout + 1), outlet of a
# receiver in recirculation elsewhere; bc_II enforces the same inlet temperature everywhere
function response_matrix(pair)
    gm = zeros(length(tm), nb, nb)
    for i in 1:nb, j in 1:nb
        gm[:, i, j] = i == j ? g_self_m .+ 1 : pair(dist[i, j])
    end
    return gm
end

# BC-IV with TRCM pairs (Rose et al., 2026)
g_field_trcm = bc_II(response_matrix(trcm_pair)) .- 1
# BC-IV with FLS pairs (first simplification of Rose et al., 2026)
g_field_fls = bc_II(response_matrix(fls_pair)) .- 1
# BC-IV with TRCM pairs up to 12 m, FLS beyond
cut_pair(r) = r <= 12.0 + 1e-6 ? trcm_pair(r) : fls_pair(r)
g_field_cut = bc_II(response_matrix(cut_pair)) .- 1
# Current hybrid approach: BC-II FLS field shifted onto the TRCM at tc = 7 d
g_field_hyb = outlet_transfer_function(tm, ex, ks, Cs, H, D, V̇, fluid; xy,
    settings = TRCMSettings(nz = nz))

for k in (12, 60, 120, 240)
    @printf("t = %4.1f yr  TRCM pairs g = %.3f | FLS pairs %+.2f %% | cut 12 m %+.2f %% | hybrid %+.2f %%\n",
        tm[k] / yr, g_field_trcm[k], 100 * (g_field_fls[k] / g_field_trcm[k] - 1),
        100 * (g_field_cut[k] / g_field_trcm[k] - 1), 100 * (g_field_hyb[k] / g_field_trcm[k] - 1))
end

fig3 = Figure(size = (1100, 450))
ax4 = Axis(fig3[1, 1], xlabel = "Time (yr)", ylabel = "g (-)",
    title = "Mixed outlet g-function, 5×5 borefield, B = 6 m")
lines!(ax4, tm ./ yr, g_field_trcm, linewidth = 3, color = :black, label = "BC-IV, TRCM pairs")
lines!(ax4, tm ./ yr, g_field_fls, linewidth = 2, color = :tomato, linestyle = :dash,
    label = "BC-IV, FLS pairs")
lines!(ax4, tm ./ yr, g_field_cut, linewidth = 2, color = :darkgreen, linestyle = :dot,
    label = "BC-IV, TRCM pairs ≤ 12 m")
lines!(ax4, tm ./ yr, g_field_hyb, linewidth = 2, color = :steelblue, linestyle = :dashdot,
    label = "Hybrid TRCM + BC-II FLS")
axislegend(ax4, position = :lt)

ax5 = Axis(fig3[1, 2], xlabel = "Time (yr)", ylabel = "Difference with TRCM pairs (%)",
    title = "Relative difference")
lines!(ax5, tm ./ yr, 100 .* (g_field_fls ./ g_field_trcm .- 1), color = :tomato,
    label = "BC-IV, FLS pairs")
lines!(ax5, tm ./ yr, 100 .* (g_field_cut ./ g_field_trcm .- 1), color = :darkgreen,
    label = "BC-IV, TRCM pairs ≤ 12 m")
lines!(ax5, tm ./ yr, 100 .* (g_field_hyb ./ g_field_trcm .- 1), color = :steelblue,
    label = "Hybrid TRCM + BC-II FLS")
hlines!(ax5, 0.0, color = :gray, linestyle = :dot)
axislegend(ax5, position = :lt)
display(fig3)
