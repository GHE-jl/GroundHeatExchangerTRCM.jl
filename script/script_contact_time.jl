# Contact time and far-field radius of the hybrid TRCM + finite line source (FLS).
#
# The TRCM of `trcm` (without a buried depth) has no vertical conduction in the ground and holds the
# ground at its initial temperature at rs; `outlet_transfer_function` of GroundHeatExchanger.jl
# therefore switches to the FLS increase after a contact time tc. The reference is `trcm` with a
# buried depth D: the same TRCM columns with vertical conduction in the ground, the buried depth,
# the ground surface and the ground below the toe (rs = 300 m), i.e. a model valid from seconds to
# decades.
#   Figure 1: outlet transfer function of the reference, of the TRCM alone and of the hybrid
#   Figure 2: maximum relative error of the hybrid (1 h → 50 yr) vs tc, for several rs
# A contact time too early misses the borehole capacities, too late the axial effects that the
# TRCM lacks; from ~3 days on the error is flat (~0.5 %, the vertical discretization of the
# TRCM), whatever rs as long as the far field is not felt at tc. Defaults: tc = 7 days
# (the horizon of the published short-term ANN) and rs = 20 m.

import Pkg; Pkg.activate(@__DIR__)

using CairoMakie
using GroundHeatExchanger: outlet_transfer_function
using Printf
using ThermalResistanceCapacitanceModel

# Ground, borehole and fluid
ks = 2.5      # Ground thermal conductivity [W/mK]
Cs = 2.2e6    # Ground volumetric heat capacity [J/m³K]
α  = ks / Cs  # Ground thermal diffusivity [m²/s]
H  = 150.0    # Borehole length [m]
D  = 4.0      # Buried depth [m]
V̇  = 0.5e-3   # Total flow rate [m³/s]
fluid = FluidProperties(0.6, 4180.0, 1000.0, 1.0e-3)

# Single U-loop
rb = 0.075    # Borehole radius [m]
ro = 0.0167   # Pipe outer radius [m]
ri = 0.0137   # Pipe inner radius [m]
s  = 0.06     # Shank spacing [m]
kg = 1.5      # Grout thermal conductivity [W/mK]
Cg = 2.0e6    # Grout volumetric heat capacity [J/m³K]
kp = 0.4      # Pipe thermal conductivity [W/mK]
Cp = 1.8e6    # Pipe volumetric heat capacity [J/m³K]
ex = SingleUTube(rb, ro, ri, s, kg, Cg, kp, Cp)

# Log-spaced time vector: 10 s → 50 yr; errors are measured after one hour
day  = 86400.0
yr   = 365.25day
t    = exp10.(range(1, log10(50yr), length = 160))
late = t .>= 3600

# Reference: TRCM with vertical conduction (buried depth D: surface, toe)
g_ref = trcm(t, ex, ks, Cs, H, V̇, fluid; D)
relative(g) = abs.(g .- g_ref) ./ g_ref

# TRCM alone (no vertical conduction, rs = 500 m: unbounded over 50 years)
g_trcm = trcm(t, ex, ks, Cs, H, V̇, fluid; settings = TRCMSettings(rs = 500.0, ns = 40))

# Hybrid with the default settings (tc = 7 days, rs = 20 m)
g_hyb = outlet_transfer_function(t, ex, ks, Cs, H, D, V̇, fluid)

println("Outlet transfer function: reference, TRCM alone and hybrid (relative to the reference)")
for (label, y) in (("1 day", day), ("7 days", 7day), ("1 year", yr), ("10 years", 10yr),
    ("50 years", 50yr))
    k = argmin(abs.(t .- y))
    @printf("   %-9s reference %.4f   TRCM %+6.2f %%   hybrid %+6.2f %%\n", label, g_ref[k],
        100 * (g_trcm[k] / g_ref[k] - 1), 100 * (g_hyb[k] / g_ref[k] - 1))
end

# Maximum error of the hybrid vs the contact time, for one far-field radius (NaN when the far field
# is felt before tc, α·tc/rs² > 0.3)
tc = [3600, 6 * 3600, day, 3day, 7day, 21day, 60day, 180day]
function hybrid_errors(rs)
    ns = ceil(Int, 24 * log(rs / rb) / log(20 / rb))   # Same shell ratio as the default
    return [α * tc_i / rs^2 > 0.3 ? NaN : 100 * maximum(relative(outlet_transfer_function(t, ex,
        ks, Cs, H, D, V̇, fluid; tc = tc_i, settings = TRCMSettings(rs = rs, ns = ns)))[late])
        for tc_i in tc]
end

# Hybrid rs = 5 m
err_5 = hybrid_errors(5.0)
# Hybrid rs = 20 m (default)
err_20 = hybrid_errors(20.0)
# Hybrid rs = 50 m
err_50 = hybrid_errors(50.0)

println("\nHybrid: maximum relative error (1 h → 50 yr) [%]")
@printf("   %-8s", "tc [d]"); foreach(x -> @printf(" %7.2f", x / day), tc); println()
@printf("   %-8s", "rs = 5"); foreach(e -> @printf(" %7.3f", e), err_5); println()
@printf("   %-8s", "rs = 20"); foreach(e -> @printf(" %7.3f", e), err_20); println()
@printf("   %-8s", "rs = 50"); foreach(e -> @printf(" %7.3f", e), err_50); println()

# Figure 1: outlet transfer functions
fig1 = Figure(size = (750, 480))
ax1 = Axis(fig1[1, 1], xlabel = "Time (s)", ylabel = "g_out (-)", xscale = log10,
    title = "Outlet transfer function, single U-loop")
lines!(ax1, t, g_ref, linewidth = 3, color = :black, label = "TRCM with vertical conduction")
lines!(ax1, t, g_trcm, linewidth = 2, color = :steelblue, linestyle = :dash,
    label = "TRCM alone (no vertical conduction)")
lines!(ax1, t, g_hyb, linewidth = 2, color = :tomato, linestyle = :dot,
    label = "Hybrid TRCM + FLS (tc = 7 d)")
vlines!(ax1, [7day], color = :gray, linestyle = :dashdot)
axislegend(ax1, position = :lt)
display(fig1)

# Figure 2: hybrid error vs contact time
fig2 = Figure(size = (750, 480))
ax2 = Axis(fig2[1, 1], xlabel = "Contact time tc (days)", ylabel = "Max relative error (%)",
    xscale = log10, yscale = log10, title = "Hybrid TRCM + FLS vs the reference")
scatterlines!(ax2, tc ./ day, err_5, label = "rs = 5 m")
scatterlines!(ax2, tc ./ day, err_20, label = "rs = 20 m (default)")
scatterlines!(ax2, tc ./ day, err_50, label = "rs = 50 m")
vlines!(ax2, [7.0], color = :gray, linestyle = :dashdot)
axislegend(ax2, position = :rt)
display(fig2)
