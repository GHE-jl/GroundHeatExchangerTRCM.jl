# Outlet transfer functions of the exchanger types of the package for the same ground, depth
# and flow rate, from 10 s to 25 years:
#   Figure 1: outlet transfer function of the TRCM alone (solid, ground held at 0 at rs = 20 m) and
#             of the hybrid TRCM + finite line source (dashed)
# Reference
#   Pasquier, P., & Marcotte, D. (2014). Joint use of quasi-3D response model and spectral method
#   to simulate borehole heat exchanger. Geothermics, 51, 281–299.

import Pkg; Pkg.activate(@__DIR__)

using CairoMakie
using GroundHeatExchanger: outlet_transfer_function
using Printf
using ThermalResistanceCapacitanceModel

# Ground, borehole and fluid
ks = 2.5      # Ground thermal conductivity [W/mK]
Cs = 2.2e6    # Ground volumetric heat capacity [J/m³K]
H  = 150.0    # Borehole length [m]
D  = 4.0      # Buried depth [m]
V̇  = 0.5e-3   # Total flow rate [m³/s]
fluid = FluidProperties(10.0, :water)   # Water at 10 °C (CoolProp)

# Grout and pipes of the U-loops
rb = 0.075    # Borehole radius [m]
ro = 0.0167   # Pipe outer radius [m]
ri = 0.0137   # Pipe inner radius [m]
kg = 1.5      # Grout thermal conductivity [W/mK]
Cg = 2.0e6    # Grout volumetric heat capacity [J/m³K]
kp = 0.4      # Pipe thermal conductivity [W/mK]
Cp = 1.8e6    # Pipe volumetric heat capacity [J/m³K]

# Log-spaced time vector: 10 s → 25 yr
yr = 365.25 * 86400
t  = exp10.(range(1, log10(25yr), length = 150))
k_day = argmin(abs.(t .- 86400))   # Index of t ≈ 1 day

println("Exchanger                    g_out(1 d)   g_out(25 yr, hybrid)")

# TRCM single U-loop
s_1U = 0.06   # Shank spacing [m]
ex_1U = SingleUTube(rb, ro, ri, s_1U, kg, Cg, kp, Cp)
g_trcm_1U = trcm(t, ex_1U, ks, Cs, H, V̇, fluid)
g_1U = outlet_transfer_function(t, ex_1U, ks, Cs, H, D, V̇, fluid)
@printf("Single U-loop                %10.4f   %10.4f\n", g_trcm_1U[k_day], g_1U[end])

# TRCM double U-loop (diagonal inlets)
s_2U = 0.09   # Shank spacing between opposite pipes [m]
ex_2U_diag = DoubleUTube(rb, ro, ri, s_2U, kg, Cg, kp, Cp, :diagonal)
g_trcm_2U_diag = trcm(t, ex_2U_diag, ks, Cs, H, V̇, fluid)
g_2U_diag = outlet_transfer_function(t, ex_2U_diag, ks, Cs, H, D, V̇, fluid)
@printf("Double U-loop (diagonal)     %10.4f   %10.4f\n", g_trcm_2U_diag[k_day], g_2U_diag[end])

# TRCM double U-loop (adjacent inlets)
ex_2U_adj = DoubleUTube(rb, ro, ri, s_2U, kg, Cg, kp, Cp, :adjacent)
g_trcm_2U_adj = trcm(t, ex_2U_adj, ks, Cs, H, V̇, fluid)
g_2U_adj = outlet_transfer_function(t, ex_2U_adj, ks, Cs, H, D, V̇, fluid)
@printf("Double U-loop (adjacent)     %10.4f   %10.4f\n", g_trcm_2U_adj[k_day], g_2U_adj[end])

# TRCM coaxial (inlet in the annulus)
rii = 0.0177  # Inner pipe inner radius [m]
rio = 0.0201  # Inner pipe outer radius [m]
roi = 0.0545  # Outer pipe inner radius [m]
roo = 0.0572  # Outer pipe outer radius [m]
ex_coax = Coaxial(rb, rii, rio, roi, roo, kg, Cg, kp, Cp, kp, Cp, :annulus)
g_trcm_coax = trcm(t, ex_coax, ks, Cs, H, V̇, fluid)
g_coax = outlet_transfer_function(t, ex_coax, ks, Cs, H, D, V̇, fluid)
@printf("Coaxial (annulus inlet)      %10.4f   %10.4f\n", g_trcm_coax[k_day], g_coax[end])

# TRCM standing column well
rb_scw = 0.076   # Borehole radius [m]
ro_scw = 0.05    # Dip tube outer radius [m]
ri_scw = 0.045   # Dip tube inner radius [m]
ex_scw = StandingColumn(rb_scw, ro_scw, ri_scw, kp, Cp)
g_trcm_scw = trcm(t, ex_scw, ks, Cs, H, V̇, fluid)
g_scw = outlet_transfer_function(t, ex_scw, ks, Cs, H, D, V̇, fluid)
@printf("Standing column well         %10.4f   %10.4f\n", g_trcm_scw[k_day], g_scw[end])

# Figure 1: outlet transfer functions, TRCM (solid) and hybrid TRCM + FLS (dashed)
fig1 = Figure(size = (800, 520))
ax1 = Axis(fig1[1, 1], xlabel = "Time (s)", ylabel = "g_out (-)", xscale = log10,
    title = "Outlet transfer function — TRCM (solid), hybrid TRCM + FLS (dashed)")
lines!(ax1, t, g_trcm_1U, linewidth = 2, color = :steelblue, label = "Single U-loop")
lines!(ax1, t, g_1U, linewidth = 2, color = :steelblue, linestyle = :dash)
lines!(ax1, t, g_trcm_2U_diag, linewidth = 2, color = :tomato,
    label = "Double U-loop (diagonal)")
lines!(ax1, t, g_2U_diag, linewidth = 2, color = :tomato, linestyle = :dash)
lines!(ax1, t, g_trcm_2U_adj, linewidth = 2, color = :darkorange,
    label = "Double U-loop (adjacent)")
lines!(ax1, t, g_2U_adj, linewidth = 2, color = :darkorange, linestyle = :dash)
lines!(ax1, t, g_trcm_coax, linewidth = 2, color = :darkgreen, label = "Coaxial (annulus inlet)")
lines!(ax1, t, g_coax, linewidth = 2, color = :darkgreen, linestyle = :dash)
lines!(ax1, t, g_trcm_scw, linewidth = 2, color = :purple, label = "Standing column well")
lines!(ax1, t, g_scw, linewidth = 2, color = :purple, linestyle = :dash)
axislegend(ax1, position = :lt)
display(fig1)

