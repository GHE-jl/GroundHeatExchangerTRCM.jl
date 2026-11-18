# Transient borehole resistance Rb(t) of the modified TRCM compared with the steady-state multipole
# value, for scenarios "a" and "e" of Table 2 of Pasquier & Marcotte (2012). As in their Fig. 8, the
# pipe inner-wall temperature is imposed (no convective film, `films = false`).
#
# Every curve must converge to the multipole Rb; only the short-time shape changes with the number
# of grout sub-resistances ng, i.e. with the attachment point x = 1/ng of the grout core (x = 0 for
# ng = 1), which replaces the x-sensitivity of the original TRCM discussed in the paper.
#
# Notes on Table 2:
#   - the extracted text gives ri = 0.065 m and ro = 0.022 m, which is impossible (ri > ro); a 32 mm
#     SDR-11 HDPE pipe is used instead, which reproduces the steady Rb of Fig. 8 (≈ 0.14 and
#     ≈ 0.075 mK/W);
#   - the "cp, cg, cs" column labelled MJ/kgK holds volumetric capacities in kJ/(m³K).
# Reference
#   Pasquier, P., & Marcotte, D. (2012). Short-term simulation of ground heat exchanger with an
#   improved TRCM. Renewable Energy, 46, 92–99.

import Pkg; Pkg.activate(@__DIR__)

using CairoMakie
using Printf
using ThermalResistanceCapacitanceModel
using ThermalResistanceCapacitanceModel: delta_circuit, transient_resistance

# Common parameters of scenarios a and e
ks = 3.0      # Ground thermal conductivity [W/mK]
Cs = 2.0e6    # Ground volumetric heat capacity [J/m³K]
kg = 1.5      # Grout thermal conductivity [W/mK]
Cg = 2.0e6    # Grout volumetric heat capacity [J/m³K]
rb = 0.08     # Borehole radius [m]
ro = 0.016    # Pipe outer radius [m]
ri = 0.0137   # Pipe inner radius [m]
Cp = 2.185e6  # Pipe volumetric heat capacity [J/m³K]
V̇  = 0.4e-3   # Flow rate [m³/s] (unused without films)
fluid = FluidProperties(0.6, 4180.0, 1000.0, 1.0e-3)   # Films are disabled: only Cf matters

# Log-spaced time vector: 10 s → 1e5 s, as in Fig. 8
t = exp10.(range(1, 5, length = 200))
k_600 = argmin(abs.(t .- 600))   # Index of t ≈ 10 min

println("Scenario  ng   Rb multipole   Rb(1e5 s)   Rb(600 s)")

# Scenario a: kp = 0.4 W/mK, half shank spacing 0.025 m
kp_a = 0.4    # Pipe thermal conductivity [W/mK]
s_a  = 0.05   # Shank spacing [m]
ex_a = SingleUTube(rb, ro, ri, s_a, kg, Cg, kp_a, Cp)
Rb_a = delta_circuit(ex_a, ks, V̇, fluid; films = false).Rb

Rb_a_ng1  = transient_resistance(t, ex_a, ks, Cs, V̇, fluid; films = false,
    settings = TRCMSettings(ng = 1))
Rb_a_ng2  = transient_resistance(t, ex_a, ks, Cs, V̇, fluid; films = false,
    settings = TRCMSettings(ng = 2))
Rb_a_ng6  = transient_resistance(t, ex_a, ks, Cs, V̇, fluid; films = false,
    settings = TRCMSettings(ng = 6))
Rb_a_ng12 = transient_resistance(t, ex_a, ks, Cs, V̇, fluid; films = false,
    settings = TRCMSettings(ng = 12))
@printf("   a       1   %10.5f   %9.5f   %9.5f\n", Rb_a, Rb_a_ng1[end], Rb_a_ng1[k_600])
@printf("   a       2   %10.5f   %9.5f   %9.5f\n", Rb_a, Rb_a_ng2[end], Rb_a_ng2[k_600])
@printf("   a       6   %10.5f   %9.5f   %9.5f\n", Rb_a, Rb_a_ng6[end], Rb_a_ng6[k_600])
@printf("   a      12   %10.5f   %9.5f   %9.5f\n", Rb_a, Rb_a_ng12[end], Rb_a_ng12[k_600])

# Scenario e: kp = 1.0 W/mK, half shank spacing 0.055 m
kp_e = 1.0    # Pipe thermal conductivity [W/mK]
s_e  = 0.11   # Shank spacing [m]
ex_e = SingleUTube(rb, ro, ri, s_e, kg, Cg, kp_e, Cp)
Rb_e = delta_circuit(ex_e, ks, V̇, fluid; films = false).Rb

Rb_e_ng1  = transient_resistance(t, ex_e, ks, Cs, V̇, fluid; films = false,
    settings = TRCMSettings(ng = 1))
Rb_e_ng2  = transient_resistance(t, ex_e, ks, Cs, V̇, fluid; films = false,
    settings = TRCMSettings(ng = 2))
Rb_e_ng6  = transient_resistance(t, ex_e, ks, Cs, V̇, fluid; films = false,
    settings = TRCMSettings(ng = 6))
Rb_e_ng12 = transient_resistance(t, ex_e, ks, Cs, V̇, fluid; films = false,
    settings = TRCMSettings(ng = 12))
@printf("   e       1   %10.5f   %9.5f   %9.5f\n", Rb_e, Rb_e_ng1[end], Rb_e_ng1[k_600])
@printf("   e       2   %10.5f   %9.5f   %9.5f\n", Rb_e, Rb_e_ng2[end], Rb_e_ng2[k_600])
@printf("   e       6   %10.5f   %9.5f   %9.5f\n", Rb_e, Rb_e_ng6[end], Rb_e_ng6[k_600])
@printf("   e      12   %10.5f   %9.5f   %9.5f\n", Rb_e, Rb_e_ng12[end], Rb_e_ng12[k_600])

# Figure 1: Rb(t) of scenarios a (red) and e (black) for several ng, multipole Rb dotted
fig = Figure(size = (800, 500))
ax = Axis(fig[1, 1], xlabel = "Time (s)", ylabel = "Rb(t) (mK/W)", xscale = log10,
    title = "Modified TRCM, Pasquier & Marcotte (2012) scenarios a and e")
hlines!(ax, [Rb_a], color = :firebrick, linestyle = :dot)
lines!(ax, t, Rb_a_ng1, linewidth = 2, color = :firebrick, linestyle = :dash, label = "a, ng = 1")
lines!(ax, t, Rb_a_ng2, linewidth = 2, color = :firebrick, linestyle = :dashdot,
    label = "a, ng = 2")
lines!(ax, t, Rb_a_ng6, linewidth = 2, color = :firebrick, label = "a, ng = 6")
lines!(ax, t, Rb_a_ng12, linewidth = 1, color = :firebrick, label = "a, ng = 12")
hlines!(ax, [Rb_e], color = :black, linestyle = :dot)
lines!(ax, t, Rb_e_ng1, linewidth = 2, color = :black, linestyle = :dash, label = "e, ng = 1")
lines!(ax, t, Rb_e_ng2, linewidth = 2, color = :black, linestyle = :dashdot, label = "e, ng = 2")
lines!(ax, t, Rb_e_ng6, linewidth = 2, color = :black, label = "e, ng = 6")
lines!(ax, t, Rb_e_ng12, linewidth = 1, color = :black, label = "e, ng = 12")
axislegend(ax, position = :lt)
display(fig)
