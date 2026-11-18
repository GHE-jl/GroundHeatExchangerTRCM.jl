# References

## Thermal resistance and capacity models

- **Bauer, D., Heidemann, W., Müller-Steinhagen, H., & Diersch, H.-J. G.** (2011). Thermal
  resistance and capacity models for borehole heat exchangers. *International Journal of Energy
  Research*, 35(4), 312–320. <https://doi.org/10.1002/er.1689>
  TRCMs of coaxial, single and double U-tube exchangers; grout resistances from `Rb` and `Ra`
  (Eqs. 21, 28-29, 34-38) and location of the grout capacity.
- **Bauer, D., Heidemann, W., & Diersch, H.-J. G.** (2011). Transient 3D analysis of borehole heat
  exchanger modeling. *Geothermics*, 40(4), 250–260.
- **Pasquier, P., & Marcotte, D.** (2012). Short-term simulation of ground heat exchanger with an
  improved TRCM. *Renewable Energy*, 46, 92–99. <https://doi.org/10.1016/j.renene.2012.03.014>
  The modified TRCM (sub-resistances and sub-capacities) used for every layer.
- **Pasquier, P., & Marcotte, D.** (2014). Joint use of quasi-3D response model and spectral method
  to simulate borehole heat exchanger. *Geothermics*, 51, 281–299.
  <https://doi.org/10.1016/j.geothermics.2014.02.001>
  The quasi-3D response model and its boundary conditions (Eqs. 13-18), [`trcm`](@ref).
- **Nguyen, A., Pasquier, P., & Marcotte, D.** (2015). Thermal resistance and capacity model for
  standing column wells operating under a bleed control. *Renewable Energy*, 76, 743–756.
  <https://doi.org/10.1016/j.renene.2014.11.080>
  The standing column well network, [`StandingColumn`](@ref).
- **Minaei, A., & Maerefat, M.** (2017). Thermal resistance capacity model for short-term borehole
  heat exchanger simulation with non-stiff ordinary differential equations. *Geothermics*, 70,
  260–270.
  An alternative to the stiff TRCM; the sparse-Jacobian stiff solver used here makes the stiffness
  harmless instead.
- **Pasquier, P., Zarrella, A., & Labib, R.** (2018). Application of artificial neural networks to
  near-instant construction of short-term g-functions. *Applied Thermal Engineering*, 143, 910–921.
  <https://doi.org/10.1016/j.applthermaleng.2018.07.137>
  The ANN surrogate of the quasi-3D model that this package aims to make unnecessary; its 7-day
  horizon sets the default contact time of `outlet_transfer_function`.

## Steady-state resistances

- **Javed, S., & Spitler, J.** (2017). Accuracy of borehole thermal resistance calculation methods
  for grouted single U-tube ground heat exchangers. *Applied Energy*, 187, 790–806.
- **Claesson, J., & Javed, S.** (2019). Explicit multipole formulas and thermal network models for
  calculating thermal resistances of double U-pipe borehole heat exchangers. *Science and
  Technology for the Built Environment*, 25(8), 980–992.
  <https://doi.org/10.1080/23744731.2019.1620565>
  The double U-tube Δ-circuit (Eq. 26), [`DoubleUTube`](@ref).
- **Lamarche, L.** (2021). Analytic models and effective resistances for coaxial ground heat
  exchangers. *Geothermics*, 97, 102224. <https://doi.org/10.1016/j.geothermics.2021.102224>
- **Lamarche, L.** (2023). *Fundamentals of Geothermal Heat Pump Systems: Design and Application.*
  Springer Nature Switzerland.

## Long-term ground response

- **Claesson, J., & Javed, S.** (2011). An analytical method to calculate borehole fluid
  temperatures for time-scales from minutes to decades. *ASHRAE Transactions*, 117(2), 279–288.
  The finite line source joined to the TRCM by `outlet_transfer_function` of
  GroundHeatExchanger.jl.
- **Cimmino, M., & Bernier, M.** (2014). A semi-analytical method to generate g-functions for
  geothermal bore fields. *International Journal of Heat and Mass Transfer*, 70, 641–650.
  <https://doi.org/10.1016/j.ijheatmasstransfer.2013.11.037>

## Borehole pairs and borefields

- **Rose, C., Pasquier, P., Nguyen, A., & Labib, R.** (2026). A unified framework for the
  construction of short- and long-term g-functions for parallel-connected ground heat exchangers.
  *International Journal of Heat and Mass Transfer*, 261, 128565.
  <https://doi.org/10.1016/j.ijheatmasstransfer.2026.128565>
  Single-borehole and borehole-to-borehole response functions, [`trcm_pairs`](@ref).
- **Dusseault, B., Pasquier, P., & Marcotte, D.** (2018). A block matrix formulation for efficient
  g-function construction. *Renewable Energy*, 121, 249–260.
  BC-IV borefield g-functions from the pair responses (`bc_II` of GroundHeatExchanger.jl).
