# Photovoltaic_Performance_Simulator
An open-source, Jupyter-based tool with graphical user interface for simulating photovoltaic performance under arbitrary illumination conditions. This tool is based on the theory and methods outlined in the following works:
<ol>
  <li>
  <i>'The Thermodynamic Limit of Indoor Photovoltaics Based on Energetically-Disordered Molecular Semiconductors'</i>
  <br>
  Kay, A.M., Fitzsimons, M.E., Burwell, G., Meredith, P., Armin, A., & Sandberg, O.J.,
  <br>
  <i>Sol. RRL</i>, 7(18), article 2300277, (2023) 
  <br>
  Available at: https://doi.org/10.1002/solr.202300277
  </li>
  <li>
Seunarine K., Haymoor, Z., Spence, M., Burwell, G., Kay, A.M., Meredith, P., Armin, A., & Carnie, M.J. (2024) _'Light Power Resource Availability for Energy Harvesting Photovoltaics for Self-Powered IoT'_
  </li>
  <li>
    Kay, A.M., Ahmed, S.N., Burridge, N., Riley, D.B., Armin, A., Sandberg, O.J., Haymoor, Z., Carnie, M.J., Meredith, P., & Burwell, G. (2025) ‘Optimising Photovoltaic Modules for Indoor Energy-Harvesting Systems’, <i>Journal of Physics: Energy</i>, 7(3), article 035019. Available at: https://doi.org/10.1088/2515-7655/ade38b
  </li>
  <li>
    Kay, A.M., Riley, D.B., Burwell, G., & Meredith, P. (2025) ‘Thermodynamic Principles for Optimizing Multi-Junction Photovoltaics—Exemplified for Perovskite-Based Indoor Photovoltaics’, <i>APL Energy</i>, 3, 036102. Available at: https://doi.org/10.1063/5.0266374
  </li>
  <li>
    Kay, A.M., Riley, D.B., Meredith, P., & Burwell, G. (2026) ‘Validating and Comparing Reported Indoor Photovoltaic Performance’, <i>EES Solar</i>, Available at: https://doi.org/10.1039/D6EL00150E
  </li>
</ol>

Other resources developed by the Sustainable Advanced Materials (Sêr SAM) research group are available at https://www.advmaterswansea.co.uk/resources.
Sêr SAM homepage: https://www.advmaterswansea.co.uk/. More information on the Sêr Cymru programme available at https://gov.wales/ser-cymru. Further information on the Application-Targeted and Integrated Photovoltaics (ATIP) programme that funded this work is available at https://www.swansea.ac.uk/science-and-engineering/research/atip/, as well as the Centre for Integrative Materials (CISM) at https://www.cism-swansea-semiconductors.co.uk/.

Updates:

<ul>
    <li>
_From Version 3.1 of the tool onwards, two additional simulation modes are available: (i) Indoor PV Validation, which accompanies the EES Solar publication by Kay et al. and can be used to validate indoor photovoltaic measurements using thermodynamic measurements, and (ii) Shunt Resistance Estimation, which accompanies the pre-print "Computing at the Edge Enabled by Indoor Photovoltaics and Ferroelectrics". In addition, more light sources have been added, and more systems stored in the database. Furthermore, in response to the rebranding of the National Renewable Energy Laboratory (NREL) to the National Laboratory of the Rockies (NLR), the API-requesting code for geo-meteorological modelling has had to be revised. Whilst revising this code, additional satellite resources and a "User-defined" deployment location option have been added, alongside a two-dimensional (temperature and irradiance) interpolation algorithm to speed up analysis of all six irradiance types._
  </li>
  <li>
_From Version 3.0 of the tool onwards, the underlying script has been entirely reworked, such that the input parameters are taken from an Excel file ("Parameter_Values.xlsx"), rather than the current values of the input widgets. This means that the "PV_Sim_3.0.py" script can be run independently of the GUI in Spyder._
  </li>
  <li>
_From Version 2.1 of the tool onwards, geographical+meteorological photovoltaic (GeoPV) modelling capabilities have been included into the tool, enabling location-dependent modelling of different photovoltaic technologies._
  </li>
  <li>
_From Version 2.1 of the tool onwards, a normal-incidence optical transfer matrix model has been included in the tool to simulate the photovoltaic external quantum efficiency of organic semiconductor-based photovoltaics using the optical constants of the consituent layers. This model (based on the work of Pettersson et al. [1999] and Peumans et al. [2002]) feeds directly into the thermodynamic limit calculation model._
  </li>
  <li>
_From Version 2.0 of the tool onwards, an alternative algorithm (Brent's method) is used to find the open-circuit voltage and maximum power point parameters, speeding up the code by around a factor of four. In addition, shunt and series resistances have been incorporated into the ideal diode model used to simulate the photovoltaic figures-of-merit._
  </li>
</ul>
