# Photovoltaic Performance Simulator (PV-Simulator)
An open-source, Jupyter-based tool with graphical user interface for simulating photovoltaic performance under arbitrary illumination conditions. This tool is based on theory, methods, and measurements described in contributions to the following works:
<ol>
  <li>
    <i>'Quantifying the Excitonic Static Disorder in Organic Semiconductors'</i>,
    <br>
    Kay, A.M., Sandberg, O.J., Zarrabi, N., Li, W., Zeiske, S., Kaiser, C., Meredith, P., & Armin, A.,
    <br>
    <i>Adv. Func. Mater.</i>, 32 (32), article 2113181, (2022) 
    <br>
    DOI: <a href="https://doi.org/10.1002/adfm.202113181">10.1002/adfm.202113181</a>
  </li>
  <li>
    <i>'The Thermodynamic Limit of Indoor Photovoltaics Based on Energetically-Disordered Molecular Semiconductors'</i>,
    <br>
    Kay, A.M., Fitzsimons, M.E., Burwell, G., Meredith, P., Armin, A., & Sandberg, O.J.,
    <br>
    <i>Sol. RRL</i>, 7 (18), article 2300277, (2023) 
    <br>
    DOI: <a href="https://doi.org/10.1002/solr.202300277">10.1002/solr.202300277</a>
  </li>
  <li>
    <i>'Light Power Resource Availability for Energy Harvesting Photovoltaics for Self-Powered IoT'</i>,
    <br>
    Seunarine K., Haymoor, Z., Spence, M., Burwell, G., Kay, A.M., Meredith, P., Armin, A., & Carnie, M.J.,
    <br>
    <i>J. Phys. Energy</i>, 6 (1), article 015018, (2024)
    <br>
    DOI: <a href="https://doi.org/10.1088/2515-7655/ad1764">10.1088/2515-7655/ad1764</a>
  </li>
  <li>
    <i>'Wide-Gap Perovskites for Indoor Photovoltaics'</i>,
    <br>
    Burwell, Zeiske, S., Caprioglio, P., Sandberg, O.J., Kay, A.M., Farrar. M.D., Kim., Y.R., Snaith, H.J., Meredith, P., & Armin, A.,
    <br>
    <i>Sol. RRL</i>, 8 (11), article 2400180, (2024)
    <br>
    DOI: <a href="https://doi.org/10.1002/solr.202400180">10.1002/solr.202400180</a>
  </li>
  <li>
    <i>'On the Performance Limits of Agrivoltaics--From Thermodynamic to Geo-Meteorological Considerations'</i>,
    <br>
    Kay, A.M., Riley, D.B., Sandberg, O.J., Burwell, G., Meredith, P., & Armin, A.,
    <br>
    <i>Sol. RRL</i>, 8 (18), article 2400456, (2024)
    <br>
    DOI: <a href="https://doi.org/10.1002/solr.202400456">10.1002/solr.202400456</a>
  </li>
  <li>
    <i>‘Optimising Photovoltaic Modules for Indoor Energy-Harvesting Systems’</i>,
    <br>
    Kay, A.M., Ahmed, S.N., Burridge, N., Riley, D.B., Armin, A., Sandberg, O.J., Haymoor, Z., Carnie, M.J., Meredith, P., & Burwell, G.,
    <br>
    <i>J. Phys. Energy</i>, 7 (3), article 035019, (2025)
    <br>
    DOI: <a href="https://doi.org/10.1088/2515-7655/ade38b">10.1088/2515-7655/ade38b</a>
  </li>
  <li>
    <i>‘Thermodynamic Principles for Optimizing Multi-Junction Photovoltaics—Exemplified for Perovskite-Based Indoor Photovoltaics’</i>,
    <br>
    Kay, A.M., Riley, D.B., Burwell, G., & Meredith, P.,
    <br>
    <i>APL Energy</i>, 3, article 036102, (2025)
    <br>
    DOI: <a href="https://doi.org/10.1063/5.0266374">10.1063/5.0266374</a>
  </li>
  <li>
    <i>‘Validating and Comparing Reported Indoor Photovoltaic Performance’</i>,
    <br>
    Kay, A.M., Riley, D.B., Meredith, P., & Burwell, G.,
    <br>
    <i>EES Solar</i> (2026)
    <br>
    DOI: <a href="https://doi.org/10.1039/d6el00150e">10.1039/d6el00150e</a>
  </li>
  <li>
    <i>‘Computing at the Edge Enabled by Indoor Photovoltaics and Ferrorelectrics’</i>,
    <br>
    Hoye, R.L.Z., Hellenbrand, M., Estrócio, N., Kay, A.M., Scardina, A., Mezzetti, A., Fina, I., Íñiguez-González, J., Marques, L.S., Matteucci, F., Jia, Q., MacManus-Driscoll, J.L., Offrein, B., Sánchez, F., Noheda, B., Koutsourakis, G., Freitag, M., Burwell, G., Grancini, G., & Silva, J.P.B.,
    <br>
    <i>arXiv Pre-print</i> (2026)
    <br>
    DOI: <a href="https://arxiv.org/abs/2610.05133">https://arxiv.org/abs/2610.05133</a>
  </li>
</ol>

<hr>

<b>Updates</b>:

<ul>
    <li>
<i>From Version 3.1 of the tool onwards, two additional simulation modes are available: (i) Indoor PV Validation, which accompanies the EES Solar publication by Kay et al. and can be used to validate indoor photovoltaic measurements using thermodynamic measurements, and (ii) Shunt Resistance Estimation, which accompanies the pre-print "Computing at the Edge Enabled by Indoor Photovoltaics and Ferroelectrics". In addition, more light sources have been added, and more systems stored in the database. Furthermore, in response to the rebranding of the National Renewable Energy Laboratory (NREL) to the National Laboratory of the Rockies (NLR), the API-requesting code for geo-meteorological modelling has had to be revised. Whilst revising this code, additional satellite resources and a "User-defined" deployment location option have been added, alongside a two-dimensional (temperature and irradiance) interpolation algorithm to speed up analysis of all six irradiance types.</i>
  </li>
  <li>
<i>From Version 3.0 of the tool onwards, the underlying script has been entirely reworked, such that the input parameters are taken from an Excel file ("Parameter_Values.xlsx"), rather than the current values of the input widgets. This means that the "PV_Sim_3.0.py" script can be run independently of the GUI in Spyder.</i>
  </li>
  <li>
<i>From Version 2.1 of the tool onwards, geographical+meteorological photovoltaic (GeoPV) modelling capabilities have been included into the tool, enabling location-dependent modelling of different photovoltaic technologies.</i>
  </li>
  <li>
<i>From Version 2.1 of the tool onwards, a normal-incidence optical transfer matrix model has been included in the tool to simulate the photovoltaic external quantum efficiency of organic semiconductor-based photovoltaics using the optical constants of the consituent layers. This model (based on the work of Pettersson et al. [1999] and Peumans et al. [2002]) feeds directly into the thermodynamic limit calculation model.</i>
  </li>
  <li>
<i>From Version 2.0 of the tool onwards, an alternative algorithm (Brent's method) is used to find the open-circuit voltage and maximum power point parameters, speeding up the code by around a factor of four. In addition, shunt and series resistances have been incorporated into the ideal diode model used to simulate the photovoltaic figures-of-merit.</i>
  </li>
</ul>

<hr>

Other resources developed by the Sustainable Advanced Materials (Sêr SAM) research group are available at https://www.advmaterswansea.co.uk/resources.
Sêr SAM homepage: https://www.advmaterswansea.co.uk/. More information on the Sêr Cymru programme available at https://gov.wales/ser-cymru. Further information on the Application-Targeted and Integrated Photovoltaics (ATIP) programme that funded this work is available at https://www.swansea.ac.uk/science-and-engineering/research/atip/, as well as the Centre for Integrative Materials (CISM) at https://www.cism-swansea-semiconductors.co.uk/.
