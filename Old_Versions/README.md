# Photovoltaic_Performance_Simulator
An open-source, Jupyter-based tool with graphical user interface for simulating photovoltaic performance under arbitrary illumination conditions. This tool is based on the theory outlined in:

Kay, A.M., Fitzsimons, M.E., Burwell, G., Meredith, P., Armin, A. and Sandberg, O.J. (2023), The Thermodynamic Limit of Indoor Photovoltaics Based on Energetically-Disordered Molecular Semiconductors. Sol. RRL. Accepted Author Manuscript. https://doi.org/10.1002/solr.202300277

Other resources developed by the Sustainable Advanced Materials (Sêr SAM) research group are available at https://www.advmaterswansea.co.uk/resources.
Sêr SAM homepage: https://www.advmaterswansea.co.uk/. More information on the Sêr Cymru programme available at https://gov.wales/ser-cymru. Further information on the Application-Targeted and Integrated Photovoltaics (ATIP) programme that funded this work is available at https://www.swansea.ac.uk/science-and-engineering/research/atip/, as well as the Centre for Integrative Materials (CISM) at https://www.cism-swansea-semiconductors.co.uk/.

Updates:

_From Version 2.1 of the tool onwards, geographical+meteorological photovoltaic (GeoPV) modelling capabailities have been included into the tool, enabling location-dependent modelling of different photovoltaic technologies._

_From Version 2.1 of the tool onwards, a normal-incidence optical transfer matrix model has been included in the tool to simulate the photovoltaic external quantum efficiency of organic semiconductor-based photovoltaics using the optical constants of the consituent layers. This model (based on the work of Pettersson et al. [1999] and Peumans et al. [2002]) feeds directly into the thermodynamic limit calculation model._

_From Version 2.0 of the tool onwards, an alternative algorithm (Brent's method) is used to find the open-circuit voltage and maximum power point parameters, speeding up the code by around a factor of four. In addition, shunt and series resistances have been incorporated into the ideal diode model used to simulate the photovoltaic figures-of-merit._
