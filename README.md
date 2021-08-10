# y2mn2o7-selectivity
Repository containing all code and data associated with manuscript "Selectivity in yttrium manganese oxide synthesis via local chemical potentials in hyperdimensional phase space" 
(2021).

### Code (y2mn2o7_selectivity)

To install the code as a Python package using pip, we recommend **creating a new conda 
environment with Python 3.9**.

After activating the conda environment and cloning this repository, first 
navigate to the repository folder in your terminal window and then type:
 
    pip install -r requirements.txt
    pip install .

A standalone copy of the necessary scripts in the reaction-network package was 
included to ensure future compatability with generated entries/data in this work. To 
install this package, type in your terminal window:

    cd y2mn2o7_selectivity/reaction-network
    pip install -e .

The plotting script notebooks should now be usable after beginning a Jupyter 
notebook server.

NOTE: Computed data was acquired from the Materials Project (http://materialsproject.org
), version 2021.05.13. A copy of this data can be accessed in the "computed_entries_C_Cl_K_Li_Mn_Na_O_Y_08_08_21.json.gz" file within the y2mn2o7_selectivity folder.

### exp_data
##### AnalysisData
- Contains .csv files for all data represented in Experimental Figures from Main
 Article and SI.

##### ControlPXRD
- Contains all raw PXRD data, Rietveld Refinements, background, and difference
 patterns for all control experiments performed.

##### INP
- Contains all Bruker TOPAS Input files used in the sequential Rietveld refinements
 for all experiments performed at 17-BM-B at Argonne National Lab in October of 2018.

##### Metadata
- Contains all .metadata for experiments performed at 17-BM-B at Argonne National Lab
 in October of 2018.

##### TIFF (not yet uploaded)
- Raw 2-D plate detector images from experiments performed at 17-BM-B at Argonne
 National Lab in October of 2018.

##### XYE
- Integrated X-ray diffraction patterns calculated from images in /TIFF using GSAS-II

