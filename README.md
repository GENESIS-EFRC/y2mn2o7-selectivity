# y2mn2o7-selectivity
Repository containing all code and data associated with manuscript "Selectivity in materials synthesis via local chemical potentials in hyperdimensional phase space".


### Code (y2mn2o7_selectivity)

To install the code as a Python package using pip, first navigate to the repository 
folder in your terminal window and then type:
 
    pip install -r requirements.txt
    pip install .
    
Computed data was acquired from the Materials Project (http://materialsproject.org
), version 2021_03_22. A copy of this data can be accessed in the "computed_entries
.json" file within the code folder.

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

