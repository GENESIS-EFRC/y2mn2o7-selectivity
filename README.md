# thermochemical-connectivity
Repository containing all code and data associated with manuscript "Selectivity in
 materials synthesis via thermochemical connectivity in hyperdimensional phase space".


### Code (thermochemical_connectivity)

To install the code, first navigate to the repository folder in your terminal window
, and then type:
 
    pip install -r requirements.txt
    pip install .
    

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


### computed_data
Computed data was acquired from the Materials Project (http://materialsproject.org
), version 2020_09_08. A copy of this data can be accessed in the computed_data
 folder within the package.
