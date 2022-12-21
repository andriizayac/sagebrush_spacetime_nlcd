#### Authors:  

Andrii Zaiats, Megan E. Cattau , David S. Pilliod, Rongsong Liu, Juan M. Requena-Mullor, T. Trevor Caughlin  

---  

#### Overview:

This repository stores the scripts that were used for the article "Forecasting natural regeneration of sagebrush after wildfires using demographic models and spatial matching".  

___  

The repository is organized into folders that contain tabular data, figures,  modeling results, and R scripts.

Folders:\
    - **data** and **outputs**: includes data used for the analysis in .rds format. These data are subsets of the databases referenced in the manuscript after the geospatial filtering was applied. Loading data directly may save geoprocessing time.\
    - **figures**: includes .png/.pdf figures and tables in .csv format.
    
R scripts:\
    - *helper_fns.R*: includes user functions defined in the global environment and called into R environment.  
    - *figures.R*: generates figures for the manuscript.
    
Data extraction and pre-processing. To avoid downloading and pre-processing raw data that could take up to several hours and will require specifying relative file paths, the output from the following scripts is stored in the **data** folder, or can be downloaded from (https://drive.google.com/drive/folders/1e52NiiFp9u8TpgnclTXbSb11DSmuJ-xd?usp=share_link). However, as a reference the pre-processing scripts are included.  
    - *gee_scrip.txt*: source file for Google Earth Engine data extraction.  
    - *site_pool_selection.R*: performs the selection of wildfire polygons using LTDL data layers (https://doi.org/10.5066/P98OBOLS), Bureau of Land Management mask (https://www.sciencebase.gov/catalog/item/4fc8f940e4b0bffa8ab25a65), and  Wildfire Database (https://doi.org/10.5066/P9Z2VVRT).  
    - *raster_processing.R*: extracts raster data from the Homer et al. 2020 dataset and biophysical covariates for each wildfire polygon.  
    - *pixelmatching.R*: applies unsupervised classification of the pixels within each wildfire.  
    
Data analysis and results.\
    - *null_fit.R*: fits regional model.\
    - *model_fit.R*: fits wildfire- and cluster-level models.\
    - *model_pred.R*: generates predictions from each modeling schema and saves as an intermediate product.\
    - *prederr_sensitivity.R*: combines predictions and generates error matrices.
    

    

