### Analysis code and data
In this part of the repository, source code, processed data, figures and notebooks are stored belonging to the analysis of TC change under SAI in the CESM simulations.

___

#### Contents:
- `figures` output figures of the analysis 
- `jobs` processed climate data. Typically, each directory within `jobs` contains source code (main script: name.py, job script: name.sh, and script for merging results from all ensemble members: merge.sh), logs and data for a specific environment variable, indicated by the directory name. The data typically comprises interannual mean monthly mean values for each experiment, but may include other relevant statistics.  
- `notebooks` jupyter lab notebooks in which the scientific figures are produced, typically ordered per variable of interest
-   - `paper.ipynb` main figures/tables quantifying TC change
    - `figures.ipynb` track statistics
    - `link_seeds_TCs.ipynb` main TC seed figures
    - `expand_tracker_data.ipynb` TC precipitation figures
- `src` collection of modules with commonly used computational routines
-   - `interpolate.py` vertical interpolation from hybrid to pressure levels
    - `load_SAIdata.py` convenience module to load data from all experiments
    - `physics.py` physical computations
    - `tracks.py` computations for TC (seed) track datasets

