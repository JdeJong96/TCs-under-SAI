## TC (seed) tracker
Main repository for tracking tropical cyclones in CESM data on the Snellius HPC (SLURM setup).

Tracking of tropical cyclones in CESM data broadly follows the following two-stage workflow. First, positions of TC candidates are stored for each day in the simulation by `src/RVmax_finder.py`. Then the TC candidates are stitched together, filtered and supplied with additional data by `src/Tracking_TC_RV.py`. In the first stage, only selection criteria acting at a specific time are applied (e.g. on current intensity or proximity to other candidates), while criteria acting over multiple timesteps are applied only in the second stage, e.g. stitching rules or requirements that must be met within the TC lifetime. 

### General
This directory contains the TC tracker source code and generated TC (seed) track data used for the SAI impact study. The directory is structured as follows:

- `additional` early attempts for improving TC precipitation data, not used in final study
- `doc` list of useful references
- `jobs` generated track datasets, a set of jobs created with the final version of the tracker is:
  - `Tracking_TC_RV.24hrext2` TC tracks with 24 hour extension
  - `Tracking_TC_RV.48hrext2` TC tracks with 48 hour extension
  - `Tracking_TC_RV.infext2` TC tracks with indefinite extension, used for the TC precipitation figure for better land estimates
  - `Tracking_TC_RV.seeds2` TC seed tracks, used for all TC seed analyses
  - `Tracking_TC_RV` TC tracks regular, used for all TC analyses
  Here, the extensions (...`ext2`) mean that if a TC track is stopped by regular criteria, the closed circulation criterium is lifted to append more candidate points to the end point of the regular track. TC seed tracks are constructed similarly to TC tracks, but have less restrictive stage 2 criteria.
- `src` source code
  - `RVmax_finder.original.py` preliminary stage 1 tracker file taken from van Westen et al. (see doc)
  - `RVmax_finder.py` stage 1 tracker file
  - `RVmax_finder.sh` jobscript for `RVmax_finder.py`
  - `Tracking_TC_RV.original.py` preliminary stage 2 tracker file taken from van Westen et al.
  - `Tracking_TC_RV.py` stage 2 tracker file
  - `Tracking_TC_RV.sh` jobscript for `Tracking_TC_RV.py`
  - `analysis` convenience pointer for imports to the `analysis/src` directory
  - `make_grid_file.py` script for generating a grid specification file
- `README.md` this readme 
- `environment.yml` environment specification 

___

### Run instructions

First, make sure the relevant packages are available. A conda environment can be created using the `environment.yml` file.

Store position of TC candidates with `src/RVmax_finder.py`
1. Open `RVmax_finder.py` in an editor, e.g. vim 
2. Set the right paths in the top part of the script. Note that arguments to the jobscript ($1, $2, ...) are passed on to the python script for convenience, where they are used to construct paths. So `sbatch RVmax_finder.sh arg1 arg2` will call `python RVmax_finder.py arg1 arg2` and store results in `/some/path.arg1.arg2/RV_Max/`. In the current study, `arg1` typically denotes the experiment (REF, RCP, SAI), and `arg2`, denotes the ensemble member (1-6), such that each simulation can be processed separately.  
3. Set `dry_run=True`
4. Save and run with: `>> python RVmax_finder.py arg1 arg2`
5. Check if the output looks OK
6. Set `dry_run=False` in `RVmax_finder.py`
7. Edit amount of cores and other parameters in `RVmax_finder.sh` and submit with sbatch.
8. After run, move logfile to the output folder
9. Optionally place a copy of `RVmax_finder.py` in the output folder (for archival purposes)
The jobscript may be easily resubmitted at any time, the program checks which files are already written and skips those. To recreate all files, the old output files first need to be deleted or a new output folder should be used.

Stitch TC candidate locations in time and apply selection criteria
1. Open `Tracking_TC_RV.py` in an editor, e.g. vim
2. Set the right paths and save. Note that the script takes two arguments (`sys.argv`) that are used in the paths.
3. Open `Tracking_TC_RV.sh` and create a line like 
`>> python -u Tracking_TC_RV.py arg1 arg2 > Tracking.log.arg1.arg2 &`   
for each combination of input arguments that are used in the paths. The tracking script can only run on one core, but using this approach multiple instances of the tracker are launched that each target a different simulation depending on the arguments and write to a unique logfile. Note that the `&` means that python runs as a background process (i.e. does not block), therefore a wait statement at the end of the jobscript is needed to prevent the job from exiting immediately.
4. Submit `Tracking_TC_RV.sh`. A netCDF output file is only written at the end of each python process, progress can be checked by viewing the separate logfiles (`Tracking.log.arg1.arg2`) for each process. If there is no job output, check the main logfile (`Tracking_TC_RV.log${jobid}`) for any errors.
