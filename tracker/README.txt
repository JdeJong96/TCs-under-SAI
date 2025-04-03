To track tropical cyclones in CESM data, follow this procedure:

First, positions of TC candidates are stored for each day in the simulation by RVmax_finder. Then the TC candidates are stitched together, filtered and supplied with additional data by Tracking_TC_RV.

============ Store position of TC candidates ===============
1) Open RVmax_finder.py in an editor, e.g. vim 
2) Set the right paths in the top part of the script. Note that arguments to the jobscript ($1, $2, ...) are passed on to the python script for convenience, where they are used to construct paths. So sbatch RVmax_finder.sh arg1 arg2 will call python RVmax_finder.py arg1 arg2 and store results in /some/path.arg1.arg2/RV_Max/.
3) Set dry_run=True
4) Save and run with: >> python RVmax_finder.py arg1 arg2
5) Check if the output looks OK
6) Set dry_run=False in RVmax_finder.py
7) Edit amount of cores and other parameters in RVmax_finder.sh and submit with sbatch.
8) After run, move logfile to the output folder
9) Optionally place a copy of RVmax_finder.py in the output folder (for future reference)
The jobscript may be easily resubmitted at any time, the program checks which files are already written and skips those. To recreate all files, the old output files first need to be deleted or a new output folder should be used.

=========== Stitch TC candidate locations in time and apply selection criteria ============
1) Open Tracking_TC_RV.py in an editor, e.g. vim
2) Set the right paths and save. Note that the script takes two arguments (sys.argv) that are used in the paths.
3) Open Tracking_TC_RV.sh and create a line like 
>> python -u Tracking_TC_RV.py arg1 arg2 > Tracking.log.arg1.arg2 & 
for each combination of input arguments that are used in the paths. The tracking script can only run on one core, but using this approach multiple instances of the tracker are launched that each target a different simulation depending on the arguments and write to a unique logfile. Note that the & means that python runs as a background process (i.e. does not block), therefore a wait statement at the end of the jobscript is needed to prevent the job from exiting immediately.
4) Submit Tracking_TC_RV.py. A netCDF output file is only written at the end of each python process, progress can be checked by viewing the separate logfiles (Tracking.log.arg1.arg2) for each process. If there is no job output, check the main logfile (Tracking_TC_RV.log${jobid} for any errors.
