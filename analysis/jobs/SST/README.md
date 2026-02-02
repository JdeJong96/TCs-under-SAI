In this job, the interannual mean monthly mean SST and monthly sum of days with SST >= 25C are computed.

Corresponding notebook: TCs-under-sai/analysis/notebooks/SST.ipynb

Workflow:
1) submit SST.sh to run computations for each ensemble member separately as an array job. 
2) submit regridder.sh to regrid to the rectangular CAM grid (and prevent some issues with subtle grid differences)
3) run merge.sh to compute ensemble averages and merge results for all experiments in one file (`data/SST.nc`)
4) submit landmask.sh to compute the land mask and regrid it to the CAM grid (`landmask.f02_t12.nc`)
