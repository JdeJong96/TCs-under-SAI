#!/bin/bash

module purge
module load 2025
module load CDO/2.5.2-gompi-2025a

tag=hres.ref.1
files=($(cat $TMPDIR/files.$tag))
outfile=$TMPDIR/outfile.$tag.nc
echo "writing to $outfile"

cdo select,name=SST,timestep=1 ${files[0]} $outfile # dummy

# write SST as monthly mean timeseries
cdo -O -v -f nc4 --chunkspec t=100,x=360,y=240 -ymonmean -selyear,2003/2007 -mergetime [ -shifttime,-12hours -selname,SST : ${files[@]} ] $outfile

cdo showchunkspec $outfile
cdo showformat $outfile
#rm $outfile
