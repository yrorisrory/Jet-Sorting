#!/bin/bash

input_dir="jet_sorting/final/inputs"
output_dir="jet_sorting/final/outputs"

# files=(4-1 5-1 5-2 6-1 6-2 7-1 7-2 7-3 8-1 8-2 8-3)
files=(7-3) 

for f in "${files[@]}"; do
    # infile="file:${input_dir}/HtHt/${f}.root"
    # outfile="${output_dir}/HtHt/clustered_${f}.root"
    # cmsRun jet_sorting/final/python/dynamic_config.py inputFiles=${infile} outputFile=${outfile}

    # infile="file:${input_dir}/ZtZt/${f}.root"
    # outfile="${output_dir}/ZtZt/clustered_${f}.root"
    # cmsRun jet_sorting/final/python/dynamic_config.py inputFiles=${infile} outputFile=${outfile}

    infile="file:${input_dir}/WbWb/${f}.root"
    outfile="${output_dir}/WbWb/clustered_${f}.root"
    cmsRun jet_sorting/final/python/dynamic_config.py inputFiles=${infile} outputFile=${outfile}
done