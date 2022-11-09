#!/bin/bash
# Creates .pdf, .png images and then .mp4 video files 
# from CellModeller simualtion .pickle output files
# New subdirectory in scratch/ is created to store
# the generated files
# 
# As an argument, a subdirectory path form data/
# must be passed to the script


DATA_DIR=$1
OUT_DIR=scratch/$(basename ${DATA_DIR})
# Create directory to store resutls
mkdir -p ${OUT_DIR}

# Moving to data directory
pushd ${DATA_DIR}

# Run modified CellModeller script
bash ${HOME}/apps/cellmodeller/Scripts/video.sh video.mp4 > /dev/null

# Move back to root
popd

# Move images adn video to the results directory
mv ${DATA_DIR}/{*.pdf,*.png,*.mp4} ${OUT_DIR}
