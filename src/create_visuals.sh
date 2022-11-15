#!/bin/bash
# Creates .pdf, .png images and then .mp4 video files 
# from CellModeller simualtion .pickle output files
# New subdirectory in scratch/ is created to store
# the generated files
# 
# As an argument, data/ subdirectory path, width and height of simulation page size
# centered areound (0,0) must be passed to the script
# $ bash src/create_visuals.sh <pickle-dir-path> <page-width> <page-heigth>
# Example: $ bash src/create_visuals.sh data/simpleGrowth2D 30 3

if [ "$#" != 3 ]; then
        echo "Check passed argumets"
        echo "Usage:"
        echo "$ bash src/create_visuals.sh <pickle-dir-path> <page-width> <page-heigth>"
        echo "Example:"
        echo "$ bash src/create_visuals.sh data/simpleGrowth2D 35 6"
        exit 0
fi

DATA_DIR=$1
PAGE_WIDTH=$2
PAGE_HEIGHT=$3

OUT_DIR=scratch/$(basename ${DATA_DIR})
FRAME_DIR=scratch/frames/$(basename ${DATA_DIR})

# Create directories to store resutls
mkdir -p ${OUT_DIR}
mkdir -p ${FRAME_DIR}

# Moving to data directory
pushd ${DATA_DIR}

# Run Draw2DPDF to generate pdf files
for FILE in $( ls *.pickle ); do
    echo Processing: $FILE
    python3 $HOME/apps/cellmodeller/Scripts/Draw2DPDF.py ${PAGE_WIDTH} ${PAGE_HEIGHT} ${FILE}> /dev/null
done

# Convert and resize etc. pdf files into png
for FILE in $( ls *.pdf ); do
    NAME=`basename $FILE .pdf`
    convert \
        -colorspace RGB \
        -verbose        \
        -density 150    \
        $NAME.pdf       \
        $NAME.png
done

# Run ffmpeg to generate video file
echo "Compiling video:"
ffmpeg -framerate 7 -i %*.png -vf scale=1920:1080 -r 24 video.mp4

# Move back to root
popd

# Move images adn video to the results directory
mv ${DATA_DIR}/{*.mp4} ${OUT_DIR}
mv ${DATA_DIR}/{*.pdf,*.png} ${FRAME_DIR}

