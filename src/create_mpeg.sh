#!/bin/bash

for IDX in "$@"; do
    DIR_BASENAME=$(basename $IDX)
    pushd scratch/frames/$DIR_BASENAME
    rm -f *.pdf
    ffmpeg -y -i %*.png -codec:v mpeg2video -qscale:v 2 -codec:a mp2 -b:a 192k -vf "setpts=10*PTS" video.mpeg
    popd
    mv scratch/frames/$DIR_BASENAME/video.mpeg scratch/$DIR_BASENAME
done