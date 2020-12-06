#!/bin/bash

mkdir ../block_annotations_gffs
for block in $(seq 1 262); do
file=$(find ../block_annotations/* -name "*_$block.gff" | head -n1)
cp $file ../block_annotations_gffs
done
