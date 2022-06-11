# Polysome_detect

MATLAB script to annotate polysomes in cellular tomograms based on ribosome coordinate and rotation information from subtomo analysis/averaging.

This script was first described in the manuscript by Xue *et al*, 2021. Preprint on [bioRxiv](https://doi.org/10.1101/2021.12.18.473270).

## Basic concepts

1. A motivelist(motl) table to combine and store all relevant information for each ribosomes, such as coordinates, rotations, class info, *etc.* 

2. Use ribosomes' mRNA exit to entry distance to define polysomes and their sequence.

3. Extract quantitative information on polysome arrangement, such as overall distribution, relative orientations, *etc.*
 
## How to run

Require MATLAB 2016b or equivalent. 

Put the two associated fromEuler_RELION.m and rotmat2eulang.m files in MATLAB search Path.

Open the polysome_detect.m, and run it step by step (select and F9).
