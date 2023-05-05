# Description

## captureAreas
This is the output of `Vistoseg`’s function `splitslide`, which splits the image of the whole Visium slide into different individual capture areas. The files of `captureAreas`  will be used as the input to the PNN segmentation algorithm.

## loupe
The folder contains the output from `Loupe` (https://www.10xgenomics.com/products/loupe-browser). 
There are two files for each Visium sample, i.e. `.png`, `.json`.
The `.png` shows the the results of manual annotation of `in-tissue` spots. The `.json` stores the coordinate information of spots. Both files will be piped into `spaceranger`.