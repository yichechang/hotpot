# Snakemake workflow: Crop Rectangular ROI by Human

This workflow allows users to crop rectangular ROIs by mouse
for images in a folder using napari.

## Usage
1. Copy config file `config.yaml` to your working directory.
2. Edit `config.yaml` to specify the path to your image folder,
   the pattern of image file names, and the image channel info.
3. Run the workflow with
   ```
   snakemake -s <path_to_snakefile> \
       --configfile <path_to_copied_config_file> \
       --use-conda \
       --cores <number_of_cores> \
       all
   ``` 
   Include `--conda-prefix <conda_prefix>` if you want to specify
   where conda environments are created. This is useful if you
   don't want to create conda environments on each and every run
   for different projects. 
4. Draw ROIs using napari viewer's shape layer. By default, the
   rectangle shape tool is selected. Press button `next` once 
   you have finished drawing ROI(s) for the current image. If
   no new image is set up after pressing `next`, it means you
   have finished drawing ROIs for all images. Simply close the
   napari window to allow post-processing to start.
5. In the working directory, the results can be found:
   - `results/roilist.csv` contains the ROI coordinates in each 
     image.
   - `results/cropped` directory contains the cropped images.

## Requirements
- `snakemake`
- `mamba` (optional). If prefer using `conda`, include 
  `--conda-frontend conda` in the snakemake command.

The workflow will automatically create a conda environment and
install other required packages for each rule if needed. 