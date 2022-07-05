# B2RIO: "Brain to Reverse Inference Ontology" package

## How to use this package

### Instalation
Clone this repository and install the application with `python setup.py install`


### Running the analysis
The only necessary parameter is the path to a nifti file with the image to be analyzed. We expect a 4D image with the region label in the voxels we want to evaluate and zero in all other voxels.

E.g: `b2rio --brain_path nifti_image.nii.gz`.

Apart from this, there are a number of other parameters that can be configured:

- `--n_folds`: Number of folds use to sample studies (default=150)
- `--resample`: Resample of the provided image (default=1, no resample)
- `--frac_sample`: Percentage of studies to use during sample (default=0.7)
- `--radius`: When matching voxels between the provided image and neurosynth activations, this radius is used to limit the maximum distance between co-activated voxels. (default=4)
- `--folder_results`: Folder to save results obtained (default='./')