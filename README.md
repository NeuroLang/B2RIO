# B2RIO: "Brain t(w)o Reverse Inference Ontology" package

## How to use this package

### Instalation
Clone this repository and install the application with `python setup.py install`


### Running the analysis
The only necessary parameters are the path to a nifti file with the image to be analyzed and the name of the output file. The program expects a 4D image with the region's labels (or probabilities in the probabilistic version) in the voxels we want to evaluate and zero in all the other voxels.

E.g: `b2rio --brain_path nifti_image.nii.gz --output_file results`.

Apart from this, there are a number of other parameters that can be configured:

- `--n_folds`: Number of folds use to sample studies (default=150)
- `--resample`: Resample of the provided image (default=1, no resample)
- `--frac_sample`: Percentage of studies to use during sample (default=0.7)
- `--radius`: When matching voxels between the provided image and neurosynth activations, this radius is used to limit the maximum distance between co-activated voxels. (default=4)
- `--tfIdf`: Threshold to limit the number of studies based on the tfIdf provided by Neurosynth.


## Examples

### Analizing the Area STS1 (STS) of the Julich's atlas.

Suppose we want to obtain the cognitive processes associated with the STS1 region of the left hemisphere of the Julich atlas.
To do this, we must first obtain the corresponding image. This can be done in a simple way with the following command:

`curl -O https://object.cscs.ch/v1/AUTH_227176556f3c4bb38df9feea4b91200c/hbp-d000001_jubrain-cytoatlas-Area-STS1_pub/5.3/Area-STS1_l_N10_nlin2ICBM152asym2009c_5.3_publicDOI_a8128c417aeeb34c81dabc225670faa5.nii.gz`

Then, it only remains to run B2RIO with this file as a parameter and defining a name for the output file.

`b2rio --brain_path Area-STS1_l_N10_nlin2ICBM152asym2009c_5.3_publicDOI_a8128c417aeeb34c81dabc225670faa5.nii.gz --output_file STS1`

When the process is complete, it will return a CSV file named `STS1.csv` with our results. The following are the first ten entries of the obtained file (out of a total of 18600, remember that the number of folds used can be modified with the `--n_folds` parameter)

|    | term   |   fold |       bf | topConcept   |
|---:|:-------|-------:|---------:|:-------------|
|  0 | action |      0 | 0.937732 | Action       |
|  1 | action |      1 | 0.83446  | Action       |
|  2 | action |      2 | 0.870331 | Action       |
|  3 | action |      3 | 0.834874 | Action       |
|  4 | action |      4 | 0.852323 | Action       |
|  5 | action |      5 | 0.86945  | Action       |
|  6 | action |      6 | 0.920112 | Action       |
|  7 | action |      7 | 0.82419  | Action       |
|  8 | action |      8 | 0.852273 | Action       |
|  9 | action |      9 | 0.885683 | Action       |

Moreover, if we organize the results by their Bayes Factor in a decreasing order, we can see that the most relevant results are related to cognitive processes associated with listening/audition (confirming in some way the knowledge available in the literature).

|      | term                   |   fold |      bf | topConcept   |
|-----:|:-----------------------|-------:|--------:|:-------------|
| 9256 | listening              |    106 | 2.5791  |              |
| 9192 | listening              |     42 | 2.56588 |              |
| 9164 | listening              |     14 | 2.53278 |              |
| 9269 | listening              |    119 | 2.52988 |              |
| 9178 | listening              |     28 | 2.5297  |              |
| 8984 | language comprehension |    134 | 2.51597 | Language     |
| 9280 | listening              |    130 | 2.51489 |              |
| 9265 | listening              |    115 | 2.50568 |              |
| 8929 | language comprehension |     79 | 2.50209 | Language     |
| 9167 | listening              |     17 | 2.49724 |              |
