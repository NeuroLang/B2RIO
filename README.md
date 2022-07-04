# B2RIO: "Brain to Reverse Inference Ontology" package

## How to use this package
### Computing data

The program must be run in two phases. The first (Solver), analyzes the data provided and generates the results. It can be run in parallel. The second (Analyzer), unifies the different values obtained and calculates the bayes factor for each term/region.

To be able to run the Solver, we need 3 parameters.
1) A list of the ids of the regions to be analyzed.
2) A list of the regions that make up the atlas.
3) The .nii.gz file of the atlas to be analyzed.

The list of regions provided must comply with the following specification:
a) First, a header with the names of the columns.
b) Then, a list of regions:

`Region number` '`Region name` `Hemisphere`'.

An example of this would be:
```
mapindex regionname
1 'Frontal-to-Temporal-I (GapMap) left'
2 'Frontal-to-Temporal-I (GapMap) right'
3 'Ch 123 (Basal Forebrain) left'
4 'Ch 123 (Basal Forebrain) right'
5 'Ch 4 (Basal Forebrain) left'
6 'Ch 4 (Basal Forebrain) right'
```

### How to run the solver in parallel
The solver requires a list of regions to be analyzed. Each of these regions will generate its own file with results, so that each region can be analyzed in parallel.

Several solvers with subsets of regions can be launched in a cluster at the same time (we recommend launching one region per instance).

### Analyzing data
After all the solvers are finished, it is only necessary to run the analyzer to generate the last file with the results.

### Simple example
Simple example of how to run this package to obtain results using Julich's atlas:

``` python
from b2rio.analysis import Analyzer
from b2rio.processing import Solver

s = Solver() #default parameters: n_folds=150, resample=1, random_state=42, frac_sample=0.7

regions_path = '../reverse_inference_histology/JULICH_BRAIN_CYTOARCHITECTONIC_MAPS_2_9_MNI152_2009C_NONL_ASYM.txt'
brain_path = '../reverse_inference_histology/JULICH_BRAIN_CYTOARCHITECTONIC_MAPS_2_9_MNI152_2009C_NONL_ASYM.pmaps.nii.gz'

# All the regions need to be compute to be able to run the analyzer
s.run([1], regions_path, brain_path, probabilistic_atlas=True, results_path='./')
# ...
s.run([296], regions_path, brain_path, probabilistic_atlas=True, results_path='./')

a = Analyzer()

