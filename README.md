# Overview

Given a set of voxel location (sample file provided), experimental design matrix, and a desired signal-to-noise ratio, this program generates a set of synthesized voxel activations. 
Each experimental covariate (which can be activated to varying degrees during each trial, as specified in the design matrix) is assigned a single radial basis function. The basis function's center is chosen uniformly within the confines of the specified brain, and the width is also chosen uniformly. (Note: you can associate multiple radial basis functions with each covariate by giving each covariate multiple columns in the design matrix.)

The synthesized brain image for each trial is a weighted combination of the basis functions for the sources active during that trial, plus zero-mean Gaussian noise. (The standard deviation of the Gaussian noise is controlled via the SNR parameter.)

For added realism, an optional flag allows the user to specify whether the synthetic data incorporate a synthetic hemodynamic response function.

The synthetic data can be saved in NIFTI format (this relies on code, included in the download, from [here](http://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)).

Note: this code includes the "join" function available [here](http://www.mathworks.com/matlabcentral/fileexchange/4872-join).

# Generating data

To get started, try loading in the `meta.mat` file:

`load meta.mat;`

This contains a set of voxel locations that look like a "real" brain.

Next, define a set of activation values that you want the images to exhibit.  The activation matrix should have number-of-images rows and number-of-ROIs columns.  Each "ROI" (region of interest) will be a spherical radial basis function (RBF) placed somewhere (at random) in the brain.  You can then define how active each ROI is in each image.  A useful starting point to set your activation matrix is to generate a `toeplitz` matrix:

```
K = 25;
w = toeplitz(1:K);
```
This will vary the activity of each ROI smoothly over the course of each of 25 images, and each ROI will exhibit a different timecourse.  For fancier (and possibly, more useful) setups, you could set `w` to follow the design matrix of your experiment, with number-of-timepoints rows, and number-of-regressors columns.  E.g. this allows you to generate data where each ROI responds to some feature in your experiment, which could be useful for testing decoders.

Having set up the basics, you can now generate synthetic fMRI data using:
```
generate_data(meta, w, 100, 1, 0, 'synthetic_data.nii')
```

(See documentation in `generate_data` for more details on how to set the different options.)  This will create `synthetic_data.nii`, a Nifti file that can be processed and/or viewed using standard fMRI tools.
