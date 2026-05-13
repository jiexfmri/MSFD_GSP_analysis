# Multiscale structure-function decoupling (MSFD)

This repository contains code created to conduct the main analyses described in the paper "**Multiscale structural connectome eigenmodes constrain human brain functional dynamics**"

Using graph signal processing (GSP) to study how the multiscale structural connectome (MSC) constrains functional dynamics acquired from functional MRI (fMRI) data.

## System Requirements

Download the repository, and you're good to go. Read the comments and documentation within each code for usage guidance.

*MATLAB* (versions of MATLAB from R2018b to R2024b: The MathWorks, Inc.) was used to write the code.  

Some of the MATLAB-based scripts depend on packages developed by others. Copies of these packages have been stored in the `Matlab_function/` folder to ensure version compatibility.

## File descriptions

1. `Data/`: Folder containing example data. We provide the demo dataset (5 subjects) for the calculations described in the paper.  Original human fMRI data from the HCP [Human Connectome Project](https://db.humanconnectome.org/). Please consult the link for detailed information on access, licensing, and usage terms and conditions.
2. `Matlab_function/`: Folder containing utility MATLAB functions for data analysis.
3. `Step01_Fusion_Network.m`: MATLAB script for constructing the multiscale structural connectome (MSC).
4. `Step02_Comparison_acc_rsfMRI.m`: MATLAB script for reconstructing spontaneous resting-state functional activity.
5. `Step03_Comparison_acc_task.m`: MATLAB script for reconstructing task-evoked functional activity.
6. `Step04_MSFD_Full_Pipline.m`: MATLAB script for quantifying multiscale structure–function decoupling (MSFD).

## Installation

Download the repository, and you're good to go. Read the comments and documentation within each code for usage guidance.

`Step04_MSFD_Full_Pipline.m`: MATLAB script to demonstrate how to use multiscale structural connectome eigenmodes to analyze fMRI data and compute multiscale structure–function decoupling (MSFD) index.

## Further details

Feel free to get in touch if you have any questions ([jiexiafmri@gmail.com](jiexiafmri@gmail.com))

