# HUVEC Morphology and MΦ Interaction Analysis
This repository contains the custom scripts used for the analysis of endothelial cell (HUVEC) morphology and macrophage (MΦ) adhesion during infection studies, as described in the thesis project.

## Overview
The analysis pipeline processes raw microscopy data to extract quantitative information about:
- MΦ Adhesion (via CD45 staining)
- HUVEC Morphology (via VE-cadherin and nuclei staining)
- Nematic Ordering of HUVECs (using [AFT - Alignment by Fourier Transform](https://github.com/OakesLab/AFT-Alignment_by_Fourier_Transform))
  
Custom Python and MATLAB scripts were developed to perform segmentation, morphology quantification, and nematic ordering analysis.
## Workflow Summary

1. **Imaging & Data Extraction**  
   - HUVECs were cultured, fixed, and fluorescently stained.
   - Imaging was performed using a Nikon SDCM microscope.
   - From `.nd2` files, VE-cadherin, nuclei, and CD45 channels were extracted.

2. **MΦ Adhesion Analysis**  
   - CD45-stained images were analyzed in FIJI to quantify macrophage numbers.

3. **Morphology Analysis**  
   - VE-cadherin and nuclei channels were used in a human-in-the-loop segmentation workflow with Cellpose.
   - The script `morphology_analysis.ipynb` was used to quantify cellular morphology based on the segmentations.

4. **Nematic Ordering Analysis**  
   - Cell edge segmentations were analyzed for nematic ordering using the script `nematic_ordering.ipynb`.
   - The scientific method AFT (Alignment by Fourier Transform) was applied.
   - The MATLAB script `AFT_batch.m` (part of the AFT software package) was used for visualization of local alignment vectors.

## Repository Contents

- `morphology_analysis.ipynb` – Quantitative analysis of cellular morphology from segmented images.
- `nematic_ordering.ipynb` – Analysis of nematic ordering using Fourier transform methods.
- `AFT_batch.m` – MATLAB script used for batch processing alignment vectors with AFT.
- (Optional) Example input and output data if included.

## Requirements

- **Python 3.x**  
  - Required packages: `numpy`, `pandas`, `matplotlib`, `scikit-image`, `cellpose`, `opencv-python`
- **MATLAB**
- **FIJI/ImageJ**  

## Usage

1. Preprocess and segment the microscopy images using FIJI and Cellpose (see provided Segmentation model).
2. Run `morphology_analysis.ipynb` to analyze cell shape parameters.
3. Run `nematic_ordering.ipynb` to calculate nematic order parameters.
4. Use `AFT_batch.m` in MATLAB to visualize alignment vectors.

# Video visualisation
This montage shows the geaomety and apearnace of a section of one microfluidc channel. The cells were fixed and stained for nuclei (blue), actin (red) and CD45, a macrophage marker (green).
![](images/20250312_Chip19_3D_channel10_1_maxres.gif)


# License
The content of this project itself is licensed under the [Creative Commons Attribution 4.0 Unported license](https://creativecommons.org/licenses/by/4.0/deed.en), and the underlying source code used to analyse and display that content is licensed under the [Apache 2.0 license](https://www.apache.org/licenses/LICENSE-2.0).
