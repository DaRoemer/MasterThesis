# HUVEC Morphology and MΦ Interaction Analysis

This repository contains custom scripts for analyzing endothelial cell (HUVEC) morphology and macrophage (MΦ) adhesion during infection studies, as described in the thesis project.

## Table of Contents
1. [Project Overview](#project-overview)
2. [Workflow Summary](#workflow-summary)
3. [Requirements](#requirements)
4. [Video Visualisation](#video-visualisation)
5. [Citation](#citation)
6. [License](#license)
7. [Acknowledgment](#acknowledgment)

## Project Overview
The analysis pipeline processes raw microscopy data to extract quantitative information about:
- **MΦ Adhesion** (via CD45 staining)
- **HUVEC Morphology** (via VE-cadherin and nuclei staining)
- **Nematic Ordering of HUVECs** (using [AFT - Alignment by Fourier Transform](https://github.com/OakesLab/AFT-Alignment_by_Fourier_Transform))

Custom Python and MATLAB scripts were developed for segmentation, morphology quantification, and nematic ordering analysis.

## Workflow Summary

### 1. Imaging & Data Extraction
- HUVECs were cultured, fixed, and fluorescently stained.
- Imaging was performed using a Nikon SDCM microscope.
- `.nd2` files were processed in FIJI to extract VE-cadherin, nuclei, and CD45 channels.

### 2. MΦ Adhesion Analysis
- **Script**: `MFadhesion_analysis.ipynb`
- **Steps**:
  1. Load CD45-stained images.
  2. Preprocess and count adherent MΦs.
  3. Visualize and analyze results.

### 3. Morphology Analysis
- **Script**: `morphology_analysis.ipynb`
- **Steps**:
  1. Segment images using Cellpose.
  2. Analyze morphology parameters (e.g., area, aspect ratio, tortuosity).
  3. Save results as images and data files.

### 4. Nematic Ordering Analysis
- **Script**: `nematic_ordering.ipynb`
- **Steps**:
  1. Analyze nematic ordering using AFT.
  2. Calculate order parameters and visualize alignment.
  3. Use MATLAB's `AFT_batch.m` for alignment vectors and heatmaps.

## Requirements

| Software/Tool   | Description                                                                 |
|------------------|-----------------------------------------------------------------------------|
| **Python 3.x**  | Required packages: `numpy`, `pandas`, `matplotlib`, `seaborn`, `scikit-image`, `cellpose`, `SciPy`, `Statsmodels`, `Scikit-posthocs`. |
| **MATLAB**      | Needed for running `AFT_batch.m` in the Nematic Ordering Analysis section.  |
| **FIJI/ImageJ** | Used for preprocessing and segmenting microscopy images.                   |
| **AFT Tool**    | Required for nematic ordering analysis using the AFT method.               |

### Conda Environments
| Environment File         | Purpose                                      |
|---------------------------|----------------------------------------------|
| `cellpose_env.yml`        | Segmentation tasks with Cellpose.            |
| `morphology_analysis.yml` | Analyzing cell shape parameters.             |
| `analysis_env.yml`        | Processing nematic ordering data.            |
| `aft_312_env.yml`         | Calculating order parameters in AFT scripts. |

# Video visualisation
This video illustrates the geometry and appearance of a microfluidic channel section. Cells were fixed and stained for:
- **Nuclei** (blue)
- **Actin** (red)
- **CD45** (green, macrophage marker)

<div style="text-align: center;">
  <img src="images/20250312_Chip19_3D_channel10_1_maxres.gif" alt="Microfluidic Channel Visualization">
</div>

# Citation

If you use this code for your own work, please cite the corresponding thesis and/or the AFT method where appropriate.

# License
The content of this project itself is licensed under the [Creative Commons Attribution 4.0 Unported license](https://creativecommons.org/licenses/by/4.0/deed.en), and the underlying source code used to analyse and display that content is licensed under the [Apache 2.0 license](https://www.apache.org/licenses/LICENSE-2.0).

# Acknowledgment

This work was conducted as part of my Master's thesis. I would like to express my gratitude to everyone who supported and guided me throughout this project. In particular, I would like to thank my group leader Effie Bastounis and my supervisor Marie Münkel for their invaluable support and mentorship.
