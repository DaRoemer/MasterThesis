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
   - From `.nd2` files, VE-cadherin, nuclei, and CD45 channels were extracted using FIJI.

2. **MΦ Adhesion Analysis**  
   - CD45-stained images were analyzed in FIJI to quantify MΦ numbers.
   - **Relevant Files**:
     - `MFadhesion_analysis.ipynb` – Python notebook for processing and analyzing MΦ adhesion.
   - **Key Steps**:
     - Load CD45-stained images.
     - Preprocess the images and count adherend MΦs in FIJI.
     - Analyse and vizualize quantified MΦ numbers.

3. **Morphology Analysis**  
   - VE-cadherin and nuclei channels were segmented in a human-in-the-loop segmentation workflow with Cellpose.
   - The script `morphology_analysis.ipynb` was used to quantify cellular morphology based on the segmentations.
   - **Relevant Files**:
     - `morphology_analysis.ipynb` – Main script for morphology analysis.
     - `Morphology_tools.py` – Contains helper functions for processing and analyzing morphology.
   - **Key Functions**:
     - `process_and_initialize`: Prepares labeled images, outlines, and actin images for analysis.
     - `create_propertie_image`: Generates property images (e.g., area, aspect ratio, tortuosity) for visualization.
     - `save_as_figure`: Saves the generated property images as figures.
   - **Key Steps**:
     - Segment images using Cellpose.
     - Analyze morphology parameters such as area, aspect ratio, and tortuosity.
     - Save results as images and data files.

4. **Nematic Ordering Analysis**  
   - Cell edge segmentations were analyzed for nematic ordering using the script `nematic_ordering.ipynb`.
   - The scientific method AFT (Alignment by Fourier Transform) was applied.
   - The MATLAB script `AFT_batch.m` (part of the AFT software package and adapted by me) was used for visualization of local alignment vectors. The functions from the AFT software are needed to successfully use this function.
   - **Relevant Files**:
     - `nematic_ordering.ipynb` – Main script for nematic ordering analysis.
     - `Analysis_tools.py` – Contains helper functions for processing nematic ordering data.
   - **Key Functions**:
     - `load_and_process_matlab_data`: Loads and processes MATLAB AFT data for analysis.
     - `process_data`: Processes nematic ordering data, including calculating order parameters and neighborhood sizes.
     - `process_python_AFT_data`: Processes Python-based AFT data for additional analysis.
   - **Key Steps**:
     - Analyse neamtic ordering in `nematic_ordering.ipynb`.
     - Calculate nematic order parameters and order decay to visualize distance dependet alignment.
     - Use MATLAB's `AFT_batch.m` script to generate alignment vectors and heatmaps of selected images.

## Requirements

- **Python 3.x**  
  - Required packages: `numpy`, `pandas`, `matplotlib`, `seaborn`, `scikit-image`, `cellpose`, `SciPy`, `Statsmodels`, `Scikit-posthocs`
  - Different environments are provided for specific tasks:
    - `cellpose_env.yml`: Used for segmentation tasks with Cellpose in the **Morphology Analysis** section.
    - `morphology_analysis.yml`: Used for analyzing cell shape parameters in the **Morphology Analysis** section.
    - `analysis_env.yml`: Used for processing and analyzing nematic ordering data in the **Nematic Ordering Analysis** section.
    - `aft_312_env.yml`: Used for calculating order parameter in AFT-related scripts in the **Nematic Ordering Analysis** section.
- **MATLAB**  
  - Required for running the `AFT_batch.m` script in the **Nematic Ordering Analysis** section.
- **FIJI/ImageJ**  
  - Required for preprocessing and segmenting microscopy images in the **Imaging & Data Extraction** and **MΦ Adhesion Analysis** sections.

# Video visualisation
This montage shows the geometry and appearance of a section of one microfluidc channel. The cells were fixed and stained for nuclei (blue), actin (red) and CD45, a macrophage marker (green).

<div style="text-align: center;">
  <img src="images/20250312_Chip19_3D_channel10_1_maxres.gif" alt="Microfluidic Channel Visualization">
</div>

# Citation

If you use this code for your own work, please cite the corresponding thesis and/or the AFT method where appropriate.

# License
The content of this project itself is licensed under the [Creative Commons Attribution 4.0 Unported license](https://creativecommons.org/licenses/by/4.0/deed.en), and the underlying source code used to analyse and display that content is licensed under the [Apache 2.0 license](https://www.apache.org/licenses/LICENSE-2.0).

# Acknowledgment

This work was conducted as part of my Master's thesis. I would like to express my gratitude to everyone who supported and guided me throughout this project. In particular, I would like to thank my group leader Effie Bastounis and my supervisor Marie Münkel for their invaluable support and mentorship.
