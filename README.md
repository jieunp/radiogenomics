# Radiogenomics

Online Supplement for 
Prediction of Core Signaling Pathway using Diffusion- and Perfusion-based MRI Radiomics and Next Generation Sequencing in IDH wild type Glioblastoma

## Feature Extraction (Radiomics-Code-AMC-Anew)

1. File for feature extraction: Run_test.m

2. Data preparation

		1)  All imaging data should be in Nifty type
		2)  Two imaging data are necessary: segmented mask (our example: reg_mni_ROI.nii) and base image (reg_mni_n4_ss_norm_reg_rs_T2.nii)
		3)  Above two imaging data needs to be co-registered before feature extraction
		4)  For CE-T1w/T1w/T2w/FLAIR image, the base image needs to be white-striped before feature extraction
		5)  Data structure: needs to be given like Figure 1					
		6)  In the root path, only patients’ folders are allowed. Any other file or folder will cause an error

3. Code preparation

		1)  Put all the codes within the same folder (name: Radiomics-Code-AMC-Anew)
		2)  You need to download matlab tools: imMinkowski
		3)  You need to download matlab tools: Tools for NIfTI and ANALYZE image (version 1.27)
		
4. How to use

		1)  Unzip the zipped folder (Radiomics-Code-AMC-Anew)
		2)  Open the Matlab software (Run_test.m)
		3)  Addpath to the subfolder (Radiomics-Code-AMC-Anew /imMinkowski)
		4)  Addpath to the subfolder (Radiomics-Code-AMC-Anew /NIfTI_20140122)
		5)  Scale the parameters in the code: modifiable number of bins in histogram analysis and number of adjacent voxels in the texture analysis.
		6)  Set the pixel width and slice thickness
		7)  Set the rootpath (above 2-4) and pathname
		8)  Put the ROI mask and base image name
		9)  Run 
		10) When it properly runs, subject 1, subject 2 .. appears in the background


## Analysis Pipeline (Analysis-Pipeline-Core Signaling Pathway.R)

1.	The R codes are embedded

2.	The codes have 4 parts 

		1)  Line 5-129: Feature selection via Student's t-test with false discovery rate correction
		2)  Line 134-233: Feature selection via LASSO penalization and calculate AUC for each genetic mutation
		3)  Line 234-329: Feature selection via Random Forest and find top 5 important features
		4)  Line 330- : Calculate diagnostic performance
		
3.	The data needs to prepared as csv format

4.	Feature data: columns – features, rows- patients 

5.	Reference data: Name of columns: genetic mutation coded as 1 (positive) or 0 (negative) and name of rows: patients (Figure 2)

### Figure 1

![gitfig](https://user-images.githubusercontent.com/26832081/63306404-56dc5f80-c325-11e9-942b-27d71ac99de2.png)

### Figure 2

![gitfig2](https://user-images.githubusercontent.com/26832081/63306405-5774f600-c325-11e9-8e50-3d0cdc370ff7.jpg)
