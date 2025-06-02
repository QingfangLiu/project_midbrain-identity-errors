
# 🧠 Guidance: Finding Stimulation Target Coordinates

This document outlines the process of identifying brain stimulation targets based on resting-state fMRI data using SPM (MATLAB).

---

## 📜 Script 1: Bias Correction

1. **Convert DICOM to NIfTI:**
   - Use `spm_dicom_headers` and `spm_dicom_convert`.
   - Output: Single NIfTI file (`s*.nii`).

2. **Bias Field Correction:**
   - Segment anatomical image.
   - Save the bias-corrected image (with prefix `m`).

---

## 📜 Script 2: Define Stimulation Coordinates

1. **MNI Coordinates of ROIs:**
   - Right OFC: `(28, 38, -16)`
   - Left OFC: `(-28, 38, -16)`
   - Right LPFC: `(48, 38, 20)`
   - Left LPFC: `(-48, 38, 20)`

2. **Create Spherical ROIs:**
   - Centered on above coordinates.
   - Constrained by gray matter probability map.
   - Defined in voxel space.

---

## 🧠 rsfMRI Preprocessing Steps

### 1. Convert rsfMRI DICOMs to NIfTI
- Again use `spm_dicom_convert` and `spm_dicom_headers`.
- 250 volumes → 250 files (`f*.nii`).

### 2. Realignment
- `Estimate & Reslice`
- Motion correction across 250 volumes.
- Prefix: `r`

### 3. Coregistration
- `Estimate`
- Reference: T1 after bias correction (`ms*.nii`)
- Source: Mean of resliced EPI
- Others: All resliced EPIs

### 4. Normalization
- Apply only to anatomical image (prefix `w`)
- rsfMRI images remain in native space.

### 5. Smoothing
- Applied to coregistered EPI images (likely `s6rf*.nii`)
- Use dependency handling or manual specification.

### 6. Deformations 
- Generate inverse deformation field: `y_inverse_deformation.nii`
- Image base: `s6rf*.nii`

---

## 🎯 ROI Transformation

### Unnormalize & Reslice ROIs:
1. **Normalise → Write:**
   - Use inverse deformation field to "un-normalize" 4 ROIs.
   - Output has original voxel size: `79 × 95 × 79`.

2. **Coregister → Reslice:**
   - Reference image: `s6rf*.nii`
   - Output voxel size: `104 × 96 × 58`.

---

## 📈 Functional Connectivity (GLM)

1. **GLM Setup:**
   - Calculate mean BOLD in seed (OFC) ROI.
   - Use it as the regressor in voxel-wise GLMs.
   - Include motion parameters as nuisance regressors.
   - Use unnormalized ROIs to extract voxel signals.

2. **Contrast:**
   - Simple T-test with weight `[1]`.

---

## 🎯 Identify Optimal Target Coordinates

- From `spm_T` map, find voxel in LPFC with **strongest functional connectivity** to OFC seed.
- This is the **final stimulation target coordinate**.

