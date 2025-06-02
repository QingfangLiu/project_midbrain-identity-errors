# How to Make the Network-Targeted TMS Illustration Figure

**Date:** July 5, 2023

## Step 1: Coordinate Transformation

Transform the LPFC stimulation coordinates from the subjects’ native space to MNI space using **MATLAB**.

- Output file: `myNode.txt`
- Remove `""` in the label column (6th column) by replacing them with empty strings.

## Software Used

- **BrainNet Viewer** (Version: `BrainNetViewer_20191031`)  
  Easy to use within MATLAB.

## Surface File

- Use `BrainMesh_ICBM152.nv`
- Note: The version with `_tal` appears visually identical.

## Display Options in BrainNet Viewer

### Surface Settings
- **Opacity:** Set to `0.5`

### Node Settings
- **Size:** Use `Value – raw` (easier to control than the auto option)
- **Color Mapping:**
  - **Modular 1:** Blue (or cyan as shown in the paper)
  - **Modular 2:** Red

## Additional Notes

- Used the "Add ROIs" function to preview actual ROI sizes.
- For illustration purposes, only the **Node** function was used to make the figure look cleaner.
- **Node Sizes:**
  - ROI nodes: `7` (visually matched to actual ROI size)
  - Individual stimulation coordinates: `1`

## Caveats

- Some LPFC stimulation coordinates may fall just outside the LPFC ROI boundaries when overlaid with actual ROIs.
