
Notes on how functional imaging data were organized.

For univariate analysis done in MNI space:
- Starting with 'f' images
- Realignment: not creating new images, saving info in the header
- Coregistration: not creating new images, saving info in the header
- Normalize: adding 'w' to 'f' images
- Spatial smoothing: adding 's6' to 'wf' images
- Spatial smoothing: adding 's2mm' to 'wf' images
- Overall, the output images: 'wf', 's6wf', 's2mmwf'


For multivariate analysis done in native space:
- Starting with 'f' images
- Realignment & reslice: creating new images 'rf' files
- coregistration: info saved in 'rf' header files
- Spatial smoothing: adding 's2' to 'rf' images
- Overall, the output images: 'rf', 's2rf'


For global connectedness analysis:
- Starting with 'f' images (with realignment & coregistration in the header)
- Rewrite with 3mm voxel size: 'w3f' files
- Spatial smoothing with 6mm: adding 's6mm'
- Overall, the output images: 'w3f', 's6mmw3f'




