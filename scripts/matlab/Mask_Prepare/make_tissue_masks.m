
clear 
hdr = spm_vol('rTPM.nii');

nhdr = hdr(1); 

% gray matter
x = spm_read_vols(hdr(1)); 
y = zeros(size(x));
y(x > 0.1) = 1; % thereshold at 10%
nhdr.fname = 'gm_0.1.nii';
spm_write_vol(nhdr,y);

% white matter
x = spm_read_vols(hdr(2)); 
y = zeros(size(x));
y(x > 0.9) = 1; % thereshold at 90%
nhdr.fname = 'wm_0.9.nii';
spm_write_vol(nhdr,y);

% CSF
x = spm_read_vols(hdr(3)); 
y = zeros(size(x));
y(x > 0.9) = 1; % thereshold at 90%
nhdr.fname = 'csf_0.9.nii';
spm_write_vol(nhdr,y);