
% This script was modified based on James' script from server
% to create anatomical masks based on known coordinates with different sizes

% where are these coordinates from?
coords = [...
    4 -16 -12; ... %peak of iPE effect in right midbrain
    -8 -14 -10; ... %peak of iPE effect in left midbrain
    ];

rs = [2 3 4];
tpmhdr = spm_vol_nifti('rTPM_point1_NoCerebellum.nii');
[tpmvol,XYZ] = spm_read_vols(tpmhdr);
sz = size(tpmvol);
[MX, MY, MZ] = ndgrid(1:sz(1),1:sz(2),1:sz(3));
ref_vox = round(sz/2);
radii = sqrt((MX-ref_vox(1)).^2 + (MY-ref_vox(2)).^2 + (MZ-ref_vox(3)).^2);

for r = rs
    
    ridx = find(radii < r) - sub2ind(sz,ref_vox(1),ref_vox(2),ref_vox(3));    
    idx1 = find(ismember(XYZ',coords(1,:),'rows')>0) + ridx;
    idx2 = find(ismember(XYZ',coords(2,:),'rows')>0) + ridx;
    idx = union(idx1,idx2); % conjunction of left and right
    
    newhdr = tpmhdr;
    newvol = zeros(sz);
    newhdr.fname = strcat('Midbrain_iPE_radius',num2str(r),'mm.nii');
    newvol(idx) = 1;
    spm_write_vol(newhdr,newvol);
    
end
