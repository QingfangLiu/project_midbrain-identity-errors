
%% load TMS coordinates in subjects' native space

clc; clear;
xlsfile = '/Users/qingfangliu/Dropbox/KahntLab/Project_SPEEDTMS/Analysis/SubjectConds.xlsx';
table = readtable(xlsfile);
table = table(table.Excluded==0,:); % exclude some subjects
xyz_left = table2array(table(:,26:28)); % 31*3 double
xyz_right = table2array(table(:,29:31)); % 31*3 double

subs = table2array(table(:,1));
nsubs = length(subs);

%% convert coordinates from subjects' native space to MNI space

for i = 1:nsubs
    
    % find the .mat file and load it
    subno = sprintf('Sub%d',subs(i));
    mat_path = fullfile('/Users/qingfangliu/Experiment',subno,'TMS/T1','*.mat');
    mat_file = dir(mat_path);
    mat_file_name = fullfile(mat_file.folder,mat_file.name);
    load(mat_file_name)

    % left and right coordinates of current subject
    sub_left = xyz_left(i,:); 
    sub_right = xyz_right(i,:);

    % convert left and right coordinates to MNI space 
    % 'Affine' comes from loading the mat file
    MNI_left(i,:) = Affine(1:3,:) * [sub_left 1]';
    MNI_right(i,:) = Affine(1:3,:) * [sub_right 1]';

end

%% Organize coordinates into a table and save as txt file

% create labels (can be shown on the plots)
labels = cellstr(num2str(subs));
labels_L = strcat(labels, 'L');
labels_R = strcat(labels, 'R');

% combine left and right & add other columns
dat = [MNI_left;MNI_right];
dat = round(dat * 100) / 100;
dat(:,4) = 2;  % code for color
dat(:,5) = 1;  % code for size
dat = array2table(dat);
dat(:,6) = [labels_L;labels_R];

% rename variable name and save
dat.Properties.VariableNames = ["x","y","z","color","size","label"];
writetable(dat,'myNode.txt','Delimiter',' ')  


