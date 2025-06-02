
% The purpose of this script is to reorganize all DCM files 
% into the OutDicom folder, within each of day1,day2,day3 subfolders

% the benefit of doing this step
    % For DCM files from disk, it organizes single DCM files into folders.
    % DCM files in the OutDicom folder are organized with similar
    % structures as in other subjects' folder.
% The problem of doing this 
    % DCM file names remain the same although they are in the organized
    % folders (meaningless if the files are from disk), but it's ok as all
    % files are renamed after transferred to NifTi.
    
    % If for some reason, one series has been run twice in the same
    % session, DCM files are merged in the outcome folder.
% Also notice
    % Some T1 DCM files might have been duplicated, because they have been
    % transferred out twice

fprintf('Sorting Dicom files\n')
dicomdir = fullfile(studydir,subno,'DICOM');
SessFolders = dir(fullfile(dicomdir,'Day*'));
nSess = length(SessFolders); % number of sessions based on subfolders (ideally 3 sessions)

for ss = 1:nSess
    fprintf('Session %d/%d\n',ss,nSess)
    SessName = SessFolders(ss).name;
    Sessdir = fullfile(dicomdir,SessName); 
    filelist = dir(fullfile(Sessdir, '**/*.IMA'));  % find all files under this directory and subdirectory
    filelist = filelist(~[filelist.isdir]); 

    fprintf('%d DCM files found\n',length(filelist))

    outdir = fullfile(studydir,subno,'OutDicom',SessName); % create a new folder to store copied files within each session
    if ~exist(outdir,'dir')
        mkdir(outdir)
    end

series_names = [{''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''} {''}];

% next should find all Dicom files under each Session folder

sctr = 0;

if ~isempty(filelist) % if there is DCM file then continue
    for j = 1:length(filelist) % loop over each dcm file of the folder
        hdr = dicominfo(fullfile(filelist(j).folder,filelist(j).name)); % need image processing toolbox
        s = hdr.SeriesDescription;

        %check if series description matches any already identified
        kctr = 0;
        for k = 1:length(series_names)

            if strcmp(s,series_names{k}) %if there is a match
                copyfile(fullfile(filelist(j).folder,filelist(j).name),fullfile(outdir,s,filelist(j).name));
                kctr = kctr + 1;
            else
            end
        end

        if ~kctr
            sctr = sctr + 1;
            series_names{sctr} = s;
            if ~exist(fullfile(outdir,s),'dir')
                mkdir(fullfile(outdir,s))
            else
            end
            copyfile(fullfile(filelist(j).folder,filelist(j).name),fullfile(outdir,s,filelist(j).name));
        end


    end % end loop of dcm files
end

end  % end of session

% here the logic is: 
% if a certain s has appeared before, then copying current
% dcm file into corresponding folder
% if a certain s has not appeared before, then adding this s into series name
% cell vector and increasing 'sctr' by one, and copying current dcm file
% into corresponding folder
        
                