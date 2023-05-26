%% SET UP

% STEP 1: Change path directories. Choose the NAS case but update the paths
% the first time you use it to make sure they are correct.
userPC = 'JAT_HOME';
switch userPC
    case 'JAT_HOME'
        codeLocation = 'C:\Users\Admin\Documents\Github\Eye-tracking\NewOldTask-analysis';
        edf2matLOC = 'C:\Users\Admin\Documents\Github\Eye-tracking\edf-converter-master';
        edfCheck = which('Edf2Mat.m');
        if isempty(edfCheck)
            addpath(genpath(edf2matLOC));
        end
        addpath(codeLocation)
        excelLocation = 'C:\Users\Admin\Documents\Github\Eye-tracking\NewOldTask-analysis';
        dataLocation = 'E:\Dropbox\Publications_Meta\InProgress\Eye-tracking-MD-memory\Code_test\EyeTrack\eyeTrack';
        savePreProcLocation = [dataLocation , filesep , 'eyeDATA'];
        saveCleanLocation = [savePreProcLocation , filesep , 'cleaned_eyeDATA'];
    case 'JAT_WORK'
        codeLocation = 'C:\Users\Admin\Documents\Github\Eye-tracking\NewOldTask-analysis';
        edf2matLOC = 'C:\Users\Admin\Documents\Github\Eye-tracking\edf-converter-master';
        edfCheck = which('Edf2Mat.m');
        if isempty(edfCheck)
            addpath(genpath(edf2matLOC));
        end
        addpath(codeLocation)
        excelLocation = 'C:\Users\Admin\Documents\Github\Eye-tracking\NewOldTask-analysis';
        dataLocation = 'D:\Dropbox\Publications_Meta\InProgress\Eye-tracking-MD-memory\Code_test\EyeTrack\eyeTrack';
        savePreProcLocation = [dataLocation , filesep , 'eyeDATA'];
        saveCleanLocation = [savePreProcLocation , filesep , 'cleaned_eyeDATA'];
    case 'MLD'
        % basePath = 'C:\Users\Admin\Documents\Github\Eye-tracking\NewOldTask-analysis';
        % excelLOC = basePath;
        % mainLOC = [excelLOC, '\eyeTrack'];
        % saveLOC = [mainLOC, '\eyeDATA'];
        % cleanedDataLOC = [saveLOC, '\cleaned_eyeDATA'];
    case 'NAS'
        % basePath = 'C:\Users\Admin\Documents\Github\Eye-tracking\NewOldTask-analysis';
        % excelLOC = basePath;
        % mainLOC = [excelLOC, '\eyeTrack'];
        % saveLOC = [mainLOC, '\eyeDATA'];
        % cleanedDataLOC = [saveLOC, '\cleaned_eyeDATA'];
    case 'MLD_test'
        % basePath = 'C:\Users\Admin\Documents\Github\Eye-tracking\NewOldTask-analysis';
        % excelLOC = basePath;
        % mainLOC = [excelLOC, '\eyeTrack'];
        % saveLOC = [mainLOC, '\eyeDATA'];
        % cleanedDataLOC = [saveLOC, '\cleaned_eyeDATA'];
end

% STEP 2: Change ptID to be specific to pt 
% ptID = 'AMC_PY21NO05';
% ptID = 'AMC_PY22NO09';
ptID = 'AMC_PY22NO12';
% ptID = 'AMC_PY22NO13';
% ptID = 'AMC_PY22NO16';

%% STEP 3 CONVERT From EDF to MAT
Extract_Eye_EDF(excelLocation , dataLocation, ptID, edf2matLOC)
clc
%% STEP 4 Run Initial_EyeAnalysis function
EyeAnalysis_DataExtract_v2(excelLocation, dataLocation, ptID, savePreProcLocation);
clc
%% STEP 5 Run eyeTrackProc funciton
% STEP 5: Run eyeTRACKproc.m f(x) 
eyeTRACKproc_v5(saveCleanLocation, savePreProcLocation, ptID);
clc
%% STEP 5a - plot quality check

eyeQUALITY_PS(saveCleanLocation, ptID);

% CHOOSE ONE EYE PER CONDITION - ADD to EXCEL FILE

%% STEP 6 - Create stats matrices




%% STEP 7a - plot data - Learning block - single subject
eyeTrackPlots_v2(saveCleanLocation, excelLocation, ptID , 1)

% STEP 6: Figures - line plot of pupil diameter for confidence ratings
% fig = pupil_confRatings(cleanedDataLOC, ptID);

%% STEP 7b - plot data - Learning block - single subject
eyeTrackPlots_v2(saveCleanLocation, excelLocation, ptID , 2)


