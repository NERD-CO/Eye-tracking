%% SET UP

% STEP 1: Change path directories. Choose the NAS case but update the paths
% the first time you use it to make sure they are correct.
userPC = 'JAT_WORK';
switch userPC
    case 'JAT_HOME'
        basePath = 'E:\Dropbox\Publications_Meta\InProgress\Eye-tracking-MD-memory\Code_test';
        excelLOC = [basePath, '\EyeTrack'];
        mainLOC = [excelLOC, '\eyeTrack'];
        saveLOC = [mainLOC, '\eyeDATA'];
        cleanedDataLOC = [saveLOC, '\cleaned_eyeDATA'];
    case 'JAT_WORK'
        basePath = 'D:\Dropbox\Publications_Meta\InProgress\Eye-tracking-MD-memory\Code_test';
        excelLOC = [basePath, '\EyeTrack'];
        mainLOC = [excelLOC, '\eyeTrack'];
        saveLOC = [mainLOC, '\eyeDATA'];
        cleanedDataLOC = [saveLOC, '\cleaned_eyeDATA'];
    case 'MLD'
        basePath = 'C:\Users\darwinm\Documents\MATLAB';
        excelLOC = [basePath, '\EyeTrack'];
        mainLOC = [excelLOC, '\eyeTrack'];
        saveLOC = [mainLOC, '\eyeDATA'];
        cleanedDataLOC = [saveLOC, '\cleaned_eyeDATA'];
    case 'NAS'
        basePath = 'Z:\EyeTrack manuscript\Code_Test';
        excelLOC = [basePath, '\EyeTrack'];
        mainLOC = [excelLOC, '\eyeTrack'];
        saveLOC = [mainLOC, '\eyeDATA'];
        cleanedDataLOC = [saveLOC, '\cleaned_eyeDATA'];
    case 'MLD_test'
        basePath = 'C:\Users\darwinm\Documents\Code_test';
        excelLOC = [basePath, '\EyeTrack'];
        mainLOC = [excelLOC, '\eyeTrack'];
        saveLOC = [mainLOC, '\eyeDATA'];
        cleanedDataLOC = [saveLOC, '\cleaned_eyeDATA'];
end

% STEP 2: Change ptID to be specific to pt 
ptID = 'AMC_PY21NO05';
% ptID = 'AMC_PY22NO09';
% ptID = 'AMC_PY22NO12';
% ptID = 'AMC_PY22NO13';
% ptID = 'AMC_PY22NO16';

%% STEP 3 CONVERT From EDF to MAT
Extract_Eye_EDF(excelLOC , mainLOC, ptID, basePath)

%% STEP 4 Run Initial_EyeAnalysis function
Initial_EyeAnalysis_MLD_v8(excelLOC, mainLOC, ptID, saveLOC);

%% STEP 5 Run eyeTrackProc funciton
% STEP 5: Run eyeTRACKproc.m f(x) 
eyeTRACKproc_v2(cleanedDataLOC, saveLOC, ptID);

%% 
% STEP 6: Figures - line plot of pupil diameter for confidence ratings
fig = pupil_confRatings(cleanedDataLOC, ptID);
  

