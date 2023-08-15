%% SET UP

% STEP 1: Change path directories. Choose the NAS case but update the paths
% the first time you use it to make sure they are correct.

userPC = 'MLD';

% userPC = 'JAT_WORK';

switch userPC
    case 'JAT_HOME'
        codeLocation = 'C:\Users\Admin\Documents\Github\Eye-tracking\NewOldTask-analysis';
        edf2matLOC = 'C:\Users\Admin\Documents\Github\Eye-tracking\edf-converter-master';
        edfCheck = which('Edf2Mat.m');
        boundLOC = 'C:\Users\Admin\Documents\Github\Eye-tracking\boundedLINE';
        cbrewLOC = 'C:\Users\Admin\Documents\Github\Eye-tracking\cbrewerALL';
        if isempty(edfCheck)
            addpath(genpath(edf2matLOC));
        end

        addpath(genpath(boundLOC))
        addpath(genpath(cbrewLOC))

        addpath(codeLocation)
        excelLocation = 'C:\Users\Admin\Documents\Github\Eye-tracking\NewOldTask-analysis';
        dataLocation = 'E:\Dropbox\Publications_Meta\InProgress\Eye-tracking-MD-memory\Code_test\EyeTrack\eyeTrack';
        savePreProcLocation = [dataLocation , filesep , 'eyeDATA'];
        saveCleanLocation = [savePreProcLocation , filesep , 'cleaned_eyeDATA'];
    case 'JAT_WORK'
        codeLocation = 'C:\Users\Admin\Documents\Github\Eye-tracking\NewOldTask-analysis';
        edf2matLOC = 'C:\Users\Admin\Documents\Github\Eye-tracking\edf-converter-master';
        edfCheck = which('Edf2Mat.m');
        boundLOC = 'C:\Users\Admin\Documents\Github\Eye-tracking\boundedLINE';
        cbrewLOC = 'C:\Users\Admin\Documents\Github\Eye-tracking\cbrewerALL';
        if isempty(edfCheck)
            addpath(genpath(edf2matLOC));
        end

        addpath(genpath(boundLOC))
        addpath(genpath(cbrewLOC))

        addpath(codeLocation)
        excelLocation = 'C:\Users\Admin\Documents\Github\Eye-tracking\NewOldTask-analysis';
        dataLocation = 'D:\Dropbox\Publications_Meta\InProgress\Eye-tracking-MD-memory\Code_test\EyeTrack\eyeTrack';
        savePreProcLocation = [dataLocation , filesep , 'eyeDATA'];
        saveCleanLocation = [savePreProcLocation , filesep , 'cleaned_eyeDATA'];
    case 'MLD'
        codeLocation = 'C:\Users\darwinm\Documents\MATLAB\EyeTrack';
        edf2matLOC = 'C:\Users\darwinm\Documents\Github\edf-converter-master';
        edfCheck = which('Edf2Mat.m');
        boundLOC = 'C:\Users\darwinm\Documents\MATLAB\EyeTrack\boundedLINE';
        cbrewLOC = 'C:\Users\darwinm\Documents\MATLAB\EyeTrack\cbrewerALL';
        if isempty(edfCheck)
            addpath(genpath(edf2matLOC));
        end

        addpath(genpath(boundLOC))
        addpath(genpath(cbrewLOC))

        addpath(codeLocation)
        excelLocation = 'C:\Users\darwinm\Documents\MATLAB\EyeTrack';
        dataLocation = 'C:\Users\darwinm\Documents\MATLAB\EyeTrack\eyeTrack';
        savePreProcLocation = [dataLocation , filesep , 'eyeDATA'];
        saveCleanLocation = [savePreProcLocation , filesep , 'cleaned_eyeDATA'];
        
    case 'NAS'
%         basePath = 'C:\Users\darwinm\Documents\MATLAB\EyeTrack';
%         excelLOC = basePath;
%         mainLOC = [excelLOC, '\eyeTrack'];
%         saveLOC = [mainLOC, '\eyeDATA'];
%         cleanedDataLOC = [saveLOC, '\cleaned_eyeDATA'];
    case 'MLD_test'
        % basePath = 'C:\Users\Admin\Documents\Github\Eye-tracking\NewOldTask-analysis';
        % excelLOC = basePath;
        % mainLOC = [excelLOC, '\eyeTrack'];
        % saveLOC = [mainLOC, '\eyeDATA'];
        % cleanedDataLOC = [saveLOC, '\cleaned_eyeDATA'];
end

% STEP 2: Change ptID to be specific to pt 
 ptID = 'AMC_PY21NO05';
% ptID = 'AMC_PY22NO09';
% ptID = 'AMC_PY22NO12';
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
eyeTRACKproc_PupilSize(saveCleanLocation, savePreProcLocation, ptID);
%%
%eyeTRACKproc_PupilLocation(saveCleanLocation, savePreProcLocation, ptID);
% ADD PUPIL LOCATION EXTRACT
% ADD GAZE EXTRACT
% ADD SACCADE EXTRACT


%clc
%% STEP 5a - plot quality check

eyeQUALITY_PS(saveCleanLocation, ptID);

% CHOOSE ONE EYE PER CONDITION - ADD to EXCEL FILE

%% STEP 6a - Create stats matrices Learn

% INPUTS: (cleanedDataLOC, exclLOC, ptID ,conditiON, subStep)
statSummaryL = eyeTrack_StatPrep_v1(saveCleanLocation, excelLocation, ptID , 'learn' , 1);

%% STEP 6b - Create stats matrices recog 1

% INPUTS: (cleanedDataLOC, exclLOC, ptID ,conditiON, subStep)
statSummaryR1 = eyeTrack_StatPrep_v1(saveCleanLocation, excelLocation, ptID , 'recog' , 1);
statSummaryR2 = eyeTrack_StatPrep_v1(saveCleanLocation, excelLocation, ptID , 'recog' , 2);
statSummaryR3 = eyeTrack_StatPrep_v1(saveCleanLocation, excelLocation, ptID , 'recog' , 3);
statSummaryR4 = eyeTrack_StatPrep_v1(saveCleanLocation, excelLocation, ptID , 'recog' , 4);

%% STEP 6c
[groupStatsL_1] = eyeTrack_StatPrep_allSubjects_v3(saveCleanLocation, 'learn', 1); % Only 1 condition
[groupStatsR_1] = eyeTrack_StatPrep_allSubjects_v3(saveCleanLocation, 'recog', 1);
[groupStatsR_2] = eyeTrack_StatPrep_allSubjects_v3(saveCleanLocation, 'recog', 2);
[groupStatsR_3] = eyeTrack_StatPrep_allSubjects_v3(saveCleanLocation, 'recog', 3);
[groupStatsR_4] = eyeTrack_StatPrep_allSubjects_v3(saveCleanLocation, 'recog', 4);

cd(dataLocation);
save('finalStat_GroupData.mat', 'groupStatsL_1','groupStatsR_1', 'groupStatsR_2',...
    'groupStatsR_3','groupStatsR_4') %continue with all



%copy and paste 120 5x change output "groupstats_Learning",
%"groupStats_Recog_case1" to get unique mat file for each test/run throguh

%% STEP 7a - plot data - Learning block - single subject
eyeTrackPlots_v3(saveCleanLocation, excelLocation, ptID , 1)

% STEP 6: Figures - line plot of pupil diameter for confidence ratings
% fig = pupil_confRatings(cleanedDataLOC, ptID);

%% STEP 7b - plot data - Learning block - single subject
eyeTrackPlots_v2(saveCleanLocation, excelLocation, ptID , 2)


