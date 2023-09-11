%% SET UP

% STEP 1: Change path directories. Choose the NAS case but update the paths
% the first time you use it to make sure they are correct.

userPC = 'JAT_WORK';

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
end

% STEP 2: Change ptID to be specific to pt
% ptID = 'AMC_PY21NO05';
ptID = 'AMC_PY21NO08';
% ptID = 'AMC_PY22NO09';
% ptID = 'AMC_PY22NO12';
% ptID = 'AMC_PY22NO13';
% ptID = 'AMC_PY22NO16';
% ptID = 'AMC_PY23NO23';

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

%% STEP 8 - to extract for Variant tables

% Variant 1
% column 37 trial number learn
% column 41 category ID
% column 42 picture number

% column 43 recog trial number
% column 47 category ID
% column 48 picture number
% column 56 ground truth recog

learn.v1 = table;
learn.v2 = table;
learn.v3 = table;

recog.v1 = table;
recog.v2 = table;
recog.v3 = table;

learn.v1.trialNumer = table2array(variantS.var1.dataTable(:,37));
learn.v1.CatID = table2cell(variantS.var1.dataTable(:,41));
learn.v1.PicNum = table2array(variantS.var1.dataTable(:,42));

learn.v2.trialNumer = table2array(variantS.var2.dataTable(:,37));
learn.v2.CatID = table2cell(variantS.var2.dataTable(:,41));
learn.v2.PicNum = table2array(variantS.var2.dataTable(:,42));

learn.v3.trialNumer = table2array(variantS.var3.dataTable(:,37));
learn.v3.CatID = table2cell(variantS.var3.dataTable(:,41));
learn.v3.PicNum = table2array(variantS.var3.dataTable(:,42));



recog.v1.trialNumer = table2array(variantS.var1.dataTable(:,43));
recog.v1.CatID = table2cell(variantS.var1.dataTable(:,47));
recog.v1.PicNum = table2array(variantS.var1.dataTable(:,48));
recog.v1.OLDNewTruth = table2array(variantS.var1.dataTable(:,56));

recog.v2.trialNumer = table2array(variantS.var2.dataTable(:,43));
recog.v2.CatID = table2cell(variantS.var2.dataTable(:,47));
recog.v2.PicNum = table2array(variantS.var2.dataTable(:,48));
recog.v2.OLDNewTruth = table2array(variantS.var2.dataTable(:,56));

recog.v3.trialNumer = table2array(variantS.var3.dataTable(:,43));
recog.v3.CatID = table2cell(variantS.var3.dataTable(:,47));
recog.v3.PicNum = table2array(variantS.var3.dataTable(:,48));
recog.v3.OLDNewTruth = table2array(variantS.var3.dataTable(:,56));











