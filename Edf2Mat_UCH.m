% ---Function 'Edf2Mat_UCH'---
% 
% Convert .edf file to .mat file 
% Uses Edf2Mat Matlab toolbox found at github.com/mzhuang666/edf-converter
%      NOTE: Need to have this toolbox downloaded and on path 
%
% Inputs: 1) file to convert 2) Patient ID, 3)date of recording (MMDDYYYY)
% Example function call: Edf2Mat_UCH('NO20221615110.edf', 'MW9', '01162022')
%
% Output: File saved with patient ID and date of recording
%
% Marielle L. Darwin | May 11 2022

function [] = Edf2Mat_UCH(edfFile, patientID, recordingDate)
% edfFile = 'NO20221615110.edf';
% patientID = 'MW9';
% recordingDate = '01162022';

% Set path structure
paths = [];
uiwait(msgbox('Navigate to and select folder that contains .edf file'))
%e.g. 'C:\Users\darwinm\Documents\Thompson Lab\Microwire\PatientData\MW9\'
paths.basePath = uigetdir;
paths.path_edf = [strcat(paths.basePath,'\',edfFile)];
cd(paths.basePath);

% Construct new file name
filename = [strcat('eyetrack_', patientID, '_', recordingDate, '.mat')];
% Convert chosen .edf file to .mat file and save
edfRAW = Edf2Mat(paths.path_edf);
save(filename,'edfRAW');


end