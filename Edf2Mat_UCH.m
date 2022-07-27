% ---Function 'Edf2Mat_UCH'---
% 
% Convert .edf file to .mat file 
% Uses Edf2Mat Matlab toolbox found at github.com/mzhuang666/edf-converter
%      NOTE: Need to have this toolbox downloaded and on path 
%
% Inputs: 1) file to convert 2) Patient ID, 3) block (learning 'L' or recog 'R'), 
% 4)date of recording (MMDDYYYY)
%
% Output: File saved with patient ID, block type, and date of recording
%
% Example function call: Edf2Mat_UCH('NO20221615110.edf', 'MW9', 'L', '01062022')
%
% Marielle L. Darwin & John A. Thompson | May 11 2022 | Last update: July 27 2022

function [eyeMatfile] = Edf2Mat_UCH(edfFile, patientID, block, recordingDate)
% edfFile = 'NO20221615110.edf';
% patientID = 'MW9';
% block = 'L';
% recordingDate = '01062022';

% Select folder for edf-converter
uiwait(msgbox('Navigate to and select edf-converter folder'))
edfLoc = uigetdir();
addpath(genpath(edfLoc));

% Set path structure
paths = [];
uiwait(msgbox('Navigate to and select folder that contains .edf file'))
%e.g. 'C:\Users\darwinm\Documents\Thompson Lab\Microwire\PatientData\MW9\'
paths.basePath = uigetdir;
paths.path_edf = [strcat(paths.basePath,'\',edfFile)];
cd(paths.basePath);

% Construct new file name
filename = [strcat('eyetrack_', patientID, '_', block, '_', recordingDate, '.mat')];

% Convert chosen .edf file to .mat file and save
edfRAW = Edf2Mat(paths.path_edf);
save(filename,'edfRAW');
eyeMatfile = filename;
end