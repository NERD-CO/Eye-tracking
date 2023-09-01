function [behavFILE] = NewOldTxttoMat_v3(excelLOC , DATADIR , subID ,...
    variantI , blockID)
% Convert .txt file from behavioral session into .mat file
% Outputs table to view TTL values and timestamps
%
% John A. Thompson & Marielle L. Darwin
% Created July 7 2021 | Updated April 24 2023
%
% SAVEDIR = 'C:\HOME\Dropbox\Cases\Case1\Behavioral-data\  %save location
% subID = MW9_NO10
% variantI = 1
% blockID = 'learn'
% textFname = 'NEWOLD00038438.txt'
%
% Events description within .txt file
% The events coorespond to the TTL markers for each trial
% For the learning trials, the TTL markers are the following:
% 55 = start of the experiment, 1 = stimulus ON, 2 = stimulus OFF, 3 = Question Screen Onset,
% 20 = Yes, 21 = NO, 6 = End of Delay after Response, 66 = End of Experiment.
% For the recognition trials, the TTL markers are the following:
% 55 = start of the experiment, 1 = stimulus ON, 2 = stimulus OFF, 3 = Question Screen Onset,
% 31:36 = Confidence in response, 66 = End of Experiment;
%
% Key for TTL_ID shorthand in the 'TTLinterpret' subfunction:
% 'NS'=no, sure; 'NLS'=no, less sure; 'NVU'=no, very unsure
% 'YS'=yes, sure; 'YLS'=yes, less sure; 'YVU'=yes, very unsure

cd(excelLOC)
% Load in excel file
allSUBtable = readtable('Eyetracking patient summary sheet.xlsx');

subROW = allSUBtable(matches(allSUBtable.patientID,subID),:);

subDATAdir = [DATADIR , filesep , subID , filesep , 'Behavioral-data' ,...
    filesep , 'Raw'];

cd(subDATAdir)

% Load text block by reading in every line of .txt file
% PROGRAMMTICALLY FIND TEXT FILE

block2use = upper(blockID(1));
var2use = ['var',num2str(variantI)];
varBl_col = [var2use , '_', block2use , '_txt'];

textFname = subROW.(varBl_col){1};

inputFile = fopen(textFname);       % Open file in Matlab
tLine = fgetl(inputFile);         % Read in 1st line of .txt file
allLines = cell(5000,1);                    % Create empty cell array
counter = 1;                      % Create a counter

while ischar(tLine)
    allLines{counter} = tLine;     % Fill cell array with txt lines
    counter = counter+1;           % Update counter
    tLine = fgetl(inputFile);      % Move through lines in txt file
end

allLines = allLines(cellfun(@(x) ~isempty(x), allLines, 'UniformOutput',true));

% Create cellfun to split elements in 'allLines' & output to a cell array
taskElements = cellfun(@(x) strsplit(x,';'),allLines,'UniformOutput',false);
% Cut 1st column of taskElements
% taskElements = taskElements_temp(2:width(taskElements_temp));
% Create cellfun to "loop" through 'taskElements' to extract 2nd element
TTL_ID = cellfun(@(x) x{2},taskElements,'UniformOutput',false);
TTL_IDn = cellfun(@(x) str2double(x),TTL_ID,'UniformOutput',true);
% Create cellfun to "loop" through 'taskElements' to extract 1st element
timestamp = cellfun(@(x) x{1},taskElements,'UniformOutput',false);
timestampn = cellfun(@(x) str2double(x),timestamp,'UniformOutput',true);
% Interpreter to translate TTL identifiers
[TTLdescription] = TTLinterpret(TTL_ID);

% Create output table
taskData = table(TTL_IDn,TTLdescription,timestampn,'VariableNames',{'TTLvalue','TTL_ID','Timestamp'});

% Store in .mat file
outData.taskinformation = taskData;
outData.patientID = subID;
outData.variant = variantI;
outData.block = blockID;
outData.txtFname = textFname;

% Save file
SAVEDIR = [DATADIR, filesep, subID , '\Behavioral-data\Processed' ];
cd(SAVEDIR);
savefilename = [subID , '_var' , num2str(variantI) , '_' , blockID , '.mat'];
save(savefilename,'outData');
behavFILE = savefilename;


end %%%% END OF MAIN FUNCTION


% Create subfunction to make an interpreter for TTL values
function [TTLdescription] = TTLinterpret(TTLnumber)
%Create table for TTL value key
TTL_ID = {'exp.start','stimON','stimOFF','OnsetQ','Yes','No','EndRespDelay','exp.end','NS','NLS','NVU','YVU','YLS','YS'};
TTLvalue = {'55','1','2','3','20','21','6','66','31','32','33','34','35','36'};

keyTable = table(transpose(TTL_ID),transpose(TTLvalue),'VariableNames',{'TTLI','TTLN'});

% Create empty array to preallocate output
TTLdescription = cell(size(TTLnumber));

% Loop through each TTL number in input argument
for i = 1:length(TTLnumber)
    tempN = TTLnumber{i};
    % Find matching ID for TTL value from keyTable
    TTLlabel = keyTable.TTLI(ismember(keyTable.TTLN,tempN));
    if isempty(TTLlabel)
        TTLdescription{i} = 'NO_label';
    else
        TTLdescription{i} = TTLlabel{1}; %1=value is true, 0=value is false
    end
end

end



