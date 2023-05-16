function [] = eyeTRACKproc_v2(cleanedDataLOC, saveLOC, ptID)

% Set data paths
if ~exist(cleanedDataLOC,'dir')
    mkdir(cleanedDataLOC)
end

% [For every file in eyeDATA, loop through to run analyses?]
% save cleaned files in new folder to only have raw data in this folder
% need to cd to new folder to save

% CD to mainPath for raw data
cd(saveLOC);

% % Get contents of eyeDATA as a variable
% test = dir('*.mat');
% fileList = {test.name};

% for loop to loop throguh fileList
%for fileS = 1:length(fileList)
eyeData_pt = append('eyeData_', ptID,'.mat');
tempFile_name = eyeData_pt;

% Set tempFile_name
saveFile_name1 = ['cl_', tempFile_name];

% Load in file
load(tempFile_name, 'variantS');

% Set tempEye- loop through every eye in every variant within variantS
% Go into variantS and determine # variants there
varSnum = length(fieldnames(variantS));

% Other option: Extract field names of variantS ahead of time
varSfieldN = fieldnames(variantS);
% Extract field names of each variant in variantS

for i = 1:varSnum
    % name of variant
    curVariant = variantS.(varSfieldN{i}).dataTable;

    % Loop through each eye in variant
    for eyE = 1:4

        switch eyE
            case 1
                % Get mean
                curVariant.Left_L_oT_pupilS_meanCl =...
                    cleanPSmean(curVariant.Left_L_oT_pupilS_mean);
                % Insert Nans in missing pupil values
                curVariant.Left_L_oT_pupilS_rawCL =...
                    cleanEyeProcBF(curVariant.Left_L_oT_pupilS_raw,0);
                % Shorten all vecs to same length
                curVariant.Left_L_oT_pupilS_rawCL =...
                    shortenVEC(curVariant.Left_L_oT_pupilS_rawCL);
            case 2
                curVariant.Left_R_oT_pupilS_meanCl =...
                    cleanPSmean(curVariant.Left_R_oT_pupilS_mean);

                curVariant.Left_R_oT_pupilS_rawCL =...
                    cleanEyeProcBF(curVariant.Left_R_oT_pupilS_raw,0);

                curVariant.Left_R_oT_pupilS_rawCL =...
                    shortenVEC(curVariant.Left_R_oT_pupilS_rawCL);
            case 3
                curVariant.Right_L_oT_pupilS_meanCl =...
                    cleanPSmean(curVariant.Right_L_oT_pupilS_mean);

                curVariant.Right_L_oT_pupilS_rawCL =...
                    cleanEyeProcBF(curVariant.Right_L_oT_pupilS_raw,0);

                curVariant.Right_L_oT_pupilS_rawCL =...
                    shortenVEC(curVariant.Right_L_oT_pupilS_rawCL);
            case 4
                curVariant.Right_R_oT_pupilS_meanCl =...
                    cleanPSmean(curVariant.Right_R_oT_pupilS_mean);

                curVariant.Right_R_oT_pupilS_rawCL =...
                    cleanEyeProcBF(curVariant.Right_R_oT_pupilS_raw,0);

                curVariant.Right_R_oT_pupilS_rawCL =...
                    shortenVEC(curVariant.Right_R_oT_pupilS_rawCL);

        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save cleaned data into variantS
        variantS.(varSfieldN{i}).dataTable = curVariant;

    end

end
cd(cleanedDataLOC);
save(saveFile_name1, 'variantS');
%end



end



function [outMean] = cleanPSmean(inMean)

cleanPupil_size = inMean < 120; %create logical vector
inMean(cleanPupil_size) = nan;
outMean = inMean;

end


function [newRAW] = cleanEyeProcBF(oldRAW,plotUi)

% STEP 1
trialIDS = zeros(height(oldRAW),1,'logical');
numCpoints = zeros(height(oldRAW),1);
locCpoints = cell(height(oldRAW),1);
for tsI = 1:height(oldRAW)

    tmpTrace = oldRAW{tsI};
    tmpIndex = tmpTrace < 100;

    if sum(tmpIndex) == 0
        continue
    else
        changeIndices = find(diff(tmpIndex));
        numCpoints(tsI) = sum(diff(tmpIndex) ~= 0);
        locCpoints{tsI} = changeIndices;
        trialIDS(tsI) = true;
        % FOR PLOT CHECKING................................................
        if plotUi
            plot(tmpTrace); hold on; xline(locCpoints{tsI})
            title(num2str(tsI))
            pause
            cla
        end
        % FOR PLOT CHECKING................................................
    end
end

% STEP 2
trialIDSf = find(trialIDS);
numCptsf = numCpoints(trialIDS);
locCptsf = locCpoints(trialIDS);

% Loop through each affected trial
forwardPts = cell(height(trialIDSf),1);
backwardPts = cell(height(trialIDSf),1);
for efi = 1:length(trialIDSf)

    tmpEFi = oldRAW{trialIDSf(efi)};

    % Loop through change points
    for chpi = 1:numCptsf(efi)

        changLoc = locCptsf{efi}(chpi);

        % forward point of change
        tmpWindowStart = changLoc - 5;
        tmpWindowStop = changLoc;
        if tmpWindowStart < 1
            forwardPOINT = 1; % consider tmpNRMS block of code to find index
        else

            chek1 = true;
            while chek1

                if tmpWindowStart < 1
                    forwardPOINT = 1;
                    chek1 = false;
                else
                    tmpVecC = tmpEFi(tmpWindowStart:tmpWindowStop);
                    tmpNRMS = abs(tmpVecC/rms(tmpEFi(tmpEFi > 100))-1);

                    if ~any(tmpNRMS < 0.1) % all changed
                        tmpWindowStop = tmpWindowStart;
                        tmpWindowStart = tmpWindowStart - 5;
                        continue
                    else
                        rmsChangpt = find(tmpNRMS < 0.1,1,'first');
                        efIindexVec = tmpWindowStart:tmpWindowStop;
                        forwardPOINT = efIindexVec(rmsChangpt);
                        chek1 = false;
                    end
                end
            end
        end % END of FORWARD IF/ELSE

        % backward point of change
        tmpWindowStart = changLoc;
        tmpWindowStop = changLoc + 5;
        if tmpWindowStop > numel(tmpEFi)
            backwardPOINT = changLoc;
        else
            chek2 = true;
            while chek2

                if tmpWindowStop > numel(tmpEFi)
                    backwardPOINT = numel(tmpEFi);
                    chek2 = false;
                else
                    tmpVecC = tmpEFi(tmpWindowStart:tmpWindowStop);
                    tmpNRMS = abs(tmpVecC/rms(tmpEFi(tmpEFi > 100))-1);

                    if ~any(tmpNRMS < 0.1) % all changed
                        tmpWindowStart = tmpWindowStop;
                        tmpWindowStop = tmpWindowStart + 5;
                        continue
                    else
                        rmsChangpt = find(tmpNRMS < 0.1,1,'first');
                        efIindexVec = tmpWindowStart:tmpWindowStop;
                        backwardPOINT = efIindexVec(rmsChangpt);
                        chek2 = false;
                    end
                end

            end % END of WHILE LOOP for BACKWARD
        end % END of BACKWARD IF/ELSE
        forwardPts{efi}(chpi) = forwardPOINT;
        backwardPts{efi}(chpi) = backwardPOINT;
    end % END of CHANGE POINT LOop

end

% STEP 3 - INSERT NANS
newRAW = oldRAW;

for naNdi = 1:height(trialIDSf)
    newVecEfi = newRAW{trialIDSf(naNdi)};
    for inNan = 1:numel(forwardPts{naNdi})
        nanSTART = forwardPts{naNdi}(inNan);
        nanSTOP = backwardPts{naNdi}(inNan);
        newVecEfi(nanSTART:nanSTOP) = nan;
    end
    newRAW{trialIDSf(naNdi)} = newVecEfi;
end


end




function [newTRIMvec] = shortenVEC(oldTRIMvec)

% Find shortest row in column 6
min_varEye = min(cellfun(@(x) numel(x), oldTRIMvec, 'UniformOutput', true));
% Trim all rows in col 6 to shortest length
newTRIMvec = oldTRIMvec;
for trI = 1:height(oldTRIMvec)
    newTRIMvec{trI} = oldTRIMvec{trI}(1:min_varEye);
end


end

