function [] = eyeTRACKproc_PupilLocation(cleanedDataLOC, saveLOC, ptID)

% Set data paths
if ~exist(cleanedDataLOC,'dir')
    mkdir(cleanedDataLOC)
end

% [For every file in eyeDATA, loop through to run analyses?]
% save cleaned files in new folder to only have raw data in this folder
% need to cd to new folder to save

% NO SHORTEN - because NOW Table INCLUDES time MARKERS

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

        % FOR ERROR CHECKING
        % disp(eyE)

        switch eyE
            case 1
                % DISTANCE VALUES
                % Insert Nans in missing pupil values
                curVariant.Left_L_oT_pupilS_rawCL =...
                    cleanEye_Posit_DIST(curVariant.Left_L_oT_posit_dist,0,...
                    curVariant.LearnTTLplus);
                % Shorten all vecs to same length
                % curVariant.Left_L_oT_pupilS_rawCL =...
                %     shortenVEC(curVariant.Left_L_oT_pupilS_rawCL);
            case 2
                curVariant.Left_R_oT_pupilS_meanCl =...
                    cleanPSmean(curVariant.Left_R_oT_pupilS_mean);

                curVariant.Left_R_oT_pupilS_rawCL =...
                    cleanEyeProcBF(curVariant.Left_R_oT_pupilS_raw,0,...
                    curVariant.LearnTTLplus);

                % curVariant.Left_R_oT_pupilS_rawCL =...
                %     shortenVEC(curVariant.Left_R_oT_pupilS_rawCL);
            case 3
                curVariant.Right_L_oT_pupilS_meanCl =...
                    cleanPSmean(curVariant.Right_L_oT_pupilS_mean);

                curVariant.Right_L_oT_pupilS_rawCL =...
                    cleanEyeProcBF(curVariant.Right_L_oT_pupilS_raw,0,...
                    curVariant.RecogTTLplus);

                % curVariant.Right_L_oT_pupilS_rawCL =...
                %     shortenVEC(curVariant.Right_L_oT_pupilS_rawCL);
            case 4
                curVariant.Right_R_oT_pupilS_meanCl =...
                    cleanPSmean(curVariant.Right_R_oT_pupilS_mean);

                curVariant.Right_R_oT_pupilS_rawCL =...
                    cleanEyeProcBF(curVariant.Right_R_oT_pupilS_raw,0,...
                    curVariant.RecogTTLplus);
                % curVariant.Right_R_oT_pupilS_rawCL =...
                %     shortenVEC(curVariant.Right_R_oT_pupilS_rawCL);

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

function [finalACSI] = cleanEye_Posit_DIST(oldRAW,plotUi,TTLtable)

test = 1;

oldRAW1 = oldRAW;
for oi = 1:height(oldRAW1)

    nanChk = sum(isnan(oldRAW1{oi}))/numel(oldRAW{oi});

    if nanChk == 1
        oldRAW1{oi} = NaN;
        continue;
    end

    % plot(oldRAW1{oi})
    diffAbOf = abs(diff(oldRAW{oi}));
    qlINE = quantile(diffAbOf,0.945);
    cutINDEX = diffAbOf > qlINE;
    cutINDEXi = find(cutINDEX) + 1;
    oldRawt = oldRAW{oi};
    oldRawt(cutINDEXi) = NaN;

    
    subplot(2,1,1)
    hold on
    plot(oldRAW{oi},'k-')
    subplot(2,1,2)
    plot(oldRawt,'r-')
    hold on
    % subplot(2,1,2)
    % plot(diffAbOf)
    % yline(qlINE)
    
    blockCinner = 0;
    blockCouter = 0;
    blockS = nan(100,1);
    blockInds = nan(100,2);
    for li = 1:numel(oldRawt)
        
        tmpSample = oldRawt(li);

        if ~isnan(tmpSample)
            blockCinner = blockCinner + 1;
        elseif isnan(tmpSample) && blockCinner ~= 0
            blockCouter = blockCouter + 1;
            blockS(blockCouter) = blockCinner;
            endLoc = li - 1;
            startLoc = endLoc - blockCinner + 1;
            blockInds(blockCouter,:) = [startLoc , endLoc];
            blockCinner = 0;
        else
            continue
        end
    end

    % Remove extra lines
    nanEend = find(isnan(blockS),1,'first');
    blockSf = blockS(1:nanEend-1);
    blockIndsf = blockInds(1:nanEend-1,:);

    % Get inner blocks
    for ibi = 1:length(blockSf)

    end

    oldRawt2 = oldRawt;
    for bi = 1:length(blockSf)

        oldRawt2(blockIndsf(bi,1):blockIndsf(bi,2)) = nan;



    end
    plot(oldRawt2,'r-')


    

    


    % plot(abs(diff(oldRAW{oi})))

end


end



function [outMean] = cleanPSmean(inMean)

cleanPupil_size = inMean < 120; %create logical vector
inMean(cleanPupil_size) = nan;
outMean = inMean;

end


function [finalACSI] = cleanEyeProcBF(oldRAW,plotUi,TTLtable)


% PRE STEP 1 - SUPER outlier TRIALS not single events
oldRAW_SO = cell2mat(oldRAW);
oldRmean = mean(oldRAW_SO,'omitnan');
oldRstd = std(oldRAW_SO,'omitnan');
% oldRrms = rms(oldRAW_SO);
oldUthreh = oldRmean + (oldRstd * 2);
oldDthreh = 25;

% STACK EVERY TRIAL within EYE / CONDITION

for tsI = 1:height(oldRAW)

    tmpTrace = oldRAW{tsI};
    tmpIndex = tmpTrace < oldDthreh | tmpTrace > oldUthreh;
    tmpTrace(tmpIndex) = nan;

    oldRAW{tsI} = tmpTrace;

end


% PRE STEP 2 - SUPER outlier TRIALS IN LENGTH
for tsI = 1:height(oldRAW)

    tmpTrace = oldRAW{tsI};
    tmpTTL = TTLtable{tsI}.ELNKint(5);

    if tmpTTL > 5000
        tmpTrace(5001:end) = nan;

    else 
        tmpTrace(tmpTTL + 500:end) = nan;
    end
    oldRAW{tsI} = tmpTrace;
end


% figure;
% for nfi = 1:length(oldRAW)
%     tmpEFi = oldRAW{nfi};
%     plot(tmpEFi,'Color',[0.3 0.3 0.3 0.3],'LineWidth',1.5)
%     hold on
%     xline(TTLtable{nfi}.ELNKint(5))
% end
% title('Remove OUTlier trials')
% % pause(10)
% close all



% STEP 1
trialIDS = zeros(height(oldRAW),1,'logical');
numCpoints = zeros(height(oldRAW),1);
locCpoints = cell(height(oldRAW),1);
for tsI = 1:height(oldRAW)

    tmpTrace = oldRAW{tsI};
    tmpTUpth = mean(tmpTrace, 'omitnan') + (std(tmpTrace, 'omitnan')*2);
    tmpTDpth = mean(tmpTrace, 'omitnan') - (std(tmpTrace, 'omitnan')*2);

    tmpIndex = tmpTrace < tmpTDpth | tmpTrace > tmpTUpth; % FIGURE OUT UPPER LIMIT??????????

    % plot(tmpTrace); yline(tmpTUpth); yline(tmpTDpth)
    % pause 
    % cla

    if sum(tmpIndex) == 0
        continue
    else
        changeIndices = find(diff(tmpIndex));
        numCpoints(tsI) = sum(diff(tmpIndex) ~= 0);
        locCpoints{tsI} = changeIndices;
        trialIDS(tsI) = true;
        % FOR PLOT CHECKING................................................
        if plotUi
            plot(tmpTrace);
            hold on
            plot(tmpIndex*100)
            plot(changeIndices,ones(numel(changeIndices),1)*500,'k.','MarkerSize',25)
            xline(locCpoints{tsI})
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

% PLOT STEP 3 
% figure;
% for nfi = 1:length(newRAW)
%     tmpEFi = newRAW{nfi};
%     plot(tmpEFi,'Color',[0.3 0.3 0.3 0.3],'LineWidth',1.5)
%     hold on
% end
% title('JAT Clean up')
% close all

% STEP 4 FIRST DERIVATIVE MASK
dernewRAW = newRAW;
for drIi = 1:length(newRAW)

    tmpNEW = newRAW{drIi};
    % First derivative of x
    dx = diff(tmpNEW);
    % Standard deviation of dx
    sigma = std(dx,'omitnan');
    % 2.5 times the standard deviation
    result = 2.5 * sigma;
    absDer = abs(dx);
    firDerSDmask = absDer > result;

    derNEW = tmpNEW;
    derNEW(firDerSDmask) = nan;

    dernewRAW{drIi} = derNEW;
end

% PLOT STEP 4 
% figure;
% for nfi = 1:length(dernewRAW)
%     tmpEFi = dernewRAW{nfi};
%     plot(tmpEFi,'Color',[0.3 0.3 0.3 0.3],'LineWidth',1.5)
%     hold on
% end
% title('1st Deriv')
% close all

% STEP 5
% 3rd ORDER ButterWORTH with 4 Hz cut off
% Define the filter order
order = 3;
% Define the cutoff frequency
fc = 4; % Hz
% Define the sampling frequency
fs = 1000; % Hz
% Normalize the cutoff frequency to the Nyquist frequency
Wn = fc / (fs/2);
% Design the filter coefficients
[b, a] = butter(order, Wn);
% Filter the signal using the filter coefficients

fourHzbutter = dernewRAW;
for buIi = 1:length(dernewRAW)
    tmpNEWt = dernewRAW{buIi};
    % Interpolate NaNs in the data
    nanLOCS = isnan(tmpNEWt);
    t = 1:length(tmpNEWt);

    % FOR ERROR CHECKING
    % disp(buIi)

    allnanCheck = sum(isnan(tmpNEWt))/numel(tmpNEWt);
    if allnanCheck > 0.85
        fourHzbutter{buIi} = transpose(tmpNEWt);
        continue
    end

    data_interp1 = interp1(t(~isnan(tmpNEWt)), tmpNEWt(~isnan(tmpNEWt)), t);
    if any(isnan(data_interp1))
        dataMEAN = mean(data_interp1,'omitnan');
        data_interp1(isnan(data_interp1)) = dataMEAN;
    end
    % Apply the filter to the interpolated data
    filtered_data = filtfilt(b, a, data_interp1);
    % filtered_signal = filtfilt(b, a, tmpNEWt);
    fbTfil = filtered_data;
    % PUT NANS back
    fbTfil(nanLOCS) = nan;
    fourHzbutter{buIi} = fbTfil;
end

finalNEWot = fourHzbutter;


% PLOT STEP 5 
figure;
for nfi = 1:length(finalNEWot)
    tmpEFi = finalNEWot{nfi};
    plot(tmpEFi,'Color',[0.3 0.3 0.3 0.3],'LineWidth',1.5)
    hold on
end
title('4Hz low cut')
close all
maxLEN = max(cellfun(@(x) numel(x), finalNEWot, 'UniformOutput',true));
nanMAT = nan(length(finalNEWot),maxLEN);
for ni = 1:length(finalNEWot)
    curLEN = numel(finalNEWot{ni});
    nanMAT(ni,1:curLEN) = finalNEWot{ni};
end


% FIND last NAN
firstAllNAN = find(sum(isnan(nanMAT),1) == length(finalNEWot),1,'first');

nanMATn = nanMAT(:,1:firstAllNAN);

fracNAN = zeros(1,width(nanMATn));
for ini = 1:width(nanMATn)

    fracNAN(ini) = sum(isnan(nanMATn(:,ini)))/height(nanMATn);

end

columnTrim = find(fracNAN >= .80, 1, 'first');

nanMATn2 = nanMATn(:,1:columnTrim);


[corVals] = computeASCI_eye_v2(nanMATn2);
% ccCUToff = quantile(corVals,0.12);
if std(corVals,'omitnan') < 0.29
    ccCUToff = mean(corVals,'omitnan') - (std(corVals,'omitnan')*1.65);
else
    ccCUToff = mean(corVals,'omitnan') - (std(corVals,'omitnan')*0.25);
end
corInds = corVals > ccCUToff;
% corInds = corVals > 0.6;

finalACSI = finalNEWot;
finalACSI(~corInds) = {nan};


% PLOT STEP 4 
% nanMATntemp = nanMATn(corInds,:);
% figure;
% for nfi = 1:height(nanMATntemp)
%     tmpEFi = nanMATntemp(nfi,:);
%     plot(tmpEFi,'Color',[0.3 0.3 0.3 0.3],'LineWidth',1.5)
%     hold on
% end
% title('1st Deriv')
% close all


% SMOOTH GUASSION or SANGOLY

% finalSMOOTH = fourHzbutter;
% for smIi = 1:length(fourHzbutter)
% 
%     tmpNEWts = fourHzbutter{smIi};
% 
%     finSMdata = smoothdata(tmpNEWts , "sgolay" , 60);
% 
%     finalSMOOTH{smIi} = finSMdata;
% end
% 
% figure;
% for nfi = 1:length(finalSMOOTH)
%     tmpEFi = finalSMOOTH{nfi};
%     plot(tmpEFi,'Color',[0.3 0.3 0.3 0.3],'LineWidth',1.5)
%     hold on
% end
% title('Final smooth')



end




% function [newTRIMvec] = shortenVEC(oldTRIMvec)
% 
% % Find shortest row in column 6
% min_varEye = min(cellfun(@(x) numel(x), oldTRIMvec, 'UniformOutput', true));
% % Trim all rows in col 6 to shortest length
% newTRIMvec = oldTRIMvec;
% for trI = 1:height(oldTRIMvec)
%     newTRIMvec{trI} = oldTRIMvec{trI}(1:min_varEye);
% end
% 
% 
% end

