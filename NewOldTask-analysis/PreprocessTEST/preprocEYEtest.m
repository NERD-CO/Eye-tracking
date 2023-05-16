close all
load("pc_test.mat","testData");

oldRAW = testData;

plotUi = 0;

% STEP 1
trialIDS = zeros(height(oldRAW),1,'logical');
numCpoints = zeros(height(oldRAW),1);
locCpoints = cell(height(oldRAW),1);
for tsI = 1:height(oldRAW)

    tmpTrace = oldRAW{tsI};
    tmpIndex = tmpTrace < 100 | tmpTrace > 1100;

    if sum(tmpIndex) == 0
        continue
    else
        changeIndices = find(diff(tmpIndex));
        numCpoints(tsI) = sum(diff(tmpIndex) ~= 0);
        locCpoints{tsI} = changeIndices;
        trialIDS(tsI) = true;
        % FOR PLOT CHECKING................................................
        if plotUi
            plot(tmpTrace)
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



figure;
for nfi = 1:length(newRAW)
    tmpEFi = newRAW{nfi};
    plot(tmpEFi,'Color',[0.3 0.3 0.3 0.3],'LineWidth',1.5)
    hold on
end
title('JAT Clean up')


% FIRST DERIVATIVE MASK
% First derivative of x
% dx = diff(x);
% % Standard deviation of dx
% sigma = std(dx);
% % Three times the standard deviation
% result = 3 * sigma;

dernewRAW = newRAW;
for drIi = 1:length(newRAW)

    tmpNEW = newRAW{drIi};
    dx = diff(tmpNEW);
    sigma = std(dx,'omitnan');
    result = 2.5 * sigma;
    absDer = abs(dx);
    firDerSDmask = absDer > result;

    derNEW = tmpNEW;
    derNEW(firDerSDmask) = nan;

    dernewRAW{drIi} = derNEW;

end

figure;
for nfi = 1:length(dernewRAW)
    tmpEFi = dernewRAW{nfi};
    plot(tmpEFi,'Color',[0.3 0.3 0.3 0.3],'LineWidth',1.5)
    hold on
end
title('1st Deriv')


% 3rd ORDER ButterWORTH with 4 Hz cut off
% Sampling frequency
% Fs = 1000; % Hz
% % Cutoff frequency
% Fc = 4; % Hz
% % Normalize the cutoff frequency to the Nyquist frequency
% Wn = Fc / (Fs/2);
% % Design the filter
% [b, a] = butter(3, Wn, 'low');
% % Load the input signal
% % Filter the signal using the filter coefficients

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
% Load the data to be smoothed
% load('data.mat'); % Replace with your own data

% Apply the filter to the data
% smoothed_data = filtfilt(b, a, data);


fourHzbutter = dernewRAW;
for buIi = 1:length(dernewRAW)
    tmpNEWt = dernewRAW{buIi};
    % Interpolate NaNs in the data
    nanLOCS = isnan(tmpNEWt);
    t = 1:length(tmpNEWt);
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


figure;
for nfi = 1:length(fourHzbutter)
    tmpEFi = fourHzbutter{nfi};
    plot(tmpEFi,'Color',[0.3 0.3 0.3 0.3],'LineWidth',1.5)
    hold on
end
title('4Hz low cut')


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









