
% close all
% figure;
% for ci = 1:height(curVariant.Left_R_oT_pupilS_raw)
% 
%     plot(curVariant.Left_R_oT_pupilS_raw{ci},'Color',[1 0 0 0.5])
% 
%     hold on
%     title(num2str(ci))
%     pause
%     cla
% 
% end

% MinThreshold=2
% MinDistance=5
% 3, 5, 11, 15 , 24 , 32 , 33, 37, 38, 41, 48, 49, 51, 54, 
% 55, 56, 57 , 58 , 62, 68 , 73 , 77 , 79, 82 , 85 , 87, 
% 89 , 93, 

% 1. find change points


% STEP 1
trialIDS = zeros(height(curVariant.Left_R_oT_pupilS_raw),1,'logical');
numCpoints = zeros(height(curVariant.Left_R_oT_pupilS_raw),1);
locCpoints = cell(height(curVariant.Left_R_oT_pupilS_raw),1);
for tsI = 1:height(curVariant.Left_R_oT_pupilS_raw)

    tmpTrace = curVariant.Left_R_oT_pupilS_raw{tsI};
    tmpIndex = tmpTrace < 100;

    if sum(tmpIndex) == 0
        continue
    else
        changeIndices = find(diff(tmpIndex));
        numCpoints(tsI) = sum(diff(tmpIndex) ~= 0);
        locCpoints{tsI} = changeIndices;
        trialIDS(tsI) = true;
        % FOR PLOT CHECKING................................................
        % plot(tmpTrace); hold on; xline(locCpoints{tsI})
        % title(num2str(tsI))
        % pause
        % cla
        % FOR PLOT CHECKING................................................
    end
end

% 2. Work way forwards and backwards to find the point at which change is
% less than 2% (plot(abs(t2/rms(t2)-1))) normalized root mean square

% STEP 2
trialIDSf = find(trialIDS);
numCptsf = numCpoints(trialIDS);
locCptsf = locCpoints(trialIDS);

% Loop through each affected trial
forwardPts = cell(height(trialIDSf),1);
backwardPts = cell(height(trialIDSf),1);
for efi = 1:length(trialIDSf)

    tmpEFi = curVariant.Left_R_oT_pupilS_raw{trialIDSf(efi)};

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

newRAW = curVariant.Left_R_oT_pupilS_raw;

for naNdi = 1:height(trialIDSf)
    newVecEfi = newRAW{trialIDSf(naNdi)};
    for inNan = 1:numel(forwardPts{naNdi})
        nanSTART = forwardPts{naNdi}(inNan);
        nanSTOP = backwardPts{naNdi}(inNan);
        newVecEfi(nanSTART:nanSTOP) = nan;
    end
    newRAW{trialIDSf(naNdi)} = newVecEfi;
end


% test = 1;
% 
% for nanTEST = 1:height(trialIDSf)
% 
%     tmpTrace = curVariant.Left_R_oT_pupilS_raw{trialIDSf(nanTEST)};
% 
%     % FOR PLOT CHECKING................................................
%     plot(tmpTrace);
%     coloRRs = 'rb';
%     for piPt = 1:length(forwardPts{nanTEST})
%         hold on; xline(forwardPts{nanTEST}(piPt),'Color',coloRRs(piPt));
%         xline(backwardPts{nanTEST}(piPt),'Color',coloRRs(piPt))
% 
%     end
%     title(num2str(trialIDSf(nanTEST)))
%     pause
%     cla
%     % FOR PLOT CHECKING................................................
% 
% end




% t = curVariant.Left_R_oT_pupilS_raw{3}
% findchangepts(t,Statistic="mean")
% [q,r] = findchangepts(t(1:50),'Statistic','linear','MinThreshold',10000);
% [TF,S1,S2] = ischange(t,'linear','Threshold',200);

% tStart = transpose([1, (50:50:numel(t)) + 1]); 
% tStop = transpose([(50:50:numel(t)) , numel(t)]);
% 
% numCpoints = zeros(height(tStop),1);
% locCpoints = cell(height(tStop),1);
% for tsI = 1:height(tStop)
% 
%     tTemp = t(tStart(tsI):tStop(tsI));
% 
%     [q,~] = findchangepts(tTemp,'Statistic','linear','MinThreshold',500);
% 
%     if ~isempty(q)
%         numCpoints(tsI) = numel(q);
%         locCpoints{tsI} = q;
%     end
% end


% [q,r] = findchangepts(t,'Statistic','linear','MinThreshold',10000)
% 
% varEye  = curVariant.Left_R_oT_pupilS_raw;
% for num = 1:height(varEye)
% 
%     Eye_clean = varEye{num,1};
% 
%     if sum(Eye_clean <= 120) ~= 0
%         eyeNaNindex = find(Eye_clean <= 120);
%         minEye_in = min(eyeNaNindex);
%         maxEye_in = max(eyeNaNindex);
%         eyeNaNindex2 = minEye_in-15:maxEye_in+15;
%         if min(eyeNaNindex2) < 1
%             eyeNaNindex2 = 1:maxEye_in+15;
%         elseif max(eyeNaNindex2) > length(Eye_clean)
%             eyeNaNindex2 = minEye_in-15:length(Eye_clean);
%         end
%         Eye_clean(eyeNaNindex2) = nan(length(eyeNaNindex2),1); %problem
%         varEye{num,1} = Eye_clean;
%     end
% end


figure;
for ci = 1:height(newRAW)

    plot(newRAW{ci},'Color',[1 0 0 0.5])

    hold on

end



rawMATn = zeros(100,1000);
for ci = 1:height(newRAW)

    rawMATn(ci,:) = newRAW{ci}(1:1000);

end

averP = mean(rawMATn,"omitnan");
upSTD = averP + (std(rawMATn,"omitnan")*2);
dnSTD = averP - (std(rawMATn,"omitnan")*2);

plot(averP,'LineWidth',2,'LineStyle','-','Color','k');
hold on
plot(upSTD,'LineWidth',1,'LineStyle','--','Color','r')
plot(dnSTD,'LineWidth',1,'LineStyle','--','Color','r')



