

%% Data location
pcname = getenv('COMPUTERNAME');

switch pcname
    case 'DESKTOP-FAGRV5G'
        mainDIR = 'E:\Dropbox\SfN_2022\dataOut_AMC';
        dataDIR = [mainDIR, filesep , 'dataOut_AMC\eyeDATA'];
end


cd(dataDIR)

%% Load the data

load("eyeData_AMC_PY22NO12.mat")

%%
eye1O_PS = variantS.var3.eye1O.oT_pupilS_raw;

% plot(eye1O_PS{1})

% 1. Trim raw length
% 2. Insert NaNs for values that exceed 1 SD +/-
allmean = mean(variantS.var3.eye1O.oT_pupilS_mean);
dropEind = variantS.var3.eye1O.oT_pupilS_mean > allmean;
eye1Okeep = variantS.var3.eye1O(dropEind,:);
mineye1K = min(cellfun(@(x) length(x), eye1Okeep.oT_pupilS_raw));
newminMean = mean(cellfun(@(x) mean(x), eye1Okeep.oT_pupilS_raw));
newminStd = std(cellfun(@(x) mean(x), eye1Okeep.oT_pupilS_raw));
newminThdd = newminMean - (newminStd*1);
newminThdu = newminMean + (newminStd*1);


%%


eyeMeandata = zeros(height(eye1Okeep),mineye1K);
for i = 1:height(eye1Okeep)

    tmpE = eye1Okeep.oT_pupilS_raw{i};
    newminMean = mean(tmpE);
    newminStd = std(tmpE);
    newminThdd = newminMean - (newminStd*1.25);
    newminThdu = newminMean + (newminStd*1.25);

    tmpEnI = tmpE > newminThdu | tmpE < newminThdd;
    tmpEnN = tmpE;
    tmpEnN(tmpEnI) = nan;

    tmpTrN = tmpEnN(1:mineye1K);
    eyeMeandata(i,:) = tmpTrN;

end

meanMean = mean(eyeMeandata,1,'omitnan');
stdMean = std(eyeMeandata,'omitnan');
stdU = meanMean + (stdMean * 1.5);
stdD = meanMean - (stdMean * 1.5);

smM = smoothdata(meanMean,'sgolay',100);
stU = smoothdata(stdU,'sgolay',100);
stD = smoothdata(stdD,'sgolay',100);

xpatch = [1:mineye1K , fliplr(1:mineye1K)];
ypatch = [stD , fliplr(stU)];

patch(xpatch,ypatch,'Red','FaceAlpha',0.5)
hold on
plot(smM,'k','LineWidth',3)
xlim([1 mineye1K])
ylim([0 550])


%% Plot types - pupil size
plottypeps = 1;

switch plottypeps
    case 1



    case 2



    case 3




end


%% Plot types
plottype = 1;

switch plottype
    case 1


        SD = chart.ScatterDensity("XData", x, ...
            "YData", y );

        SD.SizeData = 400;
        grid( SD, "off" )

    case 2



    case 3
        load fisheriris.mat;
        x = meas(:,1);
        y = meas(:,2);

        scatterhist(x,y,'Group',species,'Kernel','on')


    case 4



end