function [] = Initial_plots_EyePosition_fun()

% Data location
pcname = getenv('COMPUTERNAME');

switch pcname
    case 'DESKTOP-FAGRV5G'
        mainDIR = 'E:\Dropbox\SfN_2022\dataOut_AMC';
        dataDIR = [mainDIR, filesep , 'dataOut_AMC\eyeDATA'];
    case 'DESKTOP-EGFQKAI'
        mainDIR = 'C:\Users\johna\Dropbox\SfN_2022\dataOut_AMC';
        dataDIR = [mainDIR, filesep , 'dataOut_AMC\eyeDATA'];
end


cd(dataDIR)

% Load the data

%%%%% TO DO ADD MORE TIME BEFORE AND AFTER

load("eyeData_AMC_PY22NO12.mat")

%%
eye1O_PS = variantS.var3.eye1O.oT_pupilS_raw;
eye1O_PSm = variantS.var3.eye1O.oT_pupilS_mean;

eye1N_PS = variantS.var3.eye1N.oT_pupilS_raw;
eye1N_PSm = variantS.var3.eye1N.oT_pupilS_mean;

% 1. Trim raw length
% 2. Insert NaNs for values that exceed 1 SD +/-
[mineye_1O , eye1O_keep] = eyePrep(eye1O_PS,eye1O_PSm);
[mineye_1N , eye1N_keep] = eyePrep(eye1N_PS,eye1N_PSm);

[eyeMData1O] = getEyeMdata(eye1O_keep, mineye_1O);
[eyeMData1N] = getEyeMdata(eye1N_keep, mineye_1N);

plotEyepatch(eyeMData1O , mineye_1O , 'Red')
hold on
plotEyepatch(eyeMData1N , mineye_1N , 'Green')



end





function [mineye1K , eye1Okeep] = eyePrep(eyeRaw,eyeMean)

allmean = mean(eyeMean);
dropEind = eyeMean > allmean;
eye1Okeep = eyeRaw(dropEind,:);
mineye1K = min(cellfun(@(x) length(x), eye1Okeep));

end



function [eyeMData] = getEyeMdata(eyeKeep, minEyelen)

eyeMData = zeros(height(eyeKeep),minEyelen);
for i = 1:height(eyeKeep)

    tmpE = eyeKeep{i};
    newminMean = mean(tmpE);
    newminStd = std(tmpE);
    newminThdd = newminMean - (newminStd*1.25);
    newminThdu = newminMean + (newminStd*1.25);

    tmpEnI = tmpE > newminThdu | tmpE < newminThdd;
    tmpEnN = tmpE;
    tmpEnN(tmpEnI) = nan;

    tmpTrN = tmpEnN(1:minEyelen);
    eyeMData(i,:) = tmpTrN;

end


end





function [] = plotEyepatch(eyeMData , mineye1K, coloR)

meanMean = mean(eyeMData,1,'omitnan');
stdMean = std(eyeMData,'omitnan');
stdU = meanMean + (stdMean * 1.5);
stdD = meanMean - (stdMean * 1.5);

smM = smoothdata(meanMean,'sgolay',100);
stU = smoothdata(stdU,'sgolay',100);
stD = smoothdata(stdD,'sgolay',100);

xpatch = [1:mineye1K , fliplr(1:mineye1K)];
ypatch = [stD , fliplr(stU)];

patch(xpatch,ypatch,coloR,'FaceAlpha',0.5)
hold on
plot(smM,'k','LineWidth',3)
xlim([1 mineye1K])
ylim([0 550])


end



