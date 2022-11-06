function [] = Initial_plots_EyePosition_fun(dataID,varNN)


figure;
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

var2use = ['var',num2str(varNN)];

switch dataID
    case 1
        load("eyeData_AMC_PY22NO12.mat")
        var2use = 'var3';

    case 2
        load("eyeData_AMC_PY22NO13.mat")

    case 3
        load("eyeData_AMC_PY22NO16.mat")


end



%%
eye1O_PS = variantS.(var2use).eye1O.oT_pupilS_raw;
eye1O_PSm = variantS.(var2use).eye1O.oT_pupilS_mean;

eye1N_PS = variantS.(var2use).eye1N.oT_pupilS_raw;
eye1N_PSm = variantS.(var2use).eye1N.oT_pupilS_mean;

% 1. Trim raw length
% 2. Insert NaNs for values that exceed 1 SD +/-
[mineye_1O , eye1O_keep] = eyePrep(eye1O_PS,eye1O_PSm);
[mineye_1N , eye1N_keep] = eyePrep(eye1N_PS,eye1N_PSm);

[eyeMData1O] = getEyeMdata(eye1O_keep, mineye_1O);
[eyeMData1N] = getEyeMdata(eye1N_keep, mineye_1N);

tmpYm = mean(cellfun(@(x) mean(x) , eye1N_keep));
% [0 0.4470 0.7410]
% [0.8500 0.3250 0.0980]
[tmpYO] = plotEyepatch(eyeMData1O , mineye_1O , [0 0.4470 0.7410] , tmpYm);
hold on
plotEyepatch(eyeMData1N , mineye_1N , [0.8500 0.3250 0.0980],tmpYO);
legend('Old images','','New images')
% compare swarmchart

eyeMData1Om = mean(eyeMData1O,1,'omitnan');
eyeMData1Nm = mean(eyeMData1N,1,'omitnan');
figure;

x1 = ones(1,length(eyeMData1Om));
x2 = 2 * ones(1,length(eyeMData1Nm));
s1 = swarmchart(x1,eyeMData1Om,5,'filled');
s1.MarkerFaceAlpha = 0.5;
hold on
s2 = swarmchart(x2,eyeMData1Nm,5,'filled');
s2.MarkerFaceAlpha = 0.5;
ylabel('Pupil Size')
xticks([1 2])
xticklabels({'Old Images','New Images'})

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





function [ylimOUT] = plotEyepatch(eyeMData , mineye1K, coloR, baseYM)

meanMean = mean(eyeMData,1,'omitnan');
stdMean = std(eyeMData,'omitnan');
stdU = meanMean + (stdMean * 1.5);
stdD = meanMean - (stdMean * 1.5);

smM = smoothdata(meanMean,'sgolay',100);
stU = smoothdata(stdU,'sgolay',100);
stD = smoothdata(stdD,'sgolay',100);

xpatch = [1:mineye1K , fliplr(1:mineye1K)];
ypatch = [stD , fliplr(stU)];

p = patch(xpatch,ypatch,'r');
p.FaceColor = coloR;
p.EdgeColor = 'none';
p.FaceAlpha = 0.25;

hold on
plot(smM,'Color',coloR,'LineWidth',3)
xlim([1 mineye1K])

if round(max(stdU)) + 25 > baseYM
    ylim([0 round(max(stdU)) + 25])
    ylimOUT = round(max(stdU)) + 25;
else
    ylim([0 baseYM])
    ylimOUT = baseYM;
end

% Change yaxis
yl = ylim;
yticks([yl(1) round(yl(2)/2) yl(2)])

% Change xaxis
xl = xlim;
xticks([xl(1) round(xl(2)/2) xl(2)])

xlabel('Time in samples')
ylabel('Raw pupil size')


end



