function [] = getEYEstatsPlots()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



pcname = getenv('COMPUTERNAME');

switch pcname
    case 'DESKTOP-FAGRV5G'
        mainDIR = 'E:\Dropbox\SfN_2022\dataOut_AMC';
        dataDIR = [mainDIR, filesep , 'dataOut_AMC\eyeDATA'];
    case 'DESKTOP-EGFQKAI'
        mainDIR = 'C:\Users\johna\Dropbox\SfN_2022\dataOut_AMC';
        dataDIR = [mainDIR, filesep , 'dataOut_AMC\eyeDATA'];
    case 'DESKTOP-I5CPDO7'
        mainDIR = 'D:\Dropbox\SfN_2022\dataOut_AMC';
        dataDIR = [mainDIR, filesep , 'dataOut_AMC\eyeDATA'];
end


cd(dataDIR)


matDIR1 = dir('*.mat');
matDIR2 = {matDIR1.name};

varSS = {3,[1,3],1};


totCount = 1;
newALL = zeros(4,1);
oldALL = zeros(4,1);
pvalAll = zeros(4,1);
for mi = 1:length(matDIR2)

    load(matDIR2{mi})

    % Get var number
    for vi = 1:length(varSS{mi})

        var2use = ['var',num2str(varSS{mi}(vi))];

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

        eyeMData1Om = mean(eyeMData1O,1,'omitnan');
        eyeMData1Nm = mean(eyeMData1N,1,'omitnan');

        newALL(totCount) = mean(eyeMData1Nm,'omitnan');
        oldALL(totCount) = mean(eyeMData1Om,'omitnan');
        [~,pvalAll(totCount)] = ttest2(eyeMData1Om,eyeMData1Nm);
        totCount = totCount + 1;

    end

end

oldColor = [0 0.4470 0.7410];
newColor = [0.8500 0.3250 0.0980];

xAx = [ones(4,1) ; ones(4,1)*2];
yAx = [oldALL ; newALL];
coloAll = [repmat(oldColor,4,1) ; repmat(newColor,4,1)];

line(transpose([ones(4,1) , ones(4,1)*2]),...
     transpose([oldALL , newALL]),'LineWidth',0.5,'Color','k')
hold on
scatter(xAx,yAx,60,coloAll,'filled')
xlim([0.5 2.5])
xticks([1 2])
xticklabels({'Old images', 'New images'})
ylimC = ylim;
yticks([linspace(ylimC(1),ylimC(2),3)])
ylabel('Pupil size')

[~,pvalG,~,stats] = ttest2(oldALL,newALL)

axis square

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