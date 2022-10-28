function [] = Initial_EyeAnalysis()
% 'The events coorespond to the TTL markers for each trial. ', ...
%     'For the learning trials, the TTL markers are the following: 55 = start of the experiment, ', ...
%     '1 = stimulus ON, 2 = stimulus OFF, 3 = Question Screen Onset, ', ...
%     '20 = Yes (21 = NO) during learning, 6 = End of Delay after Response, ', ...
%     '66 = End of Experiment. For the recognition trials, the TTL markers are the following: ', ...
%     '55 = start of experiment, 1 = stimulus ON, 2 = stimulus OFF, 3 = Question Screen Onset, ' ...
%     '31:36 = Confidence (Yes vs. No) response, 66 = End of Experiment'

% Loop through files - determine if learn or retrieve


% learning = [55, 1, 2, 3 20, 21, 6, 66];
% recog = [55, 1, 2, 3, 31:36, 66];

excelLOC = 'C:\Users\Admin\Desktop\dataOut_AMC';
cd(excelLOC)
varTable = readtable('variantLIST.xlsx');
mainLOC = 'C:\Users\Admin\Desktop\dataOut_AMC\dataOut_AMC\eyeTrack';
cd(mainLOC)
[outFOLDS] = getfiles(mainLOC,1,nan);

for oi = 1:length(outFOLDS)

    tempCASEd = [mainLOC , filesep , outFOLDS{oi}];
    idTab = varTable(matches(varTable.Subject,outFOLDS{oi}),:);
    cd(tempCASEd)
%     [tmpCedf] = getfiles(tempCASEd,2,'edf');

    allvars = unique(idTab.Variant);
    for vi = 1:length(allvars)
        vartiTAB =  idTab(ismember(idTab.Variant,allvars(vi)),:);
        learnEDF = vartiTAB.EDF(matches(vartiTAB.Block,'learn'));
        recogEDF = vartiTAB.EDF(matches(vartiTAB.Block,'recog'));
        varianT = num2str(vi);

        [eyeProcL] = Edf2Mat_UCH(learnEDF{1}, outFOLDS{oi}, 'Learn',...
            varianT, tempCASEd, tempCASEd);
        [eyeProcR] = Edf2Mat_UCH(recogEDF{1}, outFOLDS{oi}, 'Recog',...
            varianT, tempCASEd, tempCASEd);

        [tsT_L, picT_L, fixT_L, ~, rawT_L] = ExtractEyeInfo_v2(eyeProcL);
        picT_L = getPICinfo(picT_L);
        tsT_L = getTRIinfo(tsT_L);
        [tsT_R, picT_R, fixT_R, ~, rawT_R] = ExtractEyeInfo_v2(eyeProcR);
        picT_R = getPICinfo(picT_R);
        tsT_R = getTRIinfo(tsT_R);

        % Find list of images new and old
        oldImages = ismember(picT_R.CatPICid, picT_L.CatPICid);
        newImages = ~ismember(picT_R.CatPICid, picT_L.CatPICid);

        % Get pupil size for oldImages during fixation Learn
        old_Trial = picT_R.TrialNum(oldImages);
        old_TS = picT_R.timeStamp(oldImages);

        for oTrial = 1:length(old_Trial)
            tmpOtr = old_Trial(oTrial);

            tmpOtrTAB = tsT_R(ismember(tsT_R.TrialID, tmpOtr),:);
            startTS = tmpOtrTAB.timeStamp(tmpOtrTAB.TTLid == 1);
            endTS = tmpOtrTAB.timeStamp(tmpOtrTAB.TTLid == 2);
            
            [tsBlk_OUT] = getTSBlock(startTS,endTS,rawT_R);

            pupilS_1 = tsBlk_OUT.PupilS(:,1);
            pupilS_2 = tsBlk_OUT.PupilS(:,2);

            pos_1 = [tsBlk_OUT.PosX(:,1) , tsBlk_OUT.PosY(:,1)];
            pos_1c = cleanUPpos(pos_1);
            pos_2 = [tsBlk_OUT.PosX(:,2) , tsBlk_OUT.PosY(:,2)];
            pos_2c = cleanUPpos(pos_2);

            % Clean up - nans and high values
            plot(pos_2c(:,1),pos_2c(:,2),'o')

            test = 1;




        end





        % Get pupil size for oldImages during fixation Recog



    end




end


end






















function [outfiles] = getfiles(dirIN,stage,ftype)

cd(dirIN)
switch stage
    case 1

        foldeS = dir();
        foldeS2 = {foldeS.name};
        foldeS3 = foldeS2(~ismember(foldeS2,{'.','..'}));
        outfiles = foldeS3;
    case 2

        filES = dir(['*.',ftype]);
        filES2 = {filES.name};
        outfiles = filES2;

end


end


% Pic info
function [picINFOtab] = getPICinfo(inTable)

picLOCa = inTable.PicLocation;
fparTs = split(picLOCa,'\');
numCAT = fparTs(:,11);
catSubn = cellfun(@(x) str2double(x(1)), numCAT, 'UniformOutput',true);
catID = cellfun(@(x) x(2:end), numCAT, 'UniformOutput',false);
picJPG = inTable.Picture;
picNUMs = cellfun(@(x) split(x,'.'), picJPG, 'UniformOutput',false);
picNUM = cellfun(@(x) str2double(x{1}), picNUMs, 'UniformOutput',true);
picINFOtab = inTable;
picINFOtab.CatNUM = catSubn;
picINFOtab.CatID = catID;
picINFOtab.PicNUM = picNUM;

tmpCombine = cell(length(catSubn),1);
for pi = 1:length(catSubn)
    tmpCombine{pi} = [num2str(catSubn(pi)),'.',num2str(picNUM(pi))];
end
picINFOtab.CatPICid = tmpCombine;

end


% Trial info
function [trialINFOtab] = getTRIinfo(inTable)

ttlIDt = inTable.TTLid;

ttltrialID = zeros(length(ttlIDt),1);
trialcount = 0;
for ti = 1:length(ttltrialID)
    tmpT = ttlIDt(ti);
    if tmpT == 1
        trialcount = trialcount + 1;
        ttltrialID(ti) = trialcount;
    else
        ttltrialID(ti) = trialcount;
    end

end

% start and end
ttltrialID(ttlIDt == 55) = 0;
ttltrialID(ttlIDt == 66) = 0;

trialINFOtab = inTable;
trialINFOtab.TrialID = ttltrialID;


end



% Get TS info
function [tsBlk_OUT] = getTSBlock(startI,endI,rawT)

[~, eyeTTL1_i] = min(abs(double(startI) - rawT.Time)); 
[~, eyeTTL2_i] = min(abs(double(endI) - rawT.Time)); 

tsBlk_OUT = rawT(eyeTTL1_i:eyeTTL2_i,:);


end


% Clean up Pos
function [cleanPOS] = cleanUPpos(posIN)

% nans
posT1 = posIN(~isnan(posIN(:,1)),:);
posT2 = posT1(~isnan(posT1(:,2)),:);
% cleanPOS = posT2;

% long values
posL1 = posT2(posT2(:,1) > 1,:);
posL2 = posL1(posL1(:,2) > 1,:);

posS = smoothdata(posL2,2,'gaussian',20);

cleanPOS = posS;

end