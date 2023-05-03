function [] = Initial_EyeAnalysis_MLD_v8(excelLOC , mainLOC, ptID, saveLOC)
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

% Location of variant.xlsx
cd(excelLOC)
varTable = readtable('variantLIST.xlsx');
% Location of eye tracking folders
cd(mainLOC)
[outFOLDS] = getfiles(mainLOC,1,nan);

% Match specific patient ID to location in outFOLDS
ptIdx = find(strcmp(outFOLDS, ptID));

tempCASEd = [mainLOC , filesep , outFOLDS{ptIdx}];
idTab = varTable(matches(varTable.Subject,outFOLDS{ptIdx}),:);
cd(tempCASEd)

allvars = unique(idTab.Variant);
variantS = struct;
for vi = 1:length(allvars)

    % Get Variant and Condition

    variantTAB = idTab(ismember(idTab.Variant,allvars(vi)),:);
    learn_matFile = variantTAB.MatFile{matches(variantTAB.Block,'learn')};

    [tsT_L, picT_L, ~, ~, rawT_L] = ExtractEyeInfo_v2(learn_matFile);

    picT_L = getPICinfo(picT_L);
    tsT_L = getTRIinfo(tsT_L);

    % load in information for image paths and whether new/old in recog
    cd(excelLOC)
    load('newOld_stimID_all.mat','stimAll');

    % add in compareAns code here
    % Extract variables for ground truth old vs new- 1=present in learn block
    groundTruthRecog_var1 = stimAll.stimNewOld_var1;
    groundTruthRecog_var2 = stimAll.stimNewOld_var2;
    groundTruthRecog_var3 = stimAll.stimNewOld_var3;

    % Patient answers - need to change per pt - can add input arg for pt file
    cd(tempCASEd); % Used to create variable 'outData' with TTL values

    % Convert behavioral txt files to mat files
    recogTXT = variantTAB.Behavior{matches(variantTAB.Block,'recog')};

    % Recog file
    block = 'recog'; %Will be either 'learn' or 'recog'
    patientID = idTab.Subject{1};

    % Run function and load in new .mat file
    [behavFILE_recog] = NewOldTxttoMat_v2(recogTXT,patientID,vi,block, tempCASEd);
    load(behavFILE_recog, 'outData') %loc and file name as input

    cd(mainLOC);

    % Extract answers from outData - 31:36 = Confidence (Yes vs. No) response
    % 31,32,33=no,new & 34,35,36=yes,old
    ttlValues = str2double(outData.taskinformation.TTLvalue);
    confRatings = ttlValues(ttlValues(:,1) >= 31 & ttlValues(:,1) <= 36,:);

    % Logical array of confRatings: GREATER than 34
    confRatings_logical = logical(confRatings >= 34); %1=yes,old & 0=no,new

    confVal_1 = logical(confRatings == 31); %1=31, 0= not 31
    confVal_2 = logical(confRatings == 32); %1=32, 0= not 32
    confVal_3 = logical(confRatings == 33); %1=33, 0= not 33
    confVal_4 = logical(confRatings == 34); %1=34, 0= not 34
    confVal_5 = logical(confRatings == 35); %1=35, 0= not 35
    confVal_6 = logical(confRatings == 36); %1=36, 0= not 36

    confVal_34 = logical(confRatings == 33 | confRatings == 34); %1=33-34, 0=not 33-34
    confVal_12 = logical(confRatings == 31 | confRatings == 32); %1=31-32, 0=not 31-32
    confVal_56 = logical(confRatings == 35 | confRatings == 36); %1=35-36, 0=not 35-36
    confVal_1256 = ~logical(confVal_34);                         %1=31+32+35+36, 0=33+34
    confVal_16 = logical(confRatings == 31 | confRatings == 36); %1=31+36, 0=32+33+34+35

    % Can compare ground truth to patient answers
    % 1=yes,old & 0=no,new
    % Combine arrays into a table
    switch allvars(vi)
        case 1
            groundTruthRecog = groundTruthRecog_var1;
        case 2
            groundTruthRecog = groundTruthRecog_var2;
        case 3
            groundTruthRecog = groundTruthRecog_var3;
    end

    % CONSIDER RE-NAMING RECOG INFO
    recogINFO = table(confRatings, confRatings_logical, groundTruthRecog, ...
        confVal_1, confVal_2, confVal_3, confVal_4, confVal_5, confVal_6,...
        confVal_34, confVal_12, confVal_56, confVal_1256, confVal_16, 'VariableNames', ...
        {'TTLs','confRatings', 'groundTruth', 'New_sure', 'New_less_sure', 'New_unsure',...
        'Old_unsure', 'Old_less_sure', 'Old_sure', 'All_unsure', ...
        'New_sureLess_sure', 'Old_sureLess sure', 'All_sureLess_sure', 'All_sure'});

    cd(tempCASEd)
    recog_matFile = variantTAB.MatFile{matches(variantTAB.Block,'recog')};

    [tsT_R, picT_R, ~, ~, rawT_R] = ExtractEyeInfo_v2(recog_matFile);

    picT_R = getPICinfo(picT_R);
    tsT_R = getTRIinfo(tsT_R);

    trialID_R = unique(tsT_R.TrialID(tsT_R.TrialID ~= 0));
    trialID_L = unique(tsT_L.TrialID(tsT_L.TrialID ~= 0));

    [leftEYE_Recog , rightEYE_Recog] = getEYErawEpoch(rawT_R , tsT_R , trialID_R);

    [leftEYE_Learn , rightEYE_Learn] = getEYErawEpoch(rawT_L , tsT_L , trialID_L);

    [outTABLE] = fixANDcombineTables(leftEYE_Learn,rightEYE_Learn,leftEYE_Recog,...
        rightEYE_Recog , picT_L, picT_R, recogINFO);

    varianTNum = ['var',num2str(allvars(vi))];

    variantS.(varianTNum).dataTable = outTABLE;

end
% subject
saveFname = ['eyeData_',outFOLDS{ptIdx},'.mat'];
cd(saveLOC)
save(saveFname,"variantS");


%end

end % END MAIN FUNCTION



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
catID = fparTs(:,11);
% catSubn = cellfun(@(x) str2double(x(1)), numCAT, 'UniformOutput',true);
% catID = cellfun(@(x) x(2:end), numCAT, 'UniformOutput',false);
picJPG = inTable.Picture;
picNUMs = cellfun(@(x) split(x,'.'), picJPG, 'UniformOutput',false);
picNUM = cellfun(@(x) str2double(x{1}), picNUMs, 'UniformOutput',true);
picINFOtab = inTable;
picINFOtab.CatID = catID;
picINFOtab.PicNUM = picNUM;

% tmpCombine = cell(length(catSubn),1);
% for pi = 1:length(catSubn)
%     tmpCombine{pi} = [num2str(catSubn(pi)),'.',num2str(picNUM(pi))];
% end
% picINFOtab.CatPICid = tmpCombine;

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
% posT1 = posIN(~isnan(posIN(:,1)),:);
% posT2 = posT1(~isnan(posT1(:,2)),:);
% cleanPOS = posT2;

% if any(isnan(posIN))
%     keyboard
% end

% long values
% posL1 = posT2(posT2(:,1) > 1,:);
% posL2 = posL1(posL1(:,2) > 1,:);

posS1 = smoothdata(posIN(:,1),'gaussian',40);
posS2 = smoothdata(posIN(:,2),'gaussian',40);

cleanPOS = [posS1 , posS2];

end



function [outEye1 , outEye2] = createEYEtable(ps1,ps2,pos1,pos2)

for ei = 1:2
    switch ei
        case 1
            eyedata = {pos1};
            eyeCen = {mean(pos1)};
            eyeSD = {std(pos1)};
            Q3pos = quantile(pos1,0.75);
            Q1pos = quantile(pos1,0.25);
            eyeCD = (Q3pos - Q1pos) / (Q3pos + Q1pos);
            eyedist = {pdist2(eyeCen{1},eyedata{1},'euclidean')};
            pupdata = {ps1};
            pupCen = mean(ps1);
            pupSD = std(ps1);
            Q3pup = quantile(ps1,0.75);
            Q1pup = quantile(ps1,0.25);
            pupCD = (Q3pup - Q1pup) / (Q3pup + Q1pup);

            outEye1 = table(eyedata,eyeCen,eyeSD,eyeCD,eyedist,...
                pupdata,pupCen,pupSD,pupCD,'VariableNames',...
                {'oT_posit_raw','oT_posit_cen','oT_posit_sd','oT_posit_cd',...
                'oT_posit_dist','oT_pupilS_raw','oT_pupilS_mean','oT_pupilS_sd',...
                'oT_pupilS_cd'});


        case 2

            eyedata = {pos2};
            eyeCen = {mean(pos2)};
            eyeSD = {std(pos2)};
            Q3pos = quantile(pos2,0.75);
            Q1pos = quantile(pos2,0.25);
            eyeCD = (Q3pos - Q1pos) / (Q3pos + Q1pos);
            eyedist = {pdist2(eyeCen{1},eyedata{1},'euclidean')};
            pupdata = {ps2};
            pupCen = mean(ps2);
            pupSD = std(ps2);
            Q3pup = quantile(ps2,0.75);
            Q1pup = quantile(ps2,0.25);
            pupCD = (Q3pup - Q1pup) / (Q3pup + Q1pup);

            outEye2 = table(eyedata,eyeCen,eyeSD,eyeCD,eyedist,...
                pupdata,pupCen,pupSD,pupCD,'VariableNames',...
                {'oT_posit_raw','oT_posit_cen','oT_posit_sd','oT_posit_cd',...
                'oT_posit_dist','oT_pupilS_raw','oT_pupilS_mean','oT_pupilS_sd',...
                'oT_pupilS_cd'});

    end


end
end





function [leftEYE , rightEYE] = getEYErawEpoch(rawTIME , ttlTABLE , trialsOfInt)

% Make N by 2 matrix of fieldname + value type
bothVarNames = [['oT_posit_raw', "cell"]; ...
			['oT_posit_cen', "cell"]; ...
			['oT_posit_sd', "cell"]; ...
			['oT_posit_cd', "double"]; ...
			['oT_posit_dist', "cell"]; ...
			['oT_pupilS_raw', "cell"];...
            ['oT_pupilS_mean', "double"];...
            ['oT_pupilS_sd', "double"];...
            ['oT_pupilS_cd', "double"]];

% Make table using fieldnames & value types from above
leftEYE = table('Size',[0,size(bothVarNames,1)],... 
	'VariableNames', bothVarNames(:,1),...
	'VariableTypes', bothVarNames(:,2));

rightEYE = table('Size',[0,size(bothVarNames,1)],... 
	'VariableNames', bothVarNames(:,1),...
	'VariableTypes', bothVarNames(:,2));


for tttrialir = 1:length(trialsOfInt)
    tmpOtr = trialsOfInt(tttrialir);

    tmpOtrTAB = ttlTABLE(ismember(ttlTABLE.TrialID, tmpOtr),:);
    startTS = tmpOtrTAB.timeStamp(tmpOtrTAB.TTLid == 1);
    endTS = tmpOtrTAB.timeStamp(tmpOtrTAB.TTLid == 2);

    [tsBlk_OUT] = getTSBlock(startTS,endTS,rawTIME);

    pupilS_1 = tsBlk_OUT.PupilS(:,1);
    pupilS_2 = tsBlk_OUT.PupilS(:,2);

    % NEED TO clean up PupilSize

    % Clean up - nans and high values
    pos_1 = [tsBlk_OUT.PosX(:,1) , tsBlk_OUT.PosY(:,1)];
    pos_1c = cleanUPpos(pos_1);
    pos_2 = [tsBlk_OUT.PosX(:,2) , tsBlk_OUT.PosY(:,2)];
    pos_2c = cleanUPpos(pos_2);

    [left1 , right2] = createEYEtable(pupilS_1,pupilS_2,pos_1c,pos_2c);
    % left1.catID = old_catID(oTrial);
    % right2.catID = old_catID(oTrial);

    leftEYE(tttrialir,:) = left1;
    rightEYE(tttrialir,:) = right2;

end



end




function [outTABLE] = fixANDcombineTables(left_learn,right_learn,left_recog,...
    right_recog, picture_learn, picture_recog, recogInfo)

% Left Learn
leftlearnVars = left_learn.Properties.VariableNames;
leftlearnVarsN = cellfun(@(x) ['Left_L_', x], leftlearnVars, 'UniformOutput',false);
left_learn = renamevars(left_learn, leftlearnVars, leftlearnVarsN);

% Left Recog
leftRecogVars = left_recog.Properties.VariableNames;
leftRecoVarsN = cellfun(@(x) ['Left_R_', x], leftRecogVars, 'UniformOutput',false);
left_recog = renamevars(left_recog, leftRecogVars, leftRecoVarsN);

% Right Learn
rightlearnVars = right_learn.Properties.VariableNames;
rightlearnVarsN = cellfun(@(x) ['Right_L_', x], rightlearnVars, 'UniformOutput',false);
right_learn = renamevars(right_learn, rightlearnVars, rightlearnVarsN);

% Right Recog
rightRecogVars = right_recog.Properties.VariableNames;
rightRecogVarsN = cellfun(@(x) ['Right_R_', x], rightRecogVars, 'UniformOutput',false);
right_recog = renamevars(right_recog, rightRecogVars, rightRecogVarsN);

% Picture learn
learnPicVars = picture_learn.Properties.VariableNames;
learnPicVarsN = cellfun(@(x) ['Learn_', x], learnPicVars, 'UniformOutput',false);
picture_learn = renamevars(picture_learn, learnPicVars, learnPicVarsN);

% Picture recog
recogPicVars = picture_recog.Properties.VariableNames;
recogPicVarsN = cellfun(@(x) ['Recog_', x], recogPicVars, 'UniformOutput',false);
picture_recog = renamevars(picture_recog, recogPicVars, recogPicVarsN);

outTABLE = [left_learn , left_recog , right_learn , right_recog ,...
    picture_learn , picture_recog , recogInfo];




end


