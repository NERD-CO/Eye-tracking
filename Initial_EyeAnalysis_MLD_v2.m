function [] = Initial_EyeAnalysis_MLD_v2()
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

% Location of the variant.xlsx
%excelLOC = 'C:\Users\Admin\Desktop\dataOut_AMC';
excelLOC = 'C:\Users\darwinm\Documents\MATLAB\EyeTrack';
cd(excelLOC)
varTable = readtable('variantLIST.xlsx');
% Location of eye tracking folders [from Teams]
%mainLOC = 'C:\Users\Admin\Desktop\dataOut_AMC\dataOut_AMC\eyeTrack';
mainLOC = 'C:\Users\darwinm\Documents\MATLAB\EyeTrack\eyeTrack';
cd(mainLOC)
[outFOLDS] = getfiles(mainLOC,1,nan);

for oi = 1:length(outFOLDS)

    tempCASEd = [mainLOC , filesep , outFOLDS{oi}];
    idTab = varTable(matches(varTable.Subject,outFOLDS{oi}),:);
    cd(tempCASEd)                                               %% Q for JAT %%%%%%%%
%     [tmpCedf] = getfiles(tempCASEd,2,'edf');

    allvars = unique(idTab.Variant);
    variantS = struct;
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
        
        % load in information for image paths and whether new/old in recog
        load('C:\Users\darwinm\Documents\MATLAB\EyeTrack\newOld_stimID_all.mat');
        
        % extract whether new or old from table
        newOld_inRecog = logical(stimAll.(['stimNewOld_var', varianT]));
        
        picT_L = getPICinfo(picT_L);
        tsT_L = getTRIinfo(tsT_L);
        [tsT_R, picT_R, fixT_R, ~, rawT_R] = ExtractEyeInfo_v2(eyeProcR);
        picT_R = getPICinfo(picT_R);
        tsT_R = getTRIinfo(tsT_R);

        % Find list of images new and old
        oldImages = newOld_inRecog;
        newImages = ~newOld_inRecog;

        % Get pupil size for oldImages during fixation Recog
        old_Trial = picT_R.TrialNum(oldImages);
        old_TS = picT_R.timeStamp(oldImages);
        old_catID = picT_R.CatID(oldImages);


        alleye_O1 = table;
        alleye_O2 = table;
        for oTrial = 1:length(old_Trial)
            tmpOtr = old_Trial(oTrial);

            tmpOtrTAB = tsT_R(ismember(tsT_R.TrialID, tmpOtr),:);
            startTS = tmpOtrTAB.timeStamp(tmpOtrTAB.TTLid == 1);
            endTS = tmpOtrTAB.timeStamp(tmpOtrTAB.TTLid == 2);
            
            [tsBlk_OUT] = getTSBlock(startTS,endTS,rawT_R);

            pupilS_1 = tsBlk_OUT.PupilS(:,1);
            pupilS_2 = tsBlk_OUT.PupilS(:,2);

            % NEED TO clean up PupilSize
            
            % Clean up - nans and high values
            pos_1 = [tsBlk_OUT.PosX(:,1) , tsBlk_OUT.PosY(:,1)];
            pos_1c = cleanUPpos(pos_1);
            pos_2 = [tsBlk_OUT.PosX(:,2) , tsBlk_OUT.PosY(:,2)];
            pos_2c = cleanUPpos(pos_2);

            [eye1 , eye2] = createEYEtable(pupilS_1,pupilS_2,pos_1c,pos_2c);
            eye1.catID = old_catID(oTrial);
            eye2.catID = old_catID(oTrial);

            alleye_O1 = [alleye_O1 ; eye1];
            alleye_O2 = [alleye_O2 ; eye2];           
            
        end
        varianTNum = ['var',num2str(allvars(vi))];
        variantS.(varianTNum).eye1O = alleye_O1;
        variantS.(varianTNum).eye2O = alleye_O2;
        
        % ADD CODE HERE ---- CREATE A NEW FOR LOOP SIMILAR to 59-84 for NEW
        % IMAGES!!!!!!!!! HAVE FUN!!!!!!!!!!!!!
        
        %%%%%%%% MLD %%%%%%%%%%
            
        % Get pupil size for newImages during fixation Recog
        new_Trial = picT_R.TrialNum(newImages);
        new_TS = picT_R.timeStamp(newImages);
        new_catID = picT_R.CatID(newImages);

        alleye_N1 = table;
        alleye_N2 = table;
        
        for nTrial = 1:length(new_Trial)
            tmpNtr = new_Trial(nTrial);

            tmpNtrTAB = tsT_R(ismember(tsT_R.TrialID, tmpNtr),:);
            startTS = tmpNtrTAB.timeStamp(tmpNtrTAB.TTLid == 1);
            endTS = tmpNtrTAB.timeStamp(tmpNtrTAB.TTLid == 2);
            
            [tsBlk_OUT] = getTSBlock(startTS,endTS,rawT_R);

            pupilS_1 = tsBlk_OUT.PupilS(:,1);
            pupilS_2 = tsBlk_OUT.PupilS(:,2);

            % NEED TO clean up PupilSize
            
            % Clean up - nans and high values
            pos_1 = [tsBlk_OUT.PosX(:,1) , tsBlk_OUT.PosY(:,1)];
            pos_1c = cleanUPpos(pos_1);
            pos_2 = [tsBlk_OUT.PosX(:,2) , tsBlk_OUT.PosY(:,2)];
            pos_2c = cleanUPpos(pos_2);

            [eye1 , eye2] = createEYEtable(pupilS_1,pupilS_2,pos_1c,pos_2c);
            eye1.catID = new_catID(nTrial);
            eye2.catID = new_catID(nTrial);

            alleye_N1 = [alleye_N1 ; eye1];
            alleye_N2 = [alleye_N2 ; eye2];
        
        
        end
        
    variantS.(varianTNum).eye1N = alleye_N1;
    variantS.(varianTNum).eye2N = alleye_N2;   
        
    % subject
    saveFname = ['eyeData_',outFOLDS{oi},'.mat'];
    %saveLOC = 'C:\Users\Admin\Desktop\dataOut_AMC\dataOut_AMC\eyeDATA';
    saveLOC = 'C:\Users\darwinm\Documents\MATLAB\EyeTrack\eyeTrack\eyeDATA';
    cd(saveLOC)
    save(saveFname,"variantS");

end


end






















% function [outfiles] = getfiles(dirIN,stage,ftype)
% 
% cd(dirIN)
% switch stage
%     case 1
% 
%         foldeS = dir();
%         foldeS2 = {foldeS.name};
%         foldeS3 = foldeS2(~ismember(foldeS2,{'.','..'}));
%         outfiles = foldeS3;
%     case 2
% 
%         filES = dir(['*.',ftype]);
%         filES2 = {filES.name};
%         outfiles = filES2;
% 
% end
% 
% 
% end


% % Pic info
% function [picINFOtab] = getPICinfo(inTable)
% 
% picLOCa = inTable.PicLocation;
% fparTs = split(picLOCa,'\');
% numCAT = fparTs(:,11);
% catSubn = cellfun(@(x) str2double(x(1)), numCAT, 'UniformOutput',true);
% catID = cellfun(@(x) x(2:end), numCAT, 'UniformOutput',false);
% picJPG = inTable.Picture;
% picNUMs = cellfun(@(x) split(x,'.'), picJPG, 'UniformOutput',false);
% picNUM = cellfun(@(x) str2double(x{1}), picNUMs, 'UniformOutput',true);
% picINFOtab = inTable;
% picINFOtab.CatNUM = catSubn;
% picINFOtab.CatID = catID;
% picINFOtab.PicNUM = picNUM;
% 
% tmpCombine = cell(length(catSubn),1);
% for pi = 1:length(catSubn)
%     tmpCombine{pi} = [num2str(catSubn(pi)),'.',num2str(picNUM(pi))];
% end
% picINFOtab.CatPICid = tmpCombine;
% 
% end


% % Trial info
% function [trialINFOtab] = getTRIinfo(inTable)
% 
% ttlIDt = inTable.TTLid;
% 
% ttltrialID = zeros(length(ttlIDt),1);
% trialcount = 0;
% for ti = 1:length(ttltrialID)
%     tmpT = ttlIDt(ti);
%     if tmpT == 1
%         trialcount = trialcount + 1;
%         ttltrialID(ti) = trialcount;
%     else
%         ttltrialID(ti) = trialcount;
%     end
% 
% end
% 
% % start and end
% ttltrialID(ttlIDt == 55) = 0;
% ttltrialID(ttlIDt == 66) = 0;
% 
% trialINFOtab = inTable;
% trialINFOtab.TrialID = ttltrialID;
% 
% 
% end
% 
% 
% 
% % Get TS info
% function [tsBlk_OUT] = getTSBlock(startI,endI,rawT)
% 
% [~, eyeTTL1_i] = min(abs(double(startI) - rawT.Time)); 
% [~, eyeTTL2_i] = min(abs(double(endI) - rawT.Time)); 
% 
% tsBlk_OUT = rawT(eyeTTL1_i:eyeTTL2_i,:);
% 
% 
% end
% 
% 
% % Clean up Pos
% function [cleanPOS] = cleanUPpos(posIN)
% 
% % nans
% posT1 = posIN(~isnan(posIN(:,1)),:);
% posT2 = posT1(~isnan(posT1(:,2)),:);
% % cleanPOS = posT2;
% 
% % long values
% posL1 = posT2(posT2(:,1) > 1,:);
% posL2 = posL1(posL1(:,2) > 1,:);
% 
% posS1 = smoothdata(posL2(:,1),'gaussian',20);
% posS2 = smoothdata(posL2(:,2),'gaussian',20);
% 
% cleanPOS = [posS1 , posS2];
% 
% end



% function [outEye1 , outEye2] = createEYEtable(ps1,ps2,pos1,pos2)
% 
% for ei = 1:2
%     switch ei
%         case 1
%             eyedata = {pos1};
%             eyeCen = {mean(pos1)};
%             eyeSD = {std(pos1)};
%             Q3pos = quantile(pos1,0.75);
%             Q1pos = quantile(pos1,0.25);
%             eyeCD = (Q3pos - Q1pos) / (Q3pos + Q1pos);
%             eyedist = {pdist2(eyeCen{1},eyedata{1},'euclidean')};
%             pupdata = {ps1};
%             pupCen = mean(ps1);
%             pupSD = std(ps1);
%             Q3pup = quantile(ps1,0.75);
%             Q1pup = quantile(ps1,0.25);
%             pupCD = (Q3pup - Q1pup) / (Q3pup + Q1pup);
% 
%             outEye1 = table(eyedata,eyeCen,eyeSD,eyeCD,eyedist,...
%                 pupdata,pupCen,pupSD,pupCD,'VariableNames',...
%                 {'oT_posit_raw','oT_posit_cen','oT_posit_sd','oT_posit_cd',...
%                 'oT_posit_dist','oT_pupilS_raw','oT_pupilS_mean','oT_pupilS_sd',...
%                 'oT_pupilS_cd'});
% 
% 
%         case 2
% 
%             eyedata = {pos2};
%             eyeCen = {mean(pos2)};
%             eyeSD = {std(pos2)};
%             Q3pos = quantile(pos2,0.75);
%             Q1pos = quantile(pos2,0.25);
%             eyeCD = (Q3pos - Q1pos) / (Q3pos + Q1pos);
%             eyedist = {pdist2(eyeCen{1},eyedata{1},'euclidean')};
%             pupdata = {ps2};
%             pupCen = mean(ps2);
%             pupSD = std(ps2);
%             Q3pup = quantile(ps2,0.75);
%             Q1pup = quantile(ps2,0.25);
%             pupCD = (Q3pup - Q1pup) / (Q3pup + Q1pup);
% 
%             outEye2 = table(eyedata,eyeCen,eyeSD,eyeCD,eyedist,...
%                 pupdata,pupCen,pupSD,pupCD,'VariableNames',...
%                 {'oT_posit_raw','oT_posit_cen','oT_posit_sd','oT_posit_cd',...
%                 'oT_posit_dist','oT_pupilS_raw','oT_pupilS_mean','oT_pupilS_sd',...
%                 'oT_pupilS_cd'});
% 
%     end
% 
% 
% end
end