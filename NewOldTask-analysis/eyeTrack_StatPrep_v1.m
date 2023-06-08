function [] = eyeTrack_StatPrep_v1(cleanedDataLOC, exclLOC, ptID ,conditiON, subStep)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



cd(exclLOC)
varTable = readtable('variantLIST.xlsx');
ptTABLE = varTable(matches(varTable.Subject,ptID),:);

cd(cleanedDataLOC)
ptFname = ['cl_eyeData_',ptID,'.mat'];
load(ptFname, 'variantS')

switch conditiON

    % '21' = YES or '20' = NO
    case 'learn'
        
        variNUM = fieldnames(variantS);
        for vi = 1:length(variNUM)

            varPtInfoLearn = ptTABLE(matches(ptTABLE.Block,'learn') &...
                ismember(ptTABLE.Variant,str2double(variNUM{vi}(end))),:);

            varDATA = variantS.(variNUM{vi}).dataTable;

            if matches(varPtInfoLearn.Eye2Use{1},'Left')
                learnDATA = varDATA.Left_L_oT_pupilS_rawCL;
            else
                learnDATA = varDATA.Right_L_oT_pupilS_rawCL;
            end

            learnRespN = varDATA.LearnResp;
            learnResp = learnRespN;
            learnResp(matches(learnRespN,'20')) = repmat({'Yes'},sum(matches(learnRespN,'20')),1);
            learnResp(matches(learnRespN,'21')) = repmat({'No'},sum(matches(learnRespN,'21')),1);

            learnGT = varDATA.Recog_CatID;

            learnTTL = varDATA.LearnTTLplus;

            switch str2double(variNUM{vi}(end))
                case 1
                    gt_learnR = 'smallAnimal';
                case 2

                case 3

            end

            switch subStep
                case 1
                    subYes_correct = matches(learnResp,'Yes') & matches(learnGT , gt_learnR);
                    subNo_correct = ~matches(learnResp,'Yes') & matches(learnGT , gt_learnR);
                    subYes_incorrect = matches(learnResp,'Yes') & ~matches(learnGT , gt_learnR);
                    subNo_incorrect  = ~matches(learnResp,'Yes') & ~matches(learnGT , gt_learnR);
            end

            variNum = variNUM{vi};
            inData = learnDATA;
            trialGroupIDs = {'Subject_Yes_Correct','Subject_No_Correct',...
                'Subject_Yes_Incorrect','Subject_No_Incorrect'}; 
            trialGroupInds = [subYes_correct , subNo_correct ,...
                subYes_incorrect , subNo_incorrect];

            [statsOutput] = getStatsTable(variNum , inData , learnTTL ,...
                trialGroupIDs, trialGroupInds);

        end % VARIANT LOOP

    case 'recog'




end


%%%%% CREATE Variable length table/struct to output trial pupil info based
%%%%% on selections above
%%%%% AM HERE 6/7/2023


















end




function [statsOutput] = getStatsTable(variantID , pupilData , ttlDATA ,...
    trialGroupIDS, trialGroupIndicies)



statsOutput = struct;
groupAll = ones(1000,1)*-1000;
epochAll = ones(1000,1)*-1000;
meanAll = ones(1000,1)*-1000;
sdAll = ones(1000,1)*-1000;
aCount = 1;

for ci = 1:length(trialGroupIDS)

    % Trial Group Indicies
    catROWSi1 = trialGroupIndicies(:,ci);
    % Missing Row Indicies
    catROWSi2 = cellfun(@(x) numel(x) , pupilData, 'UniformOutput', true) == 1;
    catROWSi3 = (catROWSi1 & ~catROWSi2);

    tmpCatRaw = pupilData(catROWSi3);
    tmpCatTTL = ttlDATA(catROWSi3);

    % sampleLEN_L = 1402; % % 1000 ms + 400 ms + 2 ms

    % CREATE Matrix average baseline
    % -400 -- stimOn
    % stimOn -- stimOff
    % stimOff -- Question
    % Delay -- +200
    tmpCatCellRaw = cell(sum(catROWSi3),4); % NEED TO CHECK
    for bi = 1:height(tmpCatTTL)
        startVALS = [1 2 3 6];
        endVALS   = [2 3 4 7];
        for tli = 1:4
            ttlEpochStart = tmpCatTTL{bi}.ELNKint(startVALS(tli));
            ttlEpochEnd = tmpCatTTL{bi}.ELNKint(endVALS(tli));
            tmpCatCellRaw{bi,tli} = tmpCatRaw{bi}(ttlEpochStart:ttlEpochEnd);
        end
    end

    % Clean up sample offsets and create matrices for RAW
    tmpCatMatRaw = cell(1,4);
    for rri = 1:width(tmpCatCellRaw)
        tmpColumn = tmpCatCellRaw(:,rri);

        tmpColLen = cellfun(@(x) numel(x), tmpColumn, 'UniformOutput',true);

        if numel(unique(tmpColLen)) == 1
            tmpCatMatRaw{rri} = cell2mat(tmpColumn);
        else
            minColLen = min(unique(tmpColLen));
            shortenCol = cellfun(@(x) x(1:minColLen), tmpColumn, 'UniformOutput',false);
            tmpCatMatRaw{rri} = cell2mat(shortenCol);
        end
    end

    

    % Leave RAW Alone 
    % Create normalized by baseline
    tmpCatMatNormB = cell(1,4); % Epochs
    for trialI = 1:sum(catROWSi3)

        baseLINE =  tmpCatMatRaw{1}(trialI,:);
        meanBASE = mean(baseLINE,'omitnan');

        for epochi = 1:width(tmpCatCellRaw)
            tmpCatMatNormB{epochi}(trialI,:) = ((tmpCatMatRaw{epochi}(trialI,:) / meanBASE) * 100)-100;
        end
    end

    %%% DO SOMETHING WITH BLOCKS -----------------------------------
 
    statsOutput.RawMat{ci} = tmpCatMatRaw;
    statsOutput.NormMat{ci} = tmpCatMatNormB;
    
    for sst = 1:width(tmpCatMatNormB)

        groupAll(aCount) = ci;
        epochAll(aCount) = sst;
        meanAll(aCount) = mean(tmpCatMatNormB{sst}, 'all', 'omitnan');
        sdAll(aCount) = std(tmpCatMatNormB{sst}, [], 'all', 'omitnan');
        aCount = aCount + 1;

    end



end

% Remove nans from summary table
remINDEXs = find(groupAll == -1000, 1, 'first');
sumArray = [groupAll(1:remINDEXs-1) , epochAll(1:remINDEXs-1) ,...
             meanAll(1:remINDEXs-1) , sdAll(1:remINDEXs-1)];
% Create summary table
summTable = array2table(sumArray,'VariableNames',{'Group','Epoch','Mean','SD'});


statsOutput.GroupIDS = trialGroupIDS;
statsOutput.VarID = variantID;

end

