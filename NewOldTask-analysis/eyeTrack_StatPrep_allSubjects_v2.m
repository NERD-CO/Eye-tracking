function [allVariants] = eyeTrack_StatPrep_allSubjects_v2(cleanedDataLOC, exclLOC,conditiON, subStep)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%%% STEP 1
% Get rid of the ptID input argument 
ptID = ["AMC_PY21NO05", "AMC_PY22NO09", "AMC_PY22NO12"];
cd(cleanedDataLOC)
allPatients = struct; % create empty struct
for allPTs = 1:length(ptID)

    % Load variantS for the patient
    ptID2 = ptID(allPTs);
    ptFname = append('cl_eyeData_',ptID2,'.mat');
    load(ptFname, 'variantS')

    % Add to allPatients struct
    allPatients.(ptID{allPTs}) = variantS;
end 

%ptFname = ['cl_eyeData_',ptID,'.mat'];
%load(ptFname, 'variantS')

%%%% STEP 2

%allPatients.MW5 = variantS;




% Create a block of code to loop through and import all subjects
% Create a struct for all patients
% patientStruct ----|
%                   |----Patient1
%                   |----Patient2
%                   |----Patient3
%                   |----Patient4
%                   |----Patient5 -----|
%                                      |-----variantS


%%%% STEP 3
% Within patient combine variants -- [for example stick together var1 and
% var2 tables into one big table

% Q: within for loop on line 9?

% A: loop through allPatients.MW5
%create variable that will have all field names from top 
patientFieldNames = fieldnames(allPatients);
allVariants = table;
for pi = 1:length(patientFieldNames)
    tempPatient = allPatients.(patientFieldNames{pi});
    tempVariants = fieldnames(tempPatient);
    for vi = 1:length(tempVariants)
        allVariants = [allVariants; tempPatient.(tempVariants{vi})];
    end
end


%%%% STEP 4
% Create one big table of all patients 





% cd(exclLOC)
% varTable = readtable('variantLIST.xlsx');
% ptTABLE = varTable(matches(varTable.Subject,ptID),:);



switch conditiON

    % '21' = YES or '20' = NO
    case 'learn'

        variNUM = fieldnames(variantS);
        allVariants = struct;
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
            allVariants.(variNUM{vi}) = statsOutput;

        end % VARIANT LOOP

    case 'recog'

   
        %%% GET RID OF FOR loop - JUST ONE TABLE
        %%% REPLACE the assignment to varDATA with the ONE BIG TABLE


        variNUM = fieldnames(variantS);
        allVariants = struct;
        for vi = 1:length(variNUM)

            varPtInfoRecog = ptTABLE(matches(ptTABLE.Block,'recog') &...
                ismember(ptTABLE.Variant,str2double(variNUM{vi}(end))),:);

            varDATA = variantS.(variNUM{vi}).dataTable;

            if matches(varPtInfoRecog.Eye2Use{1},'Left')
                RecogDATA = varDATA.Left_R_oT_pupilS_rawCL;
            else
                RecogDATA = varDATA.Right_R_oT_pupilS_rawCL;
            end

            RecogTTL = varDATA.RecogTTLplus;
            RecogResp = varDATA.RecogResp;
            RecogGT = varDATA.groundTruth;

            switch subStep
                case 1
                    subNewY1_cor = RecogResp == 31 & RecogGT;
                    subNewY2_cor = RecogResp == 32 & RecogGT;
                    subNewY3_cor = RecogResp == 33 & RecogGT;
                    subNewN4_cor = RecogResp == 34 & ~RecogGT;
                    subNewN5_cor = RecogResp == 35 & ~RecogGT;
                    subNewN6_cor = RecogResp == 36 & ~RecogGT;
                    subNewY1_ncr = RecogResp == 31 & ~RecogGT;
                    subNewY2_ncr = RecogResp == 32 & ~RecogGT;
                    subNewY3_ncr = RecogResp == 33 & ~RecogGT;
                    subNewN4_ncr = RecogResp == 34 & RecogGT;
                    subNewN5_ncr = RecogResp == 35 & RecogGT;
                    subNewN6_ncr = RecogResp == 36 & RecogGT;

                    variNum = variNUM{vi};
                    inData = RecogDATA;
                    trialGroupIDs = {'Subject_Y1_cor','Subject_Y2_cor',...
                        'Subject_Y3_cor','Subject_N4_cor','Subject_N5_cor',...
                        'Subject_N6_cor','Subject_Y1_ncr','Subject_Y2_ncr',...
                        'Subject_Y3_ncr','Subject_N4_ncr','Subject_N5_ncr'...
                        'Subject_N6_ncr'};
                    trialGroupInds = [subNewY1_cor , subNewY2_cor , subNewY3_cor,...
                        subNewN4_cor , subNewN5_cor , subNewN6_cor , subNewY1_ncr,...
                        subNewY2_ncr , subNewY3_ncr , subNewN4_ncr , subNewN5_ncr,...
                        subNewN6_ncr];

                    [statsOutput] = getStatsTable(variNum , inData , RecogTTL ,...
                        trialGroupIDs, trialGroupInds);
                    allVariants.(variNUM{vi}) = statsOutput;

                case 2
                    subNewYa_cor = ismember(RecogResp, 31:33) & RecogGT;
                    subNewNa_cor = ismember(RecogResp, 34:36) & ~RecogGT;
                    subNewYa_ncr = ismember(RecogResp, 31:33) & ~RecogGT;
                    subNewNa_ncr = ismember(RecogResp, 34:36) & RecogGT;

                    variNum = variNUM{vi};
                    inData = RecogDATA;
                    trialGroupIDs = {'Subject_YA_cor','Subject_NA_cor',...
                        'Subject_YA_ncr','Subject_NA_ncr'};
                    trialGroupInds = [subNewYa_cor , subNewNa_cor , subNewYa_ncr,...
                        subNewNa_ncr];

                    [statsOutput] = getStatsTable(variNum , inData , RecogTTL ,...
                        trialGroupIDs, trialGroupInds);
                    allVariants.(variNUM{vi}) = statsOutput;

                case 3
                    subNewYa2_cor = ismember(RecogResp, 31:32) & RecogGT;
                    subNewNa2_cor = ismember(RecogResp, 35:36) & ~RecogGT;
                    subNewUa2_cor = (ismember(RecogResp, 33) & RecogGT) |...
                        (ismember(RecogResp, 34) & ~RecogGT);
                    subNewYa2_ncr = ismember(RecogResp, 31:32) & ~RecogGT;
                    subNewNa2_ncr = ismember(RecogResp, 35:36) & RecogGT;
                    subNewUa2_ncr = (ismember(RecogResp, 33) & ~RecogGT) |...
                        (ismember(RecogResp, 34) & RecogGT);

                    variNum = variNUM{vi};
                    inData = RecogDATA;
                    trialGroupIDs = {'Subject_Y2_cor','Subject_N2_cor',...
                        'Subject_U2_cor','Subject_Y2_ncr','Subject_N2_ncr'...
                        'Subject_U2_ncr'};
                    trialGroupInds = [subNewYa2_cor , subNewNa2_cor , subNewUa2_cor,...
                        subNewYa2_ncr , subNewNa2_ncr , subNewUa2_ncr];

                    [statsOutput] = getStatsTable(variNum , inData , RecogTTL ,...
                        trialGroupIDs, trialGroupInds);
                    allVariants.(variNUM{vi}) = statsOutput;

                case 4
                    subNewYa1_cor = ismember(RecogResp, 31) & RecogGT;
                    subNewNa1_cor = ismember(RecogResp, 36) & ~RecogGT;
                    subNewMa1_cor = (ismember(RecogResp, 32:33) & RecogGT) |...
                        (ismember(RecogResp, 34:35) & ~RecogGT);
                    subNewYa1_ncr = ismember(RecogResp, 31) & ~RecogGT;
                    subNewNa1_ncr = ismember(RecogResp, 36) & RecogGT;
                    subNewMa1_ncr = (ismember(RecogResp, 32:33) & ~RecogGT) |...
                        (ismember(RecogResp, 34:35) & RecogGT);

                    variNum = variNUM{vi};
                    inData = RecogDATA;
                    trialGroupIDs = {'Subject_Y1_cor','Subject_N1_cor',...
                        'Subject_M1_cor','Subject_Y1_ncr','Subject_N1_ncr'...
                        'Subject_M1_ncr'};
                    trialGroupInds = [subNewYa1_cor , subNewNa1_cor , subNewMa1_cor,...
                        subNewYa1_ncr , subNewNa1_ncr , subNewMa1_ncr];

                    [statsOutput] = getStatsTable(variNum , inData , RecogTTL ,...
                        trialGroupIDs, trialGroupInds);
                    allVariants.(variNUM{vi}) = statsOutput;


            end


        end


end





switch conditiON

    case 'learn'

        for vi = 1:length(variNUM)

            ydata = allVariants.(variNUM{vi}).StatTable.PupilSize;
            condition2u = allVariants.(variNUM{vi}).StatTable.GroupID;
            epoch2u = allVariants.(variNUM{vi}).StatTable.EpochID;

            [~,sumStats_condition,stats_condition] = anova1(ydata,condition2u,"off");
            [~,sumStats_epcoh,stats_epoch] = anova1(ydata,epoch2u,"off");

            [results_condition,~,~,gnames_condition] = multcompare(stats_condition);
            [results_epoch,~,~,gnames_epoch] = multcompare(stats_epoch);

            % Create Table
            tbl_condition = array2table(results_condition,"VariableNames", ...
                ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
            tbl_condition.("Group A") = gnames_condition(tbl_condition.("Group A"));
            tbl_condition.("Group B") = gnames_condition(tbl_condition.("Group B"));

            tbl_epoch = array2table(results_epoch,"VariableNames", ...
                ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
            tbl_epoch.("Group A") = gnames_epoch(tbl_epoch.("Group A"));
            tbl_epoch.("Group B") = gnames_epoch(tbl_epoch.("Group B"));

            % Create stat table
            anovaSTAT_condition = cell2table(sumStats_condition(2:4,:),...
                'VariableNames',sumStats_condition(1,:));

            anovaSTAT_epoch = cell2table(sumStats_epcoh(2:4,:),...
                'VariableNames',sumStats_epcoh(1,:));

            allVariants.(variNUM{vi}).Condition.Posthoc = tbl_condition;
            allVariants.(variNUM{vi}).Condition.Anova = anovaSTAT_condition;

            allVariants.(variNUM{vi}).Epoch.Posthoc = tbl_epoch;
            allVariants.(variNUM{vi}).Epoch.Anova = anovaSTAT_epoch;

        end

    case 'recog'


        for vi = 1:length(variNUM)

            ydata = allVariants.(variNUM{vi}).StatTable.PupilSize;
            condition2u = allVariants.(variNUM{vi}).StatTable.GroupID;
            epoch2u = allVariants.(variNUM{vi}).StatTable.EpochID;

            [~,sumStats_condition,stats_condition] = anova1(ydata,condition2u,"off");
            [~,sumStats_epcoh,stats_epoch] = anova1(ydata,epoch2u,"off");

            [results_condition,~,~,gnames_condition] = multcompare(stats_condition,"Display","off");
            [results_epoch,~,~,gnames_epoch] = multcompare(stats_epoch,"Display","off");

            % Create Table
            tbl_condition = array2table(results_condition,"VariableNames", ...
                ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
            tbl_condition.("Group A") = gnames_condition(tbl_condition.("Group A"));
            tbl_condition.("Group B") = gnames_condition(tbl_condition.("Group B"));

            tbl_epoch = array2table(results_epoch,"VariableNames", ...
                ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
            tbl_epoch.("Group A") = gnames_epoch(tbl_epoch.("Group A"));
            tbl_epoch.("Group B") = gnames_epoch(tbl_epoch.("Group B"));

            % Create stat table
            anovaSTAT_condition = cell2table(sumStats_condition(2:4,:),...
                'VariableNames',sumStats_condition(1,:));

            anovaSTAT_epoch = cell2table(sumStats_epcoh(2:4,:),...
                'VariableNames',sumStats_epcoh(1,:));

            allVariants.(variNUM{vi}).Condition.Posthoc = tbl_condition;
            allVariants.(variNUM{vi}).Condition.Anova = anovaSTAT_condition;

            allVariants.(variNUM{vi}).Epoch.Posthoc = tbl_epoch;
            allVariants.(variNUM{vi}).Epoch.Anova = anovaSTAT_epoch;

        end




end








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

% CREATE Matrix average baseline
% -400 -- stimOn
% stimOn -- stimOff
% stimOff -- Question
% Delay -- +200
epochIdenti = {'M400_StiON','StiON_StiOFF','StimOFF_Q','End_P200'};


% Expanded stats table for multicompare
% ---- C1: Group ID ---- C2: Epoch ID: ---- C3: Mean Pupil Size

totalTrialn = 0;
for trCi = 1:width(statsOutput.NormMat)
    tmpTrcI = height(statsOutput.NormMat{trCi}{1});
    totalTrialn = totalTrialn + tmpTrcI;
end


allgroup = cell(totalTrialn*4,1);
allepoch = cell(totalTrialn*4,1);
allpupilM = nan(totalTrialn*4,1);
tCounti = 1;
for gi = 1:width(statsOutput.NormMat)

    tmpGroupI = statsOutput.NormMat{gi};
    for ei = 1:width(tmpGroupI)

        tmpEpochIgi = tmpGroupI{ei};

        for ti = 1:height(tmpEpochIgi)
            allgroup{tCounti} = trialGroupIDS{gi};
            allepoch{tCounti} = epochIdenti{ei};

            tmpTRiali = tmpEpochIgi(ti,:);
            allpupilM(tCounti) = mean(tmpTRiali,'omitnan');
            tCounti = tCounti + 1;
        end
    end
end

statTableFin = table(allgroup , allepoch , allpupilM, 'VariableNames',...
    {'GroupID','EpochID','PupilSize'});


statsOutput.GroupIDS = trialGroupIDS;
statsOutput.VarID = variantID;
statsOutput.SummaryTab = summTable;
statsOutput.StatTable = statTableFin;

end

