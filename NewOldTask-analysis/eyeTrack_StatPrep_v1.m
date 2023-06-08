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
        end % VARIANT LOOP

    case 'recog'




end


%%%%% CREATE Variable length table/struct to output trial pupil info based
%%%%% on selections above
%%%%% AM HERE 6/7/2023

for ci = 1:length(catIDs)

    catROWSi = matches(learnCAT,catIDs{ci});
    % check for single nans
    catROWSi2 = cellfun(@(x) numel(x) , learnDATA, 'UniformOutput', true) == 1;

    catROWSi3 = (catROWSi & ~catROWSi2);

    tmpCatRaw = learnDATA(catROWSi3);
    tmpCatTTL = learnTTL(catROWSi3);

    sampleLEN_L = 1402; % % 1000 ms + 400 ms + 2 ms

    % CREATE Matrix average baseline
    tmpCatMatBase = nan(height(tmpCatRaw),sampleLEN_L); % NEED TO CHECK
    for ti = 1:height(tmpCatTTL)
        ttlStartB = tmpCatTTL{ti}.ELNKint(1);
        ttlStopB = tmpCatTTL{ti}.ELNKint(2);
        baseLINE =  tmpCatRaw{ti}(ttlStartB:ttlStopB);
        meanBASE = mean(baseLINE,'omitnan');

        ttlStartS = tmpCatTTL{ti}.ELNKint(2);
        ttlStopS = ttlStartS + 800; % 800 ms from stimulus on
        stIMUlus =  tmpCatRaw{ti}(ttlStartS:ttlStopS);

        tmpCatMatBase(ti,:) = ((([baseLINE , stIMUlus]) / meanBASE) * 100)-100;
    end

    %%%% BASELINE to SOME PERIOD AFTER STIM -
    %%%% USE TTL around STIM ON

    %%%%% NORMALIZE

    % tmpMatRaw = zeros(height(tmpCatRaw),750);
    % % tmpMatRaw = zeros(height(tmpCatRaw),numel(tmpCatRaw{1}));
    % for tmi = 1:height(tmpCatRaw)
    %     % tmpMatRaw(tmi,:) = tmpCatRaw{tmi};
    %     tmpMatRaw(tmi,:) = tmpCatRaw{tmi}(1:750);
    % end

    bloCKs = ceil(sampleLEN_L/45);
    blstart = 1:45:bloCKs*45;
    blstop = [blstart(2:end) - 1 , sampleLEN_L];

    blMEANS = zeros(1,bloCKs);
    blUstds = zeros(1,bloCKs);
    blDstds = zeros(1,bloCKs);

    for bil = 1:bloCKs

        tCMBcol = tmpCatMatBase(:,blstart(bil):blstop(bil));

        blMEANS(bil) = mean(tCMBcol, 'all', 'omitnan');
        tmpBLstd = std(tCMBcol, [], 'all', 'omitnan');
        blUstds(bil) = blMEANS(bil) + (tmpBLstd*0.5);
        blDstds(bil) = blMEANS(bil) - (tmpBLstd*0.5);


    end

    xTICKS = blstop - round(45/2);

    % tmpCatMean = mean(tmpCatMatBase,"omitnan");
    hold on
    plot(xTICKS, blMEANS ,'Marker', 'o','MarkerFaceColor', catColorSrgb(ci,:),...
        'Color',catColorSrgb(ci,:))
    line([xTICKS ; xTICKS] , [blUstds ; blDstds],'Color',catColorSrgb(ci,:),...
        'LineWidth',0.75)


end
















end