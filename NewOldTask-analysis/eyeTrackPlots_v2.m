function [] = eyeTrackPlots_v2(cleanedDataLOC, exclLOC, ptID , plotNUM)


cd(exclLOC)
varTable = readtable('variantLIST.xlsx');
ptTABLE = varTable(matches(varTable.Subject,ptID),:);

cd(cleanedDataLOC)
ptFname = ['cl_eyeData_',ptID,'.mat'];
load(ptFname, 'variantS')

close all
switch plotNUM
    case 1 % learning categories

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

            learnTTL = varDATA.LearnTTLplus;

            learnCAT = varDATA.Learn_CatID;

            catIDs = unique(learnCAT);

            catColorS = [239 71 111 ;...
                255 209 102;...
                6 214 160; ...
                17 138 178;...
                7 59 76];

            catColorSrgb = catColorS/255;

            figure;
            tiledlayout(1,3)

            nexttile

            for ci = 1:length(catIDs)

                catROWSi = matches(learnCAT,catIDs{ci});
                % check for single nans
                catROWSi2 = cellfun(@(x) numel(x) , learnDATA, 'UniformOutput', true) == 1;

                catROWSi3 = (catROWSi & ~catROWSi2);

                tmpCatRaw = learnDATA(catROWSi3);
                tmpCatTTL = learnTTL(catROWSi3);

                % sampleLEN_L = 1402; % % 1000 ms + 400 ms + 2 ms

                sampleLEN_L = 1202; % % 800 ms + 300 ms + 2 ms


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
            % legend(catIDs)
            xlim([1 sampleLEN_L + 1])

            xticks([1 251 451 651 851 1052])
            xticklabels([-250 0 200 400 600 800])
            xline(251,'-','Stimulus on','LabelVerticalAlignment','bottom')
            xlabel('Time in ms')
            ylabel('Percent change from baseline')

            title([variNUM{vi} , ' Learn - Timeline'])


            nexttile

            for ci = 1:length(catIDs)

                catROWSi = ~matches(learnCAT,catIDs{ci});
                % check for single nans
                catROWSi2 = cellfun(@(x) numel(x) , learnDATA, 'UniformOutput', true) == 1;

                catROWSi3 = (catROWSi & ~catROWSi2);

                tmpCatRawF = learnDATA(catROWSi3);

                tmpCatTTL = learnTTL(catROWSi3);

                % CREATE Matrix average baseline
                tmpCatMatBase = nan(height(tmpCatRawF),801);
                for ti = 1:height(tmpCatTTL)
                    ttlStartB = tmpCatTTL{ti}.ELNKint(1);
                    ttlStopB = tmpCatTTL{ti}.ELNKint(2);
                    baseLINE =  tmpCatRawF{ti}(ttlStartB:ttlStopB);
                    meanBASE = mean(baseLINE,'omitnan');

                    ttlStartS = tmpCatTTL{ti}.ELNKint(2);
                    ttlStopS = ttlStartS + 800;
                    stIMUlus =  tmpCatRawF{ti}(ttlStartS:ttlStopS);

                    tmpCatMatBase(ti,:) = (((stIMUlus) / meanBASE) * 100) - 100;
                end

                tmpAllRaw = reshape(tmpCatMatBase,1,numel(tmpCatMatBase));
                tmpAllRnoN = transpose(tmpAllRaw(~isnan(tmpAllRaw)));

                pd_eye = fitdist(tmpAllRnoN,'Kernel','Kernel','epanechnikov');
                % x_values = 50:1:250;
                eyeX_values = min(tmpAllRnoN):2:max(tmpAllRnoN);
                eyeY_values = pdf(pd_eye,eyeX_values);

                hold on
                plot(eyeX_values,eyeY_values, 'Color', catColorSrgb(ci,:), 'LineWidth',3)


            end
            % legend(catIDs)
            xlim([min(tmpAllRnoN) max(tmpAllRnoN)])
            title([variNUM{vi} , ' Learn - Density'])
            ylabel('Prob. density')
            xlabel('Percent change from baseline')

            nexttile

            bxData = zeros(801*100,1);
            grData = cell(801*100,1);
            startIND = 1;
            stopIND = 801;
            for ti = 1:height(learnDATA)

                tmpCatRaw = learnDATA{ti};

                if numel(tmpCatRaw) == 1
                    continue
                end

                ttlStartB = learnTTL{ti}.ELNKint(1);
                ttlStopB = learnTTL{ti}.ELNKint(2);
                baseLINE =  tmpCatRaw(ttlStartB:ttlStopB);
                meanBASE = mean(baseLINE,'omitnan');

                ttlStartS = learnTTL{ti}.ELNKint(2);
                ttlStopS = ttlStartS + 800;
                stIMUlus =  tmpCatRaw(ttlStartS:ttlStopS);

                tmpCatMatBase = (((stIMUlus) / meanBASE) * 100) - 100;

                bxData(startIND:stopIND,1) = tmpCatMatBase;
                grData(startIND:stopIND,1) = repmat(learnCAT(ti),801,1);

                startIND = stopIND + 1;
                stopIND = startIND + 801 - 1;

            end

            % REMOVE empty brackets
            keepIND = cellfun(@(x) ~isempty(x), grData, 'UniformOutput',true);
            bxDataN = bxData(keepIND);
            grDataN = grData(keepIND);

            nanIndex = ~isnan(bxDataN);
            bxDataf = bxDataN(nanIndex);
            grDataf = grDataN(nanIndex);

            grDataf = categorical(grDataf);
            b = boxchart(double(bxDataf),'GroupByColor',grDataf);
            for bbi = 1:5
                b(bbi).BoxFaceColor = catColorSrgb(bbi,:);
                b(bbi).BoxFaceColor = catColorSrgb(bbi,:);
                b(bbi).BoxMedianLineColor = catColorSrgb(bbi,:);
                b(bbi).MarkerColor = catColorSrgb(bbi,:);
            end
            legend('Location','northeast')
            title([variNUM{vi} , ' Learn - Boxplot'])
            ylabel('Percent change from baseline')
            xticks([])




        end

    case 2

        variNUM = fieldnames(variantS);
        for vi = 1:length(variNUM)

            varPtInfoLearn = ptTABLE(matches(ptTABLE.Block,'recog') &...
                ismember(ptTABLE.Variant,str2double(variNUM{vi}(end))),:);

            varDATA = variantS.(variNUM{vi}).dataTable;

            if matches(varPtInfoLearn.Eye2Use{1},'Left')
                recogDATA = varDATA.Left_L_oT_pupilS_rawCL;
            else
                recogDATA = varDATA.Right_L_oT_pupilS_rawCL;
            end

            recogTTL = varDATA.RecogTTLplus;

            recogNvO = varDATA.groundTruth;

            % catIDs = unique(learnCAT);

            oVnColorS = [120, 0, 0 ;...
                         102, 155, 188];

            catColorSrgb = oVnColorS/255;

            sampleLEN_R = 1402;

            figure;
            tiledlayout(1,2)

            nexttile

            for ci = 1:2

                switch ci
                    case 1
                        nvoROWSi = recogNvO == 0;
                    case 2
                        nvoROWSi = recogNvO == 1;
                end
                % check for single nans
                nvoROWSi2 = cellfun(@(x) numel(x) , recogDATA, 'UniformOutput', true) == 1;

                nvoROWSi3 = (nvoROWSi & ~nvoROWSi2);

                tmpNvORaw = recogDATA(nvoROWSi3);
                tmpNvOTTL = recogTTL(nvoROWSi3);

                % CREATE Matrix average baseline
                tmpNvOMatBase = nan(height(tmpNvORaw),sampleLEN_R); % FIGURE OUT NEW LENGTH
                for ti = 1:height(tmpNvOTTL)
                    ttlStartB = tmpNvOTTL{ti}.ELNKint(1);
                    ttlStopB = tmpNvOTTL{ti}.ELNKint(2);
                    baseLINE =  tmpNvORaw{ti}(ttlStartB:ttlStopB);
                    meanBASE = mean(baseLINE,'omitnan');

                    ttlStartS = tmpNvOTTL{ti}.ELNKint(2);
                    ttlStopS = ttlStartS + 1000;
                    stIMUlus =  tmpNvORaw{ti}(ttlStartS:ttlStopS);

                    tmpNvOMatBase(ti,:) = ((([baseLINE , stIMUlus]) / meanBASE) * 100)-100;
                    plot((([baseLINE , stIMUlus] / meanBASE)*100)-100,...
                        'Color',[catColorSrgb(ci,:) 0.2],'LineWidth',0.5); 
                    xline(401)
                    hold on

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

                bloCKs = ceil(sampleLEN_R/45);
                blstart = 1:45:bloCKs*45;
                blstop = [blstart(2:end) - 1 , sampleLEN_R];

                blMEANS = zeros(1,bloCKs);
                blUstds = zeros(1,bloCKs);
                blDstds = zeros(1,bloCKs);

                for bil = 1:bloCKs
                    tCMBcol = tmpNvOMatBase(:,blstart(bil):blstop(bil));

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
                line([xTICKS ; xTICKS] , [blUstds ; blDstds],'Color',[catColorSrgb(ci,:) 0.3],...
                    'LineWidth',0.75)

            end
            % legend(catIDs)
            xlim([1 sampleLEN_R])

            xticks([1 401 601 801 1001 1201 1402])
            xticklabels([-400 0 200 400 600 800 1000])
            xline(401,'-','Stimulus on','LabelVerticalAlignment','bottom')
            xlabel('Time in ms')
            ylabel('Percent change from baseline')

            title([variNUM{vi} , ' Recog - Timeline'])
            legend('Old','','','','','','','','','','','','','','','',...
                 '','','','','','','','','','New')


            nexttile
            
            statData = struct;
            for ci = 1:2

                switch ci
                    case 1
                        nvoROWSi = recogNvO == 0;
                    case 2
                        nvoROWSi = recogNvO == 1;
                end
                % check for single nans
                nvoROWSi2 = cellfun(@(x) numel(x) , recogDATA, 'UniformOutput', true) == 1;
                nvoROWSi3 = (nvoROWSi & ~nvoROWSi2);
                tmpNvORaw = recogDATA(nvoROWSi3);
                tmpNvOTTL = recogTTL(nvoROWSi3);

                % CREATE Matrix average baseline
                tmpNvOMatBase = nan(height(tmpNvORaw),sampleLEN_R); % FIGURE OUT NEW LENGTH
                tmpNvOonlyStim = nan(height(tmpNvORaw),1001);
                tmpNvOonlyStimMtrial = nan(height(tmpNvORaw),1);
                for ti = 1:height(tmpNvOTTL)
                    ttlStartB = tmpNvOTTL{ti}.ELNKint(1);
                    ttlStopB = tmpNvOTTL{ti}.ELNKint(2);
                    baseLINE =  tmpNvORaw{ti}(ttlStartB:ttlStopB);
                    meanBASE = mean(baseLINE,'omitnan');

                    ttlStartS = tmpNvOTTL{ti}.ELNKint(2);
                    ttlStopS = ttlStartS + 1000;
                    stIMUlus =  tmpNvORaw{ti}(ttlStartS:ttlStopS);

                    tmpNvOMatBase(ti,:) = ((([baseLINE , stIMUlus]) / meanBASE) * 100)-100;

                    tmpNvOonlyStim(ti,:) = tmpNvOMatBase(ti,ttlStartS:ttlStartS+1000);

                    tmpNvOonlyStimMtrial(ti) = mean(tmpNvOonlyStim(ti,:),'omitnan');
                end

                tmpAllRaw = reshape(tmpNvOonlyStim,1,numel(tmpNvOonlyStim));
                tmpAllRnoN = transpose(tmpAllRaw(~isnan(tmpAllRaw)));
                % pd_eye = fitdist(tmpAllRnoN,'Kernel','Kernel','epanechnikov');
                % % x_values = 50:1:250;
                % eyeX_values = min(tmpAllRnoN):2:max(tmpAllRnoN);
                % eyeY_values = pdf(pd_eye,eyeX_values);
                % meanVALUE = median(tmpAllRnoN, 'omitnan');
                hold on
                % plot(eyeX_values,eyeY_values, 'Color', catColorSrgb(ci,:), 'LineWidth',3)
                % switch ci
                %     case 1
                %         meanLtxt = 'Median OLD';
                %         horzinAl = 'right';
                %     case 2
                %         meanLtxt = 'Median NEW';
                %         horzinAl = 'left';
                % end
                % xline(meanVALUE , '-', meanLtxt, 'Color', catColorSrgb(ci,:),...
                %     'LabelVerticalAlignment','bottom','LabelHorizontalAlignment',...
                %     horzinAl)

                switch ci
                    case 1
                        statData.all.olddata = tmpAllRnoN;
                        statData.sum.olddata = tmpNvOonlyStimMtrial;
                    case 2
                        statData.all.newdata = tmpAllRnoN;
                        statData.sum.newdata = tmpNvOonlyStimMtrial;
                end

                b = boxchart(zeros(numel(tmpAllRnoN),1) + ci,tmpAllRnoN);
                b.BoxEdgeColor = catColorSrgb(ci,:);
                b.BoxFaceColor = catColorSrgb(ci,:);
                b.BoxMedianLineColor = catColorSrgb(ci,:);
                b.MarkerColor = catColorSrgb(ci,:);

                if ci == 1
                    yHeight = abs(min(b.YData) - max(b.YData))*0.05;
                    yTxtH = max(b.YData) - yHeight;
                end

                % % plot(eyeX_values,eyeY_values, 'Color', catColorSrgb(ci,:), 'LineWidth',3)
                % s.MarkerFaceAlpha = 0.2;
                % s.MarkerEdgeAlpha = 0.2;

                swarmchart(zeros(numel(tmpNvOonlyStimMtrial),1) + ci,...
                    tmpNvOonlyStimMtrial, 40, catColorSrgb(ci,:));

            end
            % legend(catIDs)
            % xlim([min(tmpAllRnoN) max(tmpAllRnoN)])
            title([variNUM{vi} , ' Recog - Boxchart'])
            ylabel('Percent change from baseline')
            xticks([1 2])
            xticklabels({'old', 'new'})


            % [~,pValue,~,statINFO] = ttest2(statData.sum.olddata , statData.sum.newdata);
            [pValue,~,statINFO] = ranksum(statData.sum.olddata , statData.sum.newdata);
            text(1.2, yTxtH, ['p = ' num2str(round(pValue,3))]);




        end

end