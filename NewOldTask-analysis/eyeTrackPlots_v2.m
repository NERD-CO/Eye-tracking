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

                catROWSi = ~matches(learnCAT,catIDs{ci});
                % check for single nans
                catROWSi2 = cellfun(@(x) numel(x) , learnDATA, 'UniformOutput', true) == 1;

                catROWSi3 = ~(catROWSi | catROWSi2);

                tmpCatRaw = learnDATA(catROWSi3);
                tmpCatTTL = learnTTL(catROWSi3);

                % CREATE Matrix average baseline
                tmpCatMatBase = nan(height(tmpCatRaw),1052);
                for ti = 1:height(tmpCatTTL)
                    ttlStartB = tmpCatTTL{ti}.ELNKint(1);
                    ttlStopB = tmpCatTTL{ti}.ELNKint(2);
                    baseLINE =  tmpCatRaw{ti}(ttlStartB:ttlStopB);
                    meanBASE = mean(baseLINE);

                    ttlStartS = tmpCatTTL{ti}.ELNKint(2);
                    ttlStopS = ttlStartS + 800;
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

                bloCKs = ceil(1052/45);
                blstart = 1:45:bloCKs*45;
                blstop = [blstart(2:end) - 1 , 1052];

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
            xlim([1 1053])

            xticks([1 201 401 601 801 1001])
            xticklabels([-200 0 200 400 600 800])
            xline(201,'-','Stimulus on','LabelVerticalAlignment','bottom')
            xlabel('Time in ms')
            ylabel('Percent change from baseline')

            title([variNUM{vi} , ' Learn - Timeline'])


            nexttile

            for ci = 1:length(catIDs)

                catROWSi = ~matches(learnCAT,catIDs{ci});
                % check for single nans
                catROWSi2 = cellfun(@(x) numel(x) , learnDATA, 'UniformOutput', true) == 1;

                catROWSi3 = ~(catROWSi | catROWSi2);

                tmpCatRawF = learnDATA(catROWSi3);

                tmpCatTTL = learnTTL(catROWSi3);

                % CREATE Matrix average baseline
                tmpCatMatBase = nan(height(tmpCatRawF),801);
                for ti = 1:height(tmpCatTTL)
                    ttlStartB = tmpCatTTL{ti}.ELNKint(1);
                    ttlStopB = tmpCatTTL{ti}.ELNKint(2);
                    baseLINE =  tmpCatRawF{ti}(ttlStartB:ttlStopB);
                    meanBASE = mean(baseLINE);

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
                meanBASE = mean(baseLINE);

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

                catROWSi = ~matches(learnCAT,catIDs{ci});
                % check for single nans
                catROWSi2 = cellfun(@(x) numel(x) , learnDATA, 'UniformOutput', true) == 1;

                catROWSi3 = ~(catROWSi | catROWSi2);

                tmpCatRaw = learnDATA(catROWSi3);
                tmpCatTTL = learnTTL(catROWSi3);

                % CREATE Matrix average baseline
                tmpCatMatBase = nan(height(tmpCatRaw),1052);
                for ti = 1:height(tmpCatTTL)
                    ttlStartB = tmpCatTTL{ti}.ELNKint(1);
                    ttlStopB = tmpCatTTL{ti}.ELNKint(2);
                    baseLINE =  tmpCatRaw{ti}(ttlStartB:ttlStopB);
                    meanBASE = mean(baseLINE);

                    ttlStartS = tmpCatTTL{ti}.ELNKint(2);
                    ttlStopS = ttlStartS + 800;
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

                bloCKs = ceil(1052/45);
                blstart = 1:45:bloCKs*45;
                blstop = [blstart(2:end) - 1 , 1052];

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
            xlim([1 1053])

            xticks([1 201 401 601 801 1001])
            xticklabels([-200 0 200 400 600 800])
            xline(201,'-','Stimulus on','LabelVerticalAlignment','bottom')
            xlabel('Time in ms')
            ylabel('Percent change from baseline')

            title([variNUM{vi} , ' Learn - Timeline'])









        end

end