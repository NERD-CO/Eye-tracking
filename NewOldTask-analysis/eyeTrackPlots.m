function [] = eyeTrackPlots(cleanedDataLOC, exclLOC, ptID , plotNUM)


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

            varDATA = variantS.(variNUM{vi}).dataTable;

            learnLeft = varDATA.Left_L_oT_pupilS_rawCL;
            learnRight = varDATA.Right_L_oT_pupilS_rawCL;

            learnCAT = varDATA.Learn_CatID;

            catIDs = unique(learnCAT);

            catColorS = [239 71 111 ;...
                         255 209 102;...
                         6 214 160; ...
                         17 138 178;...
                         7 59 76];

            catColorSrgb = catColorS/255;

            figure;
            for eyeI = 1:2
                subplot(1,2,eyeI)
                switch eyeI
                    case 1
                        eyeDATA = learnLeft;
                    case 2
                        eyeDATA = learnRight;
                end
                for ci = 1:length(catIDs)

                    catROWSi = matches(learnCAT,catIDs{ci});
                    tmpCatRaw = eyeDATA(catROWSi);
                    tmpMatRaw = zeros(height(tmpCatRaw),750);
                    % tmpMatRaw = zeros(height(tmpCatRaw),numel(tmpCatRaw{1}));
                    for tmi = 1:height(tmpCatRaw)
                        % tmpMatRaw(tmi,:) = tmpCatRaw{tmi};
                        tmpMatRaw(tmi,:) = tmpCatRaw{tmi}(1:750);
                    end

                    tmpCatMean = mean(tmpMatRaw,"omitnan");
                    hold on
                    plot(tmpCatMean , 'Color', catColorSrgb(ci,:))

                end
                legend(catIDs)
                xlim([1 numel(tmpCatRaw{1})])

            end

            figure;
            for eyeI = 1:2
                subplot(1,2,eyeI)
                switch eyeI
                    case 1
                        eyeDATA = learnLeft;
                    case 2
                        eyeDATA = learnRight;
                end
                for ci = 1:length(catIDs)

                    catROWSi = matches(learnCAT,catIDs{ci});
                    tmpCatRaw = eyeDATA(catROWSi);

                    tmpCatRawF = cellfun(@(x) x(1:750), tmpCatRaw, 'UniformOutput',false);

                    tmpAllRaw = reshape([tmpCatRawF{:}],1,numel([tmpCatRawF{:}]));
                    tmpAllRnoN = transpose(tmpAllRaw(~isnan(tmpAllRaw)));
                    
                    pd_eye = fitdist(tmpAllRnoN,'Kernel','Kernel','epanechnikov');
                    % x_values = 50:1:250;
                    eyeX_values = min(tmpAllRnoN):2:max(tmpAllRnoN);
                    eyeY_values = pdf(pd_eye,eyeX_values);
                    
                    hold on
                    plot(eyeX_values,eyeY_values, 'Color', catColorSrgb(ci,:), 'LineWidth',3)
       

                end
                legend(catIDs)
                xlim([200 750])

            end

            figure;
            for eyeI = 1:2
                subplot(1,2,eyeI)
                switch eyeI
                    case 1
                        eyeDATA = learnLeft;
                    case 2
                        eyeDATA = learnRight;
                end

                bxData = zeros(750*100,1);
                grData = cell(750*100,1);
                startIND = 1;
                stopIND = 750;
                for ti = 1:height(eyeDATA)

                    tmpCatRaw = eyeDATA{ti};

                    bxData(startIND:stopIND,1) = tmpCatRaw(1:750);
                    grData(startIND:stopIND,1) = repmat(learnCAT(ti),750,1);

                    startIND = stopIND + 1;
                    stopIND = startIND + 750 - 1; 

                end

                nanIndex = ~isnan(bxData);
                bxDataf = bxData(nanIndex);
                grDataf = grData(nanIndex);
                
                grDataf = categorical(grDataf);
                b = boxchart(double(bxDataf),'GroupByColor',grDataf);
                for bbi = 1:5
                    b(bbi).BoxFaceColor = catColorSrgb(bbi,:);
                    b(bbi).BoxFaceColor = catColorSrgb(bbi,:);
                    b(bbi).BoxMedianLineColor = catColorSrgb(bbi,:);
                    b(bbi).MarkerColor = catColorSrgb(bbi,:);
                end
                legend('Location','northeast')

            end







        end

    case 2











end

end