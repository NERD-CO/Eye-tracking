%%% GT: 1=old, 0=new | Subj: 1=old, 0=new


%%%%%%%input arg as a patient instead of looping through all the patients, and add plotON 

% Loop through each variant in variantS for confidence rating plots

% % In the cleaned_eyeDATA folder...
% cleanedDataFolder = struct2table(dir('cleaned_eyeDATA'));
% 
%    % For each file...
%    tableVals = cleanedDataFolder(cleanedDataFolder.isdir < 1,:);
%    fileNames_temp = tableVals(:,1);
%    fileNames = char(fileNames_temp{:,1});
%    for iFile = 1:height(fileNames)
% 
%     % load variantS
%     ptFile = fileNames(iFile,:);

    cd(cleanedDataLOC);
    ptFile = 'cl_eyeData_AMC_PY21NO05.mat';
    
     load(ptFile, 'variantS');              %%%% ptFile = input

    % determine # of variants and variant names in variantS
    varNames = char(fieldnames(variantS));
    totVars = height(varNames);
    
        % For each variant stored within variantS...
        for iVars = 1:totVars
            
            currentVar = cellstr(convertCharsToStrings(varNames(1,:)));
            
            % Name variables according to variant 
            eyeLold = variantS.(currentVar{1}).eye1O;
            eyeRold = variantS.(currentVar{1}).eye2O;
            eyeLnew = variantS.(currentVar{1}).eye1N;
            eyeRnew = variantS.(currentVar{1}).eye2N;
                      
            % Extract raw pupil data
            eyeLold_rawPupil = eyeLold.oT_pupilS_raw;
            eyeRold_rawPupil = eyeRold.oT_pupilS_raw;
            eyeLnew_rawPupil = eyeLnew.oT_pupilS_raw;
            eyeRnew_rawPupil = eyeRnew.oT_pupilS_raw;
    
            % Calculate ground truth new and old mean/SD for each eye
                % L old
                    % Matrix to hold individual lines - 50 x 1005
                    oldLTest = zeros(height(eyeLold_rawPupil),length(eyeLold_rawPupil{1}));
                    
                    for i = 1:height(eyeLold_rawPupil)
                         eyeLold_gTruth_rawPupil2 = eyeLold_rawPupil{i,1};
                         oldLTest(i,:) = eyeLold_gTruth_rawPupil2;
                    end
                
                    % Mean & SD of pupil diameter across trials (like ERP grand average)
                    meanLold = mean(oldLTest,"omitnan"); 
                    stdLold = std(oldLTest,"omitnan");                      

                        % Std for plotting
                        upperStd_Lold = meanLold + (stdLold*1.5);
                        lowerStd_Lold = meanLold - (stdLold*1.5);

                 % R old
                    % Matrix to hold individual lines - 50 x 1005
                    oldRTest = zeros(height(eyeRold_rawPupil),length(eyeRold_rawPupil{1}));
                    
                    for i = 1:height(eyeRold_rawPupil)
                         eyeRold_gTruth_rawPupil2 = eyeRold_rawPupil{i,1};
                         oldRTest(i,:) = eyeRold_gTruth_rawPupil2;
                    end
                
                    % Mean & SD of pupil diameter across trials (like ERP grand average)
                    meanRold = mean(oldRTest,"omitnan"); 
                    stdRold = std(oldRTest,"omitnan");     

                        % Std for plotting
                        upperStd_Rold = meanRold + (stdRold*1.5);
                        lowerStd_Rold = meanRold - (stdRold*1.5);

                % L new
                    % Matrix to hold individual lines - 50 x 1005
                    newLTest = zeros(height(eyeLnew_rawPupil),length(eyeLnew_rawPupil{1}));
                    
                    for i = 1:height(eyeLnew_rawPupil)
                         eyeLnew_gTruth_rawPupil2 = eyeLnew_rawPupil{i,1};
                         newLTest(i,:) = eyeLnew_gTruth_rawPupil2;
                    end
                
                    % Mean & SD of pupil diameter across trials (like ERP grand average)
                    meanLnew = mean(newLTest,"omitnan");
                    stdLnew = std(newLTest,"omitnan");
                       
                        % Std for plotting
                        upperStd_Lnew = meanLnew + (stdLnew*1.5);
                        lowerStd_Lnew = meanLnew - (stdLnew*1.5);

                 % R new
                    % Matrix to hold individual lines - 50 x 1005
                    newRTest = zeros(height(eyeRnew_rawPupil),length(eyeRnew_rawPupil{1}));
                    
                    for i = 1:height(eyeRnew_rawPupil)
                         eyeRnew_gTruth_rawPupil2 = eyeRnew_rawPupil{i,1};
                         newRTest(i,:) = eyeRnew_gTruth_rawPupil2;
                    end
                
                    % Mean & SD of pupil diameter across trials (like ERP grand average)
                    meanRnew = mean(newRTest,"omitnan"); 
                    stdRnew = std(newRTest,"omitnan");
                       
                        % Std for plotting
                        upperStd_Rnew = meanRnew + (stdRnew*1.5);
                        lowerStd_Rnew = meanRnew - (stdRnew*1.5);
               
             % Pupil info for all 100 trials - Combine eye1N and eye1O to get 100
                eyeLall = vertcat(eyeLold, eyeLnew);
                eyeRall = vertcat(eyeRold, eyeRnew);
            
             % Tables separated by subjective "old" and "new" ratings                 
                idx_eyeL_subjOld = eyeLall.confRatings > 0;
                eyeL_subjOld = eyeLall(idx_eyeL_subjOld, :);
        
                idx_eyeL_subjNew = eyeLall.confRatings < 1;
                eyeL_subjNew = eyeLall(idx_eyeL_subjNew, :);
        
                idx_eyeR_subjOld = eyeRall.confRatings > 0;
                eyeR_subjOld = eyeRall(idx_eyeR_subjOld, :);
        
                idx_eyeR_subjNew = eyeRall.confRatings < 1;
                eyeR_subjNew = eyeRall(idx_eyeR_subjNew, :);
    
             % Extract raw pupil data from each of the 4 tables
                eyeLoldS_rawPupil = eyeL_subjOld.oT_pupilS_raw;
                eyeRoldS_rawPupil = eyeR_subjOld.oT_pupilS_raw;
                eyeLnewS_rawPupil = eyeL_subjNew.oT_pupilS_raw;
                eyeRnewS_rawPupil = eyeR_subjNew.oT_pupilS_raw;   
             
             % Calculate subj ratings new and old mean/SD for each eye
                 % L old
                    % Matrix to hold individual lines - 50 x 1005
                    oldLTestS = zeros(height(eyeLoldS_rawPupil),length(eyeLold_rawPupil{1}));
                    
                    for i = 1:height(eyeLoldS_rawPupil)
                         eyeLold_subj_rawPupil2 = eyeLoldS_rawPupil{i,1};
                         oldLTestS(i,:) = eyeLold_subj_rawPupil2;
                    end
                
                    % Mean & SD of pupil diameter across trials (like ERP grand average)
                    meanLoldS = mean(oldLTestS,"omitnan"); 
                    stdLoldS = std(oldLTestS,"omitnan");                      

                        % Std for plotting
                        upperStd_LoldS = meanLoldS + (stdLoldS*1.5);
                        lowerStd_LoldS = meanLoldS - (stdLoldS*1.5);  
                 % R old
                    % Matrix to hold individual lines - 50 x 1005
                    oldRTestS = zeros(height(eyeRoldS_rawPupil),length(eyeRold_rawPupil{1}));
                    
                    for i = 1:height(eyeRoldS_rawPupil)
                         eyeRold_subj_rawPupil2 = eyeRoldS_rawPupil{i,1};
                         oldRTestS(i,:) = eyeRold_subj_rawPupil2;
                    end
                
                    % Mean & SD of pupil diameter across trials (like ERP grand average)
                    meanRoldS = mean(oldRTestS,"omitnan"); 
                    stdRoldS = std(oldRTestS,"omitnan");                      

                        % Std for plotting
                        upperStd_RoldS = meanRoldS + (stdRoldS*1.5);
                        lowerStd_RoldS = meanRoldS - (stdRoldS*1.5);    
                 % L new
                    % Matrix to hold individual lines - 50 x 1005
                    newLTestS = zeros(height(eyeLnewS_rawPupil),length(eyeLnew_rawPupil{1}));
                    
                    for i = 1:height(eyeLnewS_rawPupil)
                         eyeLnew_subj_rawPupil2 = eyeLnewS_rawPupil{i,1};
                         newLTestS(i,:) = eyeLnew_subj_rawPupil2;
                    end
                
                    % Mean & SD of pupil diameter across trials (like ERP grand average)
                    meanLnewS = mean(newLTestS,"omitnan"); 
                    stdLnewS = std(newLTestS,"omitnan");                      

                        % Std for plotting
                        upperStd_LnewS = meanLnewS + (stdLnewS*1.5);
                        lowerStd_LnewS = meanLnewS - (stdLnewS*1.5);  
                 % R new
                    % Matrix to hold individual lines - 50 x 1005
                    newRTestS = zeros(height(eyeRnewS_rawPupil),length(eyeRnew_rawPupil{1}));
                    
                    for i = 1:height(eyeRnewS_rawPupil)
                         eyeRnew_subj_rawPupil2 = eyeRnewS_rawPupil{i,1};
                         newRTestS(i,:) = eyeRnew_subj_rawPupil2;
                    end
                
                    % Mean & SD of pupil diameter across trials (like ERP grand average)
                    meanRnewS = mean(newRTestS,"omitnan"); 
                    stdRnewS = std(newRTestS,"omitnan");                      

                        % Std for plotting
                        upperStd_RnewS = meanRnewS + (stdRnewS*1.5);
                        lowerStd_RnewS = meanRnewS - (stdRnewS*1.5); 

             % Normalize means
                % Combine new and old means first for each eye
                meanRtotGT = [meanRnew meanRold];
                meanLtotGT = [meanLnew meanLold];
                meanRtotS = [meanRnewS meanRoldS];
                meanLtotS = [meanLnewS meanLoldS];

                % Normalize
                meanRtotGT_N = normalize(meanRtotGT, 'range');
                meanLtotGT_N = normalize(meanLtotGT, 'range'); 
                meanRtotS_N = normalize(meanRtotS, 'range');
                meanLtotS_N = normalize(meanLtotS, 'range');

                % Split back into new and old
                    % Ground truth new and old for L and R eye
                    meanRnewGT_N = meanRtotGT_N(1:length(meanRnew));
                    meanRoldGT_N = meanRtotGT_N((length(meanRnew)+1):length(meanRtotGT_N));
                    meanLnewGT_N = meanLtotGT_N(1:length(meanLnew));
                    meanLoldGT_N = meanLtotGT_N((length(meanLnew)+1):length(meanLtotGT_N));  
                    
                    % Subjective rating new and old for L and R eye
                    meanRnewS_N = meanRtotS_N(1:length(meanRnew));
                    meanRoldS_N = meanRtotS_N((length(meanRnew)+1):length(meanRtotS_N));
                    meanLnewS_N = meanLtotS_N(1:length(meanLnew));
                    meanLoldS_N = meanLtotS_N((length(meanLnew)+1):length(meanLtotS_N));



        end
       

   %end        %goes with 'for iFile...'
   

%% Tiled plots - L old/new x subj/ground truth = 4 lines

% L eye - ground truth vs subjective
figure;
fig = tiledlayout(2,1);
l
tile1 = nexttile; % L eye
plot(meanLnewGT_N)
hold on
plot(meanLoldGT_N)
plot(meanLnewS_N)
plot(meanLoldS_N)
title(tile1,'L eye new and old subj vs. ground truth (GT)')
legend(tile1, 'GT new', 'GT old', 'Subj new', 'Subj old')
%legend(tile1, 'GT new', 'GT old', 'Subj old', 'Subj new')
xlim([1 1000]); % x axis
xticks([1 250 500 750 1000]); % which plot ticks to keep visible

% R eye - ground truth vs subjective
tile2 = nexttile; % R eye
plot(meanRnewGT_N)
hold on
plot(meanRoldGT_N)
plot(meanRnewS_N)
plot(meanRoldS_N)
title(tile2, 'R eye new and old subj vs. ground truth (GT)')
legend(tile2, 'GT new', 'GT old', 'Subj new', 'Subj old')
%legend(tile2, 'GT new', 'GT old','Subj old','Subj new')
xlim([1 1000]); % x axis
xticks([1 250 500 750 1000]); % which plot ticks to keep visible
