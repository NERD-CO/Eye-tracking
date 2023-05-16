% Loop through each variant in variantS for confidence rating plots

function fig = pupil_confRatings (cleanedDataLOC, ptID)

% % In the cleaned_eyeDATA folder...
cd(cleanedDataLOC);
ptFile = append('cl_eyeData_', ptID, '.mat');
load(ptFile, 'variantS');

% Determine # of variants and variant names in variantS
varNames = fieldnames(variantS);
totVars = height(varNames);

% For each variant stored within variantS...
meanS = struct;

for iVars = 1:totVars
    
    currentVar = varNames{iVars};
    
    % Name variables according to variant 
    eyeLold = variantS.(currentVar).eye1O;
    eyeRold = variantS.(currentVar).eye2O;
    eyeLnew = variantS.(currentVar).eye1N;
    eyeRnew = variantS.(currentVar).eye2N;

    % Pupil info for all 100 trials - Combine eye1N and eye1O to get 100
    eyeLall = vertcat(eyeLold, eyeLnew);
    eyeRall = vertcat(eyeRold, eyeRnew); 

        % Trim oT_pupilS_raw if not all uniform size

            % Find shortest row in column 6
                varEyeL = eyeLall{:,6};
                min_varEyeL = min(cellfun(@(x) numel(x), varEyeL, 'UniformOutput', true));
                % Trim all rows in col 6 to shortest length
                for trim = 1:height(varEyeL)
                    varEyeL{trim} = varEyeL{trim}(1:min_varEyeL);
                end
                % Input cleaned values back into table
                eyeLall.oT_pupilS_raw = varEyeL;
    
                varEyeR = eyeRall{:,6};
                min_varEyeR = min(cellfun(@(x) numel(x), varEyeR, 'UniformOutput', true));
                % Trim all rows in col 6 to shortest length
                for trim = 1:height(varEyeR)
                    varEyeR{trim} = varEyeR{trim}(1:min_varEyeR);
                end
                % Input cleaned values back into table
                eyeRall.oT_pupilS_raw = varEyeR;
        
                
    % Extract pupil data from eyeLall & eyeRall according to desired bins
    
        % Tables separated by different confidence rating bins
%           confL_1 = eyeLall(eyeLall.New_unsure, :);           
%           confL_2 = eyeLall(eyeLall.New_less_sure, :);
%           confL_3 = eyeLall(eyeLall.New_unsure, :);
%           confL_4 = eyeLall(eyeLall.Old_unsure, :);
%           confL_5 = eyeLall(eyeLall.Old_less_sure, :);
%           confL_6 = eyeLall(eyeLall.Old_sure, :);
          confL_34 = eyeLall(eyeLall.All_unsure, :);
%           confL_12 = eyeLall(eyeLall.New_sureLess_sure, :);
%           confL_56 = eyeLall(eyeLall.Old_sureLess_sure, :);
%           confL_1256 = eyeLall(eyeLall.All_sureLess_sure, :);
          confL_16 = eyeLall(eyeLall.All_sure, :);
        
%           confR_1 = eyeRall(eyeRall.New_unsure, :);  
%           confR_2 = eyeRall(eyeRall.New_less_sure, :);
%           confR_3 = eyeRall(eyeRall.New_unsure, :);
%           confR_4 = eyeRall(eyeRall.Old_unsure, :);
%           confR_5 = eyeRall(eyeRall.Old_less_sure, :);
%           confR_6 = eyeRall(eyeRall.Old_sure, :);
          confR_34 = eyeRall(eyeRall.All_unsure, :);
%           confR_12 = eyeRall(eyeRall.New_sureLess_sure, :);
%           confR_56 = eyeRall(eyeRall.Old_sureLess_sure, :);
%           confR_1256 = eyeRall(eyeRall.All_sureLess_sure, :);
          confR_16 = eyeRall(eyeRall.All_sure, :); 

        % Extract relevant raw pupil data from binned tables
        pupilL_34 = confL_34.oT_pupilS_raw;
        pupilL_16 = confL_16.oT_pupilS_raw;
        pupilR_34 = confR_34.oT_pupilS_raw;
        pupilR_16 = confR_16.oT_pupilS_raw;
          
    % Calculate mean/SD
        % L 34
            % Matrix to hold individual lines - 100 x 1005
            matrixL34 = zeros(height(pupilL_34),length(pupilL_34{1}));
            
            for i = 1:height(pupilL_34)
                 pupilL_34_2 = pupilL_34{i,1};
                 matrixL34(i,:) = pupilL_34_2;
            end
        
            % Mean & SD of pupil diameter across trials (like ERP grand average)
            meanL34 = mean(matrixL34,"omitnan"); 
            stdL34 = std(matrixL34,"omitnan");                      

                % Std for plotting
                upperStd_L34 = meanL34 + (stdL34*1.5);
                lowerStd_L34 = meanL34 - (stdL34*1.5);  
    
        % L 16
            % Matrix to hold individual lines - 100 x 1005
            matrixL16 = zeros(height(pupilL_16),length(pupilL_16{1}));
            
            for i = 1:height(pupilL_16)
                 pupilL_16_2 = pupilL_16{i,1};
                 matrixL16(i,:) = pupilL_16_2;
            end
        
            % Mean & SD of pupil diameter across trials (like ERP grand average)
            meanL16 = mean(matrixL16,"omitnan"); 
            stdL16 = std(matrixL16,"omitnan");                      

                % Std for plotting
                upperStd_L16 = meanL16 + (stdL16*1.5);
                lowerStd_L16 = meanL16 - (stdL16*1.5);

        % R 34
            % Matrix to hold individual lines - 100 x 1005
            matrixR34 = zeros(height(pupilR_34),length(pupilR_34{1}));
            
            for i = 1:height(pupilR_34)
                 pupilR_34_2 = pupilR_34{i,1};
                 matrixR34(i,:) = pupilR_34_2;
            end
        
            % Mean & SD of pupil diameter across trials (like ERP grand average)
            meanR34 = mean(matrixR34,"omitnan"); 
            stdR34 = std(matrixR34,"omitnan");                      

                % Std for plotting
                upperStd_R34 = meanR34 + (stdR34*1.5);
                lowerStd_R34 = meanR34 - (stdR34*1.5);  
    
        % R 16
            % Matrix to hold individual lines - 100 x 1005
            matrixR16 = zeros(height(pupilR_16),length(pupilR_16{1}));
            
            for i = 1:height(pupilR_16)
                 pupilR_16_2 = pupilR_16{i,1};
                 matrixR16(i,:) = pupilR_16_2;
            end
        
            % Mean & SD of pupil diameter across trials (like ERP grand average)
            meanR16 = mean(matrixR16,"omitnan"); 
            stdR16 = std(matrixR16,"omitnan");                      

                % Std for plotting
                upperStd_R16 = meanR16 + (stdR16*1.5);
                lowerStd_R16 = meanR16 - (stdR16*1.5);

     % Normalize means
            meanL34_N = normalize(meanL34, 'range');
            meanL16_N = normalize(meanL16, 'range'); 
            meanR34_N = normalize(meanR34, 'range');
            meanR16_N = normalize(meanR16, 'range');

    % Store in struct meanS
    meanS.(currentVar).meanL34_N = meanL34_N;
    meanS.(currentVar).meanL16_N = meanL16_N;
    meanS.(currentVar).meanR34_N = meanR34_N;
    meanS.(currentVar).meanR16_N = meanR16_N;

% Tiled plots for each patient- Sure vs. unsure confidence for each eye 

figure;
fig = tiledlayout(2,1);
tile1 = nexttile; % L eye
plot(meanL16_N)
hold on
plot(meanL34_N)
%title(tile1,'Pupil Size Sure vs. Unsure Recognition')
legend(tile1, 'L Sure', 'L Unsure')
xlim([1 1000]); % x axis
xticks([1 250 500 750 1000]); % which plot ticks to keep visible

tile2 = nexttile; % R eye
plot(meanR16_N)
hold on
plot(meanR34_N)
%title(tile2,'Pupil Size Sure vs. Unsure Recognition')
legend(tile2, 'R Sure', 'R Unsure')
xlim([1 1000]); % x axis
xticks([1 250 500 750 1000]); % which plot ticks to keep visible

title(fig,{[ptID, ' ', currentVar]; 'Pupil Size Sure vs. Unsure Recognition'},...
    'Interpreter', 'none');  

end

% Tiled plots for all patients together - TBD

    % Extract data from struct
    



end
