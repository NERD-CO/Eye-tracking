%function [] = eyeTRACKproc(mainPath)

% Set data paths
mainPath = 'C:\Users\darwinm\Documents\MATLAB\EyeTrack\eyeTrack\eyeDATA';
cleanedData_path = [mainPath, '\cleaned_eyeDATA'];


% [For every file in eyeDATA, loop through to run analyses?]
    % save cleaned files in new folder to only have raw data in this folder
    % need to cd to new folder to save
    
    % Get contents of eyeDATA as a variable
    test = dir('*.mat');
    fileList = {test.name};
    
% for loop to loop throguh fileList
for fileS = 1:length(fileList)

    % CD to mainPath for raw data
    cd(mainPath);
    
    tempFile_name = fileList{fileS}; 

% Set tempFile_name
    %tempFile_name = "eyeData_AMC_" + num2str(ptID) + ".mat";
    saveFile_name1 = ['cl_', tempFile_name];

% Load in file
load(tempFile_name, 'variantS'); %%%%%%%%%%%%%%%%%%%%%%%%rename specific to patient?

% Set tempEye- loop through every eye in every variant within variantS
    % Go into variantS and determine # variants there
    varSnum = length(fieldnames(variantS));
    
    % Other option: Extract field names of variantS ahead of time
    varSfieldN = fieldnames(variantS);
        % Extract field names of each variant in variantS
   
    
    for i = 1:varSnum
       % name of variant
       currentVariant = variantS.(varSfieldN{i});
       %currentVariant = fieldnames(variantS(varSnum)); %%variantS.(variable name for field{varSnum})
       %tempEye_variantOnly = char("variantS." + currentVariant); 
       
       % create # field names in var3 and then actual field names
       eyeNum = length(fieldnames(currentVariant));
        % Variable for field names of current variant
        varCurrent_fieldN = fieldnames(currentVariant);
       
       % Loop through each eye in variant
       for eyE = 1:eyeNum
            tempEye = currentVariant.(varCurrent_fieldN{eyE});
          
%%%%%%%%%%%%%% EYE PROCCESSING CODE %%%%%%%%%%%%%%%%%%%%%%
        % chop data into trials and make sure all are same length - trim back
        % to make all same size/length. find shortest trial time and make all rest the same
        cleanPupil_size = tempEye.oT_pupilS_mean < 120;
        tempEye_clean = tempEye(~cleanPupil_size,:);
      
        varEye  = tempEye_clean.oT_pupilS_raw;
        for num = 1:height(varEye)

        Eye_clean = varEye{num,1};

            if sum(Eye_clean <= 120) ~= 0
            eyeNaNindex = find(Eye_clean <= 120);
            minEye_in = min(eyeNaNindex);
            maxEye_in = max(eyeNaNindex);
            eyeNaNindex2 = minEye_in-15:maxEye_in+15;
                if min(eyeNaNindex2) < 1
                    eyeNaNindex2 = 1:maxEye_in+15;
                elseif max(eyeNaNindex2) > length(Eye_clean)
                    eyeNaNindex2 = minEye_in-15:length(Eye_clean);
                end
            Eye_clean(eyeNaNindex2) = nan(length(eyeNaNindex2),1); %problem
            varEye{num,1} = Eye_clean;
            end
        end

        % Input cleaned values in varEye into tempEye_clean.oT_pupilS_raw
        tempEye_clean.oT_pupilS_raw = varEye;

        %%% Trim values to be all the same length

        % Find shortest row in column 6
            min_varEye = min(cellfun(@(x) numel(x), varEye, 'UniformOutput', true));

            % Trim all rows in col 6 to shortest length
            for trim = 1:height(varEye)
                varEye{trim} = varEye{trim}(1:min_varEye);

            end

            % Input cleaned values in varEye into tempEye_clean.oT_pupilS_raw
            tempEye_clean.oT_pupilS_raw = varEye;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Save cleaned data into variantS
            variantS.(varSfieldN{i}).(varCurrent_fieldN{eyE}) = tempEye_clean;
                
       end
       
    end 
            cd(cleanedData_path);
            save(saveFile_name1, 'variantS'); 
end    



%end