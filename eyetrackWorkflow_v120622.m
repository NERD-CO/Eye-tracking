%% Loc files/f(x)'s for preprocessing


%% Preprocessing

% Convert EDF to MAT - Edf2Mat_UCH function
[eyeProcName] = Edf2Mat_UCH(edfFile, patientID, block, variant, fileDIR, saveDIR);
    % Example function call: Edf2Mat_UCH('NO20221615110.edf', 'MW9', 'L' , '1')

% John's f(x) for gaze, saccades, etc tables - ExtractEyeInfo_v2 function
[tsTable, picTable, fixTab, saccTab, rawTab] = ExtractEyeInfo_v2(eyeMatfile);
    % Input: eyeMatfile generated output from 'Edf2Mat_UCH.m' function

%%    
% Clean up data
    % 1. subf(x) to look at 1 eye in 1 session at a time
        %
    
    
    % plot lines to visualize data issues
    mainPath = 'C:\Users\darwinm\Documents\MATLAB\EyeTrack\eyeTrack\eyeDATA';
    tempFile_name = 'eyeData_AMC_PY22NO16.mat';
    load(tempFile_name, 'variantS');
    tempEye = variantS.var1.eye2O;
    
    % chop data into trials and make sure all are same length - trim back
    % to make all same size/length. find shortest trial time and make all rest the same
    cleanPupil_size = tempEye.oT_pupilS_mean < 120;
    tempEye_clean = tempEye(~cleanPupil_size,:);
    
    
    figure;
    % loop through raw pupil size data in tempEye and plot
    for i = 1:height(tempEye_clean)
    pupilRaw = tempEye_clean.oT_pupilS_raw{i};   
        
    plot(pupilRaw, 'k')
    %pause
    hold on
    end
    
    %column 6 - in each individual value, where does it drop below 120,
    %replace with NaN. then, make all same length by trimming off the back 
    %and make it all into a
    %f(x). input = file name, var #, eye selection
    
    %%% Replace values under 120 with NaN
    
    %A(A==yourvalue)=NewValue;
    
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
    
    figure;
    % loop through raw pupil size data in tempEye and plot
    for i = 1:height(tempEye_clean)
    pupilRaw = tempEye_clean.oT_pupilS_raw{i};   
        
    plot(pupilRaw, 'k')
    %pause
    hold on
    end  
    
 
    
    