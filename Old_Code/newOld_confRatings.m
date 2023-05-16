%%%% PUPIL SIZE OF NEW VS OLD CONFIDENCE RATINGS

%%% NEXT: Loop through all pts and variants (but also plot each separately too)

% Recognition stimuli presentation
    % Load ground truth data
    excelLOC = 'C:\Users\darwinm\Documents\MATLAB\EyeTrack';
    cd(excelLOC);
    load('newOld_stimID_all.mat','stimAll');

    % Extract variables for ground truth old vs new- 1=present in learn block
    groundTruthRecog_var1 = stimAll.stimNewOld_var1;
    groundTruthRecog_var2 = stimAll.stimNewOld_var2;
    groundTruthRecog_var3 = stimAll.stimNewOld_var3;

% Patient answers
    % Load recognition file - var 1
    behavDataLOC = 'C:\Users\darwinm\Documents\Thompson Lab\Microwire\PatientData\MW5';
    cd(behavDataLOC);
    load('recog_7.16.21_MW5.mat', 'outData')

    % Extract answers from outData - 31:36 = Confidence (Yes vs. No) response
                % 31,32,33=yes,old & 34,35,36=no,new 
    ttlValues = str2double(outData.taskinformation.TTLvalue);
    confRatings = ttlValues(ttlValues(:,1)>=31&ttlValues(:,1)<=36,:);

    % Logical array of confRatings: 1=yes,old & 0=no,new
    confRatings_logical = logical(confRatings>=31&confRatings<=33);

    % Another logical array eliminating 'not sure' (33,34) %%TO DO

% Can compare ground truth to patient answers
                % 1=yes,old & 0=no,new
    % Combine arrays into a table
    compareAns = table(confRatings_logical, groundTruthRecog_var1, ...
        'VariableNames', {'confRatings', 'groundTruth'});

%%%%% NOTE: add this in before cleaning data- specifically before chopping data

% Pupil size as a function of ground truth new vs. old

    % Load in pupil data struct 
    pupilDataLOC = 'C:\Users\darwinm\Documents\MATLAB\EyeTrack\eyeTrack\eyeDATA\cleaned_eyeDATA';
    cd(pupilDataLOC);
    load('cl_eyeData_AMC_PY21NO05.mat', 'variantS');

             %eye1O and eye1N = 100 trials for recog block

    % Mean and SD of pupil size for R OLD
        eyeRold_raw = variantS.var1.eye2O.oT_pupilS_raw;
        % Matrix to hold individual lines - 50 x 1005
        oldRTest = zeros(height(eyeRold_raw),length(eyeRold_raw{1}));
        
        for i = 1:height(eyeRold_raw)
             eyeRold_raw2 = variantS.var1.eye2O.oT_pupilS_raw{i,1};
             oldRTest(i,:) = eyeRold_raw2;
        end
    
        % Mean & SD of pupil diameter across trials (like ERP grand average)
        meanRold = mean(oldRTest,"omitnan");
        stdRold = std(oldRTest,"omitnan");

    % Mean and SD of pupil size for L OLD
        eyeLold_raw = variantS.var1.eye1O.oT_pupilS_raw;
        % Matrix to hold individual lines - 50 x 1005
        oldLTest = zeros(height(eyeLold_raw),length(eyeLold_raw{1}));
        
        for i = 1:height(eyeLold_raw)
             eyeLold_raw2 = variantS.var1.eye1O.oT_pupilS_raw{i,1};
             oldLTest(i,:) = eyeLold_raw2;
        end
    
        % Mean p& SDupil diameter across trials (like ERP grand average)
        meanLold = mean(oldLTest,"omitnan"); 
        stdLold = std(oldLTest,"omitnan");
   
    % Mean and SD of pupil size for R NEW
        eyeRnew_raw = variantS.var1.eye2N.oT_pupilS_raw;
        % Matrix to hold individual lines - 50 x 1005
        newRTest = zeros(height(eyeRnew_raw),length(eyeRnew_raw{1}));
        
        for i = 1:height(eyeRnew_raw)
             eyeRnew_raw2 = variantS.var1.eye2N.oT_pupilS_raw{i,1};
             newRTest(i,:) = eyeRnew_raw2;
        end
    
        % Mean & SD pupil diameter across trials (like ERP grand average)
        meanRnew = mean(newRTest,"omitnan");
        stdRnew = std(newRTest,"omitnan");

    % Mean and SD of pupil size for L NEW
        eyeLnew_raw = variantS.var1.eye1N.oT_pupilS_raw;
        % Matrix to hold individual lines - 50 x 1005
        newLTest = zeros(height(eyeLnew_raw),length(eyeLnew_raw{1}));
        
        for i = 1:height(eyeLnew_raw)
             eyeLnew_raw2 = variantS.var1.eye1N.oT_pupilS_raw{i,1};
             newLTest(i,:) = eyeLnew_raw2;
        end
    
        % Mean & SD pupil diameter across trials (like ERP grand average)
        meanLnew = mean(newLTest,"omitnan");
        stdLnew = std(newLTest,"omitnan");
        
    % Compare "grand average" for new vs old in a plot
    figure;
    hold on;
    plot(meanLnew, 'r')
    plot(meanRnew, 'm')
    plot(meanLold, 'b')
    plot(meanRold, 'c')

    legend('meanLnew','meanRnew', 'meanLold', 'meanRold');

    % x axis
    xlim([1 1000]);
    
    % which plot ticks to keep visible
    xticks([1 250 500 750 1000]);



% Pupil size as a function of patient ratings new vs. old

    % Get pupil info for all 100 trials
        
        % Combine eye1N and eye1O to get 100
    
        % Code for old vs new 
    
        % Get trial # info?
            
    
        % Map onto conf ratings
    
        % Plot pupil size as a f(x) of old vs new conf ratings


