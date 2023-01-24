% STEP 1: Put EDF files into a folder with patient ID 
    % Path example: C:\Users\darwinm\Documents\MATLAB\EyeTrack\eyeTrack\AMC_PY21NO05

% STEP 2: Update variantLIST excel sheet with new patient
    % Path example: C:\Users\darwinm\Documents\MATLAB\EyeTrack\variantLIST

% STEP 3: Run Initial_EyeAnalysis_MLD_v2.m f(x)

userPC = 'JAT';
switch userPC
    case 'JAT'
        excelLOC = 'C:\Users\Admin\Documents\Github\Eye-tracking';
        mainLOC = 'E:\Dropbox\SfN_2022\dataOut_AMC\dataOut_AMC\eyeTrack';
        saveLOC = 'E:\Dropbox\SfN_2022\dataOut_AMC\dataOut_AMC\eyeDATA';
    case 'MD'
%         excelLOC = 'C:\Users\Admin\Documents\Github\Eye-tracking';
%         mainLOC = 'E:\Dropbox\SfN_2022\dataOut_AMC\dataOut_AMC\eyeTrack';
end

Initial_EyeAnalysis_MLD_v2(excelLOC , mainLOC , saveLOC)


%%
% STEP 4: Run eyeTRACKproc.m f(x)  

mainPath = saveLOC;
eyeTRACKproc(mainPath);



    
   

