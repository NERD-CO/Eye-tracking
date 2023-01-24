% STEP 1: Put EDF files into a folder with patient ID 
    % Path example: C:\Users\darwinm\Documents\MATLAB\EyeTrack\eyeTrack\AMC_PY21NO05

% STEP 2: Update variantLIST excel sheet with new patient
    % Path example: C:\Users\darwinm\Documents\MATLAB\EyeTrack\variantLIST

% STEP 3: Run Initial_EyeAnalysis_MLD_v2.m f(x)
[] = Initial_EyeAnalysis_MLD_v2()
    
% STEP 4: Run eyeTRACKproc.m f(x)  
[] = eyeTRACKproc(mainPath);



    
   

