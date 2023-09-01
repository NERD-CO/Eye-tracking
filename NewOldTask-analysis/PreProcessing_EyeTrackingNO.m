%% SETUP code for NewOld delay preprocessing

% First run NewOldTxttoMat_v2.mat

curPC = getenv('COMPUTERNAME');

switch curPC
    case 'DESKTOP-I5CPDO7'   % JAT WORK desktop


        excelLOC = 'Z:\MW_JAT_Backup\EyeTrack_Manuscript';
        DATADIR = 'Z:\MW_JAT_Backup\EyeTrack_Manuscript\Patient_folders';



    case 'otherwise'



end



%% PROCESS BEHAVIORAL FILES

cd(excelLOC)
% Load in excel file
allSUBtable = readtable('Eyetracking patient summary sheet.xlsx');

for si = 1:height(allSUBtable)
% for si = 1:height(allSUBtable)

    subROW = allSUBtable(si,:);
    subID = subROW.patientID{1};

    for vi = 1:3
        tmpVar = ['variant',num2str(vi)];
        if subROW.(tmpVar)

            blmps = {'learn','recog'};
            for bi = 1:2

                NewOldTxttoMat_v3(excelLOC , DATADIR , subID , vi , blmps{bi});

            end
        else
            continue
        end
    end
end



%% PROCESS TTL FILES





