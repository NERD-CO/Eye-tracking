function [] = runSETUP_NewOldDelay(stage2run)
% SETUP code for NewOld delay preprocessing
% First run NewOldTxttoMat_v2.mat

curPC = getenv('COMPUTERNAME');

switch curPC
    case 'DESKTOP-I5CPDO7'   % JAT WORK desktop

        excelLOC = 'Z:\MW_JAT_Backup\EyeTrack_Manuscript';
        DATADIR = 'Z:\MW_JAT_Backup\EyeTrack_Manuscript\Patient_folders';

        matNWBloc = 'C:\Users\Admin\Documents\MATLAB\matnwb-2.5.0.0';
        addpath(genpath(matNWBloc));

    case 'otherwise'


end


switch stage2run
    case 1

        % PROCESS BEHAVIORAL FILES

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

    case 2


        % PROCESS TTL FILES

        % Load in RAW event files and RAW timestamps for events

        % Load in Processed Behavioral File table

        % FIGURE OUT new LABELS and TTL TIMES

        % RE-INSERT BACK INTO FILTER NWB

        cd(excelLOC)
        % Load in excel file
        allSUBtable = readtable('Eyetracking patient summary sheet.xlsx');
        for si = 1:1

            subROW = allSUBtable(si,:);
            subID = subROW.patientID{1};

            for vi = 1:3
                tmpVar = ['variant',num2str(vi)];
                if subROW.(tmpVar)

                    blmps = {'L','R'};
                    blmpsT = {'learn','recog'};
                    for bi = 1:2

                        % NWB f names
                        tmpNWBcol = ['nwbS_v' , num2str(vi) , '_' , blmps{bi}];
                        rawNWBname = subROW.(tmpNWBcol){1};
                        filtNWBname = replace(rawNWBname , 'filter' , 'raw');

                        % NWB LOC
                        nwbDATAdir = [DATADIR , filesep , subID , filesep , 'NWBprocessing-data\NWB_Data'];

                        % BEHAVIORAL LOC
                        subDATAdir = [DATADIR , filesep , subID , filesep , 'Behavioral-data' ,...
                            filesep , 'Processed'];

                        var2use = ['var',num2str(vi)];
                        textFname = [subID, '_' , var2use , '_' , blmpsT{bi} , '.mat'];

                        % NEED TO LOAD
                        cd(nwbDATAdir)
                        rawNWB = nwbRead(rawNWBname);

                        eventSTRS = rawNWB.acquisition.get('events').data.load();
                        eventTS = rawNWB.acquisition.get('events').timestamps.load();


                    end
                else
                    continue
                end
            end
        end


end  % END OF SWITCH CASE



end % END OF main function



















