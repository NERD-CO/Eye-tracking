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
                        filtNWBname = subROW.(tmpNWBcol){1};
                        rawNWBname = replace(filtNWBname , 'filter' , 'raw');

                        % NWB LOC
                        nwbDATAdir = [DATADIR , filesep , subID , filesep , 'NWBprocessing-data\NWB_Data'];

                        % BEHAVIORAL LOC
                        subDATAdir = [DATADIR , filesep , subID , filesep , 'Behavioral-data' ,...
                            filesep , 'Processed'];

                        var2use = ['var',num2str(vi)];
                        matFname = [subID, '_' , var2use , '_' , blmpsT{bi} , '.mat'];

                        % Load acquisition event info
                        cd(nwbDATAdir)
                        rawNWB = nwbRead(rawNWBname);

                        eventSTRS = rawNWB.acquisition.get('events').data.load();
                        eventTS = rawNWB.acquisition.get('events').timestamps.load();

                        % Load behavioral file table
                        cd(subDATAdir)

                        load(matFname , 'outData')

                        taskTable = outData.taskinformation;

                        % Create processing function 
                        [taskTTLtable] = alignTTLevents(taskTable,eventTS,eventSTRS);

                        % Save new Table
                        % CD to new folder and SAVE
                        subDATAttldir = [DATADIR , filesep , subID , filesep , 'Behavioral-data' ,...
                            filesep , 'ProcessedTTL'];
                        cd(subDATAttldir)
                        saveTTLname = [subID, '_' , var2use , '_' , blmpsT{bi} , '_ttl.mat'];
                        save(saveTTLname , "taskTTLtable");

                        % Add to filtered NWB
                        % CD to NWB
                        cd(nwbDATAdir)
                        % Load filter
                        filterNWB = nwbRead(filtNWBname);
                        % Add field
                        behStrEvents = types.core.AnnotationSeries('data', taskTTLtable.TTL_ID,...
                            'data_unit', 'NA', 'timestamps', taskTTLtable.Timestamp, 'description', 'TTLs');

                        cellStrVals = cellfun(@(x) num2str(x) , num2cell(taskTTLtable.TTLvalue), 'UniformOutput',false);

                        behMatEvents = types.core.AnnotationSeries('data', cellStrVals,...
                            'data_unit', 'NA', 'timestamps', taskTTLtable.Timestamp, 'description', 'TTLs');

                        behNLXEvents = types.core.AnnotationSeries('data', taskTTLtable.TTLstrings,...
                            'data_unit', 'NA', 'timestamps', taskTTLtable.TTLtsTamps, 'description', 'TTLs');

                        filterNWB.acquisition.set('MatStrEvents', behStrEvents)
                        filterNWB.acquisition.set('MatIntEvents', behMatEvents)
                        filterNWB.acquisition.set('NLXEvents', behNLXEvents)

                        % Save 
                        nwbExport(filterNWB , filtNWBname);

                        % clear both NWBs
                        clear rawNWB filterNWB


                    end
                else
                    continue
                end
            end
        end


end  % END OF SWITCH CASE



end % END OF main function






function [fullTABLE] = alignTTLevents(tskTABLE,evtTstamps,evtSTrings)


evtSTrings2 = cellstr(evtSTrings);
% Remove any START and STOP
eventPROC_1_ind = contains(evtSTrings2,'AcqSystem');
eventPROC_1 = evtSTrings2(contains(evtSTrings2,'AcqSystem'));

eventPROC_2 = extractBetween(eventPROC_1,'(',')');
eventHEX = hex2dec(eventPROC_2);

% If all Zeros - take every other zero
fracZero = sum(eventHEX == 0) / length(eventHEX);

if fracZero > 0.8

    getZeroTimes = evtTstamps(eventPROC_1_ind);
    getEvenStrings = evtSTrings2(eventPROC_1_ind);
    offSetZeros = diff(getZeroTimes)/1000000;

    % NEED TO MAKE OFFSET A FRACTION ----- 

    sharPoffsetS = find(offSetZeros > 12);
    % sharPoffset + 1 = start of new segment

    for si = 1:length(sharPoffsetS)
        tmpSlength = numel(sharPoffsetS(si) + 1:length(getZeroTimes));
        if tmpSlength == height(tskTABLE)
            startINDEX = sharPoffsetS(si) + 1;
        end
    end

    trimEvents = startINDEX:length(getZeroTimes);

    finalEvtStrings = getEvenStrings(trimEvents);
    finalEvtTimes = getZeroTimes(trimEvents);


else


    test = 1;





end


tskTABLE.TTLstrings = finalEvtStrings;
tskTABLE.TTLtsTamps = finalEvtTimes;

fullTABLE = tskTABLE;




end











