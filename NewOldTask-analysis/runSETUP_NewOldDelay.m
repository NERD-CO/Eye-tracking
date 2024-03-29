function [] = runSETUP_NewOldDelay(stage2run)
% SETUP code for NewOld delay preprocessing
% First run NewOldTxttoMat_v2.mat

curPC = getenv('COMPUTERNAME');

switch curPC
    case 'DESKTOP-I5CPDO7'   % JAT WORK desktop

        excelLOC = 'Y:\MW_JAT_Backup\EyeTrack_Manuscript';
        DATADIR = 'Y:\MW_JAT_Backup\EyeTrack_Manuscript\Patient_folders';

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

        for si = 8:8
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
        load("NLX_Indices_NOcases.mat","allCell");
        nlxINDEX = cell2table(allCell,'VariableNames',{'SubID','VarI','Block','StartI','StopI'});
        for si = 10:10

            subROW = allSUBtable(si,:);
            subID = subROW.patientID{1};
            subNLX = nlxINDEX(matches(nlxINDEX.SubID,subID),:);

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

                        % Extract Var/Block from NLX Table
                        varROW = subNLX.VarI == vi;
                        blkROW = matches(subNLX.Block,upper(blmpsT{bi}(1)));

                        subNLXrow = subNLX(varROW & blkROW,:);

                        % Create processing function
                        [taskTTLtable] = alignTTLevents(taskTable,eventTS,eventSTRS,subNLXrow);

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

    case 3 % UPDATE downsample TTL
        cd(excelLOC)
        % Load in excel file
        allSUBtable = readtable('Eyetracking patient summary sheet.xlsx');
        % load("NLX_Indices_NOcases.mat","allCell");
        % nlxINDEX = cell2table(allCell,'VariableNames',{'SubID','VarI','Block','StartI','StopI'});
        for si = 1:10

            subROW = allSUBtable(si,:);
            subID = subROW.patientID{1};
            % subNLX = nlxINDEX(matches(nlxINDEX.SubID,subID),:);

            for vi = 1:3
                tmpVar = ['variant',num2str(vi)];
                if subROW.(tmpVar)

                    blmps = {'L','R'};
                    % blmpsT = {'learn','recog'};
                    for bi = 1:2

                        % NWB f names
                        tmpNWBcol = ['nwbS_v' , num2str(vi) , '_' , blmps{bi}];
                        filtNWBname = subROW.(tmpNWBcol){1};

                        % NWB LOC
                        nwbDATAdir = [DATADIR , filesep , subID , filesep , 'NWBprocessing-data\NWB_Data'];

                        % Load acquisition event info
                        cd(nwbDATAdir)
                        filterNWB = nwbRead(filtNWBname);

                        nonDNS_ts = filterNWB.processing.get('ecephys')...
                            .nwbdatainterface.get('LFP').electricalseries ...
                            .get('MacroWireSeries').timestamps.load();

                        dnsDATA = filterNWB.processing.get('ecephys')...
                            .nwbdatainterface.get('LFP').electricalseries ...
                            .get('MacroWireSeries').data.load();

                        dnsTScheck = isequal(size(nonDNS_ts,1),size(dnsDATA,2));

                        if ~dnsTScheck

                            DNS_ts = downsample(nonDNS_ts,8);

                            % reduce data input 
                            % dnsSTEPS = round(linspace(1,length(nonDNS_ts),ceil(length(nonDNS_ts)/8)));

                            dns2check = isequal(size(DNS_ts,1),size(dnsDATA,2));

                            if dns2check

                                filterNWB.processing.get('ecephys')...
                                    .nwbdatainterface.get('LFP').electricalseries ...
                                    .get('MacroWireSeries').timestamps = DNS_ts;

                                filterNWB.processing.get('ecephys')...
                                    .nwbdatainterface.get('LFP').electricalseries ...
                                    .get('MacroWireSeries').description = 'Fs = 500 : macro timestamps have been downsampled';

                                nwbExport(filterNWB , filtNWBname)

                            end
                        end
                    end
                end
            end
        end



end  % END OF SWITCH CASE



end % END OF main function






function [fullTABLE] = alignTTLevents(tskTABLE,evtTstamps,evtSTrings,subINFO)

evtSTrings2 = cellstr(evtSTrings);
% Remove any START and STOP
eventPROC_1_ind = contains(evtSTrings2,'AcqSystem');
eventPROC_1 = evtSTrings2(contains(evtSTrings2,'AcqSystem'));

eventPROC_2 = extractBetween(eventPROC_1,'(',')');
eventHEX = hex2dec(eventPROC_2);

% If all Zeros - take every other zero
fracZero = sum(eventHEX == 0) / length(eventHEX);

if fracZero > 0.8

    eventTStamps2 = evtTstamps(eventPROC_1_ind);

    if length(subINFO.StartI{1}) == 1

        finalEvtStrings = eventPROC_2(subINFO.StartI{1}:subINFO.StopI{1});
        finalEvtTimes = eventTStamps2(subINFO.StartI{1}:subINFO.StopI{1});

        if length(finalEvtTimes) < height(tskTABLE)
            if subINFO.StartI{1} == 1
                finalEvtStrings = [finalEvtStrings ; {'NaN'}];
                finalEvtTimes = [finalEvtTimes ; NaN];
            end
        end

    elseif length(subINFO.StartI{1}) == 2

        segment1 = subINFO.StartI{1}(1):subINFO.StopI{1}(1);
        segment2 = subINFO.StartI{1}(2):subINFO.StopI{1}(2);

        finalEvtStrings = eventPROC_2([segment1 , segment2]);
        finalEvtTimes = eventTStamps2([segment1 , segment2]);

    end


else

    if subINFO.StartI{1}(1) > 0

        % eventPROC_2 = eventPROC_2
        eventPROC_3 = eventPROC_2(~ismember(eventHEX,[0 255]));
        % eventPROC_2 = eventPROC_2(eventHEX ~= 255);

        if height(tskTABLE) > 502
            tskTABLE = tskTABLE(tskTABLE.TTLvalue ~= 0,:);
        end

        eventTStamps2 = evtTstamps(eventPROC_1_ind);
        eventTStamps3 = eventTStamps2(~ismember(eventHEX,[0 255]));

        finalEvtStrings = eventPROC_3(subINFO.StartI{1}:subINFO.StopI{1});
        finalEvtTimes = eventTStamps3(subINFO.StartI{1}:subINFO.StopI{1});

        if length(finalEvtTimes) < height(tskTABLE)
            if subINFO.StartI{1} == 1
                belowOffset = height(tskTABLE) - length(finalEvtTimes);
                finalEvtStrings = [finalEvtStrings ; repmat({'NaN'},belowOffset,1)];
                finalEvtTimes = [finalEvtTimes ; nan(belowOffset,1)];
            else
                test = 1;
            end
        end

    else

        eventPROC_2 = eventPROC_2(~ismember(eventHEX,[0 255]));
        eventTStamps2 = evtTstamps(eventPROC_1_ind);
        eventTStamps3 = eventTStamps2(eventHEX ~= 0);

        if height(tskTABLE) > 502
            tskTABLE = tskTABLE(tskTABLE.TTLvalue ~= 0,:);
        end

        finalEvtStrings = [repmat({'NaN'},(subINFO.StartI{1}(1)*-1)-1,1) ; eventPROC_2];
        finalEvtTimes = [nan((subINFO.StartI{1}(1)*-1)-1,1) ; eventTStamps3];




    end


end




tskTABLE.TTLstrings = finalEvtStrings;
tskTABLE.TTLtsTamps = finalEvtTimes;

fullTABLE = tskTABLE;




end












