function [] = findSTART_STOPinds()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


curPC = getenv('COMPUTERNAME');

switch curPC
    case 'DESKTOP-I5CPDO7'   % JAT WORK desktop

        excelLOC = 'Z:\MW_JAT_Backup\EyeTrack_Manuscript';
        DATADIR = 'Z:\MW_JAT_Backup\EyeTrack_Manuscript\Patient_folders';

        matNWBloc = 'C:\Users\Admin\Documents\MATLAB\matnwb-2.5.0.0';
        addpath(genpath(matNWBloc));

    case 'otherwise'


end




cd(excelLOC)
% Load in excel file
allSUBtable = readtable('Eyetracking patient summary sheet.xlsx');
allCOUNT = 1;
allCell = cell(100,5);
for si = 6:6

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

                % CHECK NWB_LOC folder for DAYS
                tmpCheckd = dir(nwbDATAdir);
                tmpCheckd2 = {tmpCheckd.name};
                tmpCheckd3 = tmpCheckd2(~ismember(tmpCheckd2,{'..','.'}));

                if ~any(contains(tmpCheckd3,'variant'))
                     nwbDATAdir = [DATADIR , filesep , subID , filesep , 'NWBprocessing-data\NWB_Data'];
                else

                    tmpDIRchck = [nwbDATAdir , filesep , tmpVar];
                    nwbDATAdir = tmpDIRchck;

                end




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
                [startNLXttlind , stopNLXttlind] = alignTTLevents(taskTable,eventTS,eventSTRS,allCOUNT);

                allCell{allCOUNT,1} = subID;
                allCell{allCOUNT,2} = vi;
                allCell{allCOUNT,3} = blmps{bi};
                allCell{allCOUNT,4} = startNLXttlind;
                allCell{allCOUNT,5} = stopNLXttlind;

                allCOUNT = allCOUNT + 1;



            end


        end


    end

    disp(['Num ', num2str(si) ,' out of 10 done'])

end

% 
rows2keep = cellfun(@(x) ~isempty(x), allCell(:,1), 'UniformOutput', true);
allCell = allCell(rows2keep,:);

% cd('Z:\MW_JAT_Backup\EyeTrack_Manuscript')
save("NLX_Indices_NOcases.mat","allCell");




end





















function [startINDEX , stopINDEX] = alignTTLevents(tskTABLE,evtTstamps,evtSTrings,allCOUNT)

startINDEX = 0;

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
    % getEvenStrings = evtSTrings2(eventPROC_1_ind);
    offSetZeros = diff(getZeroTimes)/1000000;

    % NEED TO MAKE OFFSET A FRACTION -----

    sharPoffsetS = find(offSetZeros >= 8);
    % sharPoffset + 1 = start of new segment

    if length(getZeroTimes) == height(tskTABLE)

        startINDEX = 1;
        stopINDEX = height(tskTABLE);

    else

        for si = 1:length(sharPoffsetS)
            tmpSlength = numel(sharPoffsetS(si) + 1:length(getZeroTimes)); % spikes in front of actual data
            if tmpSlength == height(tskTABLE)
                startINDEX = sharPoffsetS(si) + 1;
                stopINDEX = length(getZeroTimes);
            end

            if sharPoffsetS(si) >= 502
                startINDEX = 1;
                stopINDEX = 502;


            end

        end

        if startINDEX == 0

            tTimes = diff(tskTABLE.Timestamp);
            zTimes = diff(getZeroTimes);
            zTimes = round(zTimes / 100);

            if length(tTimes) > length(zTimes)
                tTimes = tTimes(1:length(zTimes));
            elseif length(zTimes) > length(tTimes)
                zTimes = zTimes(1:length(tTimes));
            end
            
            % zTimes = zTimes([1:132,134:length(zTimes)] );
            % if length(zTimes) > length(tTimes)
            %     frontOffset = length(zTimes) - length(tTimes);
            %     % zTimes = zTimes(1+frontOffset:length(zTimes));
            %     zTimes = zTimes(1:length(tTimes));
            % end

            [corRel,~] = corr(zTimes,tTimes);

            if corRel > 0.88
                startINDEX = 1;
                stopINDEX = length(getZeroTimes);

            else % FIND WHERE OFFSET IS
                % trim longer
                % IF A CERAIN CASE - MW3 1 learn

                if allCOUNT == 11

                    zTimes = diff(getZeroTimes);
                    zTimes = round(zTimes / 100);
                    zTimes = zTimes([1:132,134:length(zTimes)] );
                    [corRel,~] = corr(zTimes,tTimes);
                    if corRel > 0.88
                        startINDEX = [1 , 134];
                        stopINDEX = [132, length(getZeroTimes)];
                    end

                end



            end
        end

    end

    % trimEvents = startINDEX:length(getZeroTimes);
    % 
    % finalEvtStrings = getEvenStrings(trimEvents);
    % finalEvtTimes = getZeroTimes(trimEvents);


else % HEX cases


    % eventPROC_1_ind = contains(evtSTrings2,'AcqSystem');
    % eventPROC_1 = evtSTrings2(contains(evtSTrings2,'AcqSystem'));

    % eventPROC_2 = extractBetween(eventPROC_1,'(',')');
    % eventHEX = hex2dec(eventPROC_2);

    if length(eventHEX) > 502

        % Remove zeros
        eventHEX = eventHEX(eventHEX ~= 0);
        % Find experiment starts
        expSTARTS = find(eventHEX == 55);

        if height(tskTABLE) > 502
            % remove no labels
            noLabR = ~matches(tskTABLE.TTL_ID,'NO_label');
            tskTABLE = tskTABLE(noLabR,:);
        end

        for si = 1:length(expSTARTS)
            tmpSlength = numel(expSTARTS(si):length(eventHEX)); % spikes in front of actual data
            if tmpSlength == height(tskTABLE)
                startINDEX = expSTARTS(si);
                stopINDEX = length(eventHEX);
            end

            % if sharPoffsetS(si) >= 502
            %     startINDEX = 1;
            %     stopINDEX = 502;
            % 
            % end
        end

        if startINDEX == 0
            if eventHEX(1) == 55

                startINDEX = 1;
                stopINDEX = length(eventHEX);

            end
        end

    else % TTLS ARE MISSING

        if height(tskTABLE) > 502
            % remove no labels
            noLabR = ~matches(tskTABLE.TTL_ID,'NO_label');
            tskTABLE = tskTABLE(noLabR,:);
        end

        % Remove zeros
        eventHEX = eventHEX(eventHEX ~= 0);
        % Remove 255
        eventHEX = eventHEX(eventHEX ~= 255);

        if length(eventHEX) == height(tskTABLE)

            startINDEX = 1;
            stopINDEX = 502;

        elseif eventHEX(1) ~= 55
            % CHECK IF MISSING FROM front
            frontOFFset = 502 - length(eventHEX);
            ttlValoff = tskTABLE.TTLvalue(frontOFFset+1:502);

            if isequal(ttlValoff,eventHEX)

                startINDEX = -(frontOFFset+1);
                stopINDEX = -502;

            end

        else

            startINDEX = 1;
            stopINDEX = length(eventHEX);

        end


    end





end


% tskTABLE.TTLstrings = finalEvtStrings;
% tskTABLE.TTLtsTamps = finalEvtTimes;

% fullTABLE = tskTABLE;




end
























