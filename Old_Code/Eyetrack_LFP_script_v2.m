%--------------------------------------------------------------------------
% Load in data
%--------------------------------------------------------------------------

%------- NWB -------%

% Add matnwb interface library to path
uiwait(msgbox('Navigate to and select matnwb folder'))
matnwbLoc = uigetdir();
addpath(genpath(matnwbLoc));

% Add subject folder to path
uiwait(msgbox('Navigate to and select patient folder'))
patientLoc = uigetdir();
addpath(genpath(patientLoc));

% Load NWB file
  %To do future % uiwait(msgbox('Navigate to and select nwb file to load'))
                % filename = uigetdir();
filename = 'MW9_Session_2_filter.nwb';
nw = nwbRead(filename);

% Pop up to browse contents of NWB file
util.nwbTree(nw);

%------- Eye -------%

% Convert .edf file to .mat
edfFile = 'NO20221615110.edf';
patientID = 'MW9';
block = 'L';
recordingDate = '01062022';
[eyeMatfile] = Edf2Mat_UCH(edfFile, patientID, block, recordingDate);

% Create tables from relevant file data
[tsTable, picTable, fixTab, saccTab, rawTab] = ExtractEyeInfo_v2(eyeMatfile);

%--------------------------------------------------------------------------
% Extract data of interest: TS of LFPs, TS of events
%--------------------------------------------------------------------------

%------- NWB -------%

% Navigate to MW LFP timestamps
nw.processing.get('ecephys');
elecSeries =  nw.processing.get('ecephys').nwbdatainterface.get('LFP').electricalseries;
nwSeries = elecSeries.get('MacroWireSeries');
nw_TS = nwSeries.timestamps.load;  
    
% Downsample TS by factor of 8 to match length of nw_data
nw_TSdn = downsample(nw_TS, 8); %conceptually equal to time column of rawTab in eyetracking

% Load macro data (downsampled by factor of 8 and filtered)
nw_data = nwSeries.data.load(); 

% Isolate channel 15 for now... and convert to class 'double'
nw_data_chan = 15; 
channelData = double(nw_data(nw_data_chan, :)); 

% Load event IDs
nwEventID = nw.acquisition.get('events').data.load();

% Load event TS
nwEventTS = nw.acquisition.get('events').timestamps.load();

% Isolate specific TTLs from nwEventsTS
    %extract time marker for ttl 1 in ephys by creating an index
        %time between ttl 1 and ttl 2 = time of 1st trial
        %row index (integer) of macrowire data that matches ttl1 and ttl2
    % Time markers for ttl 1 and 2
    ttl1 = nwEventTS(2);
    ttl2 = nwEventTS(3);

% Find closest matching ephys TS to ttl1 TS
[~, ttl1_index] = min(abs(ttl1 - nw_TSdn)); %a = how close they are in microsec %b = index

% Find closest matching  ephys TS to ttl2 TS
[~, ttl2_index] = min(abs(ttl2 - nw_TSdn));    

% Extract LFP data b/t ttl1_index and ttl2_index
sampleLFP = channelData(:,ttl1_index:ttl2_index);

% Extract LFP data from ttl1_index + 100ms (sr = 500 * fraction of time we want .1 = 50)
ttl1_index_100ms = ttl1_index + 50;
sampleLFP2 = channelData(:,ttl1_index:ttl1_index_100ms); 
 
%------- Eye -------%

% Extract data from rawTab: time, PosX, PosY
eye_TS = rawTab.Time; 
eye_PosX = rawTab.PosX;
eye_PosY = rawTab.PosY;

% Behavioral events (ttls) of eye tracking converted to class 'double' (double = 64 bit)
eye_TTL = double(tsTable.timeStamp);

% Isolate 1st and 2nd TTL
eyeTTL1 = eye_TTL(1);
eyeTTL2 = eye_TTL(2); 

% Find eyetracking data between TTL1 and TTL2
[~, eyeTTL1_i] = min(abs(eyeTTL1 - eye_TS)); 
[~, eyeTTL2_i] = min(abs(eyeTTL2 - eye_TS)); 

% Extract eye data b/t eyeTTL1_i and eyeTTL2_i
sample_eye = eye_PosX(eyeTTL1_i:eyeTTL2_i, :);
    