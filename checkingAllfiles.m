% Check allsession data for MW9

%% CD to directory

directory = 'Y:\';
patDir = 'PatientData_MW\MW9\NWBProcessing\';
tmpDir = [directory , patDir];


%% load in allsessions file

sessDir = [tmpDir , 'Session_Data\'];
cd(sessDir)
load('allSessionData.mat')

putBeh = zeros(length(saveSessAll),1,'logical');
durAll = zeros(length(saveSessAll),1);
dateAll = cell(length(saveSessAll),1);
sessAAL = transpose(1:length(saveSessAll));
for pi = 1:length(putBeh)

    tmpStart = saveSessAll{pi,1}.StartIndex;
    tmpStop = saveSessAll{pi,1}.StopIndex;
    tmpDur = tmpStop - tmpStart;
    tmpStrg = saveSessAll{pi,1}.SessionInfo.header(7);

    if tmpDur > 400
        putBeh(pi) = true;
        durAll(pi) = tmpDur;
        [dateCol] = getTextInfo(tmpStrg);
        dateAll{pi} = dateCol;

    end
end

sessionToView = saveSessAll(putBeh);
durAll2 = durAll(durAll ~= 0);
dateAll2 = dateAll(durAll ~= 0);
sessAAL2 = sessAAL(durAll ~= 0);

table(sessAAL2, durAll2,dateAll2,'VariableNames',{'SessNum','Duration','Date'})


%% Check for Starting Recording

stINDs = find(matches(saveSessAll{2,1}.SessionInfo.tInfo,'Starting Recording'));
firstTTls = stINDs + 1;
block1start = 40;
block1end = 537;
block2start = 542;
block2end = 1046;

% Get block ttl count n = 6
start = [40; 542; 2; 523; 3; 521];
stop = [537; 1046; 518; 1026; 516; 1015];
numEvents = stop - start;

%% Load in event files and match block numbers

rawtBloc = 'Y:\PatientData_MW\MW9\Behavioral-Data\RawData';
raweBloc =  'Y:\PatientData_MW\MW9\Eye-tracking\Raw';
proeBloc = 'Y:\PatientData_MW\MW9\Eye-tracking\Processed';

[txtList] = getFlist(rawtBloc, 1);
[eyeList] = getFlist(raweBloc, 2);
dateID = cell(length(txtList),1);
blockt = cell(length(txtList),1);
varT = cell(length(txtList),1);
texTfn = cell(length(txtList),1);
eyefn = cell(length(txtList),1);
eyePfn = cell(length(txtList),1);

for bi = 1:length(txtList)
    cd(rawtBloc)
    tmpTxt = txtList{bi};
    tmpSp = split(tmpTxt,'.');
    tmpRm = extractAfter(tmpSp,'NEWOLDDELAY5');
    tmpRm2 = tmpRm{cellfun(@(x) ~isempty(x), tmpRm)};
    dateID{bi} = tmpRm2;

    % get block id
    txtTAB = readtable(tmpTxt);
    txtLINE = table2cell(txtTAB(1,3));
    % first is variant [1, 2, 3]
    % second is block [1 = learn, 2 = recog]
    txtSp = split(txtLINE,',');
    if matches(txtSp{2},'1')
        blockt{bi} = 'Learn';
    elseif matches(txtSp{2},'2')
        blockt{bi} = 'Recog';
    end
    
    varT{bi} = ['Var_', txtSp{1}];

    % find EDF file
    cd(raweBloc)
    eyeSp = split(eyeList,'.');
    eyeSpA = eyeSp(:,:,1);
    eyeRM = extractAfter(eyeSpA,'NO');

    eyeFind = eyeList{matches(eyeRM,tmpRm2)};
    texTfn{bi} = tmpTxt;
    eyefn{bi} = eyeFind;

    eyePrName = Edf2Mat_UCH(eyeFind, 'MW9', blockt{bi}(1), txtSp{1}, raweBloc, proeBloc);

    eyePfn{bi} = eyePrName;

end





%% Align eye-track / ttl event ids





















%%











function [output] = getTextInfo(input)

slice1 = split(input,' ');
% slice2 = slice1(7);
% slice3 = split(slice2,'_');
dAY = slice1{2};
hOUR = slice1{3};
timeColl = [dAY ' ' hOUR];
output = datetime(timeColl, 'InputFormat', 'yyyy/MM/dd HH:mm:SS');


end



function [flist] = getFlist(dirIN, type)

cd(dirIN)

switch type
    case 1
        flin1 = dir('*.txt');
    case 2
        flin1 = dir('*.edf');
end
flin2 = {flin1.name};
flist = flin2;


end



