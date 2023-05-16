function [] = Extract_Eye_EDF(excelLOC , mainLOC, ptID, edf2matLOC)
% 'The events coorespond to the TTL markers for each trial. ', ...
%     'For the learning trials, the TTL markers are the following: 55 = start of the experiment, ', ...
%     '1 = stimulus ON, 2 = stimulus OFF, 3 = Question Screen Onset, ', ...
%     '20 = Yes (21 = NO) during learning, 6 = End of Delay after Response, ', ...
%     '66 = End of Experiment. For the recognition trials, the TTL markers are the following: ', ...
%     '55 = start of experiment, 1 = stimulus ON, 2 = stimulus OFF, 3 = Question Screen Onset, ' ...
%     '31:36 = Confidence (Yes vs. No) response, 66 = End of Experiment'

% Loop through files - determine if learn or retrieve

% learning = [55, 1, 2, 3 20, 21, 6, 66];
% recog = [55, 1, 2, 3, 31:36, 66];

% Location of variant.xlsx
cd(excelLOC)
varTable = readtable('variantLIST.xlsx');
% Location of eye tracking folders
cd(mainLOC)
[outFOLDS] = getfiles(mainLOC,1,nan);

% Match specific patient ID to location in outFOLDS
ptIdx = find(strcmp(outFOLDS, ptID));

tempCASEd = [mainLOC , filesep , outFOLDS{ptIdx}];
idTab = varTable(matches(varTable.Subject,outFOLDS{ptIdx}),:);
cd(tempCASEd)
allvars = unique(idTab.Variant);
for vi = 1:length(allvars)
    vartiTAB =  idTab(ismember(idTab.Variant,allvars(vi)),:);
    learnEDF = vartiTAB.EDF(matches(vartiTAB.Block,'learn'));
    recogEDF = vartiTAB.EDF(matches(vartiTAB.Block,'recog'));
    varianT = num2str(vi); % make sure iterator is equivalent to variant ID

    [ltmpName] = Edf2Mat_UCH(edf2matLOC, learnEDF{1}, outFOLDS{ptIdx}, 'Learn',...
        varianT, tempCASEd, tempCASEd);

    varTable.MatFile(matches(varTable.Subject,outFOLDS{ptIdx}) &...
        ismember(varTable.Variant,allvars(vi)) &...
        matches(varTable.Block,'learn')) = {ltmpName};

    [rtmpName] = Edf2Mat_UCH(edf2matLOC, recogEDF{1}, outFOLDS{ptIdx}, 'Recog',...
        varianT, tempCASEd, tempCASEd);

    varTable.MatFile(matches(varTable.Subject,outFOLDS{ptIdx}) &...
        ismember(varTable.Variant,allvars(vi)) &...
        matches(varTable.Block,'recog')) = {rtmpName};

end

cd(excelLOC)
writetable(varTable,'variantLIST.xlsx')


end % END MAIN FUNCTION




function [outfiles] = getfiles(dirIN,stage,ftype)

cd(dirIN)
switch stage
    case 1

        foldeS = dir();
        foldeS2 = {foldeS.name};
        foldeS3 = foldeS2(~ismember(foldeS2,{'.','..'}));
        outfiles = foldeS3;
    case 2

        filES = dir(['*.',ftype]);
        filES2 = {filES.name};
        outfiles = filES2;

end


end


