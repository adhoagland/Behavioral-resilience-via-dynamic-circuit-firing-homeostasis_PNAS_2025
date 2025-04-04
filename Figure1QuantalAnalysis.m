%% Zach's quantal analysis screen

load('All TRIP RNAi SynapGCaMP6 29C 0_1Hz 50Stim Ib vs Is Boutons_AllData.mat', 'ExperimentSetLabels')
load('All TRIP RNAi SynapGCaMP6 29C 0_1Hz 50Stim Ib vs Is Boutons_AllData.mat', 'ExperimentSet_AllNMJ_QC_Mean_um2_Norm_byNMJ')
load('All TRIP RNAi SynapGCaMP6 29C 0_1Hz 50Stim Ib vs Is Boutons_AllData.mat', 'ExperimentSet_Mean_RFProb_byNMJ')

qcData = vertcat(ExperimentSet_AllNMJ_QC_Mean_um2_Norm_byNMJ{:});
prData = [horzcat(ExperimentSet_Mean_RFProb_byNMJ{:})]';

% IbQCdata = qcData(1:2:end);
% IsQCdata = qcData(2:2:end);
% IbPrData = prData(1:2:end);
% IsPrData = prData(2:2:end);

Iblabels = {};
labelCount = 1;
for qq = 1:2:numel(ExperimentSetLabels)
Iblabels{labelCount,1} = ExperimentSetLabels{qq};
labelCount = labelCount+1;
end

Islabels = {};
labelCount = 1;
for qq = 2:2:numel(ExperimentSetLabels)
Islabels{labelCount,1} = ExperimentSetLabels{qq};
labelCount = labelCount+1;
end

IbQC = {};
dataCount = 1;
IBlabels = {};
IbQCdata = [];
IbPrData = [];
for qq = 1:2:numel(ExperimentSet_AllNMJ_QC_Mean_um2_Norm_byNMJ)
IbQC{dataCount,1} = ExperimentSet_AllNMJ_QC_Mean_um2_Norm_byNMJ{qq};
IbPr{dataCount,1} = ExperimentSet_Mean_RFProb_byNMJ{qq};
theseData = IbQC{dataCount,1};
theseDataPr = IbPr{dataCount,1};
IbQCdata = [IbQCdata; theseData];
IbPrData = [IbPrData; theseDataPr'];
IBlabels = vertcat(IBlabels,repmat({Iblabels{dataCount}},numel(theseData),1));
dataCount = dataCount+1;
end

IndexC = strfind(IBlabels,'Attp2');
Index = find(not(cellfun('isempty',IndexC)));
IBlabels(Index,1) = IBlabels(Index(1),1);


IsQC = {};
dataCount = 1;
ISlabels = {};
IsQCdata = [];
IsPrData = [];
for qq = 2:2:numel(ExperimentSet_AllNMJ_QC_Mean_um2_Norm_byNMJ)
IsQC{dataCount,1} = ExperimentSet_AllNMJ_QC_Mean_um2_Norm_byNMJ{qq};
IbPr{dataCount,1} = ExperimentSet_Mean_RFProb_byNMJ{qq};
theseData = IsQC{dataCount,1};
theseDataPr = IbPr{dataCount,1};
IsQCdata = [IsQCdata; theseData];
IsPrData = [IsPrData; theseDataPr'];
ISlabels = vertcat(ISlabels,repmat({Iblabels{dataCount}},numel(theseData),1));
dataCount = dataCount+1;
end

IndexC = strfind(ISlabels,'Attp2');
Index = find(not(cellfun('isempty',IndexC)));
ISlabels(Index,1) = ISlabels(Index(1),1);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ib QC graphs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

attp2sites = {'Attp2';'Cac';'Syt1';'rab3';'Munc13';'VGlut'};
attp2indices = [];
attp2simplifiedLabels = {};
for ww = 1:numel(attp2sites)
IndexC = strfind(IBlabels,attp2sites{ww});
Index = find(not(cellfun('isempty',IndexC)));
attp2indices = [attp2indices;Index];
attp2simplifiedLabels = [attp2simplifiedLabels; repmat(attp2sites(ww),numel(Index),1)];
end

attp40sites = {'Attp40';'RIM';'RBP';'Munc18';'LiprinA';'LAR';'Sh'};
attp40indices = [];
attp40simplifiedLabels = {};
for ww = 1:numel(attp40sites)
IndexC = strfind(IBlabels,attp40sites{ww});
Index = find(not(cellfun('isempty',IndexC)));
attp40indices = [attp40indices;Index];
attp40simplifiedLabels = [attp40simplifiedLabels; repmat(attp40sites(ww),numel(Index),1)];
end

catData = IbQCdata(attp2indices);
genotypeStr = attp2simplifiedLabels;
T = table(genotypeStr,catData);
genotype = categorical(T.genotypeStr);

controlSites = contains(T.genotypeStr,'Attp2');
theseControlLocs = T.genotypeStr(controlSites);
controlName = theseControlLocs{1};
uniqueClasses = unique(genotype);
nClasses = numel(uniqueClasses);
controlData = catData(controlSites);
pVals = {};
for classNo = 1:nClasses
    thisClass = string(uniqueClasses(classNo));
    testSites = contains(T.genotypeStr,thisClass);
    testData = catData(testSites);
    pVals{classNo} = ranksum(controlData,testData);
end

% catData = IbQCdata(attp2indices);

figure('Position',[10 100 1100 1100])

a1 = subplot(2,2,1);
boxchart(a1,genotype,catData)
hold on
% scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1)
scatter(genotype,catData,200,'.','jitter','on', 'jitterAmount',0.1,'markerfacecolor','none','markeredgecolor','k')

ylim(a1,[0 0.07]); 

a1 = axes('Position', get(a1, 'Position'),'Color', 'none');
set(a1,'Ytick',[]) % To hide y-tick
set(a1,'XAxisLocation','top')
set(a1,'XTick',[1/(nClasses*2):1/nClasses:1],'XTickLabels',pVals)
% set(a1,'YLim',[0 0.07])




% 
% pVals=[1 ranksum(attp2Quantal,unc13Quantal) ranksum(attp2Quantal,rbpQuantal) ranksum(attp2Quantal,cacQuantal)]
% 
% xlabel(num2str(pVals))

subplot(2,2,2)
catData = IbQCdata(attp40indices);
% genotypeStr = IBlabels(attp40indices);
genotypeStr = attp40simplifiedLabels;

T = table(genotypeStr,catData);
genotype = categorical(T.genotypeStr);

controlSites = contains(T.genotypeStr,'Attp40');
theseControlLocs = T.genotypeStr(controlSites);
controlName = theseControlLocs{1};
uniqueClasses = unique(genotype);
nClasses = numel(uniqueClasses);
controlData = catData(controlSites);
pVals = {};
for classNo = 1:nClasses
    thisClass = string(uniqueClasses(classNo));
    testSites = contains(T.genotypeStr,thisClass);
    testData = catData(testSites);
    pVals{classNo} = ranksum(controlData,testData);
end

a1 = subplot(2,2,2);
boxchart(a1,genotype,catData)
hold on
scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);
ylim(a1,[0 0.07]); 


a2 = axes('Position', get(a1, 'Position'),'Color', 'none');
set(a2,'Ytick',[]) % To hide y-tick
set(a2,'XAxisLocation','top')
set(a2,'XTick',[1/(nClasses*2):1/nClasses:1],'XTickLabels',pVals)


% 
% pVals=[1 ranksum(attp2Quantal,unc13Quantal) ranksum(attp2Quantal,rbpQuantal) ranksum(attp2Quantal,cacQuantal)]
% 
% xlabel(num2str(pVals))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ib Pr graphs %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% catData = IbQCdata(attp2indices);
catData = IbPrData(attp2indices);
genotypeStr = attp2simplifiedLabels;

T = table(genotypeStr,catData);
genotype = categorical(T.genotypeStr);

controlSites = contains(T.genotypeStr,'Attp2');
theseControlLocs = T.genotypeStr(controlSites);
controlName = theseControlLocs{1};
uniqueClasses = unique(genotype);
nClasses = numel(uniqueClasses);
controlData = catData(controlSites);
pVals = {};
for classNo = 1:nClasses
    thisClass = string(uniqueClasses(classNo));
    testSites = contains(T.genotypeStr,thisClass);
    testData = catData(testSites);
    pVals{classNo} = ranksum(controlData,testData);
end

a1 = subplot(2,2,3);
boxchart(a1,genotype,catData)
hold on
scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);
ylim(a1,[0 0.3]); 

a2 = axes('Position', get(a1, 'Position'),'Color', 'none');
set(a2,'Ytick',[]) % To hide y-tick
set(a2,'XAxisLocation','top')
set(a2,'XTick',[1/(nClasses*2):1/nClasses:1],'XTickLabels',pVals)
% ylim([0 0.3])

% 
% pVals=[1 ranksum(attp2Quantal,unc13Quantal) ranksum(attp2Quantal,rbpQuantal) ranksum(attp2Quantal,cacQuantal)]
% 
% xlabel(num2str(pVals))

catData = IbPrData(attp40indices);
% genotypeStr = IBlabels(attp40indices);
genotypeStr = attp40simplifiedLabels;

T = table(genotypeStr,catData);
genotype = categorical(T.genotypeStr);

controlSites = contains(T.genotypeStr,'Attp40');
theseControlLocs = T.genotypeStr(controlSites);
controlName = theseControlLocs{1};
uniqueClasses = unique(genotype);
nClasses = numel(uniqueClasses);
controlData = catData(controlSites);
pVals = {};
for classNo = 1:nClasses
    thisClass = string(uniqueClasses(classNo));
    testSites = contains(T.genotypeStr,thisClass);
    testData = catData(testSites);
    pVals{classNo} = ranksum(controlData,testData);
end

a1 = subplot(2,2,4);
boxchart(a1,genotype,catData)
hold on
scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);
ylim(a1,[0 0.3]); 

a2 = axes('Position', get(a1, 'Position'),'Color', 'none');
set(a2,'Ytick',[]) % To hide y-tick
set(a2,'XAxisLocation','top')
set(a2,'XTick',[1/(nClasses*2):1/nClasses:1],'XTickLabels',pVals)

% ylim([0 0.3])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Is QC graphs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position',[10 100 1100 1100])


attp2sites = {'Attp2';'Cac';'Syt1';'rab3';'Munc13';'VGlut'};
attp2indices = [];
attp2simplifiedLabels = {};
for ww = 1:numel(attp2sites)
IndexC = strfind(IBlabels,attp2sites{ww});
Index = find(not(cellfun('isempty',IndexC)));
attp2indices = [attp2indices;Index];
attp2simplifiedLabels = [attp2simplifiedLabels; repmat(attp2sites(ww),numel(Index),1)];
end

attp40sites = {'Attp40';'RIM';'RBP';'Munc18';'LiprinA';'LAR';'Sh'};
attp40indices = [];
attp40simplifiedLabels = {};
for ww = 1:numel(attp40sites)
IndexC = strfind(IBlabels,attp40sites{ww});
Index = find(not(cellfun('isempty',IndexC)));
attp40indices = [attp40indices;Index];
attp40simplifiedLabels = [attp40simplifiedLabels; repmat(attp40sites(ww),numel(Index),1)];
end

% catData = IbQCdata(attp2indices);
catData = IsQCdata(attp2indices);
genotypeStr = attp2simplifiedLabels;

T = table(genotypeStr,catData);
genotype = categorical(T.genotypeStr);

controlSites = contains(T.genotypeStr,'Attp2');
theseControlLocs = T.genotypeStr(controlSites);
controlName = theseControlLocs{1};
uniqueClasses = unique(genotype);
nClasses = numel(uniqueClasses);
controlData = catData(controlSites);
pVals = {};
for classNo = 1:nClasses
    thisClass = string(uniqueClasses(classNo));
    testSites = contains(T.genotypeStr,thisClass);
    testData = catData(testSites);
    pVals{classNo} = ranksum(controlData,testData);
end

a1 = subplot(2,2,1);
boxchart(a1,genotype,catData)
hold on
scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);
ylim(a1,[0 0.07]); 

a2 = axes('Position', get(a1, 'Position'),'Color', 'none');
set(a2,'Ytick',[]) % To hide y-tick
set(a2,'XAxisLocation','top')
set(a2,'XTick',[1/(nClasses*2):1/nClasses:1],'XTickLabels',pVals)
ylim([0 0.07])

% 
% pVals=[1 ranksum(attp2Quantal,unc13Quantal) ranksum(attp2Quantal,rbpQuantal) ranksum(attp2Quantal,cacQuantal)]
% 
% xlabel(num2str(pVals))

catData = IsQCdata(attp40indices);
% genotypeStr = IBlabels(attp40indices);
genotypeStr = attp40simplifiedLabels;

T = table(genotypeStr,catData);
genotype = categorical(T.genotypeStr);

controlSites = contains(T.genotypeStr,'Attp40');
theseControlLocs = T.genotypeStr(controlSites);
controlName = theseControlLocs{1};
uniqueClasses = unique(genotype);
nClasses = numel(uniqueClasses);
controlData = catData(controlSites);
pVals = {};
for classNo = 1:nClasses
    thisClass = string(uniqueClasses(classNo));
    testSites = contains(T.genotypeStr,thisClass);
    testData = catData(testSites);
    pVals{classNo} = ranksum(controlData,testData);
end

a1 = subplot(2,2,2);
boxchart(a1,genotype,catData)
hold on
scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);
ylim(a1,[0 0.07]); 

a2 = axes('Position', get(a1, 'Position'),'Color', 'none');
set(a2,'Ytick',[]) % To hide y-tick
set(a2,'XAxisLocation','top')
set(a2,'XTick',[1/(nClasses*2):1/nClasses:1],'XTickLabels',pVals)

ylim([0 0.07])

% 
% pVals=[1 ranksum(attp2Quantal,unc13Quantal) ranksum(attp2Quantal,rbpQuantal) ranksum(attp2Quantal,cacQuantal)]
% 
% xlabel(num2str(pVals))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Is Pr graphs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% catData = IbQCdata(attp2indices);
catData = IsPrData(attp2indices);
genotypeStr = attp2simplifiedLabels;

T = table(genotypeStr,catData);
genotype = categorical(T.genotypeStr);

controlSites = contains(T.genotypeStr,'Attp2');
theseControlLocs = T.genotypeStr(controlSites);
controlName = theseControlLocs{1};
uniqueClasses = unique(genotype);
nClasses = numel(uniqueClasses);
controlData = catData(controlSites);
pVals = {};
for classNo = 1:nClasses
    thisClass = string(uniqueClasses(classNo));
    testSites = contains(T.genotypeStr,thisClass);
    testData = catData(testSites);
    pVals{classNo} = ranksum(controlData,testData);
end

a1 = subplot(2,2,3);
boxchart(a1,genotype,catData)
hold on
scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);
ylim(a1,[0 0.8]); 

a2 = axes('Position', get(a1, 'Position'),'Color', 'none');
set(a2,'Ytick',[]) % To hide y-tick
set(a2,'XAxisLocation','top')
set(a2,'XTick',[1/(nClasses*2):1/nClasses:1],'XTickLabels',pVals)

ylim([0 0.8])

% 
% pVals=[1 ranksum(attp2Quantal,unc13Quantal) ranksum(attp2Quantal,rbpQuantal) ranksum(attp2Quantal,cacQuantal)]
% 
% xlabel(num2str(pVals))

catData = IsPrData(attp40indices);
% genotypeStr = IBlabels(attp40indices);
genotypeStr = attp40simplifiedLabels;

T = table(genotypeStr,catData);
genotype = categorical(T.genotypeStr);

controlSites = contains(T.genotypeStr,'Attp40');
theseControlLocs = T.genotypeStr(controlSites);
controlName = theseControlLocs{1};
uniqueClasses = unique(genotype);
nClasses = numel(uniqueClasses);
controlData = catData(controlSites);
pVals = {};
for classNo = 1:nClasses
    thisClass = string(uniqueClasses(classNo));
    testSites = contains(T.genotypeStr,thisClass);
    testData = catData(testSites);
    pVals{classNo} = ranksum(controlData,testData);
end

a1 = subplot(2,2,4);
boxchart(a1,genotype,catData)
hold on
scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);
ylim(a1,[0 0.8]); 

a2 = axes('Position', get(a1, 'Position'),'Color', 'none');
set(a2,'Ytick',[]) % To hide y-tick
set(a2,'XAxisLocation','top')
set(a2,'XTick',[1/(nClasses*2):1/nClasses:1],'XTickLabels',pVals)

% ylim([0 0.8])



