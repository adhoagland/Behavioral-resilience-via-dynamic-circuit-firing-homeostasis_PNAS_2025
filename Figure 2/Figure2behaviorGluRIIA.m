movieDuration = 600; % in sec

load('GluRIIAbehaviorailData.mat')

controlMeanCurve=[];
controldistanceStepMvals=[];
controldistanceStep=[];
controldistancesMvals=[];
controlmvals=[];
controlturningAngles = [];
controldistanceStepMean = [];
controlvelocityMean = [];
minTurningAngle = [];
for aa = 1:size(control,2)
    for subInd = 1:size(control(aa).meanCurv,1)
    controlMeanCurve = [controlMeanCurve;control(aa).meanCurv{subInd}];
    controldistanceStepMvals = [controldistanceStepMvals; control(aa).distanceStepsMvals{subInd}];
    controlmvals=[controlmvals; control(aa).mValCell{subInd}];
    controlturningAngles=[controlturningAngles; control(aa).trajectoryAngles{subInd}];
%     minTurningAngle=[minTurningAngle; min(control(aa).trajectoryAngles{subInd})];
    controldistanceStepMean = [controldistanceStepMean; nanmean(control(aa).distanceStepsMvals{subInd})];
    controlvelocityMean = [controlvelocityMean; (control(aa).distanceStepsMvals{subInd})];

    end
    
    controldistancesMvals=[controldistancesMvals;control(aa).distancesMvals];
    controlmeanVel = controldistancesMvals./movieDuration;
    
end

gluriiaMeanCurve=[];
gluriiadistanceStepMvals=[];
gluriiadistanceStep=[];
gluriiadistancesMvals=[];
gluriiamvals=[];
gluriiaturningAngles = [];
gluriiadistanceStepMean = [];
gluriiavelocityMean = [];
minTurningAngle = [];
for aa = 1:size(gluriia,2)
    for subInd = 1:size(gluriia(aa).meanCurv,1)
    gluriiaMeanCurve = [gluriiaMeanCurve;gluriia(aa).meanCurv{subInd}];
    gluriiadistanceStepMvals = [gluriiadistanceStepMvals; gluriia(aa).distanceStepsMvals{subInd}];
    gluriiamvals=[gluriiamvals; gluriia(aa).mValCell{subInd}];
    gluriiaturningAngles=[gluriiaturningAngles; gluriia(aa).trajectoryAngles{subInd}];
%     minTurningAngle=[minTurningAngle; min(gluriia(aa).trajectoryAngles{subInd})];
    gluriiadistanceStepMean = [gluriiadistanceStepMean; nanmean(gluriia(aa).distanceStepsMvals{subInd})];
    gluriiavelocityMean = [gluriiavelocityMean; (gluriia(aa).distanceStepsMvals{subInd})];

    end
    
    gluriiadistancesMvals=[gluriiadistancesMvals;gluriia(aa).distancesMvals];
    gluriiameanVel = gluriiadistancesMvals./movieDuration;
    
end

rbpMeanCurve=[];
rbpdistanceStepMvals=[];
rbpdistanceStep=[];
rbpdistancesMvals=[];
rbpmvals=[];
rbpturningAngles = [];
rbpdistanceStepMean = [];
rbpvelocityMean = [];
minTurningAngle = [];
for aa = 1:size(rbp,2)
    for subInd = 1:size(rbp(aa).meanCurv,1)
    rbpMeanCurve = [rbpMeanCurve;rbp(aa).meanCurv{subInd}];
    rbpdistanceStepMvals = [rbpdistanceStepMvals; rbp(aa).distanceStepsMvals{subInd}];
    rbpmvals=[rbpmvals; rbp(aa).mValCell{subInd}];
    rbpturningAngles=[rbpturningAngles; rbp(aa).trajectoryAngles{subInd}];
%     minTurningAngle=[minTurningAngle; min(rbp(aa).trajectoryAngles{subInd})];
    rbpdistanceStepMean = [rbpdistanceStepMean; nanmean(rbp(aa).distanceStepsMvals{subInd})];
    rbpvelocityMean = [rbpvelocityMean; (rbp(aa).distanceStepsMvals{subInd})];

    end
    
    rbpdistancesMvals=[rbpdistancesMvals;rbp(aa).distancesMvals];
    rbpmeanVel = rbpdistancesMvals./movieDuration;
    
end

unc13MeanCurve=[];
unc13distanceStepMvals=[];
unc13distanceStep=[];
unc13distancesMvals=[];
unc13mvals=[];
unc13turningAngles = [];
unc13distanceStepMean = [];
unc13velocityMean = [];
minTurningAngle = [];
for aa = 1:size(unc13,2)
    for subInd = 1:size(unc13(aa).meanCurv,1)
    unc13MeanCurve = [unc13MeanCurve;unc13(aa).meanCurv{subInd}];
    unc13distanceStepMvals = [unc13distanceStepMvals; unc13(aa).distanceStepsMvals{subInd}];
    unc13mvals=[unc13mvals; unc13(aa).mValCell{subInd}];
    unc13turningAngles=[unc13turningAngles; unc13(aa).trajectoryAngles{subInd}];
%     minTurningAngle=[minTurningAngle; min(unc13(aa).trajectoryAngles{subInd})];
    unc13distanceStepMean = [unc13distanceStepMean; nanmean(unc13(aa).distanceStepsMvals{subInd})];
    unc13velocityMean = [unc13velocityMean; (unc13(aa).distanceStepsMvals{subInd})];

    end
    
    unc13distancesMvals=[unc13distancesMvals;unc13(aa).distancesMvals];
    unc13meanVel = unc13distancesMvals./movieDuration;
    
end

%% cdf for curvatures

minThresh = 0.1

controlMeanCurve(controlMeanCurve>minThresh)=[];
gluriiaMeanCurve(gluriiaMeanCurve>minThresh)=[];
rbpMeanCurve(rbpMeanCurve>minThresh)=[];
unc13MeanCurve(unc13MeanCurve>minThresh)=[];

% % gluriia curvature
% h1 = histogram(controlMeanCurve(:));
% hold on
% h2 = histogram(gluriiaMeanCurve(:));
% h1.Normalization = 'probability';
% h1.BinWidth = 0.002;
% h2.Normalization = 'probability';
% h2.BinWidth = 0.002;
% % xlim([0 0.1])

% gluriia

figure('units','normalized','outerposition',[.1 .1 0.8 .6])

subplot(1,2,1)

binEdges = 0.001:0.002:0.1;

[ctData, ~] = histcounts((controlMeanCurve(:)),binEdges,'Normalization','probability');
[expData, ~] = histcounts((gluriiaMeanCurve(:)),binEdges,'Normalization','probability');

xAxis = binEdges;
xAxis(end)=[];
b = bar(xAxis,[ctData' expData'],1.5, 'grouped');

[pVal,h,stats]=ranksum(controlMeanCurve(:),gluriiaMeanCurve(:));
title(['gluriia-RNAi p=' num2str(pVal)]);

b(2).FaceColor = 'r';
set(gca,'fontsize',12,'FontWeight','bold')

xlim([0 0.1])
ylim([0 0.125])
% xtickangle(45)
ylabel('Probability density')
xlabel('Mean body curvature')

subplot(1,2,2)
cdfplot(controlMeanCurve)
hold on
cdfplot(gluriiaMeanCurve)

savefig('gluriia curvature.fig')
saveas(gcf,'gluriia curvature.tif')
saveas(gcf,'gluriia curvature.eps')


save('bodyCurvatureDist.mat','controlMeanCurve','gluriiaMeanCurve')


%%  0-180 deg histograms for turning angles


% gluriia

figure('units','normalized','outerposition',[.1 .1 0.8 .6])
subplot(1,2,1)

binEdges = 0:10:180;
xlabelArray = cell(numel(binEdges),1);
for yy = 1:length(binEdges)
xlabelArray{yy}=num2str(binEdges(yy));
end

[ctData, ~] = histcounts(abs(controlturningAngles),binEdges,'Normalization','probability');
[expData, ~] = histcounts(abs(gluriiaturningAngles),binEdges,'Normalization','probability');

[pVal,h,stats]=ranksum(abs(controlturningAngles),abs(gluriiaturningAngles));

xAxis = 1:numel(ctData);
b = bar(xAxis,[ctData' expData'],1.5, 'grouped');
b(2).FaceColor = 'r';
xticks([1 9 18])
xticklabels({'0','90','180'})
set(gca,'XTickLabel',{'0','90','180'},'fontsize',12,'FontWeight','bold')

xtickangle(45)
ylabel('Probability density')
title(['gluriia-RNAi p=' num2str(pVal) 'Turning Angle']);
ylim([0 0.35])

saveas(gca,'gluriia turning angles 180.tif')
saveas(gca,'gluriia turning angles 180.eps')
saveas(gca,'gluriia turning angles 180.fig')

close all


save('turningAngles.mat','controlturningAngles','gluriiaturningAngles')

%%  -180 - 180 deg histograms for turning angles

% gluriia
figure
binEdges = -180:10:180;
xlabelArray = cell(numel(binEdges),1);
for yy = 1:length(binEdges)
xlabelArray{yy}=num2str(binEdges(yy));
end

[ctData, ~] = histcounts((controlturningAngles),binEdges,'Normalization','probability');
[expData, ~] = histcounts((gluriiaturningAngles),binEdges,'Normalization','probability');

[pVal,h,stats]=ranksum((controlturningAngles),(gluriiaturningAngles));

xAxis = 1:numel(ctData);
b = bar(xAxis,[ctData' expData'],1.1, 'grouped');
title('Turning Angle');
b(2).FaceColor = 'r';
xticks([1 19 36])
xticklabels({'-180','0','180'})
set(gca,'fontsize',12,'FontWeight','bold')

xtickangle(45)
ylabel('Probability density')
title('Cac-RNAi Turning Angles');
saveas(gca,'Cac turning angles.tif')
saveas(gca,'Cac turning angles.eps')


save('turningAngles.mat','controlturningAngles','gluriiaturningAngles')

%%
%%%%%%%%%%% total distances
%%%%%%%%%%%%%%%%%%%%%%
% rbp

% controltotDist = [control(:).distancesMvals];
% rbptotDist = [rbp(:).distancesMvals];

% controltotDist=controltotDist(:);
% rbptotDist=rbptotDist(:);
% 
% rbptotDist(isoutlier(rbptotDist))=[];

% meanDistCt = nanmean(controltotDist);
% semDistCt = nanstd(controltotDist/sqrt(length(controltotDist)));
% meanDistExp = nanmean(rbptotDist);
% semDistExp = nanstd(rbptotDist/sqrt(length(rbptotDist)));
% pvalsAmp = ranksum(controltotDist,rbptotDist);

% cc = linspecer(2);
% figure('units','normalized','outerposition',[.1 .1 0.4 .9])
% 
% hold on
% superbar([meanDistCt meanDistExp],'E',[semDistCt semDistExp],'BarFaceColor', 'none', 'BarEdgeColor',cc,'PStarFontSize',35,'BarWidth',0.4)
% scatter(repmat(1.35,length(controltotDist),1),controltotDist,250,cc(1,:),'.','jitter','on', 'jitterAmount',0.1);
% scatter(repmat(2.35,length(rbptotDist),1),rbptotDist,250,cc(2,:),'.','jitter','on', 'jitterAmount',.1);
% set(gca,'fontsize',20)
% 
% set(gca,'XTick',[1 2])
% xticklabels({'Ct (Attp2)','RBP-RNAi'})    
% ylabel('Mean total distance (mm)')
% ylim([0 1200])
% 
% % xtickangle(45)
% 
% savefig('RBP distance.fig')
% saveas(gcf,'RBP distance.tif')
% saveas(gcf,'RBP distance.eps')

% gluriia

controltotDist = [control(:).distancesMvals];
gluriiatotDist = [gluriia(:).distancesMvals];

controltotDist=controltotDist(:);
gluriiatotDist=gluriiatotDist(:);


% meanDistCt = nanmean(controltotDist);
% semDistCt = nanstd(controltotDist/sqrt(length(controltotDist)));
% meanDistExp = nanmean(gluriiatotDist);
% semDistExp = nanstd(gluriiatotDist/sqrt(length(gluriiatotDist)));
% pvalsAmp = ranksum(controltotDist,gluriiatotDist);

% cc = linspecer(2);
% figure('units','normalized','outerposition',[.1 .1 0.4 .9])
% hold on
% superbar([meanDistCt meanDistExp],'E',[semDistCt semDistExp],'BarFaceColor', 'none', 'BarEdgeColor',cc,'PStarFontSize',35,'BarWidth',0.4)
% scatter(repmat(1.35,length(controltotDist),1),controltotDist,250,cc(1,:),'.','jitter','on', 'jitterAmount',0.1);
% scatter(repmat(2.35,length(gluriiatotDist),1),gluriiatotDist,250,cc(2,:),'.','jitter','on', 'jitterAmount',.1);
% set(gca,'fontsize',20)
% 
% set(gca,'XTick',[1 2])
% xticklabels({'Ct (Attp2)','Cac-RNAi'})    
% ylabel('Mean total distance (mm)')
% 
% % xtickangle(45)
% 
% savefig('Cac distance.fig')
% saveas(gcf,'Cac distance.tif')
% saveas(gcf,'Cac distance.eps')

% unc13
% 
% controltotDist = [control(:).distancesMvals];
% unc13totDist = [unc13(:).distancesMvals];
% 
% controltotDist=controltotDist(:);
% unc13totDist=unc13totDist(:);

% rbptotDist(isoutlier(rbptotDist))=[];
% 
% meanDistCt = nanmean(controltotDist);
% semDistCt = nanstd(controltotDist/sqrt(length(controltotDist)));
% meanDistExp = nanmean(unc13totDist);
% semDistExp = nanstd(unc13totDist/sqrt(length(unc13totDist)));
% pvalsAmp = ranksum(controltotDist,unc13totDist);

% cc = linspecer(2);
% figure('units','normalized','outerposition',[.1 .1 0.4 .9])
% hold on
% superbar([meanDistCt meanDistExp],'E',[semDistCt semDistExp],'BarFaceColor', 'none', 'BarEdgeColor',cc,'PStarFontSize',35,'BarWidth',0.4)
% scatter(repmat(1.35,length(controltotDist),1),controltotDist,250,cc(1,:),'.','jitter','on', 'jitterAmount',0.1);
% scatter(repmat(2.35,length(unc13totDist),1),unc13totDist,250,cc(2,:),'.','jitter','on', 'jitterAmount',.1);
% set(gca,'fontsize',20)
% 
% set(gca,'XTick',[1 2])
% xticklabels({'Ct (Attp2)','Unc13-RNAi'})    
% ylabel('Mean total distance (mm)')
% 
% % xtickangle(45)
% 
% savefig('Unc13 distance.fig')
% saveas(gcf,'Unc13 distance.tif')
% saveas(gcf,'Unc13 distance.eps')

save('totalDistances.mat','controltotDist','gluriiatotDist')

%% 


%   VIOLIN PLOT

dataSavedperLarva{1}=controltotDist;
dataSavedperLarva{2}=gluriiatotDist;

catData = [dataSavedperLarva{1}; dataSavedperLarva{2}];
genotypeStr = [repmat({'Controls'},numel(dataSavedperLarva{1}),1);...
    repmat({'GluRIIA'},numel(dataSavedperLarva{2}),1)];

T = table(genotypeStr,catData);

genotypeOrder = {'Controls','GluRIIA'};
genotype = categorical(T.genotypeStr,genotypeOrder);

subplot(2,4,1)
    
boxchart(genotype,catData)  
hold on
scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);

pVals=[1 ranksum(dataSavedperLarva{1},dataSavedperLarva{2})];

ylim([0 ceil(max(catData))+0.1*ceil(max(catData))])

ylabel('Total distance traveled (mm)')

pVals(1) = NaN;
pVals(2) = ranksum(dataSavedperLarva{1},dataSavedperLarva{2});

xlabel(num2str(pVals))
%         xlim([-0.05,5.05])

set(gca, 'FontName', 'Arial')
set(gca,'fontsize',10)


%% Total distances / velocity

analyzeRNAis = 1;
recordingLengthInMinutes = 10;

figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9])
subplot(2,3,1)
    
% velocity
    dataSavedperLarva{1}=controltotDist./(recordingLengthInMinutes*60);
    dataSavedperLarva{2}=gluriiatotDist./(recordingLengthInMinutes*60);

    boxDataOutlierRemoved = [];
    outlierVals = [];
    
    boxNames = [];
    outLierData = {};
    clear nPoints boxDataOutlierRemoved
    for ww = 1:numel(dataSavedperLarva)
        nPoints(ww,1) = numel(dataSavedperLarva{ww});
    end

    maxLength = max(nPoints);
    
    for ww = 1:numel(dataSavedperLarva)
        
        theseData = dataSavedperLarva{ww};
        
        thisOutlierInd =(isoutlier(theseData,'mean'));
        
        numOutliers = sum(thisOutlierInd);
        
        if numOutliers>0
            outLierData{ww,1}=theseData(thisOutlierInd);
        else
            outLierData{ww,1}=[];
        end
        
        theseDataOutlierRemoved = theseData;
        theseDataOutlierRemoved(thisOutlierInd)=[];

        pointBuffer = maxLength-numel(theseDataOutlierRemoved);
        nanBuffer = repmat([NaN],pointBuffer,1);
        data2cat = [theseDataOutlierRemoved; nanBuffer];
        boxDataOutlierRemoved(:,ww) = data2cat'; % for violin plot
        

    end
    
  

    boxDataOutlierRemoved = abs(boxDataOutlierRemoved);
    x = 1:numel(dataSavedperLarva);
    
   dataSavedperLarva{1}=boxDataOutlierRemoved(~isnan(boxDataOutlierRemoved(:,1)),1);
   dataSavedperLarva{2}=boxDataOutlierRemoved(~isnan(boxDataOutlierRemoved(:,2)),2);

    catData = [dataSavedperLarva{1}; dataSavedperLarva{2}];
    genotypeStr = [repmat({'Controls'},numel(dataSavedperLarva{1}),1);...
        repmat({'GluRIIA'},numel(dataSavedperLarva{2}),1)];

    T = table(genotypeStr,catData);

    genotypeOrder = {'Controls','GluRIIA'};
    genotype = categorical(T.genotypeStr,genotypeOrder);

    boxchart(genotype,catData)  
    hold on
    scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);


    set(gca, 'FontName', 'Arial')
    set(gca,'fontsize',10)
    
    ranksumResults(1) = NaN;
    ranksumResults(2) = ranksum(dataSavedperLarva{1},dataSavedperLarva{2});
   
    newNames = {};
    for zz = 1:length(ranksumResults)
        thisGene = [num2str(ranksumResults(zz))];
        newNames{zz}=thisGene;
    end
    
 mean_X = mean(dataSavedperLarva{1});       % Compute mean
        std_X = std(dataSavedperLarva{1});         % Compute standard deviation
        CV1 = (std_X / mean_X) * 100; 
        
        mean_X = mean(dataSavedperLarva{2});       % Compute mean
        std_X = std(dataSavedperLarva{2});         % Compute standard deviation
        CV2 = (std_X / mean_X) * 100; 
        
   

    
    
% xlabel(controlfns{fieldNum})
    
%     set(gca,'XTick',[1:length(newNames)])
%         set(gca,'xtick',[])
    set(gca,'XTickLabel',newNames)
    xticklabels(newNames')
%     xlim([-0.05,5.05])
%     ylim([0 180])
    hLeg = legend('example');
    set(hLeg,'visible','off')
    ylabel('Velocity (mm/s)');
    

% DISTANCE STEP LENGTH
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,2)


% gluriia

controltotDistStep = [control(:).distanceStepsMvals];
gluriiatotDistStep = [gluriia(:).distanceStepsMvals];

controltotDistStep=controltotDistStep(:);
gluriiatotDistStep=gluriiatotDistStep(:);

clear controldistanceStep
qqCnt = 1;
for qq = 1:length(controltotDistStep)
    
    if ~isempty(controltotDistStep{qq})
        controldistanceStep(qqCnt,1)= mean(controltotDistStep{qq});
            qqCnt =qqCnt+ 1;
    end
end
controldistanceStep(isoutlier(controldistanceStep))=[];

clear gluriiadistanceStep
qqCnt = 1;
for qq = 1:length(gluriiatotDistStep)
    
    if ~isempty(gluriiatotDistStep{qq})
        gluriiadistanceStep(qqCnt,1)= mean(gluriiatotDistStep{qq});
            qqCnt =qqCnt+ 1;
    end
end

gluriiadistanceStep(isoutlier(gluriiadistanceStep))=[];


save('distanceSteps.mat','controldistanceStep','gluriiadistanceStep');

% Distance Step  VIOLIN PLOT

controldistanceStep(controldistanceStep==0)=[];
gluriiadistanceStep(gluriiadistanceStep==0)=[];

dataSavedperLarva{1}=controldistanceStep;
dataSavedperLarva{2}=gluriiadistanceStep;


    boxDataOutlierRemoved = [];
    outlierVals = [];
    
    boxNames = [];
    outLierData = {};
    clear nPoints boxDataOutlierRemoved
    for ww = 1:numel(dataSavedperLarva)
        nPoints(ww,1) = numel(dataSavedperLarva{ww});
    end

    maxLength = max(nPoints);
    
    for ww = 1:numel(dataSavedperLarva)
        
        theseData = dataSavedperLarva{ww};
        
        thisOutlierInd =(isoutlier(theseData,'median'));
        
        numOutliers = sum(thisOutlierInd);
        
        if numOutliers>0
            outLierData{ww,1}=theseData(thisOutlierInd);
        else
            outLierData{ww,1}=[];
        end
        
        theseDataOutlierRemoved = theseData;
        theseDataOutlierRemoved(thisOutlierInd)=[];

        pointBuffer = maxLength-numel(theseDataOutlierRemoved);
        nanBuffer = repmat([NaN],pointBuffer,1);
        data2cat = [theseDataOutlierRemoved; nanBuffer];
        boxDataOutlierRemoved(:,ww) = data2cat'; % for violin plot
        

    end
    
        boxDataOutlierRemoved = abs(boxDataOutlierRemoved);
        x = 1:numel(dataSavedperLarva);

    dataSavedperLarva{1}=boxDataOutlierRemoved(~isnan(boxDataOutlierRemoved(:,1)),1);
    dataSavedperLarva{2}=boxDataOutlierRemoved(~isnan(boxDataOutlierRemoved(:,2)),2);

   
    catData = [dataSavedperLarva{1}; dataSavedperLarva{2}];
    genotypeStr = [repmat({'Controls'},numel(dataSavedperLarva{1}),1);...
        repmat({'GluRIIA'},numel(dataSavedperLarva{2}),1)];

    T = table(genotypeStr,catData);

    genotypeOrder = {'Controls','GluRIIA'};
    genotype = categorical(T.genotypeStr,genotypeOrder);

    boxchart(genotype,catData)  
    hold on
    scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);

    set(gca, 'FontName', 'Arial')
    set(gca,'fontsize',10)
    
    ranksumResults(1) = NaN;
    ranksumResults(2) = ranksum(dataSavedperLarva{1},dataSavedperLarva{2});
 
    newNames = {};
    for zz = 1:length(ranksumResults)
        thisGene = [num2str(ranksumResults(zz))];
        newNames{zz}=thisGene;
    end
    

% xlabel(controlfns{fieldNum})
    
%     set(gca,'XTick',[1:length(newNames)])
%         set(gca,'xtick',[])
    set(gca,'XTickLabel',newNames)
    xticklabels(newNames')
%     xlim([-0.05,5.05])
%     ylim([0 180])
    hLeg = legend('example');
    set(hLeg,'visible','off')

    ylabel('Distance step length')
    
%%%%%%%%%%% fraction reorientation time
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

% gluriia

controlfractTimeReor = [control(:).fractTimeReorientation];
gluriiafractTimeReor = [gluriia(:).fractTimeReorientation];

controlfractTimeReor=controlfractTimeReor(:);
gluriiafractTimeReor=gluriiafractTimeReor(:);

% gluriiafractTimeReor(isoutlier(gluriiafractTimeReor))=[];

save('fractionTimeReorienting.mat','controlfractTimeReor','gluriiafractTimeReor')

%   VIOLIN PLOT

subplot(2,3,3)

dataSavedperLarva{1}=controlfractTimeReor;
dataSavedperLarva{2}=gluriiafractTimeReor;

      boxDataOutlierRemoved = [];
    outlierVals = [];
    
    boxNames = [];
    outLierData = {};
    clear nPoints boxDataOutlierRemoved
    for ww = 1:numel(dataSavedperLarva)
        nPoints(ww,1) = numel(dataSavedperLarva{ww});
    end

    maxLength = max(nPoints);
    
    for ww = 1:numel(dataSavedperLarva)
        
        theseData = dataSavedperLarva{ww};
        
        thisOutlierInd =(isoutlier(theseData,'mean'));
        
        numOutliers = sum(thisOutlierInd);
        
        if numOutliers>0
            outLierData{ww,1}=theseData(thisOutlierInd);
        else
            outLierData{ww,1}=[];
        end
        
        theseDataOutlierRemoved = theseData;
        theseDataOutlierRemoved(thisOutlierInd)=[];

        pointBuffer = maxLength-numel(theseDataOutlierRemoved);
        nanBuffer = repmat([NaN],pointBuffer,1);
        data2cat = [theseDataOutlierRemoved; nanBuffer];
        boxDataOutlierRemoved(:,ww) = data2cat'; % for violin plot
        
    end
    

        boxDataOutlierRemoved = abs(boxDataOutlierRemoved);
        
        x = 1:numel(dataSavedperLarva);

    dataSavedperLarva{1}=boxDataOutlierRemoved(~isnan(boxDataOutlierRemoved(:,1)),1);
    dataSavedperLarva{2}=boxDataOutlierRemoved(~isnan(boxDataOutlierRemoved(:,2)),2);

    catData = [dataSavedperLarva{1}; dataSavedperLarva{2}];
    genotypeStr = [repmat({'Controls'},numel(dataSavedperLarva{1}),1);...
        repmat({'GluRIIA'},numel(dataSavedperLarva{2}),1)];
    
    T = table(genotypeStr,catData);
    genotypeOrder = {'Controls','GluRIIA'};
    genotype = categorical(T.genotypeStr,genotypeOrder);
    
    boxchart(genotype,catData)  
    hold on
    scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);

    set(gca, 'FontName', 'Arial')
    set(gca,'fontsize',10)
    
    ranksumResults(1) = NaN;
    ranksumResults(2) = ranksum(dataSavedperLarva{1},dataSavedperLarva{2});

    newNames = {};
    for zz = 1:length(ranksumResults)
        thisGene = [num2str(ranksumResults(zz))];
        newNames{zz}=thisGene;
    end
    

% xlabel(controlfns{fieldNum})
    
%     set(gca,'XTick',[1:length(newNames)])
%         set(gca,'xtick',[])
    set(gca,'XTickLabel',newNames)
    xticklabels(newNames')
%     xlim([-0.05,5.05])
%     ylim([0 180])
    hLeg = legend('example');
    set(hLeg,'visible','off')

    ylabel('Fraction time reorienting')

    
%%%%%%%%%%% time spent reorienting
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

% gluriia

controlreorBoutDur = [control(:).reorBoutDur];
gluriiareorBoutDur = [gluriia(:).reorBoutDur];

controlreorBoutDur=controlreorBoutDur(:);
gluriiareorBoutDur=gluriiareorBoutDur(:);

clear controlreorBoutDuration
qqCnt = 1;
for qq = 1:length(controlreorBoutDur)
    
    if ~isempty(controlreorBoutDur{qq})
        controlreorBoutDuration(qqCnt,1)= mean(controlreorBoutDur{qq});
            qqCnt =qqCnt+ 1;
    end
end

clear gluriiareorBoutDuration
qqCnt = 1;
for qq = 1:length(gluriiareorBoutDur)
    
    if ~isempty(gluriiareorBoutDur{qq})
        gluriiareorBoutDuration(qqCnt,1)= mean(gluriiareorBoutDur{qq});
            qqCnt =qqCnt+ 1;
    end
end

controlreorBoutDuration = controlreorBoutDuration*(1/fps);
gluriiareorBoutDuration = gluriiareorBoutDuration*(1/fps);

save('timeSpentReorienting.mat','controlreorBoutDuration','gluriiareorBoutDuration')


%   VIOLIN PLOT

subplot(2,3,4)

dataSavedperLarva{1}=controlreorBoutDuration;
dataSavedperLarva{2}=gluriiareorBoutDuration;

boxNames = [];
clear nPoints
for ww = 1:numel(dataSavedperLarva)
    nPoints(ww,1) = numel(dataSavedperLarva{ww});
end
        
        maxLength = max(nPoints);
        boxData = [];
        for ww = 1:numel(dataSavedperLarva)
            theseData = dataSavedperLarva{ww};
            pointBuffer = maxLength-numel(theseData);
            nanBuffer = repmat([NaN],pointBuffer,1);
            data2cat = [theseData; nanBuffer];
            boxData(:,ww) = data2cat;
        end
        
    boxData = abs(boxData);
    x = 1:numel(dataSavedperLarva);

    dataSavedperLarva{1}=boxData(~isnan(boxData(:,1)),1);
    dataSavedperLarva{2}=boxData(~isnan(boxData(:,2)),2);

    
    catData = [dataSavedperLarva{1}; dataSavedperLarva{2}];
    genotypeStr = [repmat({'Controls'},numel(dataSavedperLarva{1}),1);...
    repmat({'GluRIIA'},numel(dataSavedperLarva{2}),1)];
    
    T = table(genotypeStr,catData);
    genotypeOrder = {'Controls','GluRIIA'};
    genotype = categorical(T.genotypeStr,genotypeOrder);
    
    boxchart(genotype,catData)  
    hold on
    scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);

        
ylabel('Total time reorienting')

pVals(1) = NaN;
pVals(2) = ranksum(dataSavedperLarva{1},dataSavedperLarva{2});

     
% xlim([-0.05,5.05])
xlabel(num2str(pVals))
ylim([0 30])
set(gca, 'FontName', 'Arial')
set(gca,'fontsize',10)

ylabel('Time spent reorienting')


%%%%%%%%%%% time spent during continuous trajectory path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% gluriia

controltimeContTraj = [control(:).contPathDur];
gluriiaTimeContTraj = [gluriia(:).contPathDur];

controltimeContTraj=controltimeContTraj(:);
gluriiaTimeContTraj=gluriiaTimeContTraj(:);

clear controlcontinuousTimeOnTraj
qqCnt = 1;
for qq = 1:length(controltimeContTraj)
    
    if ~isempty(controltimeContTraj{qq})
        controlcontinuousTimeOnTraj(qqCnt,1)= mean(controltimeContTraj{qq});
            qqCnt =qqCnt+ 1;
    end
end

clear gluriiaContinuousTimeOnTraj
qqCnt = 1;
for qq = 1:length(gluriiaTimeContTraj)
    
    if ~isempty(gluriiaTimeContTraj{qq})
        gluriiaContinuousTimeOnTraj(qqCnt,1)= mean(gluriiaTimeContTraj{qq});
            qqCnt =qqCnt+ 1;
    end
end

controlcontinuousTimeOnTraj = controlcontinuousTimeOnTraj*(1/fps);
gluriiaContinuousTimeOnTraj = gluriiaContinuousTimeOnTraj*(1/fps);


save('timeOnContinuousTrajectories.mat','controlcontinuousTimeOnTraj','gluriiaContinuousTimeOnTraj')

%   VIOLIN PLOT

subplot(2,3,5)

dataSavedperLarva{1}=controlcontinuousTimeOnTraj;
dataSavedperLarva{2}=gluriiaContinuousTimeOnTraj;


boxNames = [];
clear nPoints
for ww = 1:numel(dataSavedperLarva)
    nPoints(ww,1) = numel(dataSavedperLarva{ww});
end
        
        maxLength = max(nPoints);
        boxData = [];
        for ww = 1:numel(dataSavedperLarva)
            theseData = dataSavedperLarva{ww};
            pointBuffer = maxLength-numel(theseData);
            nanBuffer = repmat([NaN],pointBuffer,1);
            data2cat = [theseData; nanBuffer];
            boxData(:,ww) = data2cat;
        end
        
        boxData = abs(boxData);
        x = 1:numel(dataSavedperLarva);
        
        
    dataSavedperLarva{1}=boxData(~isnan(boxData(:,1)),1);
    dataSavedperLarva{2}=boxData(~isnan(boxData(:,2)),2);

% figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9])
subplot(2,3,5)

    catData = [dataSavedperLarva{1}; dataSavedperLarva{2}];
    genotypeStr = [repmat({'Controls'},numel(dataSavedperLarva{1}),1);...
        repmat({'GluRIIA'},numel(dataSavedperLarva{2}),1)];
    
    T = table(genotypeStr,catData);
    genotypeOrder = {'Controls','GluRIIA'};
    genotype = categorical(T.genotypeStr,genotypeOrder);
    
    boxchart(genotype,catData)  
    hold on
    scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);


ylabel('continuous time on traj path')

pVals(1) = NaN;
pVals(2) = ranksum(dataSavedperLarva{1},dataSavedperLarva{2});

%         xlim([-0.05,5.05])

xlabel(num2str(pVals))

set(gca, 'FontName', 'Arial')
set(gca,'fontsize',10)
ylim([0 30])


%%%%%%%%%%% Velocity during continuous trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,6)

minLengthInSeconds4contTraj = 2; % a trajectory has to be at least this long in order to use it in the velocity calculation
%control
for ff = 1:length(control)
    
    
    thisPlate = control(ff).meanCurv;
    thisPlateMvals = control(ff).mValCell;
    thisPlateCoords = control(ff).coordinates;

    thisPlateDistSteps = control(ff).distanceStepsMvals;
    
    clear contTrajVelocity
    contTrajVelocity=[];
    for wellN = 1:length(thisPlate)
        thisCurvature=thisPlate{wellN};
        if isempty(thisCurvature)
            contTrajVelocity(wellN,1)=NaN;
            continue
        end
        mVals = thisPlateMvals{wellN};
        x = thisPlateCoords{wellN}(:,1);
        y = thisPlateCoords{wellN}(:,2);
        
        if x==0
            contTrajVelocity(wellN,1)=NaN;
            continue
        end

        thisCurveSmoothed = smooth(thisCurvature,20);
        [~,~,bins] = histcounts(thisCurveSmoothed, 100);
        binDivLow = bins<=75;
        boundLow = mean(thisCurveSmoothed(binDivLow));
        highCurveInds = thisCurveSmoothed>boundLow;
        highCurveIndsTemp = ~highCurveInds;
        curveProps = regionprops(highCurveIndsTemp,'Area','PixelIdxList');
        curveAreas = [curveProps(:).Area];

        inds2merge = curveAreas<periodsOfHighCurvatureMinLength*fps;

        pixInds = find(inds2merge);
        for zz = 1:numel(pixInds)
            thisSpan = curveProps(pixInds(zz)).PixelIdxList;
            highCurveInds(thisSpan) = true;
        end
    
%         plot(thisCurveSmoothed)
%         hold on
%         plot(highCurveInds*0.01)
%         
        contTrajPeriods = ~highCurveInds;
        contTrajPeriodsProps = regionprops(contTrajPeriods,'Area','PixelIdxList');
        contAreas = [contTrajPeriodsProps(:).Area];
        contAreasInSec = contAreas/fps;
        indsBelowThresh = contAreasInSec<minLengthInSeconds4contTraj;
        contTrajPeriodsProps(indsBelowThresh)=[];
        
        clear meanVel
        for trajReg = 1:length(contTrajPeriodsProps)
            pixIdx = contTrajPeriodsProps(trajReg).PixelIdxList;
            xVals = x(pixIdx);
            yVals = y(pixIdx);
%              plot(xVals,yVals)
%              hold on
            thisDist = sqrt(diff(xVals).^2+diff(yVals).^2).*sf;
            thisDist = sum(thisDist);
            trajTime = length(xVals)/fps;
            meanVel(trajReg,1) = thisDist/trajTime;
        end
        
        contTrajVelocity(wellN,1) = mean(meanVel);
        
        
        
    end
    
    control(ff).contTrajVelocity=contTrajVelocity;
    
end

%gluriia
for ff = 1:length(gluriia)
    
    
    thisPlate = gluriia(ff).meanCurv;
    thisPlateMvals = gluriia(ff).mValCell;
    thisPlateCoords = gluriia(ff).coordinates;

    thisPlateDistSteps = gluriia(ff).distanceStepsMvals;
    
    clear contTrajVelocity
    contTrajVelocity=[];
    for wellN = 1:length(thisPlate)
        thisCurvature=thisPlate{wellN};
        if isempty(thisCurvature)
            contTrajVelocity(wellN,1)=NaN;
            continue
        end
        mVals = thisPlateMvals{wellN};
        x = thisPlateCoords{wellN}(:,1);
        y = thisPlateCoords{wellN}(:,2);
        
        if x==0
            contTrajVelocity(wellN,1)=NaN;
            continue
        end

        thisCurveSmoothed = smooth(thisCurvature,20);
        [~,~,bins] = histcounts(thisCurveSmoothed, 100);
        binDivLow = bins<=75;
        boundLow = mean(thisCurveSmoothed(binDivLow));
        highCurveInds = thisCurveSmoothed>boundLow;
        highCurveIndsTemp = ~highCurveInds;
        curveProps = regionprops(highCurveIndsTemp,'Area','PixelIdxList');
        curveAreas = [curveProps(:).Area];

        inds2merge = curveAreas<periodsOfHighCurvatureMinLength*fps;

        pixInds = find(inds2merge);
        for zz = 1:numel(pixInds)
            thisSpan = curveProps(pixInds(zz)).PixelIdxList;
            highCurveInds(thisSpan) = true;
        end
    
%         plot(thisCurveSmoothed)
%         hold on
%         plot(highCurveInds*0.01)
%         
        contTrajPeriods = ~highCurveInds;
        contTrajPeriodsProps = regionprops(contTrajPeriods,'Area','PixelIdxList');
        contAreas = [contTrajPeriodsProps(:).Area];
        contAreasInSec = contAreas/fps;
        indsBelowThresh = contAreasInSec<minLengthInSeconds4contTraj;
        contTrajPeriodsProps(indsBelowThresh)=[];
        
        clear meanVel
        for trajReg = 1:length(contTrajPeriodsProps)
            pixIdx = contTrajPeriodsProps(trajReg).PixelIdxList;
            xVals = x(pixIdx);
            yVals = y(pixIdx);
%              plot(xVals,yVals)
%              hold on
            thisDist = sqrt(diff(xVals).^2+diff(yVals).^2).*sf;
            thisDist = sum(thisDist);
            trajTime = length(xVals)/fps;
            meanVel(trajReg,1) = thisDist/trajTime;
        end
        
        contTrajVelocity(wellN,1) = mean(meanVel);
        
        
        
    end
    
    gluriia(ff).contTrajVelocity=contTrajVelocity;
    
end


controltimeContTraj = [control(:).contTrajVelocity];
controltimeContTraj=controltimeContTraj(:);

gluriiaTimeContTraj = [gluriia(:).contTrajVelocity];
gluriiaTimeContTraj=gluriiaTimeContTraj(:);

%   VIOLIN PLOT

dataSavedperLarva{1}=controltimeContTraj;
dataSavedperLarva{2}=gluriiaTimeContTraj;


boxNames = [];
clear nPoints
for ww = 1:numel(dataSavedperLarva)
    nPoints(ww,1) = numel(dataSavedperLarva{ww});
end
        
        maxLength = max(nPoints);
        boxData = [];
        for ww = 1:numel(dataSavedperLarva)
            theseData = dataSavedperLarva{ww};
            pointBuffer = maxLength-numel(theseData);
            nanBuffer = repmat([NaN],pointBuffer,1);
            data2cat = [theseData; nanBuffer];
            boxData(:,ww) = data2cat;
        end
        
        boxData = abs(boxData);
        x = 1:numel(dataSavedperLarva);
        
    dataSavedperLarva{1}=boxData(~isnan(boxData(:,1)),1);
    dataSavedperLarva{2}=boxData(~isnan(boxData(:,2)),2);
 
    catData = [dataSavedperLarva{1}; dataSavedperLarva{2}];
    genotypeStr = [repmat({'Controls'},numel(dataSavedperLarva{1}),1);...
    repmat({'GluRIIA'},numel(dataSavedperLarva{2}),1)];
    
    T = table(genotypeStr,catData);
    genotypeOrder = {'Controls','GluRIIA'};
    genotype = categorical(T.genotypeStr,genotypeOrder);
    
    boxchart(genotype,catData)  
    hold on
    scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);

        
        

ylabel('continuous time on traj path')

pVals(1) = NaN;
pVals(2) = ranksum(dataSavedperLarva{1},dataSavedperLarva{2});



xlabel(num2str(pVals))

set(gca, 'FontName', 'Arial')
set(gca,'fontsize',10)
ylim([0 1])

