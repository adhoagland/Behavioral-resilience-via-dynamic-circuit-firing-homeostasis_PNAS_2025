load('RNAiBehavioralData.mat')

attp2MeanCurve=[];
attp2distanceStepMvals=[];
attp2distanceStep=[];
attp2distancesMvals=[];
attp2mvals=[];
attp2turningAngles = [];
attp2distanceStepMean = [];
attp2velocityMean = [];
minTurningAngle = [];
for aa = 1:size(attp2,2)
    for subInd = 1:size(attp2(aa).meanCurv,1)
    attp2MeanCurve = [attp2MeanCurve;attp2(aa).meanCurv{subInd}];
    attp2distanceStepMvals = [attp2distanceStepMvals; attp2(aa).distanceStepsMvals{subInd}];
    attp2mvals=[attp2mvals; attp2(aa).mValCell{subInd}];
    attp2turningAngles=[attp2turningAngles; attp2(aa).trajectoryAngles{subInd}];
%     minTurningAngle=[minTurningAngle; min(attp2(aa).trajectoryAngles{subInd})];
    attp2distanceStepMean = [attp2distanceStepMean; nanmean(attp2(aa).distanceStepsMvals{subInd})];
    attp2velocityMean = [attp2velocityMean; (attp2(aa).distanceStepsMvals{subInd})];

    end
    
    attp2distancesMvals=[attp2distancesMvals;attp2(aa).distancesMvals];
    attp2meanVel = attp2distancesMvals./movieDuration;
    
end

cacMeanCurve=[];
cacdistanceStepMvals=[];
cacdistanceStep=[];
cacdistancesMvals=[];
cacmvals=[];
cacturningAngles = [];
cacdistanceStepMean = [];
cacvelocityMean = [];
minTurningAngle = [];
for aa = 1:size(cac,2)
    for subInd = 1:size(cac(aa).meanCurv,1)
    cacMeanCurve = [cacMeanCurve;cac(aa).meanCurv{subInd}];
    cacdistanceStepMvals = [cacdistanceStepMvals; cac(aa).distanceStepsMvals{subInd}];
    cacmvals=[cacmvals; cac(aa).mValCell{subInd}];
    cacturningAngles=[cacturningAngles; cac(aa).trajectoryAngles{subInd}];
%     minTurningAngle=[minTurningAngle; min(cac(aa).trajectoryAngles{subInd})];
    cacdistanceStepMean = [cacdistanceStepMean; nanmean(cac(aa).distanceStepsMvals{subInd})];
    cacvelocityMean = [cacvelocityMean; (cac(aa).distanceStepsMvals{subInd})];

    end
    
    cacdistancesMvals=[cacdistancesMvals;cac(aa).distancesMvals];
    cacmeanVel = cacdistancesMvals./movieDuration;
    
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

attp2MeanCurve(attp2MeanCurve>minThresh)=[];
cacMeanCurve(cacMeanCurve>minThresh)=[];
rbpMeanCurve(rbpMeanCurve>minThresh)=[];
unc13MeanCurve(unc13MeanCurve>minThresh)=[];

% % cac curvature
% h1 = histogram(attp2MeanCurve(:));
% hold on
% h2 = histogram(cacMeanCurve(:));
% h1.Normalization = 'probability';
% h1.BinWidth = 0.002;
% h2.Normalization = 'probability';
% h2.BinWidth = 0.002;
% % xlim([0 0.1])

% cac

figure('units','normalized','outerposition',[.1 .1 0.8 .6])

subplot(1,2,1)

binEdges = 0.001:0.002:0.1;

[ctData, ~] = histcounts((attp2MeanCurve(:)),binEdges,'Normalization','probability');
[expData, ~] = histcounts((cacMeanCurve(:)),binEdges,'Normalization','probability');

xAxis = binEdges;
xAxis(end)=[];
b = bar(xAxis,[ctData' expData'],1.5, 'grouped');

[pVal,h,stats]=ranksum(attp2MeanCurve(:),cacMeanCurve(:));
title(['cac-RNAi p=' num2str(pVal)]);

b(2).FaceColor = 'r';
set(gca,'fontsize',12,'FontWeight','bold')

xlim([0 0.1])
ylim([0 0.125])
% xtickangle(45)
ylabel('Probability density')
xlabel('Mean body curvature')

subplot(1,2,2)
cdfplot(attp2MeanCurve)
hold on
cdfplot(cacMeanCurve)

savefig('cac curvature.fig')
saveas(gcf,'cac curvature.tif')
saveas(gcf,'cac curvature.eps')

% rbp

figure('units','normalized','outerposition',[.1 .1 0.8 .6])

subplot(1,2,1)

binEdges = 0.001:0.002:0.1;

[ctData, ~] = histcounts((attp2MeanCurve(:)),binEdges,'Normalization','probability');
[expData, ~] = histcounts((rbpMeanCurve(:)),binEdges,'Normalization','probability');

xAxis = binEdges;
xAxis(end)=[];
b = bar(xAxis,[ctData' expData'],1.5, 'grouped');

[pVal,h,stats]=ranksum(attp2MeanCurve(:),rbpMeanCurve(:));
title(['RBP-RNAi p=' num2str(pVal)]);
b(2).FaceColor = 'r';


set(gca,'fontsize',12,'FontWeight','bold')

xlim([0 0.1])
ylim([0 0.125])
% xtickangle(45)
ylabel('Probability density')
xlabel('Mean body curvature')

subplot(1,2,2)
cdfplot(attp2MeanCurve)
hold on
cdfplot(rbpMeanCurve)

savefig('rbp curvature.fig')
saveas(gcf,'rbp curvature.tif')
saveas(gcf,'rbp curvature.eps')

% unc13

figure('units','normalized','outerposition',[.1 .1 0.8 .6])

subplot(1,2,1)

binEdges = 0.001:0.002:0.1;

[ctData, ~] = histcounts((attp2MeanCurve(:)),binEdges,'Normalization','probability');
[expData, ~] = histcounts((unc13MeanCurve(:)),binEdges,'Normalization','probability');

xAxis = binEdges;
xAxis(end)=[];
b = bar(xAxis,[ctData' expData'],1.5, 'grouped');
[pVal,h,stats]=ranksum(attp2MeanCurve(:),unc13MeanCurve(:));
title(['Unc13-RNAi p=' num2str(pVal)]);
b(2).FaceColor = 'r';
set(gca,'fontsize',12,'FontWeight','bold')

xlim([0 0.1])
ylim([0 0.125])
% xtickangle(45)
ylabel('Probability density')
xlabel('Mean body curvature')

subplot(1,2,2)
cdfplot(attp2MeanCurve)
hold on
cdfplot(unc13MeanCurve)

savefig('unc13 curvature.fig')
saveas(gcf,'unc13 curvature.tif')
saveas(gcf,'unc13 curvature.eps')


save('bodyCurvatureDist.mat','attp2MeanCurve','rbpMeanCurve','unc13MeanCurve','cacMeanCurve')


%%  0-180 deg histograms for turning angles

% rbp
% figure('units','normalized','outerposition',[.1 .1 0.8 .6])


figure('units','normalized','outerposition',[.1 .1 0.8 .6])
subplot(1,2,1)

binEdges = 0:10:180;
xlabelArray = cell(numel(binEdges),1);
for yy = 1:length(binEdges)
xlabelArray{yy}=num2str(binEdges(yy));
end

[ctData, ~] = histcounts(abs(attp2turningAngles),binEdges,'Normalization','probability');
[expData, ~] = histcounts(abs(rbpturningAngles),binEdges,'Normalization','probability');

[pVal,h,stats]=ranksum(abs(attp2turningAngles),abs(rbpturningAngles));


xAxis = 1:numel(ctData);
b = bar(xAxis,[ctData' expData'],1.5, 'grouped');
b(2).FaceColor = 'r';
xticks([1 9 18])
xticklabels({'0','90','180'})
set(gca,'XTickLabel',{'0','90','180'},'fontsize',12,'FontWeight','bold')

xtickangle(45)
ylabel('Probability density')

xtickangle(45)
ylabel('Probability density')
title(['RBP-RNAi p=' num2str(pVal) 'Turning Angle']);
ylim([0 0.35])

saveas(gca,'RBP turning angles 180.tif')
saveas(gca,'RBP turning angles 180.eps')
saveas(gca,'RBP turning angles 180.fig')

close all

% cac

figure('units','normalized','outerposition',[.1 .1 0.8 .6])
subplot(1,2,1)

binEdges = 0:10:180;
xlabelArray = cell(numel(binEdges),1);
for yy = 1:length(binEdges)
xlabelArray{yy}=num2str(binEdges(yy));
end

[ctData, ~] = histcounts(abs(attp2turningAngles),binEdges,'Normalization','probability');
[expData, ~] = histcounts(abs(cacturningAngles),binEdges,'Normalization','probability');

[pVal,h,stats]=ranksum(abs(attp2turningAngles),abs(cacturningAngles));

xAxis = 1:numel(ctData);
b = bar(xAxis,[ctData' expData'],1.5, 'grouped');
b(2).FaceColor = 'r';
xticks([1 9 18])
xticklabels({'0','90','180'})
set(gca,'XTickLabel',{'0','90','180'},'fontsize',12,'FontWeight','bold')

xtickangle(45)
ylabel('Probability density')
title(['cac-RNAi p=' num2str(pVal) 'Turning Angle']);
ylim([0 0.35])

saveas(gca,'cac turning angles 180.tif')
saveas(gca,'cac turning angles 180.eps')
saveas(gca,'cac turning angles 180.fig')

close all

% unc13

figure('units','normalized','outerposition',[.1 .1 0.8 .6])
subplot(1,2,1)

binEdges = 0:10:180;
xlabelArray = cell(numel(binEdges),1);
for yy = 1:length(binEdges)
xlabelArray{yy}=num2str(binEdges(yy));
end

[ctData, ~] = histcounts(abs(attp2turningAngles),binEdges,'Normalization','probability');
[expData, ~] = histcounts(abs(unc13turningAngles),binEdges,'Normalization','probability');

[pVal,h,stats]=ranksum(abs(attp2turningAngles),abs(unc13turningAngles));

xAxis = 1:numel(ctData);
b = bar(xAxis,[ctData' expData'],1.5, 'grouped');
b(2).FaceColor = 'r';
xticks([1 9 18])
xticklabels({'0','90','180'})
set(gca,'XTickLabel',{'0','90','180'},'fontsize',12,'FontWeight','bold')

xtickangle(45)
ylabel('Probability density')
title(['unc13-RNAi p=' num2str(pVal) 'Turning Angle']);
ylim([0 0.35])

saveas(gca,'Unc13 turning angles 180.tif')
saveas(gca,'Unc13 turning angles 180.eps')
saveas(gca,'Unc13 turning angles 180.fig')

save('turningAngles.mat','attp2turningAngles','cacturningAngles','rbpturningAngles','unc13turningAngles')

%%  -180 - 180 deg histograms for turning angles

% rbp
binEdges = -180:10:180;
xlabelArray = cell(numel(binEdges),1);
for yy = 1:length(binEdges)
xlabelArray{yy}=num2str(binEdges(yy));
end

[ctData, ~] = histcounts((attp2turningAngles),binEdges,'Normalization','probability');
[expData, ~] = histcounts((rbpturningAngles),binEdges,'Normalization','probability');

[pVal,h,stats]=ranksum((attp2turningAngles),(rbpturningAngles));

xAxis = 1:numel(ctData);
b = bar(xAxis,[ctData' expData'],1.1, 'grouped');
b(2).FaceColor = 'r';
xticks([1 19 36])
xticklabels({'-180','0','180'})
set(gca,'fontsize',12,'FontWeight','bold')

xtickangle(45)
ylabel('Probability density')
title(['RBP-RNAi p=' num2str(pVal) 'Turning Angle']);
saveas(gca,'RBP turning angles.tif')
saveas(gca,'RBP turning angles.eps')

% cac
figure
binEdges = -180:10:180;
xlabelArray = cell(numel(binEdges),1);
for yy = 1:length(binEdges)
xlabelArray{yy}=num2str(binEdges(yy));
end

[ctData, ~] = histcounts((attp2turningAngles),binEdges,'Normalization','probability');
[expData, ~] = histcounts((cacturningAngles),binEdges,'Normalization','probability');

[pVal,h,stats]=ranksum((attp2turningAngles),(cacturningAngles));

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

% unc13
figure
binEdges = -180:10:180;
xlabelArray = cell(numel(binEdges),1);
for yy = 1:length(binEdges)
xlabelArray{yy}=num2str(binEdges(yy));
end

[ctData, ~] = histcounts((attp2turningAngles),binEdges,'Normalization','probability');
[expData, ~] = histcounts((unc13turningAngles),binEdges,'Normalization','probability');

xAxis = 1:numel(ctData);
b = bar(xAxis,[ctData' expData'],1.1, 'grouped');
title('Turning Angle');
b(2).FaceColor = 'r';
xticks([1 19 36])
xticklabels({'-180','0','180'})
set(gca,'fontsize',12,'FontWeight','bold')

xtickangle(45)
ylabel('Probability density')
title('Unc13-RNAi Turning Angles');
saveas(gca,'Unc13 turning angles.tif')
saveas(gca,'Unc13 turning angles.eps')

save('turningAngles.mat','attp2turningAngles','cacturningAngles','rbpturningAngles','unc13turningAngles')

%%
%%%%%%%%%%% total distances
%%%%%%%%%%%%%%%%%%%%%%
% rbp

attp2totDist = [attp2(:).distancesMvals];
rbptotDist = [rbp(:).distancesMvals];

attp2totDist=attp2totDist(:);
rbptotDist=rbptotDist(:);

rbptotDist(isoutlier(rbptotDist))=[];

% meanDistCt = nanmean(attp2totDist);
% semDistCt = nanstd(attp2totDist/sqrt(length(attp2totDist)));
% meanDistExp = nanmean(rbptotDist);
% semDistExp = nanstd(rbptotDist/sqrt(length(rbptotDist)));
% pvalsAmp = ranksum(attp2totDist,rbptotDist);

% cc = linspecer(2);
% figure('units','normalized','outerposition',[.1 .1 0.4 .9])
% 
% hold on
% superbar([meanDistCt meanDistExp],'E',[semDistCt semDistExp],'BarFaceColor', 'none', 'BarEdgeColor',cc,'PStarFontSize',35,'BarWidth',0.4)
% scatter(repmat(1.35,length(attp2totDist),1),attp2totDist,250,cc(1,:),'.','jitter','on', 'jitterAmount',0.1);
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

% cac

attp2totDist = [attp2(:).distancesMvals];
cactotDist = [cac(:).distancesMvals];

attp2totDist=attp2totDist(:);
cactotDist=cactotDist(:);


% meanDistCt = nanmean(attp2totDist);
% semDistCt = nanstd(attp2totDist/sqrt(length(attp2totDist)));
% meanDistExp = nanmean(cactotDist);
% semDistExp = nanstd(cactotDist/sqrt(length(cactotDist)));
% pvalsAmp = ranksum(attp2totDist,cactotDist);

% cc = linspecer(2);
% figure('units','normalized','outerposition',[.1 .1 0.4 .9])
% hold on
% superbar([meanDistCt meanDistExp],'E',[semDistCt semDistExp],'BarFaceColor', 'none', 'BarEdgeColor',cc,'PStarFontSize',35,'BarWidth',0.4)
% scatter(repmat(1.35,length(attp2totDist),1),attp2totDist,250,cc(1,:),'.','jitter','on', 'jitterAmount',0.1);
% scatter(repmat(2.35,length(cactotDist),1),cactotDist,250,cc(2,:),'.','jitter','on', 'jitterAmount',.1);
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

attp2totDist = [attp2(:).distancesMvals];
unc13totDist = [unc13(:).distancesMvals];

attp2totDist=attp2totDist(:);
unc13totDist=unc13totDist(:);

% rbptotDist(isoutlier(rbptotDist))=[];
% 
% meanDistCt = nanmean(attp2totDist);
% semDistCt = nanstd(attp2totDist/sqrt(length(attp2totDist)));
% meanDistExp = nanmean(unc13totDist);
% semDistExp = nanstd(unc13totDist/sqrt(length(unc13totDist)));
% pvalsAmp = ranksum(attp2totDist,unc13totDist);

% cc = linspecer(2);
% figure('units','normalized','outerposition',[.1 .1 0.4 .9])
% hold on
% superbar([meanDistCt meanDistExp],'E',[semDistCt semDistExp],'BarFaceColor', 'none', 'BarEdgeColor',cc,'PStarFontSize',35,'BarWidth',0.4)
% scatter(repmat(1.35,length(attp2totDist),1),attp2totDist,250,cc(1,:),'.','jitter','on', 'jitterAmount',0.1);
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

save('totalDistances.mat','attp2totDist','cactotDist','rbptotDist','unc13totDist')

%% 

% 
% %   Box PLOT
% 
% % 
% % attp2totDist=attp2totDist(~isnan(attp2totDist));
% % unc13totDist=unc13totDist(~isnan(unc13totDist));
% % rbptotDist=rbptotDist(~isnan(rbptotDist));
% % cactotDist=cactotDist(~isnan(cactotDist));
% 
% dataSavedperLarva{1}=attp2totDist;
% dataSavedperLarva{2}=unc13totDist;
% dataSavedperLarva{3}=rbptotDist;
% dataSavedperLarva{4}=cactotDist;
% 
%  boxNames = [];
% clear nPoints
% for ww = 1:numel(dataSavedperLarva)
%     nPoints(ww,1) = numel(dataSavedperLarva{ww});
% end
%         
%         maxLength = max(nPoints);
%         boxData = [];
%         for ww = 1:numel(dataSavedperLarva)
%             theseData = dataSavedperLarva{ww};
%             pointBuffer = maxLength-numel(theseData);
%             nanBuffer = repmat([NaN],pointBuffer,1);
%             data2cat = [theseData; nanBuffer];
%             boxData(:,ww) = data2cat;
%         end
%         
%         boxData = abs(boxData);
%         x = 1:numel(dataSavedperLarva);
%         
% figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9])
% subplot(2,3,1)
% 
%    catData = [dataSavedperLarva{1}; dataSavedperLarva{2}; dataSavedperLarva{3}; dataSavedperLarva{4}];
%     genotypeStr = [repmat({'Attp2'},numel(dataSavedperLarva{1}),1);...
%         repmat({'Unc13'},numel(dataSavedperLarva{2}),1);...
%         repmat({'RBP'},numel(dataSavedperLarva{3}),1);...
%         repmat({'Cac'},numel(dataSavedperLarva{4}),1)];
%     
%     T = table(genotypeStr,catData);
% 
%     genotypeOrder = {'Attp2','Unc13','RBP','Cac'};
%     genotype = categorical(T.genotypeStr,genotypeOrder);
% 
%     boxchart(genotype,catData)  
%     hold on
%     scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);
% 
% ylabel('Total distance traveled (mm)')
% 
% pVals(1) = NaN;
% pVals(2) = ranksum(dataSavedperLarva{1},dataSavedperLarva{2});
% pVals(3) = ranksum(dataSavedperLarva{1},dataSavedperLarva{3});
% pVals(4) = ranksum(dataSavedperLarva{1},dataSavedperLarva{4});
% 
% xlabel(num2str(pVals))
%         xlim([-0.05,5.05])
% 
% set(gca, 'FontName', 'Arial')
% set(gca,'fontsize',10)


%% Total distances / velocity
clear pVals

analyzeRNAis = 1;
recordingLengthInMinutes = 10;

figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9])
subplot(2,3,1)
    
% velocity
    dataSavedperLarva{1}=attp2totDist./(recordingLengthInMinutes*60);
    dataSavedperLarva{2}=unc13totDist./(recordingLengthInMinutes*60);
    dataSavedperLarva{3}=rbptotDist./(recordingLengthInMinutes*60);
    dataSavedperLarva{4}=cactotDist./(recordingLengthInMinutes*60);

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
        dataSavedperLarva{3}=boxDataOutlierRemoved(~isnan(boxDataOutlierRemoved(:,3)),3);
        dataSavedperLarva{4}=boxDataOutlierRemoved(~isnan(boxDataOutlierRemoved(:,4)),4);

    
        catData = [dataSavedperLarva{1}; dataSavedperLarva{2}; dataSavedperLarva{3}; dataSavedperLarva{4}];
        genotypeStr = [repmat({'Attp2'},numel(dataSavedperLarva{1}),1);...
        repmat({'Unc13'},numel(dataSavedperLarva{2}),1);...
        repmat({'RBP'},numel(dataSavedperLarva{3}),1);...
        repmat({'Cac'},numel(dataSavedperLarva{4}),1)];

        T = table(genotypeStr,catData);

        genotypeOrder = {'Attp2','Unc13','RBP','Cac'};
        genotype = categorical(T.genotypeStr,genotypeOrder);

        boxchart(genotype,catData)  
        hold on
        scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);

        mean_X = mean(dataSavedperLarva{1});       % Compute mean
        std_X = std(dataSavedperLarva{1});         % Compute standard deviation
        CV1 = (std_X / mean_X) * 100; 
        
        mean_X = mean(dataSavedperLarva{2});       % Compute mean
        std_X = std(dataSavedperLarva{2});         % Compute standard deviation
        CV2 = (std_X / mean_X) * 100; 
        
        mean_X = mean(dataSavedperLarva{3});       % Compute mean
        std_X = std(dataSavedperLarva{3});         % Compute standard deviation
        CV3 = (std_X / mean_X) * 100; 
        
        mean_X = mean(dataSavedperLarva{4});       % Compute mean
        std_X = std(dataSavedperLarva{4});         % Compute standard deviation
        CV4 = (std_X / mean_X) * 100; 

    set(gca, 'FontName', 'Arial')
    set(gca,'fontsize',10)
    
    ranksumResults(1) = NaN;
    ranksumResults(2) = ranksum(dataSavedperLarva{1},dataSavedperLarva{2});
    ranksumResults(3) = ranksum(dataSavedperLarva{1},dataSavedperLarva{3});
    ranksumResults(4) = ranksum(dataSavedperLarva{1},dataSavedperLarva{4});

    newNames = {};
    for zz = 1:length(ranksumResults)
        thisGene = [num2str(ranksumResults(zz))];
        newNames{zz}=thisGene;
    end
    
% xlabel(attp2fns{fieldNum})
    
%         set(gca,'xtick',[])
    set(gca,'XTickLabel',newNames)
    xticklabels(newNames')
%     ylim([0 180])
    hLeg = legend('example');
    set(hLeg,'visible','off')
    ylabel('Velocity (mm/s)');


% DISTANCE STEP LENGTH
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,2)

% rbp

attp2totDistStep = [attp2(:).distanceStepsMvals];
rbptotDistStep = [rbp(:).distanceStepsMvals];

attp2totDistStep=attp2totDistStep(:);
rbptotDistStep=rbptotDistStep(:);

clear attp2distanceStep
qqCnt = 1;
for qq = 1:length(attp2totDistStep)
    
    if ~isempty(attp2totDistStep{qq})
        attp2distanceStep(qqCnt,1)= mean(attp2totDistStep{qq});
            qqCnt =qqCnt+ 1;
    end
end

clear rbpdistanceStep
qqCnt = 1;
for qq = 1:length(rbptotDistStep)
    
    if ~isempty(rbptotDistStep{qq})
        rbpdistanceStep(qqCnt,1)= mean(rbptotDistStep{qq});
            qqCnt =qqCnt+ 1;
    end
end

rbpdistanceStep(isoutlier(rbpdistanceStep))=[];


% cac

attp2totDistStep = [attp2(:).distanceStepsMvals];
cactotDistStep = [cac(:).distanceStepsMvals];

attp2totDistStep=attp2totDistStep(:);
cactotDistStep=cactotDistStep(:);

clear attp2distanceStep
qqCnt = 1;
for qq = 1:length(attp2totDistStep)
    
    if ~isempty(attp2totDistStep{qq})
        attp2distanceStep(qqCnt,1)= mean(attp2totDistStep{qq});
            qqCnt =qqCnt+ 1;
    end
end

clear cacdistanceStep
qqCnt = 1;
for qq = 1:length(cactotDistStep)
    
    if ~isempty(cactotDistStep{qq})
        cacdistanceStep(qqCnt,1)= mean(cactotDistStep{qq});
            qqCnt =qqCnt+ 1;
    end
end

% cacdistanceStep(isoutlier(cacdistanceStep))=[];

% unc 13

attp2totDistStep = [attp2(:).distanceStepsMvals];
unc13totDistStep = [unc13(:).distanceStepsMvals];

attp2totDistStep=attp2totDistStep(:);
unc13totDistStep=unc13totDistStep(:);

clear attp2distanceStep
qqCnt = 1;
for qq = 1:length(attp2totDistStep)
    
    if ~isempty(attp2totDistStep{qq})
        attp2distanceStep(qqCnt,1)= mean(attp2totDistStep{qq});
            qqCnt =qqCnt+ 1;
    end
end

clear unc132stepDist
qqCnt = 1;
for qq = 1:length(unc13totDistStep)
    
    if ~isempty(unc13totDistStep{qq})
        unc13distanceStep(qqCnt,1)= mean(unc13totDistStep{qq});
            qqCnt =qqCnt+ 1;
    end
end

% unc13distanceStep(isoutlier(unc13distanceStep))=[];


save('distanceSteps.mat','attp2distanceStep','cacdistanceStep','unc13distanceStep','rbpdistanceStep');

% Distance Step  VIOLIN PLOT

attp2distanceStep(attp2distanceStep==0)=[];
unc13distanceStep(unc13distanceStep==0)=[];
rbpdistanceStep(rbpdistanceStep==0)=[];
cacdistanceStep(cacdistanceStep==0)=[];

dataSavedperLarva{1}=attp2distanceStep;
dataSavedperLarva{2}=unc13distanceStep;
dataSavedperLarva{3}=rbpdistanceStep;
dataSavedperLarva{4}=cacdistanceStep;


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
        dataSavedperLarva{3}=boxDataOutlierRemoved(~isnan(boxDataOutlierRemoved(:,3)),3);
        dataSavedperLarva{4}=boxDataOutlierRemoved(~isnan(boxDataOutlierRemoved(:,4)),4);

        catData = [dataSavedperLarva{1}; dataSavedperLarva{2}; dataSavedperLarva{3}; dataSavedperLarva{4}];
        genotypeStr = [repmat({'Attp2'},numel(dataSavedperLarva{1}),1);...
        repmat({'Unc13'},numel(dataSavedperLarva{2}),1);...
        repmat({'RBP'},numel(dataSavedperLarva{3}),1);...
        repmat({'Cac'},numel(dataSavedperLarva{4}),1)];

        T = table(genotypeStr,catData);

        genotypeOrder = {'Attp2','Unc13','RBP','Cac'};
        genotype = categorical(T.genotypeStr,genotypeOrder);

        boxchart(genotype,catData)  
        hold on
        scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);


    set(gca, 'FontName', 'Arial')
    set(gca,'fontsize',10)
    
    ranksumResults(1) = NaN;
    ranksumResults(2) = ranksum(dataSavedperLarva{1},dataSavedperLarva{2});
    ranksumResults(3) = ranksum(dataSavedperLarva{1},dataSavedperLarva{3});
    ranksumResults(4) = ranksum(dataSavedperLarva{1},dataSavedperLarva{4});

    newNames = {};
    for zz = 1:length(ranksumResults)
        thisGene = [num2str(ranksumResults(zz))];
        newNames{zz}=thisGene;
    end
    

% xlabel(attp2fns{fieldNum})
    
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
attp2fractTimeReor = [attp2(:).fractTimeReorientation];
rbpfractTimeReor = [rbp(:).fractTimeReorientation];

attp2fractTimeReor=attp2fractTimeReor(:);
rbpfractTimeReor=rbpfractTimeReor(:);

% rbpfractTimeReor(isoutlier(rbpfractTimeReor))=[];

% cac

attp2fractTimeReor = [attp2(:).fractTimeReorientation];
cacfractTimeReor = [cac(:).fractTimeReorientation];

attp2fractTimeReor=attp2fractTimeReor(:);
cacfractTimeReor=cacfractTimeReor(:);

% cacfractTimeReor(isoutlier(cacfractTimeReor))=[];


% unc13

attp2fractTimeReor = [attp2(:).fractTimeReorientation];
unc13fracTimeReor = [unc13(:).fractTimeReorientation];

attp2fractTimeReor=attp2fractTimeReor(:);
unc13fracTimeReor=unc13fracTimeReor(:);

% unc13fracTimeReor(isoutlier(unc13fracTimeReor))=[];

save('fractionTimeReorienting.mat','attp2fractTimeReor','cacfractTimeReor','unc13fracTimeReor','rbpfractTimeReor')

%   VIOLIN PLOT

subplot(2,3,3)

dataSavedperLarva{1}=attp2fractTimeReor;
dataSavedperLarva{2}=unc13fracTimeReor;
dataSavedperLarva{3}=rbpfractTimeReor;
dataSavedperLarva{4}=cacfractTimeReor;

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
        dataSavedperLarva{3}=boxDataOutlierRemoved(~isnan(boxDataOutlierRemoved(:,3)),3);
        dataSavedperLarva{4}=boxDataOutlierRemoved(~isnan(boxDataOutlierRemoved(:,4)),4);

        catData = [dataSavedperLarva{1}; dataSavedperLarva{2}; dataSavedperLarva{3}; dataSavedperLarva{4}];
        genotypeStr = [repmat({'Attp2'},numel(dataSavedperLarva{1}),1);...
        repmat({'Unc13'},numel(dataSavedperLarva{2}),1);...
        repmat({'RBP'},numel(dataSavedperLarva{3}),1);...
        repmat({'Cac'},numel(dataSavedperLarva{4}),1)];

        T = table(genotypeStr,catData);

        genotypeOrder = {'Attp2','Unc13','RBP','Cac'};
        genotype = categorical(T.genotypeStr,genotypeOrder);

        boxchart(genotype,catData)  
        hold on
        scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);


    set(gca, 'FontName', 'Arial')
    set(gca,'fontsize',10)
    
    ranksumResults(1) = NaN;
    ranksumResults(2) = ranksum(dataSavedperLarva{1},dataSavedperLarva{2});
    ranksumResults(3) = ranksum(dataSavedperLarva{1},dataSavedperLarva{3});
    ranksumResults(4) = ranksum(dataSavedperLarva{1},dataSavedperLarva{4});

    newNames = {};
    for zz = 1:length(ranksumResults)
        thisGene = [num2str(ranksumResults(zz))];
        newNames{zz}=thisGene;
    end
    

% xlabel(attp2fns{fieldNum})
    
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
% rbp

attp2reorBoutDur = [attp2(:).reorBoutDur];
rbpReorBoutDur = [rbp(:).reorBoutDur];

attp2reorBoutDur=attp2reorBoutDur(:);
rbpReorBoutDur=rbpReorBoutDur(:);

clear attp2reorBoutDuration
qqCnt = 1;
for qq = 1:length(attp2reorBoutDur)
    
    if ~isempty(attp2reorBoutDur{qq})
        attp2reorBoutDuration(qqCnt,1)= mean(attp2reorBoutDur{qq});
            qqCnt =qqCnt+ 1;
    end
end

clear rbpreorBoutDuration
qqCnt = 1;
for qq = 1:length(rbpReorBoutDur)
    
    if ~isempty(rbpReorBoutDur{qq})
        rbpreorBoutDuration(qqCnt,1)= mean(rbpReorBoutDur{qq});
            qqCnt =qqCnt+ 1;
    end
end

attp2reorBoutDuration = attp2reorBoutDuration*(1/fps);
rbpreorBoutDuration = rbpreorBoutDuration*(1/fps);
% 

% cac

attp2reorBoutDur = [attp2(:).reorBoutDur];
cacreorBoutDur = [cac(:).reorBoutDur];

attp2reorBoutDur=attp2reorBoutDur(:);
cacreorBoutDur=cacreorBoutDur(:);

clear attp2reorBoutDuration
qqCnt = 1;
for qq = 1:length(attp2reorBoutDur)
    
    if ~isempty(attp2reorBoutDur{qq})
        attp2reorBoutDuration(qqCnt,1)= mean(attp2reorBoutDur{qq});
            qqCnt =qqCnt+ 1;
    end
end

clear cacreorBoutDuration
qqCnt = 1;
for qq = 1:length(cacreorBoutDur)
    
    if ~isempty(cacreorBoutDur{qq})
        cacreorBoutDuration(qqCnt,1)= mean(cacreorBoutDur{qq});
            qqCnt =qqCnt+ 1;
    end
end

attp2reorBoutDuration = attp2reorBoutDuration*(1/fps);
cacreorBoutDuration = cacreorBoutDuration*(1/fps);



% unc13

attp2reorBoutDur = [attp2(:).reorBoutDur];
unc13reorBoutDur = [unc13(:).reorBoutDur];

attp2reorBoutDur=attp2reorBoutDur(:);
unc13reorBoutDur=unc13reorBoutDur(:);

clear attp2reorBoutDuration
qqCnt = 1;
for qq = 1:length(attp2reorBoutDur)
    
    if ~isempty(attp2reorBoutDur{qq})
        attp2reorBoutDuration(qqCnt,1)= mean(attp2reorBoutDur{qq});
            qqCnt =qqCnt+ 1;
    end
end

clear unc13reorBoutDuration
qqCnt = 1;
for qq = 1:length(unc13reorBoutDur)
    
    if ~isempty(unc13reorBoutDur{qq})
        unc13reorBoutDuration(qqCnt,1)= mean(unc13reorBoutDur{qq});
            qqCnt =qqCnt+ 1;
    end
end

attp2reorBoutDuration = attp2reorBoutDuration*(1/fps);
unc13reorBoutDuration = unc13reorBoutDuration*(1/fps);


save('timeSpentReorienting.mat','attp2reorBoutDuration','rbpreorBoutDuration','unc13reorBoutDuration','cacreorBoutDuration')


%   VIOLIN PLOT

subplot(2,3,4)

dataSavedperLarva{1}=attp2reorBoutDuration;
dataSavedperLarva{2}=unc13reorBoutDuration;
dataSavedperLarva{3}=rbpreorBoutDuration;
dataSavedperLarva{4}=cacreorBoutDuration;

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
        dataSavedperLarva{3}=boxData(~isnan(boxData(:,3)),3);
        dataSavedperLarva{4}=boxData(~isnan(boxData(:,4)),4);

        catData = [dataSavedperLarva{1}; dataSavedperLarva{2}; dataSavedperLarva{3}; dataSavedperLarva{4}];
        genotypeStr = [repmat({'Attp2'},numel(dataSavedperLarva{1}),1);...
        repmat({'Unc13'},numel(dataSavedperLarva{2}),1);...
        repmat({'RBP'},numel(dataSavedperLarva{3}),1);...
        repmat({'Cac'},numel(dataSavedperLarva{4}),1)];

        T = table(genotypeStr,catData);

        genotypeOrder = {'Attp2','Unc13','RBP','Cac'};
        genotype = categorical(T.genotypeStr,genotypeOrder);

        boxchart(genotype,catData)  
        hold on
        scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);

        
ylabel('Total time reorienting')

pVals(1) = NaN;
pVals(2) = ranksum(dataSavedperLarva{1},dataSavedperLarva{2});
pVals(3) = ranksum(dataSavedperLarva{1},dataSavedperLarva{3});
pVals(4) = ranksum(dataSavedperLarva{1},dataSavedperLarva{4});

     
% xlim([-0.05,5.05])
xlabel(num2str(pVals))
ylim([0 30])
set(gca, 'FontName', 'Arial')
set(gca,'fontsize',10)

ylabel('Time spent reorienting')


%%%%%%%%%%% time spent during continuous trajectory path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

attp2timeContTraj = [attp2(:).contPathDur];
rbpTimeContTraj = [rbp(:).contPathDur];

attp2timeContTraj=attp2timeContTraj(:);
rbpTimeContTraj=rbpTimeContTraj(:);

clear attp2continuousTimeOnTraj
qqCnt = 1;
for qq = 1:length(attp2timeContTraj)
    
    if ~isempty(attp2timeContTraj{qq})
        attp2continuousTimeOnTraj(qqCnt,1)= mean(attp2timeContTraj{qq});
            qqCnt =qqCnt+ 1;
    end
end

clear rbpContinuousTimeOnTraj
qqCnt = 1;
for qq = 1:length(rbpTimeContTraj)
    
    if ~isempty(rbpTimeContTraj{qq})
        rbpContinuousTimeOnTraj(qqCnt,1)= mean(rbpTimeContTraj{qq});
            qqCnt =qqCnt+ 1;
    end
end

attp2continuousTimeOnTraj = attp2continuousTimeOnTraj*(1/fps);
rbpContinuousTimeOnTraj = rbpContinuousTimeOnTraj*(1/fps);

% cac

attp2timeContTraj = [attp2(:).contPathDur];
cacTimeContTraj = [cac(:).contPathDur];

attp2timeContTraj=attp2timeContTraj(:);
cacTimeContTraj=cacTimeContTraj(:);

clear attp2continuousTimeOnTraj
qqCnt = 1;
for qq = 1:length(attp2timeContTraj)
    
    if ~isempty(attp2timeContTraj{qq})
        attp2continuousTimeOnTraj(qqCnt,1)= mean(attp2timeContTraj{qq});
            qqCnt =qqCnt+ 1;
    end
end

clear cacContinuousTimeOnTraj
qqCnt = 1;
for qq = 1:length(cacTimeContTraj)
    
    if ~isempty(cacTimeContTraj{qq})
        cacContinuousTimeOnTraj(qqCnt,1)= mean(cacTimeContTraj{qq});
            qqCnt =qqCnt+ 1;
    end
end

attp2continuousTimeOnTraj = attp2continuousTimeOnTraj*(1/fps);
cacContinuousTimeOnTraj = cacContinuousTimeOnTraj*(1/fps);

% unc13

attp2timeContTraj = [attp2(:).contPathDur];
unc13timeContTraj = [unc13(:).contPathDur];

attp2timeContTraj=attp2timeContTraj(:);
unc13timeContTraj=unc13timeContTraj(:);

clear attp2continuousTimeOnTraj
qqCnt = 1;
for qq = 1:length(attp2timeContTraj)
    
    if ~isempty(attp2timeContTraj{qq})
        attp2continuousTimeOnTraj(qqCnt,1)= mean(attp2timeContTraj{qq});
            qqCnt =qqCnt+ 1;
    end
end

clear unc13ContinuousTimeOnTraj
qqCnt = 1;
for qq = 1:length(unc13timeContTraj)
    
    if ~isempty(unc13timeContTraj{qq})
        unc13ContinuousTimeOnTraj(qqCnt,1)= mean(unc13timeContTraj{qq});
            qqCnt =qqCnt+ 1;
    end
end

attp2continuousTimeOnTraj = attp2continuousTimeOnTraj*(1/fps);
unc13ContinuousTimeOnTraj = unc13ContinuousTimeOnTraj*(1/fps);

save('timeOnContinuousTrajectories.mat','attp2continuousTimeOnTraj','cacContinuousTimeOnTraj','rbpContinuousTimeOnTraj','unc13ContinuousTimeOnTraj')

%   VIOLIN PLOT

subplot(2,3,5)

dataSavedperLarva{1}=attp2continuousTimeOnTraj;
dataSavedperLarva{2}=unc13ContinuousTimeOnTraj;
dataSavedperLarva{3}=rbpContinuousTimeOnTraj;
dataSavedperLarva{4}=cacContinuousTimeOnTraj;

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
        
% figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9])
subplot(2,3,5)

        dataSavedperLarva{1}=boxData(~isnan(boxData(:,1)),1);
        dataSavedperLarva{2}=boxData(~isnan(boxData(:,2)),2);
        dataSavedperLarva{3}=boxData(~isnan(boxData(:,3)),3);
        dataSavedperLarva{4}=boxData(~isnan(boxData(:,4)),4);

        catData = [dataSavedperLarva{1}; dataSavedperLarva{2}; dataSavedperLarva{3}; dataSavedperLarva{4}];
        genotypeStr = [repmat({'Attp2'},numel(dataSavedperLarva{1}),1);...
        repmat({'Unc13'},numel(dataSavedperLarva{2}),1);...
        repmat({'RBP'},numel(dataSavedperLarva{3}),1);...
        repmat({'Cac'},numel(dataSavedperLarva{4}),1)];

        T = table(genotypeStr,catData);

        genotypeOrder = {'Attp2','Unc13','RBP','Cac'};
        genotype = categorical(T.genotypeStr,genotypeOrder);

        boxchart(genotype,catData)  
        hold on
        scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);

ylabel('continuous time on traj path')

pVals(1) = NaN;
pVals(2) = ranksum(dataSavedperLarva{1},dataSavedperLarva{2});
pVals(3) = ranksum(dataSavedperLarva{1},dataSavedperLarva{3});
pVals(4) = ranksum(dataSavedperLarva{1},dataSavedperLarva{4});
%         xlim([-0.05,5.05])

xlabel(num2str(pVals))

set(gca, 'FontName', 'Arial')
set(gca,'fontsize',10)
ylim([0 30])


%%%%%%%%%%% Velocity during continuous trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,6)

minLengthInSeconds4contTraj = 2; % a trajectory has to be at least this long in order to use it in the velocity calculation
%attp2
for ff = 1:length(attp2)
    
    
    thisPlate = attp2(ff).meanCurv;
    thisPlateMvals = attp2(ff).mValCell;
    thisPlateCoords = attp2(ff).coordinates;

    thisPlateDistSteps = attp2(ff).distanceStepsMvals;
    
    clear contTrajVelocity
    for wellN = 1:length(thisPlate)
        thisCurvature=thisPlate{wellN};
        mVals = thisPlateMvals{wellN};
        x = thisPlateCoords{wellN}(:,1);
        y = thisPlateCoords{wellN}(:,2);
        

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
    
    attp2(ff).contTrajVelocity=contTrajVelocity;
    
end
%unc13
for ff = 1:length(unc13)
    
    
    thisPlate = unc13(ff).meanCurv;
    thisPlateMvals = unc13(ff).mValCell;
    thisPlateCoords = unc13(ff).coordinates;

    thisPlateDistSteps = unc13(ff).distanceStepsMvals;
    
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
    
    unc13(ff).contTrajVelocity=contTrajVelocity;
    
end
%rbp
for ff = 1:length(rbp)
    
    
    thisPlate = rbp(ff).meanCurv;
    thisPlateMvals = rbp(ff).mValCell;
    thisPlateCoords = rbp(ff).coordinates;

    thisPlateDistSteps = rbp(ff).distanceStepsMvals;
    
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
    
    rbp(ff).contTrajVelocity=contTrajVelocity;
    
end
%cac
for ff = 1:length(cac)
    
    
    thisPlate = cac(ff).meanCurv;
    thisPlateMvals = cac(ff).mValCell;
    thisPlateCoords = cac(ff).coordinates;

    thisPlateDistSteps = rbp(ff).distanceStepsMvals;
    
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
    
    cac(ff).contTrajVelocity=contTrajVelocity;
    
end


attp2timeContTraj = [attp2(:).contTrajVelocity];
attp2timeContTraj=attp2timeContTraj(:);

rbpTimeContTraj = [rbp(:).contTrajVelocity];
rbpTimeContTraj=rbpTimeContTraj(:);

cacTimeContTraj = [cac(:).contTrajVelocity];
cacTimeContTraj=cacTimeContTraj(:);


unc13timeContTraj = [unc13(:).contTrajVelocity];
unc13timeContTraj=unc13timeContTraj(:);

%   VIOLIN PLOT

dataSavedperLarva{1}=attp2timeContTraj;
dataSavedperLarva{2}=unc13timeContTraj;
dataSavedperLarva{3}=rbpTimeContTraj;
dataSavedperLarva{4}=cacTimeContTraj;

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
        
% figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9])

          dataSavedperLarva{1}=boxData(~isnan(boxData(:,1)),1);
        dataSavedperLarva{2}=boxData(~isnan(boxData(:,2)),2);
        dataSavedperLarva{3}=boxData(~isnan(boxData(:,3)),3);
        dataSavedperLarva{4}=boxData(~isnan(boxData(:,4)),4);

        catData = [dataSavedperLarva{1}; dataSavedperLarva{2}; dataSavedperLarva{3}; dataSavedperLarva{4}];
        genotypeStr = [repmat({'Attp2'},numel(dataSavedperLarva{1}),1);...
        repmat({'Unc13'},numel(dataSavedperLarva{2}),1);...
        repmat({'RBP'},numel(dataSavedperLarva{3}),1);...
        repmat({'Cac'},numel(dataSavedperLarva{4}),1)];

        T = table(genotypeStr,catData);

        genotypeOrder = {'Attp2','Unc13','RBP','Cac'};
        genotype = categorical(T.genotypeStr,genotypeOrder);

        boxchart(genotype,catData)  
        hold on
        scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);

ylabel('continuous time on traj path')

pVals(1) = NaN;
pVals(2) = ranksum(dataSavedperLarva{1},dataSavedperLarva{2});
pVals(3) = ranksum(dataSavedperLarva{1},dataSavedperLarva{3});
pVals(4) = ranksum(dataSavedperLarva{1},dataSavedperLarva{4});

% xlim([-0.05,5.05])

xlabel(num2str(pVals))

set(gca, 'FontName', 'Arial')
set(gca,'fontsize',10)
ylim([0 1])

