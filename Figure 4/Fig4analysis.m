%% Plot results 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('dataTableNew.mat')
load('compiledData.mat')

analyzeRNAis = 1;
analyzeGlurIIA = 0;
analyzeShibire = 0;

onlyCountCasesThatArentConstantlyFiring = 0; % if you want to not count cases that are 100% active, since it distorts the DFF measurment

plotBar = 0;
plotBox = 1;

scatterOffset = .45;

column2analyze = [14 27]; % Fig 3: Postsynaptic 1b RNAi

dataTableNew = [dataTableScalars];

geneNames = fieldnames(genesCompiledScalars);

save('dataTableNew.mat','dataTableNew')

% figure
figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9])

for thisCol = 1:length(column2analyze)

    if analyzeRNAis
    prepNumCol = find(strcmpi(dataTableNew.Properties.VariableNames,'prepNum'));
    dateCol = find(strcmpi(dataTableNew.Properties.VariableNames,'sampleDate'));
    attp2inds =(strfind(dataTableNew.gene,'attP2'));
%     attp2inds =(strfind(dataTableNew.gene,'attP40'));

    attp2inds = find(~cellfun(@isempty,attp2inds));
    % attp40inds =(strfind(dataTableNew.gene,'attP40'));
    % attp40inds = (~cellfun(@isempty,attp40inds));
    attp40group=[];attp2group=[];attp2groupNames=[];attp40groupNames=[];
    attp2groupPrepNums = [];attp40groupPrepNums = [];
    attp2groupDates = [];attp40groupDates = [];
    elseif analyzeShibire
    prepNumCol = find(strcmpi(dataTableNew.Properties.VariableNames,'prepNum'));
    dateCol = find(strcmpi(dataTableNew.Properties.VariableNames,'sampleDate'));
    expNumCol = find(strcmpi(dataTableNew.Properties.VariableNames,'expNum'));
    attp2inds = ones(numel(dataTableNew.gene),1);
    % attp40inds =(strfind(dataTableNew.gene,'attP40'));
    % attp40inds = (~cellfun(@isempty,attp40inds));
    attp40group=[];attp2group=[];attp2groupNames=[];attp40groupNames=[];
    att6p2groupPrepNums = [];attp40groupPrepNums = [];
    attp2groupDates = [];attp40groupDates = [];
    attp2groupExpNums=[];
    else
    prepNumCol = find(strcmpi(dataTableNew.Properties.VariableNames,'prepNum'));
    dateCol = find(strcmpi(dataTableNew.Properties.VariableNames,'sampleDate'));
    attp2inds = ones(numel(dataTableNew.gene),1);
    % attp40inds =(strfind(dataTableNew.gene,'attP40'));
    % attp40inds = (~cellfun(@isempty,attp40inds));
    attp40group=[];attp2group=[];attp2groupNames=[];attp40groupNames=[];
    attp2groupPrepNums = [];attp40groupPrepNums = [];
    attp2groupDates = [];attp40groupDates = [];
    end

    attp2cnt=1;attp40cnt=1;
    clear attp2pvals attp40pvals
    thisGraphTitle = dataTableNew.Properties.VariableNames{column2analyze(thisCol)};
    disp(['Calculating ' thisGraphTitle])
    for geneNum = 1:numel(geneNames)

        if ~isempty(strmatch(geneNames{geneNum},'attP2')) || ~isempty(strmatch(geneNames{geneNum},'attP40')) 
           continue
        else
            try
                thisMat =(strfind(dataTableNew.gene,geneNames{geneNum}));
                geneInds = (~cellfun(@isempty,thisMat));
                thisInsSite = unique([dataTableNew.insertionSitesTotal(geneInds)]);

                theseData = dataTableNew{geneInds,column2analyze(thisCol)};
                theseNames = dataTableNew.gene(geneInds);
                thesePrepNums = str2num(dataTableNew.prepNum(geneInds,2:end));
                thesePrepDates = (dataTableNew.sampleDate(geneInds));
                
                
                if onlyCountCasesThatArentConstantlyFiring
                    
                    perActTime = dataTableNew{geneInds,15};
                    inds2discount = perActTime==1;
                    theseData(inds2discount)=[];
                    theseNames(inds2discount)=[];
                    thesePrepNums(inds2discount)=[];
                    thesePrepDates(inds2discount)=[];

                end
                
                if analyzeShibire ==1
                theseExpNums = dataTableNew.expNum(geneInds);
                end

                underScoreLocs = strfind(theseNames,'_');

                for ii =1:size(strfind(theseNames,'_'),1)
                    theseNames{ii}(underScoreLocs{ii})=[];
                end

                attp2Data = dataTableNew{attp2inds,column2analyze(thisCol)};
    %             attp40Data = dataTableNew{attp40inds,column2analyze};

                attp2prepNums = str2num(dataTableNew.prepNum(attp2inds,2:end));
                attp2prepDates = dataTableNew{attp2inds,dateCol};

    %             attp40prepNums = str2num(dataTableNew.prepNum(attp40inds,2:end));
    %             attp40prepDates = dataTableNew{attp40inds,dateCol};
    %             
                if ~isempty(strmatch(class(theseData),'cell'))
                    theseData=cell2mat(theseData);
                    attp2Data=cell2mat(attp2Data);
    %                 attp40Data=cell2mat(attp40Data);
                end

                if thisInsSite==1 % attp2
                    attp2group = [attp2group; theseData];
                    attp2groupNames = [attp2groupNames; theseNames];

                    attp2groupPrepNums = [attp2groupPrepNums; thesePrepNums];
                    attp2groupDates = [attp2groupDates;thesePrepDates];
                    
                    if analyzeShibire==1
                    attp2groupExpNums = [attp2groupExpNums; theseExpNums];
                    end
                    [p,~]=ranksum(attp2Data,theseData);
                    attp2pvals(attp2cnt,1)=p;
                    attp2geneNames{attp2cnt,1}=geneNames{geneNum};
                    attp2means(attp2cnt,1)=mean(theseData);
                    attp2sems(attp2cnt,1)=std(theseData)/sqrt(size(theseData,1));

                    attp2cnt = attp2cnt+1;
                elseif thisInsSite==2 % attp40
                    attp40group = [attp40group; theseData];
                    attp40groupNames = [attp40groupNames; theseNames];

                    attp40groupPrepNums = [attp40groupPrepNums; thesePrepNums];
                    attp40groupDates = [attp40groupDates;thesePrepDates];

                    [p,~]=ranksum(attp40Data,theseData);
                    attp40pvals(attp40cnt,1)=p;
                    attp40geneNames{attp40cnt,1}=geneNames{geneNum};
                    attp40means(attp40cnt,1)=mean(theseData);
                    attp40sems(attp40cnt,1)=std(theseData)/sqrt(size(theseData,1));

                    attp40cnt = attp40cnt+1;
                end

            catch
            end
        end

    end
    % 

    if analyzeRNAis

    b=repmat({'Ct (AttP2)'},numel(attp2Data),1);
    attp2group = [attp2Data;attp2group];
    attp2groupNames = [b;attp2groupNames];
    attp2groupPrepNums = [attp2prepNums;attp2groupPrepNums];
    attp2groupDates = [attp2prepDates;attp2groupDates];
%     loopIterations = [1 2 3 4];
    loopIterations = [1:numel(unique(attp2groupNames))];

    elseif analyzeGlurIIA% for Gluriia null case

    loopIterations = [2 1];
    
    elseif analyzeShibire
  
    loopIterations = [2 1];

    end

    attp2groupNamesUnique=unique(attp2groupNames);

    clear meanAttp2 semAttp2 pvalsAttp2 meanAttp2perLarva semAttp2perLarva dataSavedperLarva geneNamesperLarva
    cnt=1;
    for ii = 1:length(loopIterations)
        
    disp(attp2groupNamesUnique{loopIterations(ii)});
% for ii = 1:length(attp2groupNamesUnique)
    kInds = strfind(attp2groupNames,attp2groupNamesUnique{loopIterations(ii)});
    kInds = find(~cellfun(@isempty,kInds));

    theseData = attp2group(kInds);
    pNums = attp2groupPrepNums(kInds);
    datesNum = attp2groupDates(kInds);
    prepInfo = [pNums datesNum];
    
    if analyzeShibire ==1
    expNumInfo = attp2groupExpNums(kInds);
    prepInfo = [pNums datesNum expNumInfo];

    uniquePreps = prepInfo(:,1).*prepInfo(:,2).*expNumInfo;
    [C,ia,prepNum]=unique(uniquePreps,'stable');
    else
    uniquePreps = prepInfo(:,1).*prepInfo(:,2);
    [C,ia,prepNum]=unique(uniquePreps,'stable');
    end
    
   
%     meanAttp2(cnt,1) = nanmean(theseData);
%     semAttp2(cnt,1) = nanstd(theseData)/sqrt(size(theseData,1));
%     [p,h]=ranksum(attp2Data,theseData);
%     dataSaved{cnt,1}=theseData;
%     geneNames{cnt,1}=attp2groupNamesUnique{ii};
%     pvalsAttp2(cnt,1)=p;
    
    clear larvData 
    for larvNum = 1:max((prepNum))
    thisPrepInds = prepNum==larvNum;
    larvData(larvNum,1) = nanmean(theseData(thisPrepInds));
    end
    
    meanAttp2perLarva(cnt,1) = nanmean(larvData);
    semAttp2perLarva(cnt,1) = nanstd(larvData)/sqrt(size(larvData,1));
    dataSavedperLarva{cnt,1}=larvData;
    geneNamesperLarva{cnt,1}=attp2groupNamesUnique{loopIterations(ii)};
%     pvalsAttp2(cnt,1)=p;
    
    cnt=cnt+1;
    end

%     try
    
    clear pVals
    pVals=[];
    if analyzeRNAis
    pVals(1) = NaN;
    pVals(2) = ranksum(dataSavedperLarva{1},dataSavedperLarva{2});
    pVals(3) = ranksum(dataSavedperLarva{1},dataSavedperLarva{3});
    pVals(4) = ranksum(dataSavedperLarva{1},dataSavedperLarva{4});
    else
    pVals(1) = NaN;
    pVals(2) = ranksum(dataSavedperLarva{2},dataSavedperLarva{1});
    end
    
    pVals
    
        subplot(2,3,thisCol)
        cc = linspecer(length(meanAttp2perLarva));
        pvalsAttp2(1)=NaN';
        hold on
        
    if plotBar
       
        superbar(meanAttp2perLarva,'E',semAttp2perLarva,'P',pVals,'BarFaceColor', 'none', 'BarEdgeColor',cc,'PStarFontSize',10,'BarWidth',0.5,'BarLineWidth',.7,'ErrorbarLineWidth',.7)
        ylabel(thisGraphTitle);

        plotCnt = 1;
        for plotNum = 1:length(meanAttp2perLarva)

        % superbar(meanAttp2,'E',semAttp2,'BarFaceColor', 'none', 'BarEdgeColor',cc,'PStarFontSize',35,'BarWidth',0.4)
        scatter(repmat(plotCnt+scatterOffset,length(dataSavedperLarva{plotCnt}),1),dataSavedperLarva{plotCnt},50,'.','k','jitter','on', 'jitterAmount',0.1);

        plotCnt=plotCnt+1;

        end


    %     set(gca,'XTick',[1:length(geneNamesperLarva)])
        set(gca,'xtick',[])
        set(gca,'XTickLabel',geneNamesperLarva')
        set(gca,'XTickLabel',geneNamesperLarva')
        xlim([-0.05,5.05])
        
    elseif plotBox
%         y = randn(50,3,3);
%               x = [1 2 3.5];
%               y(1:25) = NaN;
        
      
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

        
        if size(boxData,2)==4
        catData = [boxData(:,1); boxData(:,2); boxData(:,3); boxData(:,4)];
        genotypeStr = [repmat(geneNamesperLarva(1),numel(boxData(:,1)),1);...
        repmat(geneNamesperLarva(2),numel(boxData(:,2)),1);...
        repmat(geneNamesperLarva(3),numel(boxData(:,3)),1);...
        repmat(geneNamesperLarva(4),numel(boxData(:,4)),1)];
            T = table(genotypeStr,catData);
        genotypeOrder = {'Ct (AttP2)','Munc13','RBP','cac'};
        pVals=[1 ranksum(dataSavedperLarva{1},dataSavedperLarva{2}) ranksum(dataSavedperLarva{1},dataSavedperLarva{3}) ranksum(dataSavedperLarva{1},dataSavedperLarva{4})]

        elseif size(boxData,2)==2 && analyzeRNAis
             
         catData = [boxData(:,1); boxData(:,2);];
        genotypeStr = [repmat(geneNamesperLarva(1),numel(boxData(:,1)),1);...
        repmat(geneNamesperLarva(2),numel(boxData(:,2)),1)];
            T = table(genotypeStr,catData);
%         genotypeOrder = {'Ct (AttP2)','sema2b'};
        genotypeOrder = {'Ct (AttP2)','GluRIIA'};
        pVals=[1 ranksum(dataSavedperLarva{1},dataSavedperLarva{2})]
        

        else
            
        catData = [boxData(:,1); boxData(:,2);];
        genotypeStr = [repmat(geneNamesperLarva(1),numel(boxData(:,1)),1);...
        repmat(geneNamesperLarva(2),numel(boxData(:,2)),1)];
            T = table(genotypeStr,catData);
%         genotypeOrder = {'Ct (AttP2)','sema2b'};
        genotypeOrder = {'WTGC6m','GluRIIA'};
        pVals=[1 ranksum(dataSavedperLarva{1},dataSavedperLarva{2})]

        end
        

        
        genotype = categorical(T.genotypeStr,genotypeOrder);



%         ind2remove = catData<0;
%         genotype(ind2remove)=[];
%         catData(ind2remove)=[];
%         
%         ind2remove = catData==Inf;
%         genotype(ind2remove)=[];
%         catData(ind2remove)=[];
%         
%         ind2remove = isnan(catData);
%         genotype(ind2remove)=[];
%         catData(ind2remove)=[];
        
        boxchart(genotype,catData)

        hold on
        scatter(genotype,catData,200,'.','k','jitter','on', 'jitterAmount',0.1);

%         pVals=[1 ranksum(dataSavedperLarva{1},dataSavedperLarva{2})];
        
        xlabel(num2str(pVals))
        
%         ylim([0 ceil(max(catData))+1])
% 
%             vs=violinplot(boxData,1:size(boxData,2),'Width',0.3,...
%                 'ViolinAlpha',1,'ShowData',true,'ShowMean',true);
% 
%             if analyzeRNAis
%                 vs(1).ViolinColor=[1 1 1];
%                 vs(2).ViolinColor=[18/255 97/255 160/255];
%                 vs(3).ViolinColor=[56/255 149/255 211/255];
%                 vs(4).ViolinColor=[88/255 204/255 237/255];
%                 xlim([-0.05,5.05])
% 
%             elseif analyzeGlurIIA
%                 vs(1).ViolinColor=[1 1 1];
%                 vs(2).ViolinColor=[1 0.18 0];
%                 xlim([-0.05,3.05])
% 
%             elseif analyzeShibire
%                 vs(1).ViolinColor=[56/255 149/255 211/255];
%                 vs(2).ViolinColor=[1 0.18 0];
%                 xlim([-0.05,3.05])
% 
%             end
%             
%             for xx = 1:size(boxData,2)
%                 vs(xx).ScatterPlot.MarkerFaceAlpha=1;
%                 vs(xx).ScatterPlot.SizeData = 16;
%                 vs(xx).EdgeColor=[0 0 0];
%                 vs(xx).ScatterPlot.MarkerFaceColor=[0 0 0];
%                 vs(xx).MeanPlot.Color = [0 0 0];
%             end
%             
    set(gca, 'FontName', 'Arial')
    set(gca,'fontsize',10)
    if strcmp(thisGraphTitle,'percentActiveTime')
        ylim([0 1.02])
    elseif strcmp(thisGraphTitle,'BAboutDuration_IncludingAlwaysOnAsBout')
        ylim([0 max(boxData(:))*1.02])
    end

%         hold on
%        
%         xScatterData = repmat(x,size(boxData,1),1);
%         for scatNum = 1:size(xScatterData,2)
%             xTemp = xScatterData(:,scatNum);
%             scatter(xTemp,boxData(:,scatNum),12,'k','filled')
% 
%             
%         end
%         semVals = nanstd(boxData)./[sqrt(nPoints)]';
%         eb = errorbar(x,nanmean(boxData),semVals,'.k')
%         eb.MarkerSize = 20;
%         eb.CapSize = 20;
%             ylim([0 1.1])
         
        set(gca,'fontsize',10)
        ylabel(thisGraphTitle);

%         set(gca,'XTick',[1:length(geneNamesperLarva)])
%         set(gca,'xtick',[])
        set(gca,'XTickLabel',geneNamesperLarva')
        hLeg = legend('example');
        set(hLeg,'visible','off')
    else
        disp('Neither plotting option was chosen')
    end
        
    
%     catch
%        disp('Failed because NMJs didnt reach conditions to make particular measurement') 
%         
%     end
end

