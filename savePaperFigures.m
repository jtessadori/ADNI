%% Load fig. 1 data
clear;
clc;
close all;
fig1data=readtable("fig1dataTable.txt");
tbl1=readtable("fig1PvalsState1.txt");
tbl2=readtable("fig1PvalsState2.txt");
tbl3=readtable("fig1PvalsState15.txt");
pVals=[tbl1.P_value,tbl2.P_value,tbl3.P_value];

% Plot fig. 1
close all;
clc;
fig1=figure;
fig1.Units='centimeters';
fig1.Position=[0 0 18 18];
ax=axes('OuterPosition',[0,.6,1,.4]);
hold on;
cMap=colormap('lines');
for currFeat=1:3
    boxplot(fig1data.(fig1data.Properties.VariableNames{currFeat}),fig1data.Group,'PlotStyle','compact',...
        'Color',cMap(1:5,:),'Symbol','.','Position',(1:5)+(currFeat-1)*6,'Labels',{'','','','',''})
    % Plot mean values
    for currGroup=1:5
        Xpos=(currFeat-1)*6+currGroup;
        Ypos=mean(fig1data.(fig1data.Properties.VariableNames{currFeat})(fig1data.Group==currGroup));
        plot([Xpos-.2,Xpos+.2],[Ypos,Ypos],'Color',cMap(currGroup,:),'LineWidth',2);
    end
end
for currFake=1:5 % Only used to generate legend
    plot([0,0],[0,0],'Color',cMap(currFake,:));
end

% Recover relevant data and add significance bars
tests=tbl1(:,1:2);
for currFeat=1:3
    nLines=0;
    for currTest=1:size(tbl1,1)
        if pVals(currTest,currFeat)<.05
            plot([tests{currTest,2},tests{currTest,1}]+(currFeat-1)*6,[1.05+nLines*.01,1.05+nLines*.01],'k');
            nLines=nLines+1;
        end
    end
end
legend({'HC','SMC','eMCI','lMCI','AD'});
stateNames={'dFCs 1','dFCs 2','dFCs 15'};
set(gca,'XLim',[0,22],'YLim',[-.05,1.12],'TickLength',[0,0],'XTick',3+[0,6,12],...
    'XTickLabel',stateNames);
ylabel('Fraction Time')

% Load fig. 2 data
figScore{1}=readtable("fig2dataScore1.txt");
figScore{2}=readtable("fig2dataScore2.txt");
figScore{3}=readtable("fig2dataScore3.txt");
scorePvals=readtable("fig2Pvals.txt");
stateNames={'dFCs 1','dFCs 2','dFCs 15'};

% Compute N
for currScore=1:3
    for currState=1:3
        for currGroup=1:2
            N_B(currScore,currState,currGroup)=sum((figScore{currScore}.State==currState)&(figScore{currScore}.Subpopulation==(currGroup-1))); %#ok<SAGROW> 
        end
    end
end

% Plot fig. 2
plotYlims=[-.5,32.5;19.5,31.5;-.5,6.9];
for currScore=1:3
    ax=axes('OuterPosition',[(currScore-1)/3,.3,1/3,.3]);
    hold on;
    for currFeat=1:3
        boxplot(figScore{currScore}.(figScore{currScore}.Properties.VariableNames{1})(figScore{currScore}.State==currFeat),...
            figScore{currScore}.Subpopulation(figScore{currScore}.State==currFeat),'PlotStyle','compact',...
            'Color',cMap(6:7,:),'Symbol','.','Position',(1:2)+(currFeat-1)*3,'Labels',{'',''})
        % Plot mean values
        for currGroup=1:2
            Xpos=(currFeat-1)*3+currGroup;
            Ypos=mean(figScore{currScore}.(figScore{currScore}.Properties.VariableNames{1})...
                ((figScore{currScore}.State==currFeat)&(figScore{currScore}.Subpopulation==(currGroup-1))));
            plot([Xpos-.3,Xpos+.3],[Ypos,Ypos],'Color',cMap(currGroup+5,:),'LineWidth',2);
        end

        % Recover relevant data and add significance bars
        if scorePvals{currFeat,currScore}<.05
            text(1.5+(currFeat-1)*3,plotYlims(currScore,2)-diff(plotYlims(currScore,:))*.06,'*','FontSize',16,'HorizontalAlignment','center')
%             plot((1:2)+(currFeat-1)*3,[plotYlims(currScore,2),plotYlims(currScore,2)]-diff(plotYlims(currScore,:))*.025,'k');
        end
    end
    set(gca,'TickLength',[0,0],'XTick',.5+1:3:8,...
    'XTickLabel',stateNames,'YLim',plotYlims(currScore,:),'XTickLabelRotation',0);
    title(figScore{currScore}.Properties.VariableNames{1});
end

% Load fig. 3 data
figScore{1}=readtable("fig3dataScore1.txt");
figScore{2}=readtable("fig3dataScore2.txt");
figScore{3}=readtable("fig3dataScore3.txt");
scorePvals=readtable("fig3Pvals.txt");

% Compute N
for currScore=1:3
    for currState=1:3
        for currGroup=1:2
            N_C(currScore,currState,currGroup)=sum((figScore{currScore}.State==currState)&(figScore{currScore}.Subpopulation==(currGroup-1))); %#ok<SAGROW> 
        end
    end
end

% Plot fig. 3
plotYlims=[1,3.8;.3,1.7;15,50]*1e-3;
for currScore=1:3
    ax=axes('OuterPosition',[(currScore-1)/3,0,1/3,.3]);
    hold on;
    for currFeat=1:3
        boxplot(figScore{currScore}.(figScore{currScore}.Properties.VariableNames{1})(figScore{currScore}.State==currFeat),...
            figScore{currScore}.Subpopulation(figScore{currScore}.State==currFeat),'PlotStyle','compact',...
            'Color',cMap(6:7,:),'Symbol','.','Position',(1:2)+(currFeat-1)*3,'Labels',{'',''})
        % Plot mean values
        for currGroup=1:2
            Xpos=(currFeat-1)*3+currGroup;
            Ypos=mean(figScore{currScore}.(figScore{currScore}.Properties.VariableNames{1})...
                ((figScore{currScore}.State==currFeat)&(figScore{currScore}.Subpopulation==(currGroup-1))));
            h=plot([Xpos-.3,Xpos+.3],[Ypos,Ypos],'Color',cMap(currGroup+5,:),'LineWidth',2);
            uistack(h,"down",1)
        end

        % Recover relevant data and add significance bars
        if scorePvals{currFeat,currScore}<.05
            text(1.5+(currFeat-1)*3,plotYlims(currScore,2)-diff(plotYlims(currScore,:))*.06,'*','FontSize',16,'HorizontalAlignment','center')
%             plot((1:2)+(currFeat-1)*3,[plotYlims(currScore,2),plotYlims(currScore,2)]-diff(plotYlims(currScore,:))*.025,'k');
        end
    end
    set(gca,'TickLength',[0,0],'XTick',.5+1:3:8,...
    'XTickLabel',stateNames,'YLim',plotYlims(currScore,:),'XTickLabelRotation',0);
    title(figScore{currScore}.Properties.VariableNames{1});
end

% Add panel letters
annotation('textbox','String','A','FontSize',14,'Position',[0 1 0 0],'LineStyle','none');
annotation('textbox','String','B','FontSize',14,'Position',[0 .62 0 0],'LineStyle','none');
annotation('textbox','String','C','FontSize',14,'Position',[0 .32 0 0],'LineStyle','none');
exportgraphics(fig1,'dFC_boxplots.png','Resolution',300)

%% HS-AD fraction time boxplots
% Load fig. 5 data
clear;
clc;
close all;
clear figData

fig5=figure;
fig5.Units='centimeters';
fig5.Position=[0 0 18.3 7.5];
for currExp=1:6
    figData{currExp}=readtable(sprintf("fig5dataTableExp%d.txt",currExp),'ReadVariableNames',true); %#ok<SAGROW> 
end
Pvals=readtable('fig5Pvals.txt');

% Plot fig. 5
cMap=colormap('lines');
clrs=cMap([1,5],:);
for currExp=1:size(figData,2)
    ax1=axes('OuterPosition',[mod(currExp,2)*.5,.95/3*(round(currExp/2)-1)+.05,.5,.95/3]);
    hold on;
    for currGroup=0:1
        boxplot(figData{currExp}{figData{currExp}.Group==currGroup,1:end-1},'PlotStyle','compact','Color',clrs(currGroup+1,:), ...
            'Symbol','.','Position',(1:size(Pvals,1))*3+currGroup,'Labels',{'','','','','','','','','','',''});
        for currState=1:size(Pvals,1)
            % Plot mean values
            Xpos=currState*3+currGroup;
            Ypos=mean(figData{currExp}{figData{currExp}.Group==currGroup,currState});
            h=plot([Xpos-.5,Xpos+.5],[Ypos,Ypos],'Color',clrs(currGroup+1,:),'LineWidth',2);
            uistack(h,"down",1)
        end
    end
    set(gca,'TickLength',[0,0],'XTick',(1:10)*3+.5,...
    'XTickLabel','','YLim',[-.1,1.2],'XLim',[2,35]);    
    if ismember(currExp,[1,2])
        xLbls=regexp(figData{currExp}.Properties.VariableNames(1:end-1),'\d*','Match');
        for currState=1:length(xLbls)
            text(currState*3,-.3,xLbls{currState});
        end
    end

    % Recover relevant data and add significance marks
    for currState=1:size(Pvals,1)
        if isfinite(Pvals{currState,currExp})&&Pvals{currState,currExp}<.05
            text(currState*3+.5,1.05,'*','FontSize',16,'HorizontalAlignment','center')
        end
    end
end

exportgraphics(fig5,'HShistsAllExpSorted.png','Resolution',300)
%% Load fig. 6 data (dFCs)
clear;
clc;
close all;
relStates=[1,2,15];
FC=cell(length(relStates),1);
avgdFC=cell(length(relStates),1);
for currState=1:length(relStates)
    FC{currState}=readmatrix(sprintf('fig6medianState%d.txt',relStates(currState)));
    avgdFC{currState}=readmatrix(sprintf('fig6avgdMedianState%d.txt',relStates(currState)));
end
netLbls=readmatrix("fig6networkLbls.txt");
netBegin=[0;find((netLbls(2:end)-netLbls(1:end-1))~=0)];
netEnd=[find((netLbls(2:end)-netLbls(1:end-1))~=0);length(netLbls)];

%% Plot fig. 6 (dFCs)
close all;
clc;

% Define custom colormap
cMap=ones(256,3);
cMap(1:128,[1,2])=repmat(linspace(0,1,128),2,1)';
cMap(129:256,[2,3])=repmat(linspace(1,0,128),2,1)';
loLim1=-1;
loLim2=-1;
hiLim1=1;
hiLim2=1;
imwidth=4;
imheigth=4;
fig6_1=figure('Units','centimeters','Position',[0 0 14.6 17.7]);

% Rescale cMap and set highest value to black
pp=mkpp([0 .15 .6 0],[.1 0; 1.5 .15*.1; (1-(.15*.1+1.5*(.6-.15)))/(1-.6) .15*.1+1.5*(.6-.15)]);
cMap=ppval(pp,cMap);
cMap(end,:)=0;

currRow=2;
for currState=1:length(relStates)

    % Main panel
    ax1=axes('Units','centimeters','Position',[1+(currState-1)*(imwidth+0.3),8.75+(currRow-1)*(imheigth+0.3),imwidth,imheigth]);
    imagesc(FC{currState},[loLim1,hiLim1])
    colormap(ax1,cMap)
    hold on
    for currNet=1:length(netBegin)
        plot([netBegin(currNet)+.5,netEnd(currNet)+.5],[netBegin(currNet)+.5,netBegin(currNet)+.5],'k')
        plot([netBegin(currNet)+.5,netEnd(currNet)+.5],[netEnd(currNet)+.5,netEnd(currNet)+.5],'k')
        plot([netBegin(currNet)+.5,netBegin(currNet)+.5],[netBegin(currNet)+.5,netEnd(currNet)+.5],'k')
        plot([netEnd(currNet)+.5,netEnd(currNet)+.5],[netBegin(currNet)+.5,netEnd(currNet)+.5],'k')
    end
    set(gca,'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[])
    title(sprintf('dFCs %d',relStates(currState)))

    % Left labels
    ax2=axes('Units','centimeters','InnerPosition',[0.25,8.75+(currRow-1)*(imheigth+0.3),.5,4]);
    imagesc(netLbls);
    set(gca,'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[],'Visible','Off')
    colormap(ax2,lines(7))

    % Bottom labels
    ax3=axes('Units','centimeters','InnerPosition',[1+(currState-1)*(imwidth+0.3),8+(currRow-1)*(imheigth+0.3),imwidth,.5]);
    imagesc(netLbls');
    set(gca,'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[],'Visible','Off')
    colormap(ax3,lines(7))

    % Colorbar
    if currState==3
        ax4=axes('Units','centimeters','InnerPosition',[5.25+(currState-1)*(imwidth+0.3),9+(currRow-1)*(imheigth+0.3),.5,3.5]);
        imagesc(linspace(loLim1,hiLim1,100)');
        set(gca,'XTick',[],'YTick',[],'YTickLabel',[],'XTickLabel',[],'YDir','normal')
        colormap(ax4,cMap)
        text(1,-5,num2str(loLim1),'HorizontalAlignment','center');
        text(1,107,num2str(hiLim1),'HorizontalAlignment','center');
    end
end

currRow=1;
for currState=1:length(relStates)
    % Main panel
    ax1=axes('Units','centimeters','Position',[1+(currState-1)*(imwidth+0.3),8+(currRow-1)*(imheigth+0.3),imwidth,imheigth]);
    imagesc(avgdFC{currState},[loLim2,hiLim2])
    colormap(ax1,cMap)
    hold on
    for currNet=1:length(netBegin)
        plot([netBegin(currNet)+.5,netEnd(currNet)+.5],[netBegin(currNet)+.5,netBegin(currNet)+.5],'k')
        plot([netBegin(currNet)+.5,netEnd(currNet)+.5],[netEnd(currNet)+.5,netEnd(currNet)+.5],'k')
        plot([netBegin(currNet)+.5,netBegin(currNet)+.5],[netBegin(currNet)+.5,netEnd(currNet)+.5],'k')
        plot([netEnd(currNet)+.5,netEnd(currNet)+.5],[netBegin(currNet)+.5,netEnd(currNet)+.5],'k')
    end
    set(gca,'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[])

    % Left labels
    ax2=axes('Units','centimeters','InnerPosition',[0.25,8+(currRow-1)*(imheigth+0.3),.5,4]);
    imagesc(netLbls);
    set(gca,'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[],'Visible','Off')
    colormap(ax2,lines(7))

    % Colorbar
    if currState==3
        ax4=axes('Units','centimeters','InnerPosition',[5.25+(currState-1)*(imwidth+0.3),8.25+(currRow-1)*(imheigth+0.3),.5,3.5]);
        imagesc(linspace(loLim2,hiLim2,100)');
        set(gca,'XTick',[],'YTick',[],'YTickLabel',[],'XTickLabel',[],'YDir','normal')
        colormap(ax4,cMap)
        text(1,-5,num2str(loLim2),'HorizontalAlignment','center');
        text(1,107,num2str(hiLim2),'HorizontalAlignment','center');
    end
end

% Labels

% Color bar
barXstart=1+(imwidth+0.3)-2;
barWidth=imwidth+4;
barYstart=6.75+(currRow-1)*(imheigth+0.3);
barHeight=.75;
lblAx=axes('Units','centimeters','InnerPosition',[barXstart,barYstart,barWidth,barHeight]);
imagesc(1:7);
set(gca,'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[],'Visible','Off')
colormap(lblAx,lines(7))

% Text
textboxStart=linspace(barXstart,barXstart+barWidth*6/7,7);
lblsText={'VIS','SMN','DAN','VAN','LIM','FPN','DMN'};
for currLbl=1:7
    annotation('textbox','String',lblsText{currLbl},'Units','centimeters',...
        'Position',[textboxStart(currLbl),barYstart,barWidth/7,barHeight],...
        'HorizontalAlignment','center','VerticalAlignment','middle',...
        'LineStyle','none','FontSize',12);
end


% Connectomes

% Define common variables
netGroups=cell(length(netBegin),1);
for currGroup=1:length(netBegin)
    netGroups{currGroup}=netBegin(currGroup)+1:netEnd(currGroup);
end
groupColor=[lines(7);lines(7)];
groupLbls={'L-VIS','L-SMN','L-DAN','L-VAN','L-LIM','L-FPN','L-DMN','R-VIS','R-SMN','R-DAN','R-VAN','R-LIM','R-FPN','R-DMN'};

% Recover only most significant links
nLinks=24;
for currConn=1:3
    ax5=axes('Units','centimeters','InnerPosition',[4.783*(currConn-1)+.1,0,4.783,6]);
    th=sort(reshape(FC{currConn},[],1),'descend');
    th=th(102+2*nLinks);
    currFC=FC{currConn}-th;
    circGraph(FC{currConn},'Groups',netGroups,'GroupColor',groupColor,'Colormap',repmat([.2,.2,.2],100,1),'GroupLbls',groupLbls,'Nlinks',24,'Th',.85);
%     circGraph(FC{currConn},'Groups',netGroups,'GroupColor',groupColor,'Colormap',cMap,'GroupLbls',groupLbls,'Nlinks',24);
%     exportgraphics(figConn,sprintf('conn%d.png',currConn),'Resolution',300);
end

% Add panel letters
A=annotation('textbox','String','A','FontSize',14,'Units','centimeters','Position',[0 16.8 1 1],'LineStyle','none');
B=annotation('textbox','String','B','FontSize',14,'Units','centimeters','Position',[0 5 1 1],'LineStyle','none');

% Save panel
exportgraphics(fig6_1,'dFCs.png','Resolution',300);

%% Summary figure, mini-panels
for currMini=1:3
    figMini=figure;
    figMini.Units='centimeters';
    figMini.Position=[0.1 0.1 1.2 0.6];
    plot(cumsum(randn(1e4,1)));
    set(gca,'XTick',[],'YTick',[]);
    exportgraphics(figMini,sprintf('mini%d.png',currMini),'Resolution',300);
end

%% Summary figure, DFCs sample sequence
close all;
clc;
figSeq=figure;
figSeq.Units='centimeters';
figSeq.Position=[0.1 0.1 10 1];
X=medfilt1(discretize(cumsum(randn(20,1)),5),3);
% imagesc(imresize(X,50,'nearest')')
axes('Units','centimeters','InnerPosition',[.2,.1,9.6,0.7]);
imagesc(X')
cMap=colormap('lines');
xline(.5:length(X)+.5)
set(gca,'XTick',[],'YTick',[],'Clipping','off');
X1=X==X(1);
Xstart=find([1;and(~X1(1:end-1),X1(2:end))]);
Xend=find(and(X1(1:end-1),~X1(2:end)));
hold on;
for currSeq=1:length(Xstart)
    plot([Xstart(currSeq)-.5,Xend(currSeq)+.5],[0.35,0.35],'Color',cMap(X(1),:),'LineWidth',4);
end
exportgraphics(figSeq,'dFCseqSample.png','Resolution',300);

%% Summary fig dFCs
clear;
clc;
close all;
relStates=[1,2,15];
FC=cell(length(relStates),1);
avgdFC=cell(length(relStates),1);
for currState=1:length(relStates)
    FC{currState}=readmatrix(sprintf('fig6medianState%d.txt',relStates(currState)));
    avgdFC{currState}=readmatrix(sprintf('fig6avgdMedianState%d.txt',relStates(currState)));
end
netLbls=readmatrix("fig6networkLbls.txt");
netBegin=[0;find((netLbls(2:end)-netLbls(1:end-1))~=0)];
netEnd=[find((netLbls(2:end)-netLbls(1:end-1))~=0);length(netLbls)];

% Define custom colormap
cMap=ones(256,3);
cMap(1:128,[1,2])=repmat(linspace(0,1,128),2,1)';
cMap(129:256,[2,3])=repmat(linspace(1,0,128),2,1)';
cMap=sqrt(cMap);

% Main panel
figdFCsample=figure;
figdFCsample.Units='centimeters';
figdFCsample.Position=[0 0 2 2];
seqAx=axes('Units','centimeters','InnerPosition',[.1,.1,1.9,1.9]);
imagesc(FC{1},[-.5,1]);
colormap(cMap)
hold on
for currNet=1:length(netBegin)
    plot([netBegin(currNet)+.5,netEnd(currNet)+.5],[netBegin(currNet)+.5,netBegin(currNet)+.5],'k')
    plot([netBegin(currNet)+.5,netEnd(currNet)+.5],[netEnd(currNet)+.5,netEnd(currNet)+.5],'k')
    plot([netBegin(currNet)+.5,netBegin(currNet)+.5],[netBegin(currNet)+.5,netEnd(currNet)+.5],'k')
    plot([netEnd(currNet)+.5,netEnd(currNet)+.5],[netBegin(currNet)+.5,netEnd(currNet)+.5],'k')
end
set(gca,'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[])
axis('square')

exportgraphics(gcf,'dFCsSmall.png','Resolution',300);

%% Summary fig fraction time sample
clear;
clc;
close all;
fig1data=readtable("fig1dataTable.txt");
tbl1=readtable("fig1PvalsState1.txt");
tbl2=readtable("fig1PvalsState2.txt");
tbl3=readtable("fig1PvalsState15.txt");
pVals=[tbl1.P_value,tbl2.P_value,tbl3.P_value];

close all;
clc;
figFTsample=figure;
figFTsample.Units='centimeters';
figFTsample.Position=[0.55 0.55 9 3];
hold on;
cMap=colormap('lines');
for currFeat=1:3
    boxplot(fig1data.(fig1data.Properties.VariableNames{currFeat}),fig1data.Group,'PlotStyle','compact',...
        'Color',cMap(1:5,:),'Symbol','.','Position',(1:5)+(currFeat-1)*6,'Labels',{'','','','',''})
    % Plot mean values
    for currGroup=1:5
        Xpos=(currFeat-1)*6+currGroup;
        Ypos=mean(fig1data.(fig1data.Properties.VariableNames{currFeat})(fig1data.Group==currGroup));
        plot([Xpos-.2,Xpos+.2],[Ypos,Ypos],'Color',cMap(currGroup,:),'LineWidth',2);
    end
end
set(gca,'XLim',[0,18],'YLim',[-.05,1.12],'TickLength',[0,0],'XTick',3+[0,6,12]);
exportgraphics(figFTsample,'relFreqsSmall.png','Resolution',300)

