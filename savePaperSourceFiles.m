clc;
clear;
addpath C:\Code\Common;

%% Load unsupervised sections of experiments
% currentTask={'CN','AD','EMCI','LMCI','SMC'};
% this=prepareTrajObject(currentTask);
% expDir='C:\Code\2022_01_ADNI\22_10_25_genHMM\';
% this.recoverExpResults(expDir)
load preProcData_5class.mat

% Reorder labels
clc;
cLbls=zeros(size(this.lbls));
lblVals=unique(this.lbls);
nClasses=length(lblVals);
sampleRep=1;
lblOrdr=[1,5,3,4,2];
for currClass=1:nClasses
    cLbls(this.lbls==lblOrdr(currClass))=currClass;
end

%% Multiple comparison for different medians, for all experiments
c=cell(3,1);
for currTest=1:3
    t=this.statePall{sampleRep,currTest};
    [~,~,stats]=kruskalwallis(t,cLbls,'off');
    c{currTest} = multcompare(stats,'Display','off','CriticalValueType','dunn-sidak');
%     c{currTest} = multcompare(stats,'Display','off');
end

% Produce three tables.
tbl1 = array2table(c{1}(:,[1:2,end]),"VariableNames", ...
    ["Group A","Group B","P-value"]);
tbl2 = array2table(c{2}(:,[1:2,end]),"VariableNames", ...
    ["Group A","Group B","P-value"]);
tbl3 = array2table(c{3}(:,[1:2,end]),"VariableNames", ...
    ["Group A","Group B","P-value"]);

%% Test for age in different groups
% Multiple comparison for different medians, for all experiments
for currTest=1:3
    [~,~,ageStats]=kruskalwallis(this.subjAge,cLbls,'off');
    ageMultCompare = multcompare(stats,'Display','off');
end

% Produce result table
ageTestTable = array2table(ageMultCompare(:,[1:2,end]),"VariableNames", ...
    ["Group A","Group B","P-value"]);
writetable(ageTestTable,'ageTestPvals.txt');

%% There seems to be a small number of stable clusters across experiments. 
% Recover data of fraction times by condition for a sample exp
close all
clc;
sampleRep=1;
lblVals=unique(this.lbls);
lblVals=lblVals(lblOrdr);
nClasses=length(lblVals);
stateFT=zeros(this.nSubjs,3);
condGroup=[];
this.recoverExpFeats;
for currState=1:3
    relStateFeats=[];
    for currCond=1:nClasses
        relIdx=this.lbls==lblVals(currCond);
        relStateFeats=cat(1,relStateFeats,this.stateFeats(relIdx,sampleRep,currState));
        if currState==1
            condGroup=cat(1,condGroup,ones(sum(relIdx),1)*currCond);
        end
    end
    stateFT(:,currState)=relStateFeats;
end

% Put data in a table and save them
fig1dataTable=array2table([stateFT,condGroup],'VariableNames',{'State1','State2','State15','Group'});
writetable(fig1dataTable);
writetable(tbl1,'fig1PvalsState1.txt');
writetable(tbl2,'fig1PvalsState2.txt');
writetable(tbl3,'fig1PvalsState15.txt');

%% Perform analysis on SMC + MCI block, clinical scores
relClasses={'SMC','MCI'};
relScores={'ADAS11','MMSE','CDRSB'};
ageCorrection='empirical';

[Ptable_clinScore,N_clinScore,scoreLog,lblLog,~,~,PsexTable1]=this.testClinInfo(relClasses,relScores,ageCorrection,sampleRep);
writetable(PsexTable1,'sexScoresPvals.txt');

% % Perform Bonferroni correction
% Ptable_BonCorrected=Ptable_clinScore;
% Ptable_BonCorrected{:,:}=Ptable_BonCorrected{:,:}*3*3;

% Perform Benjamini-Hochberg correction
Ptable_BH_Corrected=Ptable_clinScore;
Ptable_BH_Corrected{:,:}=reshape(fdr_BH(Ptable_BH_Corrected{:,:},.05),size(Ptable_BH_Corrected));

%% Recover clin scores in subjects split by fraction time
close all;
clc;
relScore=cell(3,1);
relLbl=cell(3,1);
featGroup=cell(3,1);
for currScore=1:3
    relScore{currScore}=[];
    relLbl{currScore}=[];
    featGroup{currScore}=[];
    for currFeat=1:3
        relScore{currScore}=cat(1,relScore{currScore},scoreLog{sampleRep,currScore,currFeat,1});
        relLbl{currScore}=cat(1,relLbl{currScore},lblLog{sampleRep,currScore,currFeat,1});
        featGroup{currScore}=cat(1,featGroup{currScore},ones(length(lblLog{sampleRep,currScore,currFeat,1}),1)*currFeat);
    end
end

% Put data in tables and save them
for currScore=1:3
    fig2dataTable=array2table([relScore{currScore},featGroup{currScore},relLbl{currScore}],'VariableNames',{relScores{currScore},'State','Subpopulation'});
    writetable(fig2dataTable,sprintf('fig2dataScore%d.txt',currScore));
end
writetable(Ptable_BH_Corrected,'fig2Pvals.txt');

%% Perform analysis on SMC + MCI block, anatomical data
relScores={'Hippocampus','Amygdala','CerebellumCortex'};
ageCorrection='empirical';

[Ptable_relVolume,N_relVolume,scoreLog,lblLog,pLarge,nLarge,PsexTable2]=this.testClinInfo(relClasses,relScores,ageCorrection,sampleRep);
writetable(PsexTable2,'sexAnatPvals.txt');

% % Perform Bonferroni correction
% Ptable_BonCorrected2=Ptable_relVolume;
% Ptable_BonCorrected2{:,:}=Ptable_BonCorrected2{:,:}*3*3;

% Perform Benjamini-Hochberg correction
Ptable_BH_Corrected2=Ptable_relVolume;
Ptable_BH_Corrected2{:,:}=reshape(fdr_BH(Ptable_BH_Corrected2{:,:},.05),size(Ptable_BH_Corrected2));

%% Recover relative volumes in subjects split by fraction time
close all;
clc;
relVol=cell(3,1);
relLbl=cell(3,1);
featGroup=cell(3,1);
for currScore=1:3
    relVol{currScore}=[];
    relLbl{currScore}=[];
    featGroup{currScore}=[];
    for currFeat=1:3
        relVol{currScore}=cat(1,relVol{currScore},scoreLog{sampleRep,currScore,currFeat,1});
        relLbl{currScore}=cat(1,relLbl{currScore},lblLog{sampleRep,currScore,currFeat,1});
        featGroup{currScore}=cat(1,featGroup{currScore},ones(length(lblLog{sampleRep,currScore,currFeat,1}),1)*currFeat);
    end
end

% Put data in tables and save them
for currScore=1:3
    fig3dataTable=array2table([relVol{currScore},featGroup{currScore},relLbl{currScore}],'VariableNames',{relScores{currScore},'State','Subpopulation'});
    writetable(fig3dataTable,sprintf('fig3dataScore%d.txt',currScore));
end
writetable(Ptable_BH_Corrected2,'fig3Pvals.txt');

%% Response to reviewer 3: repeat of the analysis above, including all populations

% Perform analysis on full pop, clinical scores
relClasses={'CN','SMC','MCI','AD'};
relScores={'ADAS11','MMSE','CDRSB'};
ageCorrection='empirical';

[Ptable_clinScore,N_clinScore,scoreLog,lblLog,~,~,PsexTable1]=this.testClinInfo(relClasses,relScores,ageCorrection,sampleRep);

% % Perform Bonferroni correction
% Ptable_BonCorrected=Ptable_clinScore;
% Ptable_BonCorrected{:,:}=Ptable_BonCorrected{:,:}*3*3;

% Perform Benjamini-Hochberg correction
Ptable_BH_Corrected=Ptable_clinScore;
Ptable_BH_Corrected{:,:}=reshape(fdr_BH(Ptable_BH_Corrected{:,:},.05),size(Ptable_BH_Corrected));

% Recover clin scores in subjects split by fraction time
close all;
clc;
relScore=cell(3,1);
relLbl=cell(3,1);
featGroup=cell(3,1);
for currScore=1:3
    relScore{currScore}=[];
    relLbl{currScore}=[];
    featGroup{currScore}=[];
    for currFeat=1:3
        relScore{currScore}=cat(1,relScore{currScore},scoreLog{sampleRep,currScore,currFeat,1});
        relLbl{currScore}=cat(1,relLbl{currScore},lblLog{sampleRep,currScore,currFeat,1});
        featGroup{currScore}=cat(1,featGroup{currScore},ones(length(lblLog{sampleRep,currScore,currFeat,1}),1)*currFeat);
    end
end

% Put data in tables and save them
for currScore=1:3
    fig2dataTable=array2table([relScore{currScore},featGroup{currScore},relLbl{currScore}],'VariableNames',{relScores{currScore},'State','Subpopulation'});
    writetable(fig2dataTable,sprintf('fig2dataScore_fullPop%d.txt',currScore));
end
writetable(Ptable_BH_Corrected,'fig2Pvals_fullPop.txt');

% Perform analysis on full pop, anatomical data
relScores={'Hippocampus','Amygdala','CerebellumCortex'};
ageCorrection='empirical';

[Ptable_relVolume,N_relVolume,scoreLog,lblLog,pLarge,nLarge,PsexTable2]=this.testClinInfo(relClasses,relScores,ageCorrection,sampleRep);

% % Perform Bonferroni correction
% Ptable_BonCorrected2=Ptable_relVolume;
% Ptable_BonCorrected2{:,:}=Ptable_BonCorrected2{:,:}*3*3;

% Perform Benjamini-Hochberg correction
Ptable_BH_Corrected2=Ptable_relVolume;
Ptable_BH_Corrected2{:,:}=reshape(fdr_BH(Ptable_BH_Corrected2{:,:},.05),size(Ptable_BH_Corrected2));

% Recover relative volumes in subjects split by fraction time
close all;
clc;
relVol=cell(3,1);
relLbl=cell(3,1);
featGroup=cell(3,1);
for currScore=1:3
    relVol{currScore}=[];
    relLbl{currScore}=[];
    featGroup{currScore}=[];
    for currFeat=1:3
        relVol{currScore}=cat(1,relVol{currScore},scoreLog{sampleRep,currScore,currFeat,1});
        relLbl{currScore}=cat(1,relLbl{currScore},lblLog{sampleRep,currScore,currFeat,1});
        featGroup{currScore}=cat(1,featGroup{currScore},ones(length(lblLog{sampleRep,currScore,currFeat,1}),1)*currFeat);
    end
end

% Put data in tables and save them
for currScore=1:3
    fig3dataTable=array2table([relVol{currScore},featGroup{currScore},relLbl{currScore}],'VariableNames',{relScores{currScore},'State','Subpopulation'});
    writetable(fig3dataTable,sprintf('fig3dataScore_fullPop%d.txt',currScore));
end
writetable(Ptable_BH_Corrected2,'fig3Pvals_fullPop.txt');

%% Response to reviewer 2: analysis of brain regions extended to all of them

relScores=this.clinInfo.Properties.VariableNames(22:62);
ageCorrection='empirical';

[Ptable_relVolume,N_relVolume,scoreLog,lblLog,pLarge,nLarge]=this.testClinInfo(relClasses,relScores,ageCorrection,sampleRep);

% % Perform Bonferroni correction
% Ptable_BonCorrected2=Ptable_relVolume;
% Ptable_BonCorrected2{:,:}=Ptable_BonCorrected2{:,:}*3*3;

% Perform Benjamini-Hochberg correction
Ptable_BH_extended2=Ptable_relVolume;
Ptable_BH_extended2{:,:}=reshape(fdr_BH(Ptable_BH_extended2{:,:},.05),size(Ptable_BH_extended2));

% Load Data
brainRegions = Ptable_BH_extended2.Properties.VariableNames; % Extract brain region names
qValues = Ptable_BH_extended2{:,:}; % Extract q-values (3x41 matrix)

% Convert q-values to -log10 scale for better visualization
logQValues = -log10(qValues);

% Define split index to divide the brain regions into two figures
nRegions = length(brainRegions);
splitIdx = ceil(nRegions / 2);

brainRegions_1 = brainRegions(1:splitIdx);
brainRegions_2 = brainRegions(splitIdx+1:end);

logQValues_1 = logQValues(:, 1:splitIdx);
logQValues_2 = logQValues(:, splitIdx+1:end);

% Plot First Half of the Brain Regions
subplot(2,1,1);
bar(logQValues_1', 'grouped'); % Transpose for grouped bar plot
ax = gca;
ax.XTick = 1:length(brainRegions_1); % Manually set all X ticks
ax.XTickLabel = brainRegions_1; % Ensure all labels are displayed
ax.XTickLabelRotation = 45; % Rotate for readability
ax.FontSize = 10;
ylabel('-log10(q-value)');
title('Brain Region Associations with dFC States (Part 1)');
yline(-log10(0.05), 'r--'); % Add significance threshold line
legend({'State 1', 'State 2', 'State 15','FDR threshold'}, 'Location', 'northeastoutside');
grid on;
set(gca, 'XTickLabelMode', 'manual'); % Prevent dynamic updating

% Plot Second Half of the Brain Regions
subplot(2,1,2);
bar(logQValues_2', 'grouped'); % Transpose for grouped bar plot
ax = gca;
ax.XTick = 1:length(brainRegions_2); % Manually set all X ticks
ax.XTickLabel = brainRegions_2; % Ensure all labels are displayed
ax.XTickLabelRotation = 45; % Rotate for readability
ax.FontSize = 10;
ylabel('-log10(q-value)');
title('Brain Region Associations with dFC States (Part 2)');
yline(-log10(0.05), 'r--'); % Add significance threshold line
legend({'State 1', 'State 2', 'State 15','FDR threshold'}, 'Location', 'northeastoutside');
grid on;
set(gca, 'XTickLabelMode', 'manual'); % Prevent dynamic updating

exportgraphics(gcf,'extendedAnatomicalComparison.png','Resolution',300);
disp('Bar plots generated successfully.');


%% Load data relative to AD-HC subjects
clear;
clc;
currDir=pwd;
currentTask={'CN','AD'};
this=prepareTrajObject(currentTask);
relExps={'C:\Code\2022_01_ADNI\22_10_25_genHMM\221027_022338_netTraining.mat'...,
    'C:\Code\2022_01_ADNI\22_07_05_genHMM_02\exp_2\exp_2.mat',...
    'C:\Code\2022_01_ADNI\22_07_05_genHMM_02\exp_3_partial\220718_110511_netTraining.mat',...
    'C:\Code\2022_01_ADNI\22_07_05_genHMM_02\exp_4\220719_010815_netTraining.mat',...
    'C:\Code\2022_01_ADNI\22_07_05_genHMM_02\exp5_partial\220719_081858_netTraining.mat',...
    'C:\Code\2022_01_ADNI\22_07_05_genHMM_02\exp6\220719_223048_netTraining.mat',...
    'C:\Code\2022_01_ADNI\22_07_05_genHMM_02\exp7\220720_054805_netTraining.mat'};

netOutCell=cell(length(relExps),1);
RBFoutCell=cell(length(relExps),1);

% currExp=1;
% fprintf('%d/%d\n',currExp,length(relExps));
% cd(fileparts(relExps{1}))
% load(relExps{currExp});
% clear netOutTrain
% clear RBFout
% blockID=zeros(size(this.lbls));
% blockID(1:25:end)=1;
% blockID=cumsum(blockID);
% for currBlock=1:max(blockID)
%     [netOutTrain(:,blockID==currBlock,:),RBFout(:,blockID==currBlock,:)]=predict(preTrndNet,dlarray(cat(4,this.feats{blockID==currBlock}),'STCB'),'Outputs',{'expLayer/out1','RBFlayer/out2'}); %#ok<SAGROW>
% end
% clear preTrndNet
% netOutCell{currExp}=extractdata(netOutTrain);
% RBFoutCell{currExp}=extractdata(RBFout);
% save("fullPopExpOuts.mat","netOutCell","RBFoutCell");

% Warning: these two blocks have to be run in separate sessions, or there
% will be an issue with class definitions for the networks

load("C:\Code\2022_01_ADNI\22_10_25_genHMM\fullPopExpOuts.mat")
for currExp=2:length(relExps)        
    fprintf('%d/%d\n',currExp,length(relExps));
    cd('C:\Code\2022_01_ADNI\22_07_05_genHMM_02')
    load(relExps{currExp});
    clear netOutTrain
    clear RBFout
    blockID=zeros(size(this.lbls));
    blockID(1:25:end)=1;
    blockID=cumsum(blockID);
    for currBlock=1:max(blockID)
        [netOutTrain(:,blockID==currBlock,:),RBFout(:,blockID==currBlock,:)]=predict(preTrndNet,dlarray(cat(4,this.feats{blockID==currBlock}),'STCB'),'Outputs',{'expLayer','RBFlayer/out2'}); %#ok<SAGROW>
    end
    clear preTrndNet
    netOutCell{currExp}=extractdata(netOutTrain);
    RBFoutCell{currExp}=extractdata(RBFout);
end

%% Recode HS clusters so that they have the same order (from more likely in
% AD to more likely in HC) for all experiments

nExps=size(relExps,2);
cAllSorted=zeros(this.nStates,nExps);
HSall=[];
cDiffAllSorted=[];
ordrAll=[];
for currExp=1:nExps
    %     HS=zeros(1,size(T1{1},2),size(T1{1},3));
    %     [~,HS(:,:,size(T1{1},3))]=max(T1{currExp}(:,:,size(T1{1},3)));
    %     for currT=size(T1{currExp},3)-1:-1:1
    %         for currSubj=1:size(HS,2)
    %             HS(:,currSubj,currT)=T2{currExp}(HS(:,currSubj,currT+1),currSubj,currT+1);
    %         end
    %     end

    [~,HS]=max(netOutCell{currExp});
    HS=squeeze(HS);
    HS1=HS(this.lbls==1,:);
    HS2=HS(this.lbls==2,:);
    c1=histcounts(reshape(HS1,[],1),.5:this.nStates+.5,'Normalization','probability');
    c2=histcounts(reshape(HS2,[],1),.5:this.nStates+.5,'Normalization','probability');
    c=histcounts(reshape(HS,[],1),.5:this.nStates+.5,'Normalization','probability');
    cDiff=c1-c2;
    [~,stateOrdr]=sort(cDiff);
    ordrAll=cat(1,stateOrdr,ordrAll);
    HSnew=zeros(size(HS));
    for currState=1:this.nStates
        HSnew(HS==stateOrdr(currState))=currState;
    end
    HSall=cat(3,HSall,HSnew);
    cAllSorted(:,currExp)=c(stateOrdr);
    cDiffAllSorted=cat(1,cDiffAllSorted,cDiff(stateOrdr));
end

%% Split subject state frequency in different cells, sorted by class
cd(currDir)
currentTask={'CN','AD'};
this=prepareTrajObject(currentTask);
this.statePall=cell(nExps,3);
for currExp=1:nExps
    this.statePall{currExp,1}=squeeze(mean(HSall(:,:,currExp)==1,2));
    this.statePall{currExp,2}=squeeze(mean(HSall(:,:,currExp)==2,2));
    this.statePall{currExp,3}=squeeze(mean(HSall(:,:,currExp)==15,2));
end

%% Compute fraction times for AD-HC subjects, for all experimtes
% HSall contains the re-ordered (accoring to AD-HC gap) clusters

% Recover numbers of all represented states
stateID=unique(HSall);
nObsStates=length(stateID);
binEdges=[stateID;stateID(end)+1]-.5;

% Now, I can compute fraction times and put them in tables
nExps=size(HSall,3);
fig5dataTable=cell(nExps,1);
P=zeros(nObsStates,nExps);
stateCounts_HC=[];
stateCounts_AD=[];
for currExp=1:nExps
    FT=zeros(size(HSall,1),nObsStates);
    for currSubj=1:size(HSall,1)
        FT(currSubj,:)=histcounts(HSall(currSubj,:,currExp),binEdges,'Normalization','probability');
    end

    % Perform significance tests. Use Benjamini-Hochberg multiple
    % comparison correction
    group=this.lbls==2;
    for currState=1:nObsStates
        P(currState,currExp)=ranksum(FT(group,currState),FT(~group,currState));
    end
    [c_pvalues, c_alpha]=fdr_BH(P,.05,false);
    C_P=reshape(c_pvalues,size(P));

    fig5dataTable{currExp}=array2table([FT,group],'VariableNames',[string(stateID);{'Group'}]);
    writetable(fig5dataTable{currExp},sprintf('fig5dataTableExp%d.txt',currExp));
    writetable(array2table(C_P),'fig5Pvals.txt');
    stateCounts_HC=cat(1,stateCounts_HC,histcounts(HSall(this.lbls==1,:,currExp)));
    stateCounts_AD=cat(1,stateCounts_AD,histcounts(HSall(this.lbls==2,:,currExp)));
end
writetable(array2table(stateCounts_HC),'fig5countsHC.txt')
writetable(array2table(stateCounts_AD),'fig5countsAD.txt')

%% Perform permutation test
% Test if cluster attributions are significantly better than chance, across
% different training runs
obsCounts=zeros(nObsStates,1);
for currState=1:nObsStates
    obsCounts(currState)=sum(prod(HSall==stateID(currState),3),'all');
end

% Perform 10K permutations, log results
nPerms=1e4;
permCounts=zeros(nObsStates,nPerms);
for currPerm=1:nPerms
    HSperm=HSall(randperm(numel(HSall)));
    HSperm=reshape(HSperm,size(HSall));
    for currState=1:nObsStates
        permCounts(currState,currPerm)=sum(prod(HSperm==stateID(currState),3),'all');
    end
    if mod(currPerm,100)==0
        fprintf('%d/%d\n',currPerm,nPerms)
    end
end

% Build table
permDataTable=array2table([obsCounts';mean(permCounts,2)';std(permCounts,[],2)'],'RowNames',...
    {'Observed counts','Mean perm. counts','SD perm. counts'},'VariableNames',...
    string(stateID));

save("permTestData_7exps.mat","obsCounts","permCounts");

%% Response to reviewer 3: permutation histogram
clear
close all
clc

load("permTestData_7exps.mat","obsCounts","permCounts");

% Loop over each target state
figure;
targetStates=[1,2,10];
for i = 1:length(targetStates)
    stateIdx = targetStates(i);
    subplot(3, 1, i)
    
    % Extract permutation data for this state
    thisPermData = permCounts(stateIdx, :);
    
    % Plot histogram of permuted counts
    histogram(thisPermData, 'FaceColor', [0.6 0.8 1], 'EdgeColor', 'black', 'Normalization','probability');
    hold on;
    
    % Plot observed count
    xline(obsCounts(stateIdx), 'r', 'LineWidth', 2);
    
    set(gca, 'XScale', 'log','XLim',[1 1e5],'YLim',[0 1]);  % Set x-axis to logarithmic scale

%     set(gca,'YLim',[0 1]);  % Set x-axis to logarithmic scale
    title(sprintf('State %d', stateIdx));
    xlabel('Permutation Counts');
    ylabel('Frequency');
    grid on;
end

sgtitle('Permutation Test Results with Observed Count Overlay');
exportgraphics(gcf,'permTestData_7exps.png','Resolution',300);

%% Response to reviewer 2: count-matched permutation test
% Perform permutation test with matched subgroup sizes

% Identify group indices
HC_indices = find(this.lbls == 1);
AD_indices = find(this.lbls == 2);

% Get AD group size (smaller group)
nAD = numel(AD_indices);

% Observed counts using matched subgroups
nObsStates = numel(stateID);
obsCounts = zeros(nObsStates, nPerms);

% Perform permutations
nPerms = 1e4;
permCounts = zeros(nObsStates, nPerms);

for currPerm = 1:nPerms
    % Randomly subsample HC group for each permutation
    randHC_subset = randsample(HC_indices, nAD, false);
    sub_HSall_HC_perm = HSall(randHC_subset, :, :);
    sub_HSall_AD_perm = HSall(AD_indices, :, :);

    % Concatenate HC and AD subgroups
    sub_HSall_perm = cat(1, sub_HSall_AD_perm, sub_HSall_HC_perm);

    % Shuffle the state assignments
    HSperm = sub_HSall_perm(randperm(numel(sub_HSall_perm)));
    HSperm = reshape(HSperm, size(sub_HSall_perm));

    for currState = 1:nObsStates
        obsCounts(currState, currPerm) = sum(prod(sub_HSall_perm == stateID(currState), 3), 'all');
        permCounts(currState, currPerm) = sum(prod(HSperm == stateID(currState), 3), 'all');
    end

    if mod(currPerm, 100) == 0
        fprintf('%d/%d permutations done\n', currPerm, nPerms);
    end
end

% Build result table
permDataTableMatched = array2table([mean(obsCounts, 2)'; std(obsCounts, [], 2)'; mean(permCounts, 2)'; std(permCounts, [], 2)'], ...
    'RowNames', {'Mean observed counts (matched)', 'SD observed counts (matched)', 'Mean perm. counts (matched)', 'SD perm. counts (matched)'}, ...
    'VariableNames', string(stateID));
writetable(permDataTableMatched,'permDataCountMatched.txt');

%% Respone to reviewer 3, transition matrices
clear
close all
clc

load("HS_fullPop.mat","HS_fullPop");

% Reorder labels
cLbls=zeros(size(this.lbls));
lblVals=unique(this.lbls);
nClasses=length(lblVals);
sampleRep=1;
lblOrdr=[1,5,3,4,2];
for currClass=1:nClasses
    cLbls(this.lbls==lblOrdr(currClass))=currClass;
end

% Compute full transition matrices
transition_counts=zeros(max(HS_fullPop,[],'all'),max(HS_fullPop,[],'all'),length(lblOrdr));
for currGroup=1:length(lblOrdr)
    groupStates=HS_fullPop(cLbls==currGroup,:);
    for currSubj = 1:size(groupStates,1)
        for currTime = 1:size(groupStates,2)-1
            transition_counts(groupStates(currSubj,currTime),groupStates(currSubj,currTime+1),currGroup)=transition_counts(groupStates(currSubj,currTime),groupStates(currSubj,currTime+1),currGroup)+1;
        end
    end
end

% Define state mapping: map 1 -> 1, 2 -> 2, 15 -> 4, others -> 3
stateMap = zeros(15,1);
stateMap([1,2,15]) = [1,2,4];  % Map key states to 1-4
stateMap(stateMap == 0) = 3;   % All others become 2

nGroups = size(transition_counts, 3);
collapsed_counts = zeros(4, 4, nGroups);  % New 4x4 matrices

for g = 1:nGroups
    % Get group matrix
    groupMatrix = transition_counts(:,:,g);
    
    % Prepare flattened transitions
    [from, to, count] = find(groupMatrix);
    
    % Map old states to new ones
    newFrom = stateMap(from);
    newTo = stateMap(to);
    
    % Accumulate counts in the reduced space
    collapsed_counts(:,:,g) = accumarray([newFrom, newTo], count, [4,4]);
end

labels = {'HC','SMC','eMCI','lMCI','AD'};
stateLabels = {'From S1','From S2','From Rest','From S15'};
colLabels = {'To S1','To S2','To Rest','To S15'};

figure;
for currState1 = 1:4
    for currState2 = 1:4
        idx = (currState1 - 1) * 4 + currState2;
        subplot(4,4,idx);
        
        % Plot groupwise proportions for this transition
        plot(squeeze(collapsed_counts(currState1,currState2,:) ./ ...
            sum(collapsed_counts(:,currState2,:),1)),'o-')
        
        set(gca, 'XTick', 1:5, 'XTickLabel', labels, 'FontSize', 8);
        
        % Add row and column labels
        if currState1 == 1
            title(colLabels{currState2}, 'FontWeight', 'bold');
        end
        if currState2 == 1
            ylabel(stateLabels{currState1}, 'FontWeight', 'bold');
        end
    end
end