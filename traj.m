classdef traj < handle
    % January 2022, Jacopo Tessadori
    
    properties
        fileFolderPath;
        descrPath1;
        descrPath2;
        clinInfo;
        timeSeries;
        FC;
        subjID;
        lbls;
        nROIs;
        nSubjs;
        nDims;
        nStates;
        RS;
        feats;
        subjAge;
    end
    
    methods
        function this=traj(FFP,DP1,DP2)
            % Constructor for traj class.
            % Required arguments are paths to data folder and to
            % description file
            this.fileFolderPath=FFP;
            this.descrPath1=DP1;
            this.descrPath2=DP2;
        end
        
        function loadRawData(this)
            %% Import database files (importfile and importfile2 are wizard-generated)
            DB1=importfile2(this.descrPath1, "ADNI3", [2, 1045]);
            DB2=importfile(this.descrPath2,[2 Inf]);

            % Generate correct labels and load patient data
            D=dir(this.fileFolderPath);
            fileNames=cat(1,{D.name});
            this.timeSeries=cell(0);
            this.lbls=cell(0);
            this.subjID=cell(0);
            this.subjAge=[];
            relIdxs={'ADAS11','MMSE','FAQ','RAVLT_immediate','RAVLT_learning','RAVLT_forgetting','CDRSB'};
            this.clinInfo=cell(length(relIdxs),1);
            for currSubj=1:height(DB1)
                patId=DB1.ID{currSubj};
                subjCond=DB1.ADNIMERGE(currSubj);
                if sum(contains(fileNames,patId))&&any(ismember({'CN','AD','MCI','EMCI','LMCI','SMC'},subjCond))
                    patIdx=fileNames{contains(fileNames,patId)};

                    % Recover entries for current subject
                    DB1curr=DB1(DB1.ID==patId,:);
                    DB2curr=DB2(DB2.PTID==patId,:);

                    % Recover classification in DB1
                    class1=DB1curr.ADNIMERGE;

                    % Recover first entry in ADNI3 for DB2
                    relEntry=find(DB2curr.COLPROT=="ADNI3",1,'first');

                    % Recover classification in first ADNI3 entry
                    class2=DB2curr.DX_bl(relEntry);

                    % Resolve conflicts
                    if isempty(class2)
                        subjCond=class1;
                    elseif isundefined(class2)
                        subjCond=class1;
                    elseif class2~=class1
                        subjCond=class2;
                    else 
                        subjCond=class2;
                    end

                    % Update time series, labels and additional clinical
                    % info
                    this.lbls=cat(1,this.lbls,subjCond);
                    subjData=load(fullfile(this.fileFolderPath,patIdx));
                    this.timeSeries=cat(1,this.timeSeries,{subjData.ms});
                    this.subjID=cat(1,this.subjID,DB1.ID(currSubj));
                    this.subjAge=cat(1,this.subjAge,DB1.AGE(currSubj));
                    for currIdx=1:length(relIdxs)
                        if isempty(DB2curr)
                            this.clinInfo{currIdx}=cat(1,this.clinInfo{currIdx},nan);
                        else
                            this.clinInfo{currIdx}=cat(1,this.clinInfo{currIdx},DB2curr.(relIdxs{currIdx})(relEntry));
                        end
                    end
                end
            end
        end

        function untrndNet=designNet(this)
            % Define network
            this.nStates=15;
            seqLength=min(cellfun(@(x)size(x,2),this.feats));
            N=100;
            M1=16;
            M2=32;
            L1=20;
            L2=1;
            %% Create Layer Graph
            % Create the layer graph variable to contain the network layers.
            lgraph = layerGraph();
            %% Add Layer Branches
            % Add the branches of the network to the layer graph. Each branch is a linear
            % array of layers.

            tempLayers=[sequenceInputLayer([size(this.feats{1},1),1],"Name","featureinput","MinLength",seqLength)
                customConvLayer([1,L1],M1,1,"STCU","timeCorrConv",[1,1])
                leakyReluLayer(.5,"Name","relu1")
                customConvLayer([N.^2,L2],M2,M1,"STCU","spaceConv1",[1,1])
                leakyReluLayer(.5,"Name","relu2")
                fullyConnectedLayer(3,"Name","fc_1")
                batchNormalizationLayer("Name","bNorm")
                GRBFlayer(this.nStates,3,"RBFlayer")
                expectationLayer("expLayer")];
%             tempLayers(end-1) = setLearnRateFactor(tempLayers(end-1),"c",10);
%             tempLayers(end-1) = setLearnRateFactor(tempLayers(end-1),"sigma",10);

            lgraph = addLayers(lgraph,tempLayers);
            untrndNet=lgraph;
        end

        function preTrndNet=preTrainNet(this)
            addpath C:\Code\Common
            % Add some degree of noise to constant time series, then
            % multiply channels
            if isempty(this.feats)||size(this.feats{1},1)~=1e4
                for currSubj=1:length(this.lbls)
                    timeSeriesSD=std(this.timeSeries{currSubj});
                    this.timeSeries{currSubj}(:,timeSeriesSD==0)=randn(size(this.timeSeries{currSubj}(:,timeSeriesSD==0)))*1e-6;
                end
                this.feats=cellfun(@(x)normalize(x)',this.timeSeries,'UniformOutput',false);
                this.feats=cellfun(@(x)reshape(permute(x,[2,1,3]).*permute(x,[2,3,1]),size(x,2),[])',this.feats,'UniformOutput',false);

                % Compute acceptable time length, remove shorter sequences and
                % trim the rest
                seqLength=prctile(cellfun(@(x)size(x,2),this.feats),5);
                toBeRemoved=cellfun(@(x)size(x,2)<seqLength,this.feats);
                this.lbls(toBeRemoved)=[];
                this.feats(toBeRemoved)=[];
                this.timeSeries(toBeRemoved)=[];
                this.feats=cellfun(@(x)x(:,1:seqLength),this.feats,'UniformOutput',false);
            end

            %% Pre-train network on all available data
            % Initialie network
            untrndNet=this.designNet;
            batchSize=ceil(this.nSubjs/36);
            maxEpochs=15;
            iteration = 0;
            learnRate = 1e-3;
            gradDecay = 0.75;
            sqGradDecay = 0.95;
            dlnet = dlnetwork(untrndNet);

            start=tic;
            for epoch = 1:maxEpochs
                % Shuffle data.
                rndOrdr=randperm(this.nSubjs);

                % Intialize useful stuff
                lastEntry=1;

                % Loop over mini-batches.
                while lastEntry+batchSize<=length(rndOrdr)
                    iteration = iteration + 1;

                    % Read mini-batch of data and add repetition-pad to
                    % ensure same-length sequences
                    relDataIdx=rndOrdr(lastEntry+1:lastEntry+batchSize);
                    tempTS=cat(4,this.feats{relDataIdx});
                    dlX=dlarray(tempTS,'STCB');
                    lastEntry=lastEntry+batchSize;

                    % Evaluate the model gradients, state, and loss using dlfeval and the
                    % modelGradients function and update the network state.
                    [grad,lossData,lgnd,state] = dlfeval(@traj.prototypeGradients,dlnet,dlX,[],1);
                    dlnet.State=state;
                    if ~exist('preTrainLog','var')
                        preTrainLog=lossData;
                    else
                        preTrainLog=cat(1,preTrainLog,lossData);
                    end

                    % Update the network parameters using the ADAM optimizer.
                    if ~exist('avg_g','var')
                        [dlnet,avg_g,avg_sqg] = adamupdate(dlnet,grad,[],[],1,learnRate,gradDecay,sqGradDecay);
                    else
                        [dlnet,avg_g,avg_sqg] = adamupdate(dlnet,grad,avg_g,avg_sqg,iteration,learnRate,gradDecay,sqGradDecay);
                    end

                    % Display the training progress.
                    if iteration>1
                        figure(1);
                        hold off;
                        D = duration(0,0,toc(start),'Format','hh:mm:ss');
                        h=plot(preTrainLog);
                        h(1).LineWidth=2;
                        xlabel("Iteration")
                        ylabel("Loss")
                        grid on
                        title("Epoch: " + epoch + ", Elapsed: " + string(D))
                        legend(lgnd,'Location','southwest');
                        drawnow
                    end
                end
            end

            % Save training graph
            DT=datetime('now','Format','yyMMdd_HHmmss');
            savefig(gcf,sprintf('%s_netTraining_fold.fig',DT));

            lgraph=dlnet.layerGraph;
            preTrndNet=dlnetwork(lgraph);

            % Save networks
            save(sprintf('%s_netTraining.mat',DT),'preTrndNet');
        end

        function n=get.nROIs(this)
            n=size(this.timeSeries{1},2);
        end
        
        function n=get.nSubjs(this)
            n=length(this.lbls);
        end
                
        function nD=get.nDims(this)
            nD=size(this.feats{1},1);
        end
    end
end