function this=prepareTrajObject(currentTask)
    % Load database
    FFP='C:\Data\2022_01_ADNI_timeSeries\TimeSeries';
    DP1='C:\Data\2022_01_ADNI_timeSeries\ADNI_data_summary_02-05-2022.xlsx';
    DP2='C:\Data\2022_01_ADNI_timeSeries\ADNIMERGE_Sept2022.csv';
    DP3='C:\Data\2022_01_ADNI_timeSeries\ADNI_ATN_info.xlsx';
    DP4='C:\Data\2022_01_ADNI_timeSeries\freesurfer_features_02112022.xlsx';
    this=traj(FFP,DP1,DP2,DP3,DP4);
    this.loadRawData;
    % this.clinInfo=cat(2,this.clinInfo{:});
    this.importATNdata;
    this.importFreesurferData;
    
    % Define task
    relIdx=ismember(this.lbls,currentTask);
    this.classNames=currentTask;
    this.subjID=this.subjID(relIdx);
    this.lbls=this.lbls(relIdx);
    this.clinInfo=this.clinInfo(relIdx,:);
    this.subjAge=this.subjAge(relIdx);
    this.subjSex=this.subjSex(relIdx);
    this.APOE4=this.APOE4(relIdx);
    lblsTemp=zeros(size(this.lbls));
    for currClass=1:length(currentTask)
        lblsTemp(ismember(this.lbls,currentTask{currClass}))=currClass;
    end
    this.lbls=lblsTemp;
    this.timeSeries=this.timeSeries(relIdx);
    
    % Add some degree of noise to constant time series, then
    % multiply channels
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
    this.clinInfo(toBeRemoved,:)=[];
    this.subjID(toBeRemoved)=[];
    this.subjAge(toBeRemoved)=[];
    this.subjSex(toBeRemoved)=[];
    this.APOE4(toBeRemoved)=[];
    this.feats=cellfun(@(x)x(:,1:seqLength),this.feats,'UniformOutput',false);
    
    % Define number of HS
    this.nStates=15;