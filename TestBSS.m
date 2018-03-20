function ResultsStruc=TestBSS(BSSfunc,BssArgs,BssOutArgs,DatasetPath,varargin)
%% Function to test BSS methods on multiple datassets. 
% Input arguments:
% BSSfunc, fanction handle ti the BSS to test
% BssArgs, cell array with all additinal arguments for BSSfunc
% NDatasets, Number of testing datasets
% DatasetPath, path to a folder with all the testing datasets
% varargin, additional parameters defined below.
%
% 'PCA', true or [false] performs dimensionality reduction with PCA before
%                           running the algorithm
% 'ExplainedVar', int between 0 and 100. [99] % Amount of explained 
%                  variance at which the PCA reduces the dataset.
% 'ShiftDataset', true or [false] creates a 3D matrix containing the
%                   dataset and shifted copies of himself. Used in IVA.
% 'Nshift', integer is the number of shifted copies. See bove.
% 'Lag', integer. Every shifted copy is shifted by this number of samples
%
%
% Output argument:
% ResultsStruc, structure with following fields:
% 
%     'CorrScore',      mean of all corrscores
%     'CorrScoreSTD'    std of all corrscores
%     'rrmse'           mean of all relative root mean square errors
%     'rrmseSTD'        std of all relative root mean square errors
%     'NumShift'        Nummber of shifted datasets
%     'SNR'             Signal Noise Ratio
%     'time'            mean of execution time in seconds
%     'timeSTD'         std of execution time in seconds
%     'W'               Cell array eith all W
% 
% Usage example:
%
% BSSfunc = @fastica;
% BssArgs = {'verbose', 'on',...
%           'displayMode', 'off',...
%           'approach', 'symm',...
%           'g', 'tanh'};
%
% NDatasets=10
% DatasetNameAndPath= '/User/Pippo/DatasetsDirectory';
% ResultsStruc = TestBSS(BSSfunc,BssArgs,NDatasets,...
%               DatasetPath,ResultsStrucName,...
%               'PCA',true,'ShiftDataset',true,'ExplainedVar',99.9);
%
% Copyright (c) 2017, Federico Barban



Params=struct('PCA',false,...
    'ExplainedVar',99.9,...
    'ShiftDataset',false,...
    'Nshift',7,...
    'Lag',4,...
   'Filter',false);
   
Params = getopt(Params,varargin{:});

fileNames=dir([DatasetPath '/*.mat']);
NDatasets=length(fileNames);
lag=Params.Lag;
SNRs=zeros(1,NDatasets);
CorrScore = zeros(1,NDatasets);
CSM = zeros(1,NDatasets);
RmsE = zeros(1,NDatasets);
time = zeros(1,NDatasets);
NArtifactComp = zeros(1,NDatasets);
Wall=cell(1,NDatasets);
Windex=strcmp(BssOutArgs,'W');
Aindex=strcmp(BssOutArgs,'A');

rms=@(x) sqrt(sum(sum(x.^2))/numel(x));


%%

    for k=1:NDatasets % datasets, here shifts of the same dataset
        %%
        OutArgs=cell(1,length(BssOutArgs));
        FName=[fileNames(k).folder filesep fileNames(k).name];
        D=load(FName);      
        EEG=double(D.data);
        Artifacts=D.artifact;
        Nsamples=length(EEG);
        
        % plotchanels(data_seg_all.artifact)
       
        
        %% Dimensionality reduction
        if Params.PCA
            mu=mean(EEG,2)';
            [wcoeff,~,latent]=pca(EEG');
            ExpVar=100*cumsum(latent)./max(cumsum(latent));
            NComponent=find(ExpVar>Params.ExplainedVar ,1);
            EEG_rid=(EEG'*wcoeff(:,1:NComponent))';
        else
            NComponent=size(EEG,1);
            wcoeff=eye(NComponent);
            mu=zeros(1,NComponent);
            EEG_rid=EEG;
        end
        
        if Params.Filter
            [EEG_fil,~]=process_emg(EEG_rid,'range',[60,350]);
            EEG_fil=EEG_fil';
        else
            EEG_fil=EEG_rid;
        end
        
        %% Shift
        if Params.ShiftDataset
            L=Params.Nshift;
            Xmat=zeros(NComponent,floor((Nsamples-lag*L)),L);
            % build shifted dataset
            for h=1:L
                lagI=lag*(h-1);
                lagE=lag*(L-h+1);
                Xmat(:,:,h)=EEG_fil(:,1+lagI:end-lagE);
                W_init(:,:,h)=eye(size(Xmat,1));
            end %h
        else
            Xmat=EEG_fil;
            W_init=eye(size(Xmat,1));
%              clear('EEG_rid');
        end

%%
        BssArgs_=[BssArgs {W_init}];
        tic;
        [OutArgs{:}]=BSSfunc(Xmat,BssArgs_{:});
        time(k)=toc;
        W=[OutArgs{Windex} ];
        A=[OutArgs{Aindex}];
        
        Wall{k}=W;
        %% find independent sources IVA G
        S=real(W(:,:,1)*EEG_rid);

        %% find artefacts and evaluate performance
        
        %%IVAG
        NIC=size(S,1);
        
            C=abs(corr(real(S)',Artifacts'));
            [~,pos]=sort(C(:),'descend');
            pos=mod(pos,NIC);
            pos(pos==0)=NIC;
            AllRms=zeros(1,NIC);
            AllCorrscore=zeros(1,NIC);
            AllCsm=zeros(1,NIC);
            for i=1:NIC
                pos=unique(pos,'stable'); 
                if isempty(A)
                    A=inv((W(:,:,1)));
                end
                Cdata = D.cleandata;
                % reconstruct clean eeg
                eeg_ricostr = (EEG - ...
                    (wcoeff(:,1:NComponent) * real(A(:,pos(1:i))*S(pos(1:i),:,1))) + ...
                    repmat(mu,Nsamples,1)');
                
                AllRms(i)=rms(eeg_ricostr-Cdata)...
                    /rms(Cdata);
                AllCorrscore(i)=mean(diag(corr(eeg_ricostr',Cdata')));
                AllCsm(i)=mean(csm(eeg_ricostr,Cdata));
            end
            
           [~,NArtifactComp(k)] = max(AllCorrscore); 

            
           [~,ind]= min(AllRms);
           RmsE(k) = AllRms(ind);
            % performance
            CorrScore(k) = AllCorrscore(NArtifactComp(k));
            CSM(k)=AllCsm(ind);
           SNRs(k)=D.snr;         
    end %k
    %%    
    ResultsStruc.AllCorrScore=CorrScore;
    ResultsStruc.Allrrmse=RmsE;
    ResultsStruc.AllTime=time;
    ResultsStruc.NumShift=Params.Nshift;
    ResultsStruc.SNR=SNRs;
    ResultsStruc.W=Wall;
    ResultsStruc.NArtifactComp=NArtifactComp;
    ResultsStruc.Lag=Params.Lag;
    ResultsStruc.BSSfun=BSSfunc;
    ResultsStruc.Csm=CSM;
    ResultsStruc.Filter=Params.Filter;
    
    

%%

% save([ResultsStrucName...
%     num2str(length(dir([ResultsStrucName '*']))+1)...
%     '.mat'],'ResultsStruc')

function properties = getopt(properties,varargin)
%GETOPT - Process paired optional arguments as 'prop1',val1,'prop2',val2,...
%
%   getopt(properties,varargin) returns a modified properties structure,
%   given an initial properties structure, and a list of paired arguments.
%   Each argumnet pair should be of the form property_name,val where
%   property_name is the name of one of the field in properties, and val is
%   the value to be assigned to that structure field.
%
%   No validation of the values is performed.
%%
% EXAMPLE:
%   properties = struct('zoom',1.0,'aspect',1.0,'gamma',1.0,'file',[],'bg',[]);
%   properties = getopt(properties,'aspect',0.76,'file','mydata.dat')
% would return:
%   properties =
%         zoom: 1
%       aspect: 0.7600
%        gamma: 1
%         file: 'mydata.dat'
%           bg: []
%
% Typical usage in a function:
%   properties = getopt(properties,varargin{:})

% Function from
% http://mathforum.org/epigone/comp.soft-sys.matlab/sloasmirsmon/bp0ndp$crq5@cui1.lmms.lmco.com

% dgleich
% 2003-11-19
% Added ability to pass a cell array of properties

if ~isempty(varargin) && (iscell(varargin{1}))
   varargin = varargin{1};
end

% Process the properties (optional input arguments)
prop_names = fieldnames(properties);
TargetField = [];
for ii=1:length(varargin)
   arg = varargin{ii};
   if isempty(TargetField)
      if ~ischar(arg)
         error('Property names must be character strings');
      end
      %f = find(strcmp(prop_names, arg));
      if isempty(find(strcmp(prop_names, arg),1)) %length(f) == 0
         error('%s ',['invalid property ''',arg,'''; must be one of:'],prop_names{:});
      end
      TargetField = arg;
   else
      properties.(TargetField) = arg;
      TargetField = '';
   end
end
if ~isempty(TargetField)
   error('Property names and values must be specified in pairs.');
end