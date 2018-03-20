%% Compare Bss
% clear
clc
% close all
NDatasets=9;
DatasetPath= '\New';

%% FastIca
BSSfunc = @fastica;
BssInArgs = {'verbose', 'off',...
          'displayMode', 'off',...
          'approach', 'symm',...
          'g', 'tanh',...
          'initGuess',};
BssOutArgs={'A', 'W'};
      
ResultsStruc1 = TestBSS(BSSfunc,BssInArgs,BssOutArgs,DatasetPath);
ResultsStruc1.BSSfun='FastICA symm';
save('\Results\ResultsStruc1.mat','-struct','ResultsStruc1')

          %% CCA
BSSfunc = @ccabss;
BssInArgs = {'W_init'};
BssOutArgs={'~','W'};
      
ResultsStruc2 = TestBSS(BSSfunc,BssInArgs,BssOutArgs,DatasetPath);
ResultsStruc2.BSSfun='CCA';
save('\Results\ResultsStruc2.mat','-struct','ResultsStruc2')

       %% InfoMax   
BSSfunc = @runica;
BssInArgs = {'weights'};
BssOutArgs={'W','~'};
      
ResultsStruc3 = TestBSS(BSSfunc,BssInArgs,BssOutArgs,DatasetPath);
ResultsStruc3.BSSfun='InfoMax';
save('\Results\ResultsStruc3.mat','-struct','ResultsStruc3')

                  %% JADE   
% BSSfunc = @jadeR;
% BssInArgs = {};
% BssOutArgs={'W'};
%       
% ResultsStruc5 = TestBSS(BSSfunc,BssInArgs,BssOutArgs,DatasetPath,'PCA',true); 

%% IVA
BSSfunc = @iva_second_orderOriginal;
BssInArgs = {'verbose',true,...
    'opt_approach','quasi',...
    'maxIter',5000,...
    'WDiffStop',1e-6,...
    'W_init'};

BssOutArgs={'W'};
      
ResultsStruc4 = TestBSS(BSSfunc,BssInArgs,BssOutArgs,DatasetPath,...
    'PCA',true,'Nshift',5,'Lag',4,'ExplainedVar',99.9,'ShiftDataset',true,'Filter',true);

ResultsStruc4.BSSfun='IVA Gauss';
save('\Results\ResultsStruc4.mat','-struct','ResultsStruc4')

%% FastIca
BSSfunc = @fastica;
BssInArgs = {'verbose', 'off',...
          'displayMode', 'off',...
          'approach', 'defl',...
          'g', 'tanh',...
          'initGuess'};
BssOutArgs={'A', 'W'};
      
ResultsStruc5 = TestBSS(BSSfunc,BssInArgs,BssOutArgs,DatasetPath);
ResultsStruc5.BSSfun='FastICA defl';
save('\Results\ResultsStruc5.mat','-struct','ResultsStruc5')

%% IVA
BSSfunc = @iva_second_orderOriginal;
BssInArgs = {'verbose',false,...
    'opt_approach','quasi',...
    'maxIter',5000,...
    'WDiffStop',1e-6,...
    'W_init'};

BssOutArgs={'W'};
      
ResultsStruc6 = TestBSS(BSSfunc,BssInArgs,BssOutArgs,DatasetPath,...
    'PCA',true,'Nshift',5,'Lag',4,'ExplainedVar',99.999,'ShiftDataset',true,'Filter',false);

ResultsStruc6.BSSfun='IVA Gauss no filter';
save('\Results\ResultsStruc6.mat','-struct','ResultsStruc6')
