%%%%%%%%%%%%%%%% Display ERP and TF Data using fieldtrip %%%%%%%%%%%%%%%%%%
oNum=3; refType = 'unipolar'; % unipolar', 'avg','bipolar','csd'
[data,layout] = getData(oNum, refType);

%%%%%%%%% Data segmentation into baseline or stimulus epoch %%%%%%%%%%%%%%%

cfg        = [];
%cfg.toilim = [-0.5 0];  % Baseline time period -0.5 to 0
cfg.toilim = [0.25 0.75]; % Stimulus time period 0.25 to 0.75
data_short = ft_redefinetrial(cfg, data);

%%%%%%%%%%%%%%%%%%%%%%%% Frequency analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg                 = [];
cfg.output          = 'powandcsd';
cfg.method          = 'mtmfft';
cfg.taper           = 'dpss';
cfg.tapsmofrq       = 4;  % Smoothing over W = +- 4 Hz. For T=0.5, W=4, TW=2, 3 tapers are needed 
cfg.foilim          = [0 50];
cfg.keeptrials      = 'yes';

data_freq           = ft_freqanalysis(cfg, data_short);

% Display
%cfg = [];
%cfg.layout        = layout;
%cfg.showlabels    = 'yes';
%cfg.showoutline   = 'yes';
%figure; ft_multiplotER(cfg,data_freq);

%%%%%%%%%%%%%%%%%%%%%%%% Connectivity Measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%
method = 'coh'; % choose between coh, imc, wpli or ppc

cfgA = []; % Configuration file for analysis
cgfD = []; % Configuration file for display
cfgD.refchannel    = 'PO3';
cfgD.layout        = layout;
cfgD.showlabels    = 'yes';
cfgD.showoutline   = 'yes';

if strcmp(method,'coh')
    cfgA.method = 'coh';
    cfgA.complex = 'abs';
    cfgD.parameter = 'cohspctrm';
    
elseif strcmp(method,'imc')
    cfgA.method = 'coh';
    cfgA.complex = 'imag';
    cfgD.parameter = 'cohspctrm';
    
elseif strcmp(method,'wpli')
    cfgA.method = 'wpli';
    cfgD.parameter = 'wplispctrm';
    
elseif strcmp(method,'ppc')
    cfgA.method = 'ppc';
    cfgD.parameter = 'ppcspctrm';
end

con_result = ft_connectivityanalysis(cfgA,data_freq);
ft_multiplotER(cfgD,con_result);