%%%%%%%%%%%%%%% Estimation of GC, PDC and DTF using Parametric Measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oNum=0; refType = 'unipolar'; % unipolar', 'avg','bipolar','csd'
[data,layout] = getData(oNum, refType);

%%%%%%%%%%% Data segmentation into baseline or stimulus epoch %%%%%%%%%%%%%
cfg          = [];
% cfg.toilim   = [-0.5 0];  % Baseline time period -0.5 to 0
cfg.toilim   = [0.25 0.75]; % Stimulus time period 0.25 to 0.75
data_short = ft_redefinetrial(cfg, data);

%%%%%%%% Multivariate Autoregressive Modeling of the Epoched Data %%%%%%%%%
cfg       = [];
cfg.order = 40; % Fs = 250 Hz, selecting 150ms duration for autoregressive modeling (Order = 0.150/(1/250) = 37.5)
data_mvar  = ft_mvaranalysis(cfg,data_short);

%%%%%%%%%%%%%%%%%%%%%%%%%% Frequency Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg        = [];
cfg.method = 'mvar';
cfg.foi    = 0:50;

data_freq  = ft_freqanalysis(cfg,data_mvar);

%%%%%%%%%%%%%%%%%%%% Display Directed Connectivity Measures %%%%%%%%%%%%%%%
method = 'granger'; % for granger, dtf, pdc

cfgA = []; % Configuration file for analysis
cgfD = []; % Configuration file for display
cfgD.refchannel     = 'PO3';
cfgD.layout         = layout;
cfgD.showlabels     = 'yes';
cfgD.showoutline    = 'yes';
cfgD.directionality = 'outflow'; %inflow or outflow; Default is inflow
cfgD.newfigure      = 'no';

if strcmp(method,'granger') % For Granger Causality
    cfgA.method     = 'granger';
    cfgD.parameter  = 'grangerspctrm';
    
elseif strcmp(method,'dtf')  % For Directed Transfer Function
    cfgA.method     = 'dtf';
    cfgD.parameter  = 'dtfspctrm';

elseif strcmp(method,'pdc') % For Partial Directed Coherence
    cfgA.method     = 'pdc';
    cfgD.parameter  = 'pdcspctrm';
end

con_result = ft_connectivityanalysis(cfgA,data_freq);
ft_multiplotER(cfgD,con_result);