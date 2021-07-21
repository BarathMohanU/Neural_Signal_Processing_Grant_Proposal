%%%%%%%%%%%%%%%% Display ERP and TF Data using fieldtrip %%%%%%%%%%%%%%%%%%
oNum=0; refType = 'unipolar'; % 'avg','bipolar','csd'
[data,layout] = getData(oNum,refType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ERP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg = [];
cfg.layout        = layout;
cfg.showlabels    = 'yes';
cfg.showoutline   = 'yes';
figure; ft_multiplotER(cfg,data);

%%%%%%%%%%%%%%%%%%%%%%% Time-frequency representation %%%%%%%%%%%%%%%%%%%%%
windowLenS       = 0.25; % Seconds
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'all';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = (1/windowLenS):(1/windowLenS):100; 
cfg.t_ftimwin    = ones(length(cfg.foi),1).* windowLenS;   % length of time window
cfg.toi          = -0.5:0.05:0.75;

TFR = ft_freqanalysis(cfg,data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Visualizing TFR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg = [];
cfg.baseline     = [-0.5 0];
cfg.baselinetype = 'db';
cfg.showlabels   = 'yes';
cfg.layout       = layout;
cfg.box          = 'yes';
cfg.colormap     = 'jet';
figure; ft_multiplotTFR(cfg,TFR);