% Source localization methods using fieldtrip

%% Step 1: Get standard MRI
% We can use data already available in fieldtrip toolbox in the folder
% template/anatomy. However, this is not in standard ctf format and needs
% to be realigned. For this realignment procedure, see the tutorial:
% http://www.fieldtriptoolbox.org/workshop/baci2017/forwardproblem/  

% Since this realignment is subjective, we instead download a standard MRI dataset from
% ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip/tutorial/Subject01.zip,
% which is in the standard ctf format. Since the raw MRI data is large, we instead save the mri data in mri.mat file

% mri = ft_read_mri('Subject01\Subject01.mri');
x = load('mri_orig.mat');
mri_orig = x.mri;

% Reslice to get in standard coordinates
cfg = [];
mri_resliced_orig = ft_volumereslice(cfg, mri_orig);
ft_sourceplot(cfg, mri_resliced_orig);

%% Segment the MRI into compartments
% cfg = [];
% cfg.output    = {'brain','skull','scalp'};
% mri_segmented  = ft_volumesegment(cfg, mri_orig);
% save mri_segmented.mat mri_segmented
load mri_segmented.mat

% Reslice for visualization
cfg = [];
mri_resliced_segmented = ft_volumereslice(cfg, mri_segmented);

% View the segmented MRI
seg_i = ft_datatype_segmentation(mri_resliced_segmented,'segmentationstyle','indexed');

cfg              = [];
cfg.funparameter = 'seg';
cfg.funcolormap  = gray(4); % distinct color per tissue
cfg.location     = 'center';
cfg.atlas        = seg_i;
ft_sourceplot(cfg, seg_i);

%% Step 2: create a mesh that will be used to get a headmodel
cfg=[];
cfg.tissue={'brain','skull','scalp'};
cfg.numvertices = [1000 1000 1000];
bnd=ft_prepare_mesh(cfg,mri_segmented);

% Display the different meshes
% figure; 
% subplot(221); ft_plot_mesh(bnd(1));
% subplot(222); ft_plot_mesh(bnd(2));
% subplot(223); ft_plot_mesh(bnd(3));

mesh_bem=bnd;
load electrode; %load the electrodes
figure, ft_plot_mesh(mesh_bem(1),'surfaceonly','yes','vertexcolor','none','facecolor',...
           'r','facealpha',0.5,'edgealpha',0.1)
ft_plot_mesh(mesh_bem(2),'surfaceonly','yes','vertexcolor','none','facecolor',...
           'g','facealpha',0.5,'edgealpha',0.1)
ft_plot_mesh(mesh_bem(3),'surfaceonly','yes','vertexcolor','none','facecolor',...
           'b','facealpha',0.5,'edgealpha',0.1)
hold on, ft_plot_sens(elec);

%% Step 3: create a headmodel and a sourcemodel and Compute lead-field

cfg=[];

if ispc % use openMEEG
    thisDir = pwd;
    openMEEGLocation = fileparts(which('om_check_geom.exe'));
    cd(openMEEGLocation);
    cfg.method = 'openmeeg'; % for EEG as mentioned in help of function
else
    cfg.method = 'dipoli';
end

cfg.tissue = {'brain','skull','scalp'};
headmodel_bem = ft_prepare_headmodel(cfg, bnd);

% Sourcemodel
cfg = [];
cfg.grid.resolution = 10;
cfg.grid.unit = 'mm';
cfg.threshold = 0.1;
cfg.smooth = 5;
cfg.headmodel = headmodel_bem;
cfg.inwardshift = 1; %shifts dipoles away from surfaces
sourcemodel = ft_prepare_sourcemodel(cfg);

figure, ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:),'vertexmarker','o','vertexsize',5);
hold on, ft_plot_mesh(mesh_bem(1),'surfaceonly','yes','vertexcolor','none','facecolor',...
           'skin','facealpha',0.5,'edgealpha',0.1);

% Lead field
cfg = [];
cfg.grid = sourcemodel;
cfg.headmodel= headmodel_bem;
cfg.elec = elec;
cfg.channel = {'EEG'};
cfg.reducerank = 3;
cfg.normalize='no';
leadfield_bem = ft_prepare_leadfield(cfg);

if ispc
    cd(thisDir);
end

%% Step 4 - Work on the EEG data       
oNum=0; refType = 'unipolar'; % 'unipolar', 'avg','bipolar','csd'
[data,layout] = getData(oNum, refType);
data.elec = elec;

%%%%%%%%% Data segmentation into baseline or stimulus epoch %%%%%%%%%%%%%%%

cfg        = [];
cfg.toilim = [-0.5 0];  % Baseline time period -0.5 to 0
data_pre   = ft_redefinetrial(cfg, data);

cfg.toilim = [0.25 0.75]; % Stimulus time period 0.25 to 0.75
data_post  = ft_redefinetrial(cfg, data);

% compute psd (power spectral density) and csd (cross spectral density) for DICS    

cfg                 = [];
cfg.output          = 'powandcsd';
cfg.method          = 'mtmfft';
cfg.taper           = 'dpss';
cfg.tapsmofrq       = 4;  % Smoothing over W = +- 4 Hz. For T=0.5, W=4, TW=2, 3 tapers are needed 
cfg.foi             = 28; % Perform analysis at this frequency
cfg.keeptrials      = 'yes';

freq_pre           = ft_freqanalysis(cfg, data_pre);
freq_post          = ft_freqanalysis(cfg, data_post);

% Source reconstruction using DICS

cfg              = [];
cfg.method       = 'dics';
cfg.grid         = leadfield_bem;
cfg.headmodel    = headmodel_bem;
cfg.frequency    = freq_post.freq;
cfg.dics.projectnoise = 'yes';
cfg.dics.lambda       = '2%';
cfg.dics.keepfilter   = 'yes';
cfg.dics.realfilter   = 'yes';

freq_source_post = ft_sourceanalysis(cfg, freq_post);
cfg.sourcemodel.filter = freq_source_post.avg.filter;
freq_source_pre = ft_sourceanalysis(cfg, freq_pre);

data_sourceNorm = freq_source_post;
data_sourceNorm.avg.pow = freq_source_post.avg.pow ./ freq_source_pre.avg.pow;

% Visualization

% step 1: Interpolate
cfg            = [];
cfg.downsample = 2;
cfg.parameter  = 'pow';
data_source_interp  = ft_sourceinterpolate(cfg, data_sourceNorm, mri_resliced_orig);

% step 2: Plot
cfg              = [];
cfg.method       = 'ortho';
cfg.funparameter = 'pow';
cfg.funcolormap    = 'jet';
ft_sourceplot(cfg, data_source_interp);