% Supplemtal figures (except dPCA which is separate script)

%                  Casper Kerren      [kerren@cbs.mpg.de]


clear
restoredefaultpath
addpath('/Users/kerrenadmin/Desktop/Postdoc/Project_1/Analyses_matlab/general_scripts_matlab/fieldtrip-20230422')
ft_defaults

[~,ftpath]=ft_version;
settings                    = [];
settings.subjects           = char('CF', 'JM', 'SO', 'AH','FC', 'HW', 'AM', 'MH','FS', 'AS', 'CB', 'KK');
settings.SubjectIDs         = char('01_CF', '02_JM', '03_SO', '06_AH','07_FC', '08_HW', '09_AM', '10_MH','11_FS', '12_AS', '13_CB', '14_KK');

settings.scalp_channels     = {'C3' 'C4'  'Cz' 'T3' 'T4' 'T5' 'T6' 'O1' 'O2' 'Oz' 'F3' 'F4' 'Fz' 'Cb1' 'Cb2'};
settings.healthyhemi        = {'R' 'LR' 'R' 'L' 'L' 'R' 'R' 'R' 'R' 'R' 'R' 'LR'};


cwd = pwd;

addpath(genpath([pwd,'/help_functions']))



settings.colour_scheme_1 = brewermap(30,'RdBu');
settings.colour_scheme_1 = settings.colour_scheme_1;

settings.nu_perm = 4096;


%%
% ---------------------------------------------------------------
% ---------------------------------------------------------------
% ---------------- SUPPLEMENTAL FIGURE 1 ------------------------
% ---------------------------------------------------------------
% ---------------------------------------------------------------
%  - plot all channels per participant


load all_channels

load mean_ripple_per_participant
load std_ripple_per_participant
load nu_ripples_per_participant

load TF_ripple_aligned

TF = TF_ripple_aligned;

timeaxis    = TF{1,1}.time;
freqaxis     = TF{1,1}.freq;

toi_idx = -1:1/1000:1;

channelsize = 20;
transp      = .05;
col         = 'r';

tmp                             = load('channels_to_exclude_all_hipp_both_hem.mat');
channels_to_exclude_all_hipp    = tmp.channels_to_exclude_all_hipp;
tmp                             = load('channels_hipp_ripples.mat');
channels_hipp_ripples           = tmp.channels_hipp_ripples;


channelsize = 15;
transp      = .05;
col         = 'r';

this_view   = 'horizontal'; % 'horizontal' 'sagittal' 'coronal'

this_brain  = 'surface_inflated_both';

ld          = load(fullfile(ftpath,'template/anatomy/', [this_brain '.mat']));
brain       = ld.mesh;


epos_all        = [];
subject_colors  = [];  % To store the color for each contact
count_hipp_chan = [];

% Generate a colormap with distinct colors for each subject
num_subjects = size(settings.subjects,1);
cmap = lines(num_subjects);  % Use 'lines' colormap for distinct colors

for isubject = 1:num_subjects
    cfg.channel     = [];
    chan_to_keep    = [];
    idx             = [];
    coordinates     = [];

    cfg.channel                 = {all_channels{isubject}.names};
    chan_to_keep                =  cfg.channel;
    idx                         = ismember(cfg.channel,chan_to_keep);
    count_hipp_chan(isubject,:) = sum(idx);
    coordinates                 = {all_channels{isubject}.coords};
    coordinates                 = coordinates(idx);
    
    epos_all                    = [epos_all; cell2mat(coordinates')];
    
    % Assign the subject's color to the corresponding contacts
    subject_colors = [subject_colors; repmat(cmap(isubject, :), sum(idx), 1)];
end


figure('units','normalized','outerposition',[0 0 1 1]);
views = [180 -90];


for i = 1:12

    % Initialize variables for each iteration
    subject_colors  = [];
    count_hipp_chan = [];
    epos_all        = [];
    cfg.channel     = [];
    chan_to_keep    = [];
    idx             = [];
    coordinates     = [];

    % Select relevant channels
    cfg.channel                 = {all_channels{i}.names};
    chan_to_keep                = intersect(cellstr(setdiff([cfg.channel], settings.scalp_channels)), char(channels_hipp_ripples{i, :}));
%     chan_to_keep                = cfg.channel;
    idx                         = ismember(cfg.channel, chan_to_keep);
    count_hipp_chan(i, :)       = sum(idx);
    coordinates                 = {all_channels{i}.coords};
    coordinates                 = coordinates(idx);
    
    epos_all                    = [epos_all; cell2mat(coordinates')];
    subject_colors = [subject_colors; repmat(cmap(i, :), sum(idx), 1)];

    % Plotting the brain mesh at the top
    if i < 5
        h1 = subplot(6, 4, i);
        pos_h    = get(h1, 'Position');
        pos_h(1) = pos_h(1)-.045;
        pos_h(2) = pos_h(2)-.01;
        pos_h(3) = pos_h(3)+.09;
        pos_h(4) = pos_h(4)+.09;
        set(h1, 'Position', pos_h);


    elseif i >= 5 && i < 9
        h1 = subplot(6, 4, i + 4);
        pos_h    = get(h1, 'Position');
        pos_h(1) = pos_h(1)-.045;
        pos_h(2) = pos_h(2)-.02;
        pos_h(3) = pos_h(3)+.09;
        pos_h(4) = pos_h(4)+.09;
        set(h1, 'Position', pos_h);

    else
        h1 = subplot(6, 4, i + 8);
        pos_h    = get(h1, 'Position');
        pos_h(1) = pos_h(1)-.045;
        pos_h(2) = pos_h(2)-.02;
        pos_h(3) = pos_h(3)+.09;
        pos_h(4) = pos_h(4)+.09;
        set(h1, 'Position', pos_h);

    end

    ft_plot_mesh(brain, 'facealpha', transp, 'facecolor', 'cortex')
    
    for ii = 1:size(epos_all, 1)
        this_MNI = epos_all(ii, :);
        ft_plot_mesh(this_MNI, 'vertexcolor', subject_colors(ii, :), 'vertexsize', channelsize);
    end

    camlight headlight;
    set(gcf, 'color', 'w', 'alpha', 0);
    view(views)

    % Determine subplot position for mean_ripple_per_participant
    if i < 5
        h = subplot(6, 4, i + 4);
        pos = get(h, 'Position'); 
        pos(2) = pos(2)+.01;
        set(h, 'Position', pos);
    elseif i >= 5 && i < 9
        h = subplot(6, 4, i + 8);
        pos = get(h, 'Position'); 
         pos(2) = pos(2)-.03;
        set(h, 'Position', pos);

        pos_h(2) = pos_h(2)-.04;
        set(h1, 'Position', pos_h);

    else
        h = subplot(6, 4, i + 12);
        pos = get(h, 'Position'); 
        pos(2) = pos(2)-.06;
        set(h, 'Position', pos);

        pos_h(2) = pos_h(2)-.08;
        set(h1, 'Position', pos_h);

    end

    % Plot mean_ripple_per_participant with bounded line
    boundedline(toi_idx, mean_ripple_per_participant(i, :), std_ripple_per_participant(i, :), 'k', 'alpha')
    xlabel('Time (sec)')
    ylabel('\muV')
    title(sprintf('patient %d, %d events', i, nu_ripples_per_participant(i, :)), 'interpreter', 'none')
    set(gca, 'tickdir', 'out')
    set(gca, 'FontSize', 16)
    axis tight


    h2 = subplot('Position',[.01, .01, .01, .01]);
    pos2    = get(h2, 'Position');
    pos2(1) = pos(1)+.12;
    pos2(2) = pos(2)+.06;
    pos2(3) = .04;
    pos2(4) = .04;
    set(h2, 'Position', pos2);
    to_plot_ripple = squeeze(mean(TF{i}.all));


    to_plot_ripple = [];
    to_plot_ripple = squeeze(mean(TF{i}.all));

    imagesc(...
        timeaxis,...
        freqaxis,...
        to_plot_ripple);
    colormap(hot)
    caxis([0 .7])
    hline(90,'w')
    axis off


    axis xy
    set(gca,'TickDir','out')
   
end

clearvars -except settings

%% dimensionality (100ms)

% ---------------------------------------------------------------
% ---------------------------------------------------------------
% ---------------- SUPPLEMENTAL FIGURE 3 ------------------------
% ---------------------------------------------------------------
% ---------------------------------------------------------------


figure('units','normalized','outerposition',[0 0 1 1]);

% dimensionality 100ms sliding window
load perf_dimensionality_100ms

subplot(3,3,1)

perf = perf_dimensionality_100ms;


correct_incorrect = {};

for isubject = 1:numel(perf)

    correct_incorrect{1}.label               = {'Channels'};
    correct_incorrect{1}.time                = perf{1,1}.time_train;

    correct_incorrect{1}.individual(isubject,1,:)   = perf{isubject}.correct.accuracy;
    correct_incorrect{1}.dimord              = 'subj_chan_time';

end
correct_incorrect{1}.avg            = squeeze(correct_incorrect{1}.individual);
correct_incorrect{2}                = correct_incorrect{1};

for isubject = 1:numel(perf)

    correct_incorrect{2}.individual(isubject,1,:)           = perf{isubject}.incorrect.accuracy;
  
end

correct_incorrect{2}.avg            = squeeze(correct_incorrect{2}.individual);


xlimits = nearest(correct_incorrect{1, 1}.time, -1):nearest(correct_incorrect{1, 1}.time, 1);
xlimits = correct_incorrect{1, 1}.time(xlimits);

cfg                     = [];
cfg.latency             = [xlimits(1) xlimits(end)];

cfg.channel             = 'all';
cfg.statistic           = 'depsamplesT';
cfg.method              = 'montecarlo'; % 'montecarlo' 'analytic';
cfg.correctm            = 'cluster'; % 'no', cluster;
cfg.alpha               = .05;
cfg.clusteralpha        = .05;
cfg.tail                = 0;
cfg.correcttail         = 'alpha'; % alpha prob no
cfg.neighbours          = [];
cfg.minnbchan           = 0;
cfg.avgovertime         = 'no';    % 'no' 'yes'
cfg.avgoverchan         = 'no';
cfg.computecritval      = 'yes';

cfg.numrandomization    = 1000;%'all';

cfg.clusterstatistic    = 'maxsum'; % 'maxsum', 'maxsize', 'wcm'
cfg.clustertail         = cfg.tail;
cfg.parameter           = 'individual';

nSub = size(correct_incorrect{1, 1}.individual  ,1);
% set up design matrix
design = zeros(2,2*nSub);
for i = 1:nSub
    design(1,i) = i;
end
for i = 1:nSub
    design(1,nSub+i) = i;
end
design(2,1:nSub)        = 1;
design(2,nSub+1:2*nSub) = 2;

cfg.design  = design;
cfg.uvar    = 1;
cfg.ivar    = 2;

% run stats
[Fieldtripstats] = ft_timelockstatistics(cfg, correct_incorrect{:});
length(find(Fieldtripstats.mask))


%% plot significant vals

d = squeeze(correct_incorrect{1}.individual-correct_incorrect{2}.individual);

m = nanmean(d);
s = nanstd(d)./sqrt(size(d,1));

boundedline(perf{1,1}.time_train,m,s,'k','cmap', flipud(settings.colour_scheme_1(23,:)));
plot(perf{1,1}.time_train,m,'k','linewidth',3);
hold on

stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline   = nan(1,numel(correct_incorrect{1}.time));
sigline(stats_time(Fieldtripstats.mask==1)) = 0;

plot(correct_incorrect{1}.time,sigline,'k','linewidth',4);

set(gca,'FontSize',20)
xlabel('ripple time (sec)')
ylabel('dim. difference')
set(gca,'TickDir','out')
xticks([-.9 0 .9])
xticklabels({'-1', '0', '1'})

axis tight
vline(0)
hline(0)

xlim([cfg.latency(1), cfg.latency(end)])

%% dimensionality 200ms sliding window
load perf_dimensionality_200ms

subplot(3,3,2)

perf = perf_dimensionality_200ms;



correct_incorrect = {};

for isubject = 1:numel(perf)

    correct_incorrect{1}.label               = {'Channels'};
    correct_incorrect{1}.time                = perf{1,1}.time_train;

    correct_incorrect{1}.individual(isubject,1,:)   = perf{isubject}.correct.accuracy;
    correct_incorrect{1}.dimord              = 'subj_chan_time';

end
correct_incorrect{1}.avg            = squeeze(correct_incorrect{1}.individual);
correct_incorrect{2}                = correct_incorrect{1};

for isubject = 1:numel(perf)

    correct_incorrect{2}.individual(isubject,1,:)           = perf{isubject}.incorrect.accuracy;
  
end

correct_incorrect{2}.avg            = squeeze(correct_incorrect{2}.individual);


xlimits = nearest(correct_incorrect{1, 1}.time, -1):nearest(correct_incorrect{1, 1}.time, 1);
xlimits = correct_incorrect{1, 1}.time(xlimits);

cfg                     = [];
cfg.latency             = [xlimits(1) xlimits(end)];

cfg.channel             = 'all';
cfg.statistic           = 'depsamplesT';
cfg.method              = 'montecarlo'; % 'montecarlo' 'analytic';
cfg.correctm            = 'cluster'; % 'no', cluster;
cfg.alpha               = .05;
cfg.clusteralpha        = .05;
cfg.tail                = 0;
cfg.correcttail         = 'alpha'; % alpha prob no
cfg.neighbours          = [];
cfg.minnbchan           = 0;
cfg.avgovertime         = 'no';    % 'no' 'yes'
cfg.avgoverchan         = 'no';
cfg.computecritval      = 'yes';

cfg.numrandomization    = 1000;%'all';

cfg.clusterstatistic    = 'maxsum'; % 'maxsum', 'maxsize', 'wcm'
cfg.clustertail         = cfg.tail;
cfg.parameter           = 'individual';

nSub = size(correct_incorrect{1, 1}.individual  ,1);
% set up design matrix
design = zeros(2,2*nSub);
for i = 1:nSub
    design(1,i) = i;
end
for i = 1:nSub
    design(1,nSub+i) = i;
end
design(2,1:nSub)        = 1;
design(2,nSub+1:2*nSub) = 2;

cfg.design  = design;
cfg.uvar    = 1;
cfg.ivar    = 2;

% run stats
[Fieldtripstats] = ft_timelockstatistics(cfg, correct_incorrect{:});
length(find(Fieldtripstats.mask))


%% plot significant vals

d = squeeze(correct_incorrect{1}.individual-correct_incorrect{2}.individual);

m = nanmean(d);
s = nanstd(d)./sqrt(size(d,1));

boundedline(perf{1,1}.time_train,m,s,'k','cmap', flipud(settings.colour_scheme_1(23,:)));
plot(perf{1,1}.time_train,m,'k','linewidth',3);
hold on

stats_time = nearest(correct_incorrect{1}.time,cfg.latency(1)):nearest(correct_incorrect{1}.time,cfg.latency(2));

sigline   = nan(1,numel(correct_incorrect{1}.time));
sigline(stats_time(Fieldtripstats.mask==1)) = 0;

plot(correct_incorrect{1}.time,sigline,'k','linewidth',4);

set(gca,'FontSize',20)
xlabel('ripple time (sec)')
ylabel('dim. difference')
set(gca,'TickDir','out')
xticks([-.9 0 .9])
xticklabels({'-1', '0', '1'})

axis tight
vline(0)
hline(0)

xlim([cfg.latency(1), cfg.latency(end)])




clearvars -except settings




%% reconstructing using components



load perf_PCA_LDA

perf = perf_PCA_LDA;

freqaxis = perf{1, 1}.time_train;
timeaxis = perf{1, 1}.time_test;

    correct                 = cell(1,size(settings.subjects,1));
    incorrect               = cell(1,size(settings.subjects,1));
    correct_dec     = [];
    incorrect_dec   = [];
sigma = .5; % STD of 2D gaussian smoothing
    for isubject = 1:size(settings.subjects,1)

        % correct
        correct{1,isubject}                                 = struct;
        correct{1,isubject}.label                           = {'chan'};
        correct{1,isubject}.dimord                          = 'chan_freq_time';
        correct{1,isubject}.freq                            = perf{1, 1}.time_train;
        correct{1,isubject}.time                            = perf{1, 1}.time_test;
        correct{1,isubject}.powspctrm(1,:,:)                = perf{isubject}.correct.accuracy;


        % incorrect
        incorrect{1,isubject}                               = correct{1,isubject};
        incorrect{1,isubject}.powspctrm(1,:,:)              = perf{isubject}.incorrect.accuracy;

        baseline_all{1,isubject}                            = correct{1,isubject};
        baseline_all{1,isubject}.powspctrm(1,:,:)           = .5*ones(size(perf{isubject}.correct.accuracy));

        correct_dec(isubject,:,:)                           = imgaussfilt(perf{isubject}.correct.accuracy,sigma);
        incorrect_dec(isubject,:,:)                         = imgaussfilt(perf{isubject}.incorrect.accuracy,sigma);
    end


%% FT stats (decoding)

cfg                     = [];
cfg.latency             = [-1 1];
cfg.frequency           = [-.2 3]; % encoding time
cfg.channel             = 'all';
cfg.statistic           = 'depsamplesT';
cfg.method              = 'montecarlo'; % 'montecarlo' 'analytic';
cfg.correctm            = 'cluster'; % 'no', cluster;
cfg.alpha               = .05;
cfg.clusteralpha        = .05;
cfg.tail                = 0;
cfg.correcttail         = 'alpha'; % alpha prob no
cfg.neighbours          = [];
cfg.minnbchan           = 0;
cfg.computecritval      = 'yes';

cfg.numrandomization    = 1000;%1000;%'all';

cfg.clusterstatistic    = 'maxsum'; % 'maxsum', 'maxsize', 'wcm'
cfg.clustertail         = cfg.tail;
cfg.parameter           = 'powspctrm';

nSub = size(settings.subjects,1);
% set up design matrix
design      = zeros(2,2*nSub);
design(1,:) = repmat(1:nSub,1,2);
design(2,:) = [1*ones(1,nSub) 2*ones(1,nSub)];

cfg.design  = design;
cfg.uvar    = 1;
cfg.ivar    = 2;

[Fieldtripstats] = ft_freqstatistics(cfg, correct{:}, incorrect{:});

length(find(Fieldtripstats.mask==1))


stats_time = nearest(correct{1, 1}.time  ,cfg.latency(1)):nearest(correct{1, 1}.time  ,cfg.latency(2));
stats_freq = nearest(freqaxis,cfg.frequency(1)):nearest(freqaxis,cfg.frequency(2));

tvals = squeeze(Fieldtripstats.stat);


%% Correlate with original data 

% get size of the different decoding analyses

xlimits = nearest(correct{1, 1}.time, -1):nearest(correct{1, 1}.time, 1);

% limit decoding from reconstructed data
correct_dec_rec     = correct_dec(:, :, xlimits);
incorrect_dec_rec   = incorrect_dec(:, :, xlimits);
% load original data
load('tvals_coarse_enc_ripple.mat')
load('mask_t_vals_coarse_enc_ripple.mat')

tvals_single_ripple_coarse  = tvals_coarse_enc_ripple;
mask_t_vals_coarse          = squeeze(mask_t_vals_coarse_enc_ripple);

% Load data from decoding analysis
load('correct_dec_all.mat')
load('incorrect_dec_all.mat')

correct_dec_all = correct_dec_all;
incorrect_dec_all = incorrect_dec_all;

[rows_large, cols_large] = size(tvals_single_ripple_coarse);

% Get the size of the smaller matrix
[rows_small, cols_small] = size(tvals);

% Calculate the number of columns to add
cols_to_add = cols_large - cols_small;

% Pad the smaller matrix with zeros to match the size of the larger one
% (only a few data points)
tvals_reconstructed_resized                 = [tvals, zeros(rows_small, cols_to_add)];
correct_dec_rec(:,:,end:end+cols_to_add)    = 0;
incorrect_dec_rec(:,:,end:end+cols_to_add)  = 0;
TOI_to_plot                                 = linspace(-1,1, size(tvals_single_ripple_coarse,2));
FOI_to_plot                                 = correct{1, 1}.freq;


for isubject = 1:size(settings.subjects,1)
    tmp_sub = [];
    tmp_sub = squeeze(correct_dec_rec(isubject, :,:));

    correct_dec_rec(isubject, :,:) = tmp_sub;

    tmp_sub = [];
    tmp_sub = squeeze(correct_dec_all(isubject, :,:));

    correct_dec_all(isubject, :,:) = tmp_sub;

    tmp_sub = [];
    tmp_sub = squeeze(incorrect_dec_rec(isubject, :,:));

    incorrect_dec_rec(isubject, :,:) = tmp_sub;

    tmp_sub = [];
    tmp_sub = squeeze(incorrect_dec_all(isubject, :,:));

    incorrect_dec_all(isubject, :,:) = tmp_sub;
end

correct_dec_mean    = squeeze(mean(correct_dec_all,2));
incorrect_dec_mean  = squeeze(mean(incorrect_dec_all,2));

correct_dec_rec_mean    = squeeze(mean(correct_dec_rec,2));
incorrect_dec_rec_mean  = squeeze(mean(incorrect_dec_rec,2));

[rho_ppp, p_ppp] = corr(correct_dec_mean-incorrect_dec_mean,correct_dec_rec_mean-incorrect_dec_rec_mean,'type', 'Spearman');
subplot(3,3,3)
 imagesc(...
    timeaxis(stats_time),...
    timeaxis(stats_time),...
    rho_ppp(1:end-5,1:end-5))
colormap(flipud(settings.colour_scheme_1))
caxis([0 1])
vline(0)


axis xy
set(gca,'FontSize',18)
ylabel('ripple time (sec) original')
xlabel('ripple time (sec) reconstructed')
set(gca,'TickDir','out')

hcb = colorbar('Location','EastOutside');

ylabel(hcb, 'Correlation', 'FontSize', 20); 

correct_dec_mean    = mean(mean(correct_dec_all,3),2);
incorrect_dec_mean  = mean(mean(incorrect_dec_all,3),2);

correct_dec_rec_mean    = mean(mean(correct_dec_rec,3),2);
incorrect_dec_rec_mean  = mean(mean(incorrect_dec_rec,3),2);

[rho_ppp, p_ppp] = corr(correct_dec_mean-incorrect_dec_mean,correct_dec_rec_mean-incorrect_dec_rec_mean,'type', 'Spearman');


