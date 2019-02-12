function [bad_chans_new, bad_chans_old, eeg_raw_proc, glob_sig_all, glob_sig_clean, mvar_ref_Z_last, mamp_ref_Z_last, line_rel_amp, chan_names]=eeg_noise_metrics_CZ(subj,sess,t_offset,t_length,T,overlap,varargin)
%EEG_NOISE_METRICS
%Written by: Julio I. Chapeton
%
%Computes the variance and range of each channels in small time windows. The averages of these quantities across all windows, z-scored across all channels, is used to flag...
%potentially bad channels as judged by variance and range being abnormally large or small (+/- 2 sigma is the default). Also allows the user to run EEG_STFT in order to look...
%for line noise or for abnormal spectra.
%
%Usage:
%EEG_NOISE_METRICS(subj,sess,t_offset,t_length,T,overlap,varargin)
%EEG_NOISE_METRICS(...,'PropertyName',PropertyValue,...)
%
%Inputs:
%    Required:
%        subj       = string  (patient ID, e.g.'NIH###', 'SZ###', or 'TJ###')
%        sess       = string/integer (date info or session number, if loading local files this must be the date info in the following format:'YYMMDD_HHMM', for ...
%                             server files this can be the date info or a session number, to use a session number you must provide a task name using the 'task' option)
%                     NOTE: this assumes that local files were extracted to '/.../subj/eeg.noreref/[dirname]/subj_YYMMDD_HHMM.###'
%        t_offset   = scalar  (how much time to skip at the beginning of a file in seconds, if input is empty the default is 60 seconds)
%        t_length   = scalar  (length of data to load in seconds, if there is not enough extra data before or after, the buffer for detrending and filtering will...
%                              be made as large as possible)
%                     NOTE: For a full clip use t_offset = 0 and t_length = Inf, there will be no buffering in this case.
%        T          = scalar  (size of the time window (in seconds). Input 'Inf' to use the whole time series with no overlap)
%        overlap    = scalar  (size of overlap for adjacent windows, as a fraction of T)
%
%    Optional:
%        (name-value pairs, not case sensitve)
%
%        specific to this code:
%        visual     = logical (provide a value of 1 to plot the time series, variances, amplitudes, and to run EEG_STFT in visual mode)
%        disp       = logical (provide a value of 1 to display the figures, or 0 to hide them (but still create/save them)...
%                              This is only relevant if the 'visual' and 'saving' options are turned on.)
%        outputdir  = string  (full path to output directory, the default is '.../.../data.noreref/session/') NOT WORKING YET
%        z_thresh   = scalar  (Threshold in units of sigma above/below which channels are flagged as potentially bad. The default is 2 sigma)
%        saving     = logical (provide a value of 1 to save the bad channels list to the bad_channs file. You will be given the option to abort before you save)
%
%        from load_eegfile:
%        loc_detrend= logical/vector (provide a value of 1 to locally detrend the data, the default parameters [1 0.5] are such that frequencies below ~1-2Hz are filtered out...
%                              Alternatively, can provide a two element vector of the form [window_size, step_size], smaller windows = higher cutoff frequencies)
%        rmline     = logical (provide a value of 1 to use a regression based method to remove 60Hz and 120Hz noise)
%        rem_sat    = scalar  (remove segments where the time series is constant for at least n samles, provide n as the value, 10 works well)
%        filter     = string  (filter data with an FIR filter, available options are:'lowpass','delta','theta','alpha','beta','lowgamma','highgamma', 'gaussian', 'ictal' or...
%                              'lowpass_sleep'. There is also the 'notch' option which implements an IIR bandstop filter at 60Hz)
%      hardwareType = string  (Options are 'depth' and 'subdural'. Default is 'subdural'. Prepending '~' is set negation.)
%        location   = string (Option for loading local or server files, options are 'local', 'server', or a full path to the 'eeg' directory...
%                             e.g. '/Volumes/shares/FRNU/dataWorking/eeg/'. If using the 'server' option you must have access to FRNU)
%        task       = string  (If files are being loaded from the server, the 'sess' input was given as a session number, then this option must be used to define which...
%                              task to load, in all other instances it will be ignored)
%                   Allowed property values:
%                   'attentionTask'
%                   'auditoryLexicalDecision'
%                   'auditoryVowels'
%                   'palRam'
%                   'palRamStim'
%                   'paRemap'
%                   'pa3'
%                   'restingState'
%                   'stimMapping'
%                   'seizure'
%                   'baseline'
%
%Outputs:
%    Required:
%        bad_chans_new  = cell  (list of channels flagged as potentially bad during the current run of EEG_NOISE_METRICS)
%        bad_chans_old  = cell  (list of channels previously rejected from previous runs of EEG_NOISE_METRICS)
%        eeg_raw_proc   = matrix  (raw data for all channels after local detrending and line noise removal)
%        glob_sig_all   = vector  (Global average signal before cleaning, includes all channels)
%        glob_sig_clean = vector  (Global average signal after cleaning, includes only 'good' channels)
%        mvar_ref_Z_last= matrix  (mean variance for each channel, z-scored across channels using sigma and mu from the final run (i.e. all chans z-scored to mu and sig of good chans))
%        mamp_ref_Z_last= matrix  (mean rms for each channel, z-scored across channels using sigma and mu from the final run (i.e. all chans z-scored to mu and sig of good chans))
%        line_rel_amp   = matrix  (Difference (in dB) between the expected signal amplitude at 60Hz and the actual amplitude of the signal at 60Hz. A value of 3 corresponds to 2x as much power as expected)
%        chan_names     = cell    (List of channel names for the channels corresponding to channelMap. For bipolar it will be of the form 'elec1-elec2')
%
%Functions Called:
%    load_eegfile
%    ts_info
%    tightfig
%    eeg_stft
%    plot_ts_gui
%
% NOTES:
% 1. If this becomes the top level channel rejection (as opposed to individual user level) then it should be ran with the same options every time. In this case...
%    many of the optional inputs should be hard coded instead.
% 2. May want to give an option such that z_thresh can be defined differently for + and -.
% 3. The recommended usage is EEG_NOISE_METRICS('NIH###',timestamp,60*5,Inf,Inf,0,'visual',1,'rmline',1,'loc_detrend',1)...
%    or EEG_NOISE_METRICS('NIH###',sess_num,60*5,Inf,Inf,0,'visual',1,'rmline',1,'task',task_flag,'loc_detrend',1) for task sessions. The default lowpass seems sharp enough to...
%    remove the huge 180Hz line noise, but maybe it should be implemented here just to make sure that it does.
%% parse inputs
inp_pars = inputParser;

% specific to this code
def_visual = 0;
def_disp = 0;
def_z_thresh = 2;
def_outputdir = fullfile('/Users/chapetonji/Desktop/temp_data',num2str(sess)); %%% MACHINE DEPENDENT
def_saving = 0;

% carry over options from 'load_eegfile'
def_loc_detrend = 0;
def_rmline = 0;
def_rem_sat = 0;
def_filter=[];
% def_ref = 'global_avg';
% expected_ref={'raw','global_avg','dev_avg','bipolar'};
def_location='server';
def_task=[];
def_hwType = 'subdural';

%
addParameter(inp_pars,'visual',def_visual,@(x) (x==1 || x==0));
addParameter(inp_pars,'disp',def_disp,@(x) (x==1 || x==0));
addParameter(inp_pars,'z_thresh',def_z_thresh,@(x) isscalar(x));
addParameter(inp_pars,'outputdir',def_outputdir,@(x) (ischar(x) || isempty(x)));
addParameter(inp_pars,'saving',def_saving,@(x) (x==1 || x==0));

% carry over options from 'load_eegfile', validation of these inputs is handled there
addParameter(inp_pars,'loc_detrend',def_loc_detrend);
addParameter(inp_pars,'rmline',def_rmline);
addParameter(inp_pars,'rem_sat',def_rem_sat);
addParameter(inp_pars,'filter',def_filter);
% carry over options from 'GetFilesandChannels', validation of these inputs is handled there
addParameter(inp_pars,'location',def_location);
addParameter(inp_pars,'task',def_task);
addParameter(inp_pars, 'hardwareType',def_hwType);

addParameter(inp_pars, 'biowulf', '0', @ischar);
addParameter(inp_pars, 'biowulf_eeg_raw_proc', '', @ismatrix);
addParameter(inp_pars, 'biowulf_Fs', 1000, @isnumeric);
addParameter(inp_pars, 'biowulf_chan_inds', 0);
addParameter(inp_pars, 'biowulf_chan_names', {});

disp(varargin);

parse(inp_pars,varargin{:})
rmline=inp_pars.Results.rmline;
loc_detrend=inp_pars.Results.loc_detrend;
z_thresh=inp_pars.Results.z_thresh;
rem_sat=inp_pars.Results.rem_sat;
visual=inp_pars.Results.visual;
disp_ind=inp_pars.Results.disp;
filt_flag=inp_pars.Results.filter;
%ref_flag=inp_pars.Results.ref;
loc_flag=inp_pars.Results.location;
task_flag=inp_pars.Results.task;
save_flag=inp_pars.Results.saving;
outputdir=inp_pars.Results.outputdir;
hwType_flag = inp_pars.Results.hardwareType;

biowulf = logical(eval(inp_pars.Results.biowulf));
 
if ~save_flag
    disp_ind=1;
end

if ~exist(outputdir, 'dir') && save_flag
    mkdir(outputdir);
end

% output
bad_chans_new = {};
bad_chans_old = {};
eeg_raw_proc = [];
glob_sig_all = [];
glob_sig_clean = [];
mvar_ref_Z_last = [];
mamp_ref_Z_last = [];
line_rel_amp = [];

%% load data
if ~biowulf

    [eeg_raw_proc, Fs, chan_inds, chan_names, sess_info]=load_eegfile(subj,sess,0,Inf,'loc_detrend',loc_detrend,'rmline',rmline,'rem_sat',rem_sat,'visual',0,...
        'filter',filt_flag,'ref','raw','location',loc_flag,'task',task_flag,'hardwareType',hwType_flag);
    if all(isnan(eeg_raw_proc(:)))
        fprintf('Error: very bad - all NaN'); keyboard;
    end

else
    % when using biowulf make sure to include varargin options:
    %
    % biowulf
    % biowulf_eeg_raw_proc
    % biowulf_Fs
    % biowulf_chan_inds
    % biowulf_chan_names
    % outputdir
    
    eeg_raw_proc = inp_pars.Results.biowulf_eeg_raw_proc;
    Fs = inp_pars.Results.biowulf_Fs;
    chan_inds = inp_pars.Results.biowulf_chan_inds;
    chan_names = inp_pars.Results.biowulf_chan_names;
    
    t_offset = 60;
    t_length = Inf;
    T = Inf;
    overlap = 0;
    
    visual = 1;
    save_flag = 1;

end

glob_sig_all=nanmean(eeg_raw_proc,2);
eeg_avg_reref=bsxfun(@minus,eeg_raw_proc,glob_sig_all); % global re-ref

% bad_chans_dir=fullfile('/Users/chapetonji/Desktop/SZ_JC/',subj) %*** DEFINED BY THE USER, SHOULD MATCH THE BADCHANSDIR FROM GetFilesandChannels.m. Maybe useful to force it by using 'assignin' or something like that?
bad_chans_dir=outputdir;
nchans=length(chan_inds);
% chan_inds=(1:numel(chan_inds))';
N_skip=round(t_offset*Fs)+1;

if isinf(t_length)
    N_eeg=size(eeg_avg_reref,1)-N_skip+1;
else
    N_eeg=round(t_length*Fs);
end
%% calculate moments & plot sorted mean variance and range for reref
Nsamples=round(T*Fs);
ovrlp_samples=round(overlap*Nsamples);
nwins=floor((N_eeg-Nsamples)/(Nsamples-ovrlp_samples))+1;
%[m1_raw, m2_raw, m3_raw, m4_raw, amp_raw,~,times_raw] = ts_info(eeg_raw,Fs,T,overlap); %%% may add an option to do ref or re-ref

[~, m2_ref, ~, ~, amp_ref, times_ref] = ts_info(eeg_avg_reref(N_skip:N_skip+N_eeg-1,:),Fs,T,overlap);
m_vars=mean(m2_ref,1);
[vars_ord,indvar]=sort(m_vars);
m_amps=mean(amp_ref,1);
[amps_ord,indamp]=sort(m_amps);

muvar_run1=mean(mean(m2_ref,1));
stdvar_run1=std(mean(m2_ref,1));
muamp_run1=mean(mean(amp_ref,1));
stdamp_run1=std(mean(amp_ref,1));

mvar_ref_Z_run1=(((mean(m2_ref,1))-mean(mean(m2_ref,1)))./std(mean(m2_ref,1)))';
mamp_ref_Z_run1=(((mean(amp_ref,1))-mean(mean(amp_ref,1)))./std(mean(amp_ref,1)))';
% rej_chans=sort(find(abs(mvar_ref_Z_run1)>=z_thresh | abs(mamp_ref_Z_run1)>=z_thresh));

ind_rej=(abs(mvar_ref_Z_run1)>=z_thresh | abs(mamp_ref_Z_run1)>=z_thresh);
%ind_rej=(abs(mvar_ref_Z_run1)>=z_thresh);

rej_chans=find(ind_rej);
rej_chans_accum=rej_chans;
chans_keep=find(~ind_rej);
if isempty(rej_chans)
    glob_sig_clean=glob_sig_all;
    bad_chans_new=[];
    bad_chans_old=[];
    if exist(fullfile(bad_chans_dir,'bad_chans.mat'))==2
        bad_chans_old=load(fullfile(bad_chans_dir,'bad_chans'),'bad_chans');
        bad_chans_old=bad_chans_old.bad_chans;
    end
    disp('No new bad channels have been found, nothing will be saved')
else
    printmat_JC([sum(ind_rej),max(mvar_ref_Z_run1),min(mvar_ref_Z_run1),max(mamp_ref_Z_run1),min(mamp_ref_Z_run1)],[],sprintf('Iter#%d',0),strjoin({'N_rejected Max_var Min_var Max_amp Min_amp'}))
end

%%% these are useful figures, however, they take a long time to produce
% subsample at Fs/4 to speed up by downsize figure
ds_k = 6;
downsample = @(x) x(1:ds_k:size(x,1), :);

if visual
    filename1 = fullfile(outputdir,'full_ts_ref');
    ds_eeg_avg_reref = downsample(eeg_avg_reref);
    if ~exist(filename1,'file') || ~save_flag
        full_ts_fig_ref=plot_ts_gui(ds_eeg_avg_reref,'scroll',0,'Fs',Fs/ds_k,'npoints',size(ds_eeg_avg_reref,1),'y_labels',chan_names,'disp',disp_ind);
    end
    
    filename2 = fullfile(outputdir,'full_ts_raw');
    ds_eeg_raw_proc = downsample(eeg_raw_proc);
    if ~exist(filename2, 'file') || ~save_flag
        full_ts_fig_raw=plot_ts_gui(ds_eeg_raw_proc,'scroll',0,'Fs',Fs/ds_k,'npoints',size(ds_eeg_raw_proc,1),'y_labels',chan_names,'disp',disp_ind);
    end
    
    if save_flag
        print(full_ts_fig_ref, filename1, '-dpng', '-r0');
        delete(full_ts_fig_ref)

        print(full_ts_fig_raw, filename2, '-dpng', '-r0');
        delete(full_ts_fig_raw)
    end
end

tic;
loop_cnt = 0;
while sum(ind_rej) || loop_cnt == 0
    loop_cnt = loop_cnt + 1;
    glob_sig_clean=nanmean(eeg_raw_proc(:,chans_keep),2);
    temp_eeg=bsxfun(@minus,eeg_raw_proc(:,chans_keep),glob_sig_clean); % global re-ref
    
    [~, m2_ref, ~, ~, amp_ref, times_ref] = ts_info(temp_eeg(N_skip:N_skip+N_eeg-1,:),Fs,T,overlap);
    muvar_last=mean(mean(m2_ref,1));
    stdvar_last=std(mean(m2_ref,1));
    muamp_last=mean(mean(amp_ref,1));
    stdamp_last=std(mean(amp_ref,1));
    
    mvar_ref_Z=(((mean(m2_ref,1))-mean(mean(m2_ref,1)))./std(mean(m2_ref,1)))';
    mamp_ref_Z=(((mean(amp_ref,1))-mean(mean(amp_ref,1)))./std(mean(amp_ref,1)))';
    
    ind_rej=(abs(mvar_ref_Z)>=z_thresh | abs(mamp_ref_Z)>=z_thresh);
    %     ind_rej=(abs(mvar_ref_Z)>=z_thresh);
    
    chans_keep(ind_rej)=[];
    rej_chans_accum=setdiff(chan_inds,chans_keep);
    printmat_JC([sum(ind_rej),max(mvar_ref_Z),min(mvar_ref_Z),max(mamp_ref_Z),min(mamp_ref_Z)],[],sprintf('Iter#%d',loop_cnt),strjoin({'N_rejected Max_var Min_var Max_amp Min_amp'}))
end

clear temp_eeg
mvar_ref_Z_last=(m_vars-muvar_last)./stdvar_last;
mamp_ref_Z_last=(m_amps-muamp_last)./stdamp_last;

if visual
    %%% vars/amps
    if disp_ind
        h1=figure;
    else
        h1=figure('visible','off');
    end
    subplot(2,2,1:2)
    stem(mvar_ref_Z_last(indvar))
    hold on
    plot([1:nchans],z_thresh*ones(nchans,1),'r')
    plot([1:nchans],-z_thresh*ones(nchans,1),'r')
    set(gca,'xtick',[1:numel(chan_names)],'xticklabel',chan_names(indvar),'Ticklength', [0 0],'XTickLabelRotation',90)
    %set(gca,'xtick',[1:numel(chan_names)],'xticklabel',[1:numel(chan_names)],'Ticklength', [0 0])
    set(gca,'LooseInset',get(gca,'TightInset'))
    xlim([0 nchans+1])
    title('Average reference variances')
    
    subplot(2,2,3:4)
    stem(mamp_ref_Z_last(indamp))
    hold on
    plot([1:nchans],z_thresh*ones(nchans,1),'r')
    plot([1:nchans],-z_thresh*ones(nchans,1),'r')
    set(gca,'xtick',[1:numel(chan_names)],'xticklabel',chan_names(indamp),'Ticklength', [0 0],'XTickLabelRotation',90)
    %set(gca,'xtick',[1:numel(chan_names)],'xticklabel',[1:numel(chan_names)],'Ticklength', [0 0])
    set(gca,'LooseInset',get(gca,'TightInset'))
    xlim([0 nchans+1])
    title('Average reference amplitudes')
    
    % set(h1,'PaperOrientation','landscape');
    % set(h1,'PaperUnits','normalized');
    % set(h1,'PaperPosition', [0 0 1 1]);
    set(h1,'units','normalized','outerposition',[0 0 1 1])
    set(gca,'LooseInset',get(gca,'TightInset'))
    
    tightfig;
    set(h1,'units','normalized','outerposition',[0 0 1 1])
    
    %%% global signals
    if disp_ind
        h2=figure;
    else
        h2=figure('visible','off');
    end
    subplot(3,2,1:2)
    plot([0:size(eeg_avg_reref,1)-1]./Fs,glob_sig_all)
    set(gca,'LooseInset',get(gca,'TightInset'))
    xlim([[0 size(eeg_avg_reref,1)-1]./Fs])
    title('Raw global signal')
    
    subplot(3,2,3:4)
    plot([0:size(eeg_avg_reref,1)-1]./Fs,glob_sig_clean)
    set(gca,'LooseInset',get(gca,'TightInset'))
    xlim([[0 size(eeg_avg_reref,1)-1]./Fs])
    title('Cleaned global signal')
    
    subplot(3,2,5:6)
    plot([0:size(eeg_avg_reref,1)-1]./Fs,abs(glob_sig_all-glob_sig_clean))
    set(gca,'LooseInset',get(gca,'TightInset'))
    xlim([[0 size(eeg_avg_reref,1)-1]./Fs])
    title('Absolute difference')
    % set(h1,'PaperOrientation','landscape');
    % set(h1,'PaperUnits','normalized');
    % set(h1,'PaperPosition', [0 0 1 1]);
    set(h2,'units','normalized','outerposition',[0 0 1 1])
    set(gca,'LooseInset',get(gca,'TightInset'))
    tightfig;
    set(h2,'units','normalized','outerposition',[0 0 1 1])
end
if save_flag
    print(h1, fullfile(outputdir,'vars_amps'), '-dpng', '-r0');
    clear h1
    
    print(h2, fullfile(outputdir,'global_sigs'), '-dpng', '-r0');
    clear h2
    
    close all
    %return
end

%%%ind_rej_ref=union(find(abs(mvar_ref_Z)>=z_thresh),find(abs(mamp_ref_Z)>=z_thresh));

%% get channels with abnormal amplitude or variance
%ind_rej=sort(find(abs(mamp_ref_Z)>=z_thresh));
% rej_chans=sort(find(abs(mvar_ref_Z)>=z_thresh | abs(mamp_ref_Z)>=z_thresh));
if isempty(rej_chans_accum)
    fprintf('no bad channels were automatically found based on amplitude or variance criteria\n');
else
    
    if any(mvar_ref_Z_last(rej_chans_accum)>=0 | mamp_ref_Z_last(rej_chans_accum)>=0)
        fprintf('the following channels have been flagged as possibly bad due to large average amplitude or variance after re-referencing\n')
        tmp=rej_chans_accum(sort(unique(find(mvar_ref_Z_last(rej_chans_accum)>=0 | mamp_ref_Z_last(rej_chans_accum)>=0))));
        printmat_JC(uint16(chan_inds(tmp)),[],strjoin(chan_names(tmp)),['index'])
        clear tmp
    end
    
    if any(mvar_ref_Z_last(rej_chans_accum)<=0 | mamp_ref_Z_last(rej_chans_accum)<=0)
        fprintf('the following channels have been flagged as possibly bad due to small average amplitude or variance after re-referencing\n')
        tmp=rej_chans_accum(sort(unique(find(mvar_ref_Z_last(rej_chans_accum)<=0 | mamp_ref_Z_last(rej_chans_accum)<=0))));
        printmat_JC(uint16(chan_inds(tmp)),[],strjoin(chan_names(tmp)),['index'])
        clear tmp
    end
end

% if ~isempty(rej_chans_accum)
%     sprintf('the following channels have been flagged as possibly bad due to small or large average amplitude or variance before re-referencing\n')
%     printmat_JC([channelMap(rej_chans_accum),rej_chans_accum],[],strjoin(chan_names(rej_chans_accum)),['number index'])
% end

%ind_rej=sort(ind_rej);
% ind_rej=union(ind_rej_raw,ind_rej_ref,'sorted');
%% check spectra for unusually large line noise and high-frequency artifacts

%inp=input('Input 1 to run EEG_STFT in visual mode, or press return to skip\nYou can access the indices of the flagged channels from EEG_STFT by using: evalin(''caller'',''ind_rej'')\n');

% T_fft = 4;
% overlap_fft = 0.25;
T_fft = 5;
overlap_fft = 0.2;
freq_range = [0, min(Fs/2, 200)];

if 1% inp == 1 %strcmpi(inp,'yes')
    
    if ~biowulf
        [AS,~,f]=eeg_stft(subj,sess,t_offset,t_length,T_fft,overlap_fft,'loc_detrend',loc_detrend,'rmline',rmline,'rem_sat',rem_sat,'visual',visual,'disp',disp_ind,'freq_range', freq_range, ...
            'filter',filt_flag,'ref','global_avg','location',loc_flag,'task',task_flag,'loading',eeg_avg_reref(N_skip:N_skip+N_eeg-1,:),'saving',save_flag,'outputdir',outputdir);
    else
        
        loading_struct = struct;
        loading_struct.eegfile_re_ref = eeg_avg_reref(N_skip:N_skip+N_eeg-1,:);
        loading_struct.Fs = Fs;
        loading_struct.chan_inds = chan_inds;
        loading_struct.chan_names = chan_names;

        [AS,~,f]=eeg_stft_CZ(subj,sess,t_offset,t_length,T_fft,overlap_fft,'visual',visual,'disp',disp_ind,'freq_range', freq_range, 'loading', loading_struct, 'saving', save_flag, 'outputdir', outputdir, 'parallel', 0);
    end
end

[~,ind55]=min(abs(f-55));
f55=f(ind55);
[~,ind60]=min(abs(f-60));
f60=f(ind60);
[~,ind65]=min(abs(f-65));
f65=f(ind65);

mu_55_65=mean(AS(:,[ind55,ind65]),2);
amp_60=AS(:,ind60);

line_thresh=3;
line_rel_amp=amp_60-mu_55_65;
line_rej=find(line_rel_amp > line_thresh);

if ~isempty(line_rej)
    fprintf('These are all the channels that have been flagged as possibly bad due to 60Hz line noise\n')
    printmat_JC(uint16(chan_inds(line_rej)),[],strjoin(chan_names(line_rej)),['index'])
end

%rej_chans_accum=union(rej_chans_accum,line_rej); % *MST 06/17 do not include line_noise rejected channels

if ~isempty(rej_chans_accum)
    fprintf('These are all the channels that have been flagged as possibly bad due to abnormal amplitude or variance:\n')
    printmat_JC(uint16(chan_inds(rej_chans_accum)),[],strjoin(chan_names(rej_chans_accum)),['index'])
end
toc
%% plot time series in small pieces
if visual
    if numel(setdiff(chan_inds,rej_chans_accum))>=numel(rej_chans_accum)
        chan_ind=[rej_chans_accum;randsample(setdiff(chan_inds,rej_chans_accum),numel(rej_chans_accum))];
    else
        chan_ind=[rej_chans_accum;randsample(setdiff(chan_inds,rej_chans_accum),numel(setdiff(chan_inds,rej_chans_accum)))];
    end
    
    if save_flag && ~isempty(chan_ind)
        % ds=downsample by a factor of ds_k (4) to increase plotting speed
        % maybe should move this ds stuff after eeg_ts_ assignment
        ds_N_eeg = floor(N_eeg / ds_k);
        ds_Fs = Fs / ds_k;
        
        npoints=min(ds_N_eeg/10,60*ds_Fs); % the minimum of 1/10 the clip length, or 1 minute
        npoints=floor(npoints);
        ind=floor(ds_N_eeg/2);
        
        eeg_ts_ref=ds_eeg_avg_reref(ind:ind+npoints,chan_ind)';
        eeg_ts_raw=ds_eeg_raw_proc(ind:ind+npoints,chan_ind)';
        T_secs = [ind-1:(ind-1)+(size(eeg_ts_ref,2)-1)] ./ ds_Fs;
        
        % RAW
        ub = max(eeg_ts_raw,[],2); lb = min(eeg_ts_raw,[],2);
        Sig_raw=nanstd(eeg_ts_raw,0,2);
        %spacing=cumsum([0; abs(ub(1:end-1))+abs(lb(2:end))]+max(std(eeg_ts,0,2)));
        spacing=cumsum([0; abs(ub(1:end-1))+abs(lb(2:end))]);%+[0;Sig(1:end-1)]);
        spacing = repmat(spacing,1,numel(T_secs));
        eegspaced_raw=eeg_ts_raw+spacing;
        
        if disp_ind
            h1=figure;
        else
            h1=figure('visible','off');
        end
        h1_ax=gca;
        %         set(h,'PaperOrientation','landscape'); % these commands are usefull if pdfs are going to be saved or if plots are going to be printed
        %         set(h,'PaperUnits','normalized');
        %         set(h,'PaperPosition', [0 0 1 1]);
        set(h1_ax,'LooseInset',get(h1_ax,'TightInset'))
        h11=plot([0,1],ones(length(chan_ind),2));
        xlabel('time in seconds')%,ylabel('electrode');
        set(h1,'units','normalized','outerposition',[1 0 1 0.475]);
        set(h1_ax,'LooseInset',get(h1_ax,'TightInset'));
        set(h1_ax,'YGrid','on')
        set(h1_ax,'XGrid','on')
        title(h1_ax,'eeg time series: no re-ref')
        h11(length(rej_chans_accum)+1).LineWidth=1.5;
        h11(length(rej_chans_accum)+1).Color='k';
        
        set(h11,'xdata',T_secs,{'ydata'},num2cell(eegspaced_raw,2))
        temp=nanmean(eegspaced_raw,2);
        % set ticks/labels
        %ylim(h1_ax,[min(eegspaced_raw(:))-Sig_raw(1) max(eegspaced_raw(:))+Sig_raw(end)])
        %xlim([min(T_secs) max(T_secs)])
        set(h1_ax,'yaxislocation','right','yticklabel',chan_names(chan_ind),'ylim',[min(eegspaced_raw(:))-Sig_raw(1) max(eegspaced_raw(:))+Sig_raw(end)])
        set(h1_ax,'xtick',[0:5:max(T_secs)],'xticklabel',[0:5:max(T_secs)],'Ticklength', [0 0],'xlim',[min(T_secs) max(T_secs)])
        %grid on; box on
        drawnow;
        
        % REREF
        ub = max(eeg_ts_ref,[],2); lb = min(eeg_ts_ref,[],2);
        Sig_ref=nanstd(eeg_ts_ref,0,2);
        %spacing=cumsum([0; abs(ub(1:end-1))+abs(lb(2:end))]+max(std(eeg_ts,0,2)));
        spacing=cumsum([0; abs(ub(1:end-1))+abs(lb(2:end))]);%+[0;Sig(1:end-1)]);
        spacing = repmat(spacing,1,numel(T_secs));
        eegspaced_ref=eeg_ts_ref+spacing;
        
        if disp_ind
            h2=figure;
        else
            h2=figure('visible','off');
        end
        h2_ax=gca;
        %         set(h1,'PaperOrientation','landscape'); % these commands are usefull if pdfs are going to be saved or if plots are going to be printed
        %         set(h1,'PaperUnits','normalized');
        %         set(h1,'PaperPosition', [0 0 1 1]);
        set(h2_ax,'LooseInset',get(h2_ax,'TightInset'))
        h21=plot([0,1],ones(length(chan_ind),2));
        %xlabel('time in seconds')%,ylabel('electrode');
        set(h2,'units','normalized','outerposition',[1 1 1 0.5]);
        set(h2_ax,'LooseInset',get(h2_ax,'TightInset'));
        set(h2_ax,'YGrid','on')
        set(h2_ax,'XGrid','on')
        title(h2_ax,'eeg time series: common avg ref')
        h21(length(rej_chans_accum)+1).LineWidth=1.5;
        h21(length(rej_chans_accum)+1).Color='k';
        
        set(h21,'xdata',T_secs,{'ydata'},num2cell(eegspaced_ref,2))
        temp=nanmean(eegspaced_ref,2);
        % set ticks/labels
        %ylim([min(eegspaced_ref(:))-Sig_ref(1) max(eegspaced_ref(:))+Sig_ref(end)])
        %xlim([min(T_secs) max(T_secs)])
        set(h2_ax,'yaxislocation','right','ytick',temp,'yticklabel',chan_names(chan_ind),'ylim',[min(eegspaced_ref(:))-Sig_ref(1) max(eegspaced_ref(:))+Sig_ref(end)])
        set(h2_ax,'xtick',[0:5:max(T_secs)],'xticklabel',[0:5:max(T_secs)],'Ticklength', [0 0],'xlim',[min(T_secs) max(T_secs)])
        %grid on; box on
        drawnow;
        
        %%% saving
        %         hgsave(h1,fullfile(outputdir,'raw_compare'))
        % saveas(h1,fullfile(outputdir,'raw_compare'),'jpg')
        print(h1, fullfile(outputdir,'raw_compare'), '-dpng', '-r0');
        
        clear h1
        %         hgsave(h2,fullfile(outputdir,'reref_compare'))
        % saveas(h2,fullfile(outputdir,'reref_compare'),'jpg')
        print(h2, fullfile(outputdir,'reref_compare'), '-dpng', '-r0');
        
        clear h2
        close all
%         
%     else % interactive mode
%         chan_ind=input('\npress return to view all channels, enter 0 to skip visualization, enter a vector of channel INDICES to view specific channels, enter ''auto'' to\nview only channels automatically detected as possibly bad, or ''compare'' to compare the bad channels to randomly selected channels\n');
%         if isempty(chan_ind)
%             chan_ind=channelInd;
%         elseif strcmp(chan_ind,'auto')
%             chan_ind=rej_chans_accum;
%         elseif strcmp(chan_ind,'compare')
%             chan_ind=[rej_chans_accum;randsample(setdiff(channelInd,rej_chans_accum),numel(rej_chans_accum))];
%         end
%         npoints=min(N_eeg/10,60*Fs); % the minimum of 1/10 the clip length, or 1 minute
%         
%         if chan_ind~=0
%             
%             h1=figure;
%             h1_ax=gca;
%             %         set(h,'PaperOrientation','landscape'); % these commands are usefull if pdfs are going to be saved or if plots are going to be printed
%             %         set(h,'PaperUnits','normalized');
%             %         set(h,'PaperPosition', [0 0 1 1]);
%             set(h1_ax,'LooseInset',get(h1_ax,'TightInset'))
%             h11=plot([0,1],ones(length(chan_ind),2));
%             xlabel('time in seconds')%,ylabel('electrode');
%             set(h1,'units','normalized','outerposition',[1 0 1 0.475]);
%             set(h1_ax,'LooseInset',get(h1_ax,'TightInset'));
%             set(h1_ax,'YGrid','on')
%             set(h1_ax,'XGrid','on')
%             title(h1_ax,'eeg time series: no re-ref')
%             h11(length(rej_chans_accum)+1).LineWidth=1.5;
%             h11(length(rej_chans_accum)+1).Color='k';
%             drawnow;
%             
%             h2=figure;
%             h2_ax=gca;
%             %         set(h1,'PaperOrientation','landscape'); % these commands are usefull if pdfs are going to be saved or if plots are going to be printed
%             %         set(h1,'PaperUnits','normalized');
%             %         set(h1,'PaperPosition', [0 0 1 1]);
%             set(h2_ax,'LooseInset',get(h2_ax,'TightInset'))
%             h21=plot([0,1],ones(length(chan_ind),2));
%             %xlabel('time in seconds')%,ylabel('electrode');
%             set(h2,'units','normalized','outerposition',[1 1 1 0.5]);
%             set(h2_ax,'LooseInset',get(h2_ax,'TightInset'));
%             set(h2_ax,'YGrid','on')
%             set(h2_ax,'XGrid','on')
%             title(h2_ax,'eeg time series: common avg ref')
%             h21(length(rej_chans_accum)+1).LineWidth=1.5;
%             h21(length(rej_chans_accum)+1).Color='k';
%             drawnow;
%             
%             btn1 = uicontrol('Style', 'pushbutton', 'String', 'next','Callback', @next_callback,'units','normalized','position',[0.95 0.95 0.025 0.025],'Fontsize',16,'Interruptible','on','UserData',0);
%             btn2 = uicontrol('Style', 'pushbutton', 'String', 'previous','Callback',@previous_callback ,'units','normalized','position',[0.925 0.95 0.025 0.025],'Fontsize',16,'Interruptible','on','UserData',0);
%             btn3 = uicontrol('Style', 'pushbutton', 'String', 'exit','Callback',@close_exit ,'units','normalized','position',[0.9 0.95 0.025 0.025],'Fontsize',16,'Interruptible','on','UserData',0);
%             drawnow;
%             
%             ind=1;
%             while ind < N_eeg
%                 if (ind+npoints)>=N_eeg
%                     npoints=N_eeg-ind;
%                 end
%                 eeg_ts_ref=eeg_avg_reref(ind:ind+npoints,chan_ind)';
%                 eeg_ts_raw=eeg_raw_proc(ind:ind+npoints,chan_ind)';
%                 T_secs = [ind-1:(ind-1)+(size(eeg_ts_ref,2)-1)]./Fs;
%                 
%                 % RAW
%                 ub = max(eeg_ts_raw,[],2); lb = min(eeg_ts_raw,[],2);
%                 Sig_raw=nanstd(eeg_ts_raw,0,2);
%                 %spacing=cumsum([0; abs(ub(1:end-1))+abs(lb(2:end))]+max(std(eeg_ts,0,2)));
%                 spacing=cumsum([0; abs(ub(1:end-1))+abs(lb(2:end))]);%+[0;Sig(1:end-1)]);
%                 spacing = repmat(spacing,1,numel(T_secs));
%                 eegspaced_raw=eeg_ts_raw+spacing;
%                 
%                 set(h11,'xdata',T_secs,{'ydata'},num2cell(eegspaced_raw,2))
%                 temp=nanmean(eegspaced_raw,2);
%                 % set ticks/labels
%                 set(h1_ax,'yaxislocation','right','ytick',temp,'yticklabel',chan_names(chan_ind),'ylim',[min(eegspaced_raw(:))-Sig_raw(1) max(eegspaced_raw(:))+Sig_raw(end)])
%                 set(h1_ax,'xtick',[0:5:max(T_secs)],'xticklabel',[0:5:max(T_secs)],'Ticklength', [0 0],'xlim',[min(T_secs) max(T_secs)])
%                 %grid on; box on
%                 drawnow;
%                 
%                 % REREF
%                 ub = max(eeg_ts_ref,[],2); lb = min(eeg_ts_ref,[],2);
%                 Sig_ref=nanstd(eeg_ts_ref,0,2);
%                 %spacing=cumsum([0; abs(ub(1:end-1))+abs(lb(2:end))]+max(std(eeg_ts,0,2)));
%                 spacing=cumsum([0; abs(ub(1:end-1))+abs(lb(2:end))]);%+[0;Sig(1:end-1)]);
%                 spacing = repmat(spacing,1,numel(T_secs));
%                 eegspaced_ref=eeg_ts_ref+spacing;
%                 
%                 set(h21,'xdata',T_secs,{'ydata'},num2cell(eegspaced_ref,2))
%                 temp=nanmean(eegspaced_ref,2);
%                 % set ticks/labels
%                 set(h2_ax,'yaxislocation','right','ytick',temp,'yticklabel',chan_names(chan_ind),'ylim',[min(eegspaced_ref(:))-Sig_ref(1) max(eegspaced_ref(:))+Sig_ref(end)])
%                 set(h2_ax,'xtick',[0:5:max(T_secs)],'xticklabel',[0:5:max(T_secs)],'Ticklength', [0 0],'xlim',[min(T_secs) max(T_secs)])
%                 %grid on; box on
%                 drawnow;
%                 
%                 % enable moving through plots using buttons
%                 uiwait
%                 b1=get(btn1,'UserData');
%                 b2=get(btn2,'UserData');
%                 bb=b1+b2; % one should always be zero
%                 
%                 set(btn1,'UserData',0); % reset buttons
%                 set(btn2,'UserData',0);
%                 
%                 if get(btn3,'UserData')
%                     %close all
%                     break
%                 elseif (ind+npoints*bb) >= 1 && (ind+npoints*bb)<N_eeg
%                     ind=ind+npoints*bb;
%                 elseif (ind+npoints*bb) < 1
%                     %ind=1;
%                     beep
%                     warning('You are at the beginning of the clip')
%                 elseif (ind+npoints*bb)>=N_eeg
%                     %ind=ind-npoints*bb;
%                     beep
%                     warning('You have reached the end of the clip')
%                 end % **
%             end %while
%         end
    end
end

%% turn channel indices to channel numbers and save the bad channels to the subject's directory if saving is turned on
%% Mike edit: turn channel indices to channel *names*

% name lookup
rej_chans_accum = chan_names(rej_chans_accum); % index --> name

if save_flag
    %save_button = questdlg('Would you like to add to the bad channel list?','see saving options?','yes','no','no'); % just being paranoid here
    if 1%strcmp(save_button,'yes')
        if ~isempty(rej_chans_accum)
            % *MST 06/17 skip user input and always use 'auto'
            bad_chans_new = rej_chans_accum;
        else
            % *MST 06/17 skip user input and always use 'auto'
            bad_chans_new=[]; %input('No bad channels were found automatically, press return to keep all or input a vector of channel INDICES to reject manually\n');
        end
        
        bad_chans_old=[];
        if exist(fullfile(bad_chans_dir,'bad_chans.mat'),'file')==2
            bad_chans_old=load(fullfile(bad_chans_dir,'bad_chans'),'bad_chans');
            bad_chans_old=bad_chans_old.bad_chans;
            if isnumeric(bad_chans_old)
                bad_chans_old = [];
            end
        end
        bad_chans=union(bad_chans_old,bad_chans_new);
        save(fullfile(bad_chans_dir,'bad_chans'),'bad_chans');
    end
else
    bad_chans_new=rej_chans_accum;
    bad_chans_old=[];
    if exist(fullfile(bad_chans_dir,'bad_chans.mat'))==2
        bad_chans_old=load(fullfile(bad_chans_dir,'bad_chans'),'bad_chans');
        bad_chans_old=bad_chans_old.bad_chans;
    end
end
close all
%clear global ind_rej % using global variables is never recommended ***
end
%% callback functions for ui objects
function close_exit(src,event)
drawnow
uiresume
set(src,'UserData',1)
end

function next_callback(src,event)
drawnow
uiresume
set(src,'UserData',1)
end

function previous_callback(src,event)
drawnow
uiresume
set(src,'UserData',-1)
end