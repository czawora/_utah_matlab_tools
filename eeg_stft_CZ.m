function [AS, tf, f, centers_pow, outputdir] = eeg_stft_CZ(subj,sess,t_offset,t_length,T,overlap,varargin)
%function [AS, tf, f, centers_pow, outputdir] = eeg_stft(subj,sess,t_offset,t_length,T,overlap,varargin)
%EEG_STFT
%Written by: Julio I. Chapeton
%
%Calculates spectra and time-frequency representation using tapered-overlapped short-time Fourier transforms.
%
%Usage:
%EEG_STFT(subj,sess,t_offset,t_length,T,overlap,varargin)
%EEG_STFT(...,'PropertyName',PropertyValue,...)
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
%        T          = scalar  (size of the time window (in seconds))
%        overlap    = scalar  (size of overlap for adjacent windows, as a fraction of T)
%
%    Optional:
%        (name-value pairs, not case sensitve)
%
%        specific to this code:
%        loading    = matrix or struct (If load_eegfile was previously ran, you can use eeg_ts as an input to skip re-loading data. The input parameters have to match those of the...
%                              original call to load_eegfile. This is usefull when calling this function inside a script that already calls load_eegfile)
%                       OR
%                     if you are running on biowulf pass in struct with the fields :
%                       * eegfile_re_ref = sample x channel matrix of data
%                       * Fs = (scalar) sampling rate of data, in hertz
%                       * chan_inds = (column scalar vector) channel index numbers for data
%                       * chan_names = (cell vector) channel names for data
%
%        visual     = logical (provide a value of 1 to enable visualization , or 0 to disable (default))
%        disp       = logical (provide a value of 1 to display the figures, or 0 to hide them (but still create/save them)...
%                              This is only relevant if the 'visual' and 'saving' options are turned on.)
%        zscore     = logical (provide a value of 1 to rescale by zscoring spectrograms for each frequency across time, will also z-score the time series in the visualization)
%        baseline   = string  (use this option to provide a baseline session (e.g. baseline = 'YYMMDD_HHMM') in order to give results in decibels with the baseline as the reference)
%        smoothing  = logical (provide a value of 1 to smooth frequencies with a gaussian (makes the results look more like wavelets))
%        saving     = logical (provide a value of 1 to save the spectrograms and spectra to pdf)
%        outputdir  = string  (full path to output directory where figures will be saved, the default is '')
%        freq_range = vector  (Range of frequencies to use for visualization, should be of the form [min_freq max_freq]. This is usefull for band-passed data. The default is [0 Fs/4])
%
%        from load_eegfile:
%        loc_detrend= logical/vector (provide a value of 1 to locally detrend the data, the default parameters are such that frequencies below ~1Hz are filtered out...
%                              Alternatively, can provide a two element vector of the form [window_size, step_size], smaller windows = higher cutoff frequencies)
%        rmline     = scalar  (implements a regression based method to remove 60Hz and 120Hz noise, this can be done iteratively, and the input should be the number or iterations)
%        rem_sat    = scalar  (remove segments where the time series is constant for at least n samles, provide n as the value, 10 works well)
%        filter     = string  (filter data with an FIR filter, available options are:'lowpass','delta','theta','alpha','beta','lowgamma','highgamma', 'gaussian', 'ictal'...
%                              'highpass',or 'lowpass_sleep'. There is also the 'notch' option which implements an IIR bandstop filter at 60Hz)
%        ref        = string  (referencing scheme, options are 'raw','global_avg', 'dev_avg', or 'bipolar', the default is 'global_avg'. The 'dev_avg' and 'bipolar' options...
%                              are pre-defined referencing schemes that use a device weighted global average, or the difference between neighboring electrodes respectively)
%                     NOTE: 1. The 'raw' option can be good for checking bad channels and referencing, but generally should not be used for analysis.
%                           2. The 'dev_avg' has no bad channel removal and has not been implemented for local files, however, this option might get removed altogether.
%        location   = string (Option for loading local or server files, options are 'local', 'server', or a full path to the 'eeg' directory...
%                             e.g. '/Volumes/shares/FRNU/dataWorking/eeg/'. If using the 'server' option you must have access to FRNU)
%      hardwareType = string (Options are 'depth' and 'subdural'. Default is 'subdural'. Prepending '~' is set negation.)
%        task       = string  (If files are being loaded from the server, and a 'sess' was given as a session number, then this option must be used to define which...
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
%        channels   = cell    (Cell of channel names to load, the default is to load all channels. For bipolar referencing any pair will be loaded where at least one of the pairs is in...
%                              the cell. NOTE: this is different from 'bad_chans', e.g. this option still uses all channels for the global signal, but only returns the requested ones)
%        bad_chans  = string (Ignores channels that have been previously flagged as bad. Options are 'session', to ignore channels that were flagged as bad for THIS...
%                             session, 'any' to ignore channels that were flagged as bad during ANY session, or a cell of channel names to reject.
%                             The default is [] (empty) which returns all channels)
%                       NOTE: The 'session' and 'any' options require that you've ran EEG_NOISE_METRICS_LOCAL
% bad_chans_local_dir = string  (full path to where the bad_chans structure were saved, GetFilesandChannels will look here for bad channels)
%        parallel   =   scalar   1 (def) or 0 indicating whether fft calculations will be done using parpool
%
%Outputs:
%    Required:
%        AS         = matrix  (amplitude spectrum up to Fnyquist (channel x frequency). Units are in dB with the reference voltage being 1 microvolt.)
%        tf         = array   (time frequency representation (channel x frequency x time). Units are in dB with the reference power being 1 microvolt^2 )
%        f          = vector  (vector of frequencies in Hz where the fft was calculated)
%        centers_pow= vector  (time (in seconds) at the center of each STFT window)
%        outputdir  = string  (directory where figures were saved)
%                     NOTE: If called with no outputs visual mode will be automatically enabled
%                           If the baseline option is enabled the units of AS and tf will be such that the references are defined by the baseline
%
%Functions Called:
%    load_eegfile
%    stft_mt
%    tightfig
%
%NOTES:
% 1. T=4 and overlap=0.5 works well for a general overview of long clips
% 2. Have not tested baseline functionality
% 3. Auto saving is only enabled for the summary figure
% 4. The figure is currently being saved as a png for space, but may add a flag for saving as .fig

%% look for user defined default params for the toolbox
usr_params_flag=false;
if exist('DT_GetUserParams.m','file')
    usr_params_flag=true;
    user_params=DT_GetUserParams;
end

%% parse inputs
inp_pars = inputParser;

% specific to this code
def_loading = [];
def_visual = 0;
def_disp = 1;
def_zscore = 0;
def_baseline = [];
def_smoothing = 0;
def_saving = 0;
def_outputdir = [];
def_freq_range = [];
def_parallel = 1;

% carry over options from 'load_eegfile'
def_loc_detrend = 0;
def_rmline = 0;
def_rem_sat = 0;
def_filter=[];
def_chans=[];
% carry over options from 'GetFilesandChannels'
def_ref = 'global_avg';
def_location='server';
def_task=[];
def_hwType = 'subdural';
def_bad_chans_local_dir = fullfile(filesep,'Users','chapetonji','Desktop','JC_Extract','bad_chans_JC',subj);
def_bad_chans=[];

%
addParameter(inp_pars,'loading',def_loading,@(x) (ismatrix(x) || isstruct(x)));
addParameter(inp_pars,'visual',def_visual,@(x) (x==1 || x==0));
addParameter(inp_pars,'disp',def_disp,@(x) (x==1 || x==0));
addParameter(inp_pars,'zscore',def_zscore,@(x) (x==1 || x==0));
addParameter(inp_pars,'baseline',def_baseline,@(x) (ischar(x) || isempty(x)));
addParameter(inp_pars,'smoothing',def_smoothing,@(x) (x==1 || x==0));
addParameter(inp_pars,'saving',def_saving,@(x) (x==1 || x==0));
addParameter(inp_pars,'outputdir',def_outputdir,@(x) (ischar(x) || isempty(x)));
addParameter(inp_pars,'freq_range',def_freq_range,@(x) (isvector(x) && numel(x)==2) || isempty(x));
addParameter(inp_pars,'parallel',def_parallel,@(x) (x==1 || x==0));

% carry over options from 'load_eegfile', validation of these inputs is handled there
addParameter(inp_pars,'loc_detrend',def_loc_detrend);
addParameter(inp_pars,'rmline',def_rmline);
addParameter(inp_pars,'rem_sat',def_rem_sat);
addParameter(inp_pars,'filter',def_filter);
addParameter(inp_pars,'channels',def_chans);
% carry over options from 'GetFilesandChannels', validation of these inputs is handled there
addParameter(inp_pars,'ref',def_ref);
addParameter(inp_pars,'location',def_location);
addParameter(inp_pars,'task',def_task);
addParameter(inp_pars, 'hardwareType',def_hwType);
addParameter(inp_pars,'bad_chans_local_dir',def_bad_chans_local_dir);
addParameter(inp_pars,'bad_chans',def_bad_chans);

parse(inp_pars,varargin{:})
if usr_params_flag && ~isempty(user_params.bad_chans_local_dir)% to keep this consistent it is best to define this through 'DT_GetUserParams', see DT_README for more info on this
    bad_chans_local_dir=fullfile(user_params.bad_chans_local_dir,subj);
else
    bad_chans_local_dir=fullfile(inp_pars.Results.bad_chans_local_dir,subj);
end

%
zscore_ind=inp_pars.Results.zscore;
saving=inp_pars.Results.saving;
BL_sess=inp_pars.Results.baseline;
smoothing=inp_pars.Results.smoothing;
outputdir=inp_pars.Results.outputdir;
%load_flag=~isempty(inp_pars.Results.loading);
freq_range=inp_pars.Results.freq_range;
visual=inp_pars.Results.visual;
disp_ind=inp_pars.Results.disp;
if ~nargout
    visual=1;
end
if ~saving
    disp_ind=1;
end

%
loc_detrend=inp_pars.Results.loc_detrend;
rmline=inp_pars.Results.rmline;
rem_sat=inp_pars.Results.rem_sat;
filt_flag=inp_pars.Results.filter;
filt_ind=~isempty(filt_flag);
ref_flag=inp_pars.Results.ref;
loc_flag=inp_pars.Results.location;
task_flag=inp_pars.Results.task;
channels=inp_pars.Results.channels;
hwType_flag = inp_pars.Results.hardwareType;
bad_chans_flag=inp_pars.Results.bad_chans;
parallel = inp_pars.Results.parallel;

if ~isempty(BL_sess) && zscore_ind
    flag=input('\nWarning: you have chosen to baseline subtract and then zscore, press return to continue, or input anything to quit\n');
    if ~isempty(flag)
        return
    end
end

% if saving
%     flag=input('Warning: saving of figures is turned on, press return to continue, or input anything to quit');
%     if ~isempty(flag)
%         return
%     end
% end

%% find baseline parameters

if ~isempty(BL_sess)
    [eeg_re_ref_BL, Fs, chan_inds]=load_eegfile(subj,BL_sess,0,Inf,'loc_detrend',loc_detrend,'rmline',rmline,'rem_sat',rem_sat,'visual',0,...
        'filter',filt_flag,'ref',ref_flag,'location',loc_flag,'task',task_flag,'channels',channels,'hardwareType',hwType_flag,'bad_chans',bad_chans_flag,'bad_chans_local_dir',bad_chans_local_dir);% right now baseline is calculated from entire recording, change '0' and 'Inf' to shorten
    
    alldata_BL = aux_func1(eeg_re_ref_BL,Fs,T,overlap,parallel);
    
    tf_BL=alldata_BL.time_frequency;
    mu_tf_BL=squeeze(mean(10*log10(tf_BL),3));
    %std_tf_BL=squeeze(std(10*log10(temp1),0,3));
    mu_amp_BL=squeeze(mean(20*log10(alldata_BL.amplitudes),1));
    %std_amp_BL=squeeze(std(20*log10(alldata_BL.amplitudes),0,2));
    freqs_FFT_BL=unique(alldata_BL.frequencies);
    
    clear alldata_BL tf_BL channelMap
end

%% load and process data
if ~isempty(inp_pars.Results.loading)
    
    if isstruct(inp_pars.Results.loading)
        
        fprintf('"loading" parameter passed in as a struct, this usage is reserved for processing on biowulf and expects certain struct fields: eegfile_re_ref, Fs, chan_inds, chan_names\n');
        
        eegfile_re_ref = inp_pars.Results.loading.eegfile_re_ref;
        Fs = inp_pars.Results.loading.Fs;
        chan_inds = inp_pars.Results.loading.chan_inds;
        chan_names = inp_pars.Results.loading.chan_names;

    else
        
        [~, Fs, chan_inds, chan_names, sess_info, eegdir]=load_eegfile(subj,sess,t_offset,t_length,'loc_detrend',loc_detrend,'rmline',rmline,'rem_sat',rem_sat,'visual',0,...
            'filter',filt_flag,'ref',ref_flag,'location',loc_flag,'task',task_flag,'channels',channels,'info_only',1,'hardwareType',hwType_flag,'bad_chans',bad_chans_flag,'bad_chans_local_dir',bad_chans_local_dir);
        eegfile_re_ref=inp_pars.Results.loading;
        % The data is in memory twice because inp_parse.Results is read only and cannot be cleared. There are probably better ways to do this
        assert( (size(eegfile_re_ref,1)==sess_info.eeg_seg_length) && (size(eegfile_re_ref,2)==length(chan_inds)),'The input parameters do not match those of the original call to ''load_eegfile''')
        %%% in order to be REALLY sure that the loaded data matches the passed parameters it may be usefull to include eegdir as part of the 'loading' value so that it can be directly checked here

    end
    
else
    [eegfile_re_ref, Fs, chan_inds, chan_names]=load_eegfile(subj,sess,t_offset,t_length,'loc_detrend',loc_detrend,'rmline',rmline,'rem_sat',rem_sat,'visual',0,...
        'filter',filt_flag,'ref',ref_flag,'location',loc_flag,'task',task_flag,'channels',channels,'hardwareType',hwType_flag,'bad_chans',bad_chans_flag,'bad_chans_local_dir',bad_chans_local_dir);
end
alldata = aux_func1(eegfile_re_ref,Fs,T,overlap,parallel);

mean_sig=sum(eegfile_re_ref,2);
T_ts=(0:size(eegfile_re_ref,1)-1)./Fs;

% data
f=unique(alldata.frequencies);
times_pow=alldata.fft_time_windows;
centers_pow = mean(times_pow,3);
centers_pow=unique(centers_pow);
clear times_pow

% dB, may be better to output the raw values since the baseline option is not enabled yet
AS=20*log10(alldata.amplitudes);
%AS=(alldata.amplitudes);
tf=10*log10(alldata.time_frequency);

if ~isempty(BL_sess)
    assert(isequal(freqs_FFT_BL,f),'the calculated frequencies for the baseline and the seizure clip do not match')
    AS=bsxfun(@minus,AS,mu_amp_BL);
    tf=bsxfun(@minus,tf,mu_tf_BL);
    clear freqs_FFT
end

str_scale=' (dB)';
if zscore_ind
    str_scale=' (zscore)';
    %     AS=bsxfun(@minus,AS,nanmean(AS)); % only z-scoring spectrograms
    %     AS=bsxfun(@rdivide,AS,nanstd(AS));
    tf=bsxfun(@minus,tf,nanmean(tf,3));
    tf=bsxfun(@rdivide,tf,nanstd(tf,[],3));
    eegfile_re_ref=zscoren(eegfile_re_ref);
end

% filter tf with gaussian
if smoothing
    filt_win=15;
    sig=3;
    gfilt=normpdf(-filt_win:filt_win,0,sig);
    
    if ~biowulf
    
%         parfor j=1:size(tf,1)
%             tf(j,:,:)=filtfilt(gfilt,1,squeeze(tf(j,:,:)));
%         end
    
    else
        
        for j=1:size(tf,1)
            tf(j,:,:)=filtfilt(gfilt,1,squeeze(tf(j,:,:)));
        end
        
    end
end
% clear alldata
% alldata = struct('frequencies',f,'amplitudes',AS,'time_frequency',tf,'fft_time_windows',times_pow); % may be better to output the variables separetly
%% spectra and time-frequency summary
if visual
    %define default font size
    f_size=14;
    %     set(0,'defaultAxesFontSize', f_size)
    
    mean_spec_raw=nanmean(AS,1);
    Mp = mean(AS,1);
    SDp = std(AS,0,1);
    sig=nanstd(tf(:));
    
    % spectrum for all channels (avg over time), usefull for looking at line noise
    %set(0,'DefaultFigureColormap',autumn); % ********* sets default colormap
    
    %     max_freq=input('Input a maximum frequency, or press enter to use the default of Fs/2 \n');
    %     if isempty(max_freq) || max_freq==0
    %         max_freq=floor(Fs/2);
    %     end
    
    if isempty(freq_range)
        min_freq=0;
        max_freq=Fs/4;
        [~,fl_ind]=min(abs(f-min_freq)); %******** minimum frequency
        [~,fh_ind]=min(abs(f-max_freq)); %******** maximum frequency
    else
        [~,fl_ind]=min(abs(f-freq_range(1))); %******** minimum frequency
        [~,fh_ind]=min(abs(f-freq_range(2))); %******** maximum frequency
    end
    
    if disp_ind
        h=figure;
    else
        h=figure('visible','off');
    end
    %     h_summary=gcf;
    set(h,'units','normalized','outerposition',[0 0 1 1])
    s1=subplot(2,4,4,'align');
    if filt_ind
        image([f(fl_ind) f(fh_ind)],[1 length(chan_inds)],AS(:,fl_ind:fh_ind))
    else
        imagesc([f(fl_ind) f(fh_ind)],[1 length(chan_inds)],AS(:,fl_ind:fh_ind))
    end
    set(gca, 'YDir', 'normal');
    set(gca,'ytick',[0.5:1:length(chan_inds)],'yticklabel',chan_names,'fontsize',10)
    set(gca,'xtick',[f(fl_ind):50:f(fh_ind)]);
    xlim([f(fl_ind) f(fh_ind)])
    %xlabel('Frequency (Hz)');
    grid on
    title(gca,'spectrum for all channels','fontsize',f_size);
    %     colorbar('location','eastoutside')
    hc1=colorbar('location','southoutside','fontsize',f_size);
    ylabel(hc1,'Amplitude (dB)','fontsize',f_size)
    set(hc1,'ydir','reverse')
    
    % mean spectrum (avg over channels), usefull for determining if the baseline periods correspond to awake or asleep, awake recordings have a peak at ~11Hz
    col=jet;
    %figure(50)
    s2=subplot(2,4,8,'align');
    %hold on
    plot(f,mean_spec_raw,'color',col(randi(64),:))
    xlim([f(fl_ind) f(fh_ind)])
    xlabel('Frequency (Hz)')
    ylabel('Amplitude (dB)')
    %     ylabel('Power \muV^2')
    %set(gca,'LooseInset',get(gca,'TightInset'))
    title('mean raw spectrum across all channels')
    
    s1Pos = get(s1,'position');
    s2Pos = get(s2,'position');
    s1Pos(3:4) = [s2Pos(3:4)];
    set(s1,'position',s1Pos);
    
    % mean spectrogram (avg over channels), usefull for global artifacts
    %figure;
    s3=subplot(2,4,[1:3],'align');
    if filt_ind
        image([min(centers_pow(:)) max(centers_pow(:))],[f(fl_ind) f(fh_ind)],squeeze(nanmean(tf(:,fl_ind:fh_ind,:))))
    else
        imagesc([min(centers_pow(:)) max(centers_pow(:))],[f(fl_ind) f(fh_ind)],squeeze(nanmean(tf(:,fl_ind:fh_ind,:))))
    end
    %caxis([-4*sig 4*sig])
    set(gca, 'YDir', 'normal');
    set(gca,'ytick',[f(fl_ind):10:f(fh_ind)])
    set(gca,'xtick',[centers_pow(1:int64(length(centers_pow)/10):end)+(1/Fs/2)],'xticklabel',[centers_pow(1:int64(length(centers_pow)/10):end)+(1/Fs/2)])
    ylim([f(fl_ind) f(fh_ind)]);ylabel('Frequency (Hz)')
    xlim([(min(centers_pow(:))-T*(1-overlap))+(1/Fs/2) (max(centers_pow(:))+T*(1-overlap))+(1/Fs/2)]);xlabel('time (s)')
    title(['mean time-frequency representation across all channels'])
    hc2=colorbar('location','eastoutside');
    ylabel(hc2,['Amplitude',str_scale])
    %xlim auto
    
    % mean signal across all electrodes
    %figure;
    s4=subplot(2,4,[5:7],'align');
    plot(T_ts,mean_sig,'k')
    set(gca,'xtick',[centers_pow(1:int64(length(centers_pow)/10):end)+(1/Fs/2)],'xticklabel',[centers_pow(1:int64(length(centers_pow)/10):end)+(1/Fs/2)])
    %set(gca,'xtick',[T_ts(1:int64(length(T_ts)/10):end)],'xticklabel',[T_ts(1:int64(length(T_ts)/10):end)])
    xlim([(min(centers_pow(:))-T*(1-overlap))+(1/Fs/2) (max(centers_pow(:))+T*(1-overlap))+(1/Fs/2)])
    ylabel('\muV')
    xlabel('Time(s)')
    title(['Mean signal across all electrodes'])
    
    s3Pos = get(s3,'position');
    s4Pos = get(s4,'position');
    s3Pos(3:4) = [s4Pos(3:4)];
    set(s3,'position',s3Pos);
    
    drawnow;
    tightfig;
    %     set(h,'PaperOrientation','landscape');
    %     set(h,'PaperUnits','normalized');
    %     set(h,'PaperPosition', [0 0 1 1]);
    %     set(gca,'LooseInset',get(gca,'TightInset'))
    set(h,'units','normalized','outerposition',[0 0 1 1])
    drawnow;
    %     pause(0.25)
    
    batch_flag=0;
    if saving
        if ~exist(outputdir, 'dir')
            mkdir(outputdir);
        end
        print(h, fullfile(outputdir,'mean_spectral'), '-dpng', '-r0');
        % if this is being called as part of processAndReref/eeg_noise_metrics the function should return here
        % and not enter interactive mode
        parent_funcs=dbstack;
        parent_funcs={parent_funcs.name};
        batch_funcs={'eeg_noise_metrics','processAndReref'};
        if any(ismember(batch_funcs,parent_funcs)) || ~disp_ind || isstruct(inp_pars.Results.loading)
            batch_flag=1;
            close(h)
            clear h
            return % right now only allowing for auto-save of the summary figure
        end
    end
    %%%%%%%%%%%% per channel plots %%%%%%%%%%%
    
    if ~batch_flag 
        % set pushbuttons for viewing per channel plots or exiting
        btn1 = uicontrol(h,'Style', 'pushbutton', 'String', 'continue','Callback', @btn1_callback,'units','normalized','position',[0.355 0.5 0.05 0.05],'Fontsize',16,'Interruptible','on');
        btn2 = uicontrol(h,'Style', 'pushbutton', 'String', 'exit','Callback',@next_callback ,'units','normalized','position',[0.355 0.45 0.05 0.05],'Fontsize',16,'Interruptible','on','UserData',0);
        set(h,'units','normalized','outerposition',[0 0 0.9 0.9]) % this is a hack so that the buttons appear in the right place
        pause(0.25)
        set(h,'units','normalized','outerposition',[0 0 1 1])
        drawnow;
        uiwait(h)
        try
            if get(btn2,'UserData')
                %close all
                delete(btn1)
                delete(btn2)
                return
            end
            delete(btn1)
            delete(btn2)
        catch
            warning('Figure was closed without the exit button')
            %return
        end
    end
    
    fprintf('\nThese are all the channel names and indices')
    %     channelInd=(1:numel(channelMap))';
    printmat_JC(uint16(chan_inds),[],strjoin(chan_names'),['index'])
    
    %close all
    plot_button = questdlg('Click ''all'' to view all channels, ''input'' to input a vector of channel indices, or ''exit''','which channels would you like to view?','all','input','exit','all');
    if strcmp(plot_button,'exit')
        return
    elseif strcmp(plot_button,'input')
        commandwindow
        plot_channels = input('\n Provide a vector with the INDEX of the channels you want to visualize\n');
    else
        plot_channels=1:length(chan_inds);
    end
    
    k=plot_channels(1);
    
    % spectra and spectrograms
    h1=figure;
    set(h1,'units','normalized','outerposition',[1 1 1 1]) %% comment for ref_testing
    
    temp=squeeze(tf(k,fl_ind:fh_ind,:));
    h1sub=subplot(3,1,1);
    if filt_ind
        h11=image([min(centers_pow(:)) max(centers_pow(:))],[f(fl_ind) f(fh_ind)],temp);
    else
        h11=imagesc([min(centers_pow(:)) max(centers_pow(:))],[f(fl_ind) f(fh_ind)],temp);
    end
    set(gca, 'YDir', 'normal');
    set(gca,'ytick',[f(fl_ind):10:f(fh_ind)])
    set(gca,'xtick',[centers_pow(1:int64(length(centers_pow)/10):end)+(1/Fs/2)],'xticklabel',[centers_pow(1:int64(length(centers_pow)/10):end)+(1/Fs/2)])
    ylim([f(fl_ind) f(fh_ind)]);ylabel('Frequency (Hz)') % *********
    xlim([(min(centers_pow(:))-T*(1-overlap))+(1/Fs/2) (max(centers_pow(:))+T*(1-overlap))+(1/Fs/2)]);xlabel('time (s)')
    title('time-frequency representation')
    c_temp=colorbar('peer',h1sub,'location','eastoutside');
    %colorbar('location','eastoutside')
    
    % vlotage trace
    %     figure
    %     h2=gcf;
    h2sub=subplot(3,1,2);
    h_ts=plot((0:size(eegfile_re_ref,1)-1)./Fs,eegfile_re_ref(:,k));
    set(gca,'xtick',[centers_pow(1:int64(length(centers_pow)/10):end)+(1/Fs/2)],'xticklabel',[centers_pow(1:int64(length(centers_pow)/10):end)+(1/Fs/2)])
    xlabel('Time(s)')
    ylabel('\muV')
    xlim([(min(centers_pow(:))-T*(1-overlap))+(1/Fs/2) (max(centers_pow(:))+T*(1-overlap))+(1/Fs/2)])
    %set(gcf,'units','normalized','outerposition',[0 0 1 1])
    h2_ax=gca;
    title(h2_ax,'Voltage trace')
    
    h3sub=subplot(3,1,3);
    hold on
    h13p=plot(f,Mp+2*SDp,'r');
    h13m=plot(f,Mp-2*SDp,'b');
    h13=plot(f,(AS(k,:)),'k');%index
    set(h3sub,'xtick',[0:10:max(f)])
    xlim([f(fl_ind) f(fh_ind)]);xlabel('frequency (Hz)');ylabel('amplitude (dB)')
    title('mean spectrum +/- 2 sigma across all channels, and spectrum for this channel')
    
    tightfig; % there are a properties that have to be set before and after tightfig for it to work properly
    
    s1Pos = get(h1sub,'position');
    s2Pos = get(h2sub,'position');
    s3Pos = get(h3sub,'position');
    s2Pos(3:4) = [s1Pos(3:4)];
    s3Pos(3:4) = [s1Pos(3:4)];
    set(h2sub,'position',s2Pos);
    set(h3sub,'position',s3Pos);
    
    set(h1,'units','normalized','outerposition',[1 1 1 1]) %% comment for ref_testing
    str_sup=sprintf('Channel %s, Index = %i out of %i',chan_names{k},1,numel(plot_channels)); % channel name and index
    %suptitle(str_sup) % suptitle can conflict with a lot of plotting stuff, consider removing it and using set(gcf,'name',...)
    set(h1,'name',str_sup)
    legend([h13p h13 h13m],'mean + std','spectrum','mean - std')
    delete(c_temp) % fix for tightfig
    colorbar('peer',h1sub,'location','eastoutside')
    
    linkaxes([h1sub,h2sub,s3,s4],'x')
    linkaxes([h1sub,s3],'y')
    hz = zoom;
    set(hz,'Enable','on');
    %     hp = pan;
    %     set(hp,'Motion','horizontal','Enable','on');
    
    btn1 = uicontrol(h1,'Style', 'pushbutton', 'String', 'next','Callback', @next_callback,'units','normalized','position',[0.5 0.65 0.05 0.05],'Fontsize',16,'Interruptible','on','UserData',0);
    btn2 = uicontrol(h1,'Style', 'pushbutton', 'String', 'previous','Callback',@previous_callback ,'units','normalized','position',[0.45 0.65 0.05 0.05],'Fontsize',16,'Interruptible','on','UserData',0);
    btn3 = uicontrol(h1,'Style', 'pushbutton', 'String', 'exit','Callback',@next_callback ,'units','normalized','position',[0.35 0.65 0.05 0.05],'Fontsize',16,'Interruptible','on','UserData',0);
    if saving
        btn4 = uicontrol(h1,'Style', 'pushbutton', 'String', 'save','Callback',@save_callback ,'units','normalized','position',[0.6 0.65 0.05 0.05],'Fontsize',16,'Interruptible','on','UserData',0);
    end
    drawnow
    %     set(0,'DefaultAxesFontSize',get(0,'FactoryAxesFontSize')); % set back to default
    
    kk=1;
    while ~isempty(ismember(kk,(1:length(plot_channels))))
        uiwait(h1)
        
        if saving
            try
                b4=get(btn4,'UserData');
            catch
                %warning('Figure was closed without the exit button')
                return
            end
            
            if b4
                str_chan=sprintf('Channel_%s_spectral',chan_names{k});
                print(h1, fullfile(outputdir,str_chan), '-dpng', '-r0');
                set(btn4,'UserData',0); % reset buttons
                uiwait(h1) % re-lock figure
            end
        end
        
        try
            b1=get(btn1,'UserData');
            b2=get(btn2,'UserData');
        catch
            %warning('Figure was closed without the exit button')
            return
        end
        bb=b1+b2; % one should always be zero
        
        set(btn1,'UserData',0); % reset buttons
        set(btn2,'UserData',0);
        
        if get(btn3,'UserData')
            close all
            return
        elseif (kk+bb) >=1 && (kk+bb)<=length(plot_channels)
            kk=kk+bb;
            %         end
        elseif (kk+bb) < 1
            %ind=1;
            beep
            warning('This is the first channel')
        elseif (kk+bb)>length(plot_channels)
            %ind=ind-npoints*bb;
            beep
            warning('This is the last channel')
        end % **
        k=plot_channels(kk);
        
        %update figures
        temp=squeeze(tf(k,fl_ind:fh_ind,:));
        set(h11,'CData',temp)
        %title(h1sub,['time-frequency representation for channel ',chan_names{k}])
        %         hz = zoom;
        %         set(hz,'Motion','vertical','Enable','on');
        drawnow
        
        set(h13p,'ydata',Mp+SDp);
        set(h13m,'ydata',Mp-SDp);
        set(h13,'ydata',AS(k,:));
        %title(h3sub,['mean spectrum +/- std across all channels, and PS for channel ',chan_names{k}])
        drawnow
        
        set(h_ts,'ydata',eegfile_re_ref(:,k))
        %         title_str=sprintf('Channel %s, Index = %i out of %i',chan_names{k},k,numel(plot_channels)); % channel name and index
        %         title(h2_ax,title_str)
        str_sup=sprintf('Channel %s, Index = %i out of %i',chan_names{k},kk,numel(plot_channels));
        %suptitle(str_sup)
        set(h1,'name',str_sup)
        legend([h13p h13 h13m],'mean + std','spectrum','mean - std')
        drawnow
    end
end
end
%% Local functions
% callback functions for ui objects
function save_callback(src,event)
drawnow
uiresume
set(src,'UserData',1)
end

function btn1_callback(src,event)
drawnow
uiresume
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

% auxiliary function for ffts and moments
function alldata_struct = aux_func1(eegfile,Fs,T,overlap, parallel)
% should calculate nwins and initialize these variables

% if parallel
% 
%     parfor m = 1:size(eegfile,2) %%% channelMap(m) is the actual channel
%         [~, amps(m,:), tf_rep(m,:,:), freqs_FFT(m,:) , ~ , t_pow(m,:,:)] = stft_mt(eegfile(:,m),Fs,T,overlap,'window',@hann);
%         %[~, amps(m,:), tf_rep(m,:,:), freqs_FFT(m,:) , ~ , t_pow(m,:,:)] = stft_mt(eegfile(:,m),Fs,T,overlap,'multitapers',3);
% 
%         %[mu1(m,:), mu2(m,:), mu3(m,:), mu4(m,:), Vamp(m,:), H(m,:), t_mu(m,:,:)] = ts_info(eegfile(:,m),Fs,T,overlap,'ztrans_full'); % use 'ztrans_full' to standardize the full data or 'ztrans_seg' to do each segment individually before finding moments
%     end
% 
% else
    
    for m = 1:size(eegfile,2) %%% channelMap(m) is the actual channel
        
        fprintf('calculating FFT for channel %d\n', m);
        [~, amps(m,:), tf_rep(m,:,:), freqs_FFT(m,:) , ~ , t_pow(m,:,:)] = stft_mt(eegfile(:,m),Fs,T,overlap,'window',@hann);
    end
    
% end

alldata_struct = struct('frequencies',freqs_FFT,'amplitudes',amps,'time_frequency',tf_rep,'fft_time_windows',t_pow);
clear amps tf_rep freqs_FFT t_pow
end