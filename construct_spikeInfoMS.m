function construct_spikeInfoMS(varargin)

fprintf('****************************************************\n');
fprintf('* %s\n', mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse inputs

p = inputParser;

p.addParameter('session_path', '', @ischar);
p.addParameter('split_path', '', @ischar);
p.addParameter('analog_pulse_fpath', '', @ischar);
p.addParameter('nev_fpath', '', @ischar);
p.addParameter('saveRoot', '', @ischar);
p.addParameter('used_jacksheet_fpath', '', @ischar);

p.addParameter('skip_plots','0', @ischar);

p.addParameter('split_fname_suffix', 'mda_chan', @ischar);
p.addParameter('firings_fname', 'firings.mda', @ischar);
p.addParameter('metrics_fname', 'metrics.json', @ischar);
p.addParameter('isol_metrics_fname', 'isol_metrics.json', @ischar);
p.addParameter('isol_pair_metrics_fname', 'isol_pair_metrics.json', @ischar);
p.addParameter('clips_raw_fname', 'clips_raw.mda', @ischar);
p.addParameter('clips_spike_fname', 'clips_spike.mda', @ischar);
p.addParameter('clips_whiten_fname', 'clips_whiten.mda', @ischar);
p.addParameter('clip_features_fname', 'clip_features.mda', @ischar);

p.addParameter('noise_overlap_max', '0.1', @ischar);
p.addParameter('snr_min', '1', @ischar);
p.addParameter('removeLargeAmpUnits', '1', @ischar);

parse(p, varargin{:});

disp(p.Results);

session_path = p.Results.session_path;
split_path = p.Results.split_path;
analog_pulse_fpath = p.Results.analog_pulse_fpath;
nev_fpath = p.Results.nev_fpath;
saveRoot = p.Results.saveRoot;
used_jacksheet_fpath = p.Results.used_jacksheet_fpath;

removeLargeAmpUnits = eval(p.Results.removeLargeAmpUnits);
skip_plots = eval(p.Results.skip_plots);

noise_overlap_max = eval(p.Results.noise_overlap_max);
snr_min = eval(p.Results.snr_min);

if ~exist(saveRoot, 'dir')
    mkdir(saveRoot);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load nsx_postProc.mat

nsx_postProc = load([session_path '/nsx_postProc.mat']);
nsx_postProc = nsx_postProc.nsx_postProc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path contents

split_path_ls = dir(split_path);
split_path_ls_names = {split_path_ls.name};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load full jacksheet + used jacksheets

split_jacksheet = readtable(used_jacksheet_fpath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% how many channels were split for processing

num_orig_split_chan = size(split_jacksheet, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get session start time + session name

session_name = split_jacksheet{1, 'RawDir'}{1};
session_name_splits = strsplit(session_name, '_');
session_start_time_str = strjoin(session_name_splits(1:2), '_');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get list of session channels in this split_fpath

% is it a directory -- AND -- does it have a done.log file?
filtered_channel_ls = cellfun( @(f) ...
                                length(f) >= 3 ... %make sure its not a hidden file
                                && exist([split_path '/' f], 'dir') ... %make sure its a directory
                                && exist([split_path '/' f '/done.log'], 'file') ... %make sure there is a ms output file present, although this does not gauruntee that ms ran successfully
                                , split_path_ls_names);
               
%sort channel directories                            
channel_dirs = sort(split_path_ls_names(filtered_channel_ls));
num_channels = sum( filtered_channel_ls );


if num_orig_split_chan ~= num_channels
   error('num_orig_split_chan ~= num_channels --- the number of split files with results does not equal the number of channels originaly split for sorting');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load sorted channel data into table

fprintf('loading channels\n');
loaded_chan = load_chan_dirs(channel_dirs, split_path, p);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the channel fig

if ~skip_plots

    fprintf('plot figs\n');

    for iChan = 1:num_channels

        chan_name = split_jacksheet{iChan, 'ChanName'}{1};
        sortfigs_saveRoot = [saveRoot '/sortFigs'];

        plotChannelSpikes('session_name', session_name, ...
                          'channel_name', chan_name, ...
                          'clip_features_fpath', loaded_chan(iChan).clips_whiten_features_fpath, ...
                          'clips', loaded_chan(iChan).clips_raw, ...
                          'firings', loaded_chan(iChan).firings, ...
                          'isol_metrics', loaded_chan(iChan).isol_metrics, ...
                          'isol_pair_metrics', loaded_chan(iChan).isol_pair_metrics, ...
                          'metrics', loaded_chan(iChan).metrics, ...
                          'mda_fpath', loaded_chan(iChan).hp_reref_whiten_fpath, ...
                          'saveDir', sortfigs_saveRoot);

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filter the units on the loaded channels

[sessUnitSummary, sessUniqueUnitID, timeStamp, jackTableUsed, extractInfoStr, waveForm, waveForm_all, waveForm_raw, waveForm_raw_all, waveForm_sort, waveForm_sort_all, metrics, aux ] = filter_units(loaded_chan, noise_overlap_max, snr_min, removeLargeAmpUnits, split_jacksheet);

if aux.multi_unit_channel_present           
    aux_fid = fopen([saveRoot '/multi_channels.txt'], 'w');    
    fclose(aux_fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the pulses
fprintf('getting the sync pulses from analog: %s and nev: %s\n', analog_pulse_fpath, nev_fpath);
pulses = getBlackrockPulses_DC_AN('ns3_fpath', analog_pulse_fpath, 'nev_fpath', nev_fpath);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract start time from jacksheet

dateInfo_column = 'FileName';
old_fmt_regex_match = regexp( split_jacksheet{1, dateInfo_column}{1} , '\d\d\d\d\d\d_\d\d\d\d', 'match');

new_fmt_regex_match = regexp( split_jacksheet{1, dateInfo_column}{1} , '\d\d\d\d\d\d\d\d-\d\d\d\d\d\d', 'match');

if isempty(old_fmt_regex_match) && isempty(new_fmt_regex_match)

    error('no old or new format datestring found in used_jacksheet{1, "FileName"}');
elseif ~isempty(old_fmt_regex_match)

    startTime_datenum = datenum(old_fmt_regex_match{1}, 'yymmdd_HHMM');
elseif ~isempty(new_fmt_regex_match)

    startTime_datenum = datenum(new_fmt_regex_match{1}, 'yyyymmdd-HHMMSS');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now we can allocate the spikeInfo struct

spikeInfo_readme = '';

spikeInfo = struct;
spikeInfo.readme = spikeInfo_readme;
spikeInfo.createdDate = datestr(datetime('now'), 'mmm.dd,yyyy HH:MM:SS');
spikeInfo.sessUniqueUnitID = sessUniqueUnitID;
spikeInfo.timeStamp = timeStamp;
spikeInfo.waveForm = waveForm;
spikeInfo.waveForm_raw = waveForm_raw;
spikeInfo.waveForm_sort = waveForm_sort;
spikeInfo.metrics = metrics;
spikeInfo.sessStr = session_start_time_str;
spikeInfo.sessDurMin = split_jacksheet{1, 'DurationMin'};
spikeInfo.startTime_datenum = startTime_datenum;
spikeInfo.extractInfoStr = extractInfoStr;
spikeInfo.pulses = {split_jacksheet{1, 'NSPsuffix'}{1} split_jacksheet{1, 'FileName'}{1} pulses};
spikeInfo.alignedTo = '';
spikeInfo.alignmentChan = '';
spikeInfo.sort_dirname = 'reref_mountainsort';
spikeInfo.sort_filenames = {};
spikeInfo.jackTableSplit = split_jacksheet;
spikeInfo.jackTableUsed = jackTableUsed;
spikeInfo.sortNoteTable = {};
spikeInfo.physio_nsx_postProc = {split_jacksheet{1, 'NSPsuffix'}{1} split_jacksheet{1, 'FileName'}{1} nsx_postProc};



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now we can allocate the spikeWaveform struct

spikeWaveform_readme = '';

spikeWaveform = struct;
spikeWaveform.readme = spikeWaveform_readme;
spikeWaveform.waveForm_all = waveForm_all;
spikeWaveform.waveForm_raw_all = waveForm_raw_all;
spikeWaveform.waveForm_sort_all = waveForm_sort_all;




% 
% traceFigs_saveRoot = [saveRoot '/traceFigs'];
% if ~exist(traceFigs_saveRoot, 'dir')
%     mkdir(traceFigs_saveRoot);
% end
% 
% %plot the quantile stat figs
% for iRefset = 1:num_refsets
%     
% %    raw_quantile_mat_fpath = [ saveRoot sprintf('/../raw_quantiles.refset%d.mat', iRefset) ];
% %    bandpass_quantile_mat_fpath = [ saveRoot sprintf('/../bandpass_quantiles.refset%d.mat', iRefset) ];
%    reref_quantile_mat_fpath = [ saveRoot sprintf('/../reref_quantiles.refset%d.mat', iRefset) ];
% %    whiten_quantile_mat_fpath = [ saveRoot sprintf('/../whiten_quantiles.refset%d.mat', iRefset) ];
% 
%    quantile_save_root = [traceFigs_saveRoot '/' sess];  
%     
%    plotWindowStats('', '', reref_quantile_mat_fpath, '', quantile_save_root);
% end


% save 
fprintf('writing out mat files \n');
saveDir_spikeInfo = [saveRoot '/' session_start_time_str '_' split_jacksheet{1, 'NSPsuffix'}{1} '_spikeInfo.mat'];
save(saveDir_spikeInfo, '-v7.3', 'spikeInfo');

saveDir_spikeWaves = [saveRoot '/' session_start_time_str '_' split_jacksheet{1, 'NSPsuffix'}{1} '_spikeWaveform.mat'];
save(saveDir_spikeWaves, '-v7.3', 'spikeWaveform');

saveDir_sortSummary = [saveRoot '/' session_start_time_str '_' split_jacksheet{1, 'NSPsuffix'}{1} '_sortSummary.csv'];
writetable(sessUnitSummary, saveDir_sortSummary);

fprintf('deleting split files\n');

rmdir(split_path, 's');
rmdir([session_path '/splits_raw'], 's');
rmdir([session_path '/splits_sort'], 's');


fprintf('construct_spikeInfoMS -- done\n');


end


function [sessUnitSummary, sessUniqueUnitID, timeStamp, jackTableUsed, extractInfoStr, waveForm, waveForm_all, waveForm_raw, waveForm_raw_all, waveForm_sort, waveForm_sort_all, metrics , aux] = filter_units(loaded_chan, noise_overlap_max, snr_min, removeLargeAmpUnits, split_jacksheet)
    
    if removeLargeAmpUnits 
        amp_thresh = 5000;
    else
        amp_thresh = Inf;
    end

    clip_timeMS_leftBound = -2;
    clip_timeMS_rightBound = 5;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define outputs 
    
    sessUnitSummary_varNames = {'NSxChanName' 'PhysChanNum' 'UnitNum' 'ChanUnitName' 'DeviceNum' 'snr' 'noise_overlap', 'num_spikes'};
    sessUnitSummary = table([],[],[],[],[],[],[],[], 'VariableNames', sessUnitSummary_varNames);
    
    sessUniqueUnitID_varNames = {'PhysChanNum' 'UnitNum' 'DeviceNum' 'CombinedNum' 'ChanUnitName' 'NSxChanName' 'NSxFileName' 'NSPsuffix' 'ChanNameNew'};
    sessUniqueUnitID = table([], [], [], [], [], [], [], [], [], 'VariableNames', sessUniqueUnitID_varNames);

    jackTableUsed = table([],[],[],[],[],[],[],[],[],[],[],[],[],[],[], 'VariableNames', split_jacksheet.Properties.VariableNames);

    timeStamp = {};
    extractInfoStr = {};

    waveForm = struct;
    waveForm_all = {};
    
    waveForm_raw = struct;
    waveForm_raw_all = {};
    
    waveForm_sort = struct;
    waveForm_sort_all = {};
    
    waveForm.infoStr = 'mean and sd in uV, time in ms; global reref, spike-band filtered (600-6000Hz)';
    waveForm.timeMS = [];
    waveForm.mean = [];
    waveForm.sd = [];

    waveForm_raw.infoStr = 'mean and sd in uV, time in ms; global reref, unfiltered (1-5000Hz)';
    waveForm_raw.timeMS = [];
    waveForm_raw.mean = [];
    waveForm_raw.sd = [];

    waveForm_sort.infoStr = 'mean and sd in uV, time in ms; global reref, spike-band filtered (600-6000Hz), whitened (the actual sorted data)';
    waveForm_sort.timeMS = [];
    waveForm_sort.mean = [];
    waveForm_sort.sd = [];
    
    
    metrics_readme = [ 'spike Info metrics contains: ' newline ...
                       'mountainsort is a struct containing isolation, noise_overlap, peak_amp, peak_noise, peak_snr, pair_overlap . See original paper for all details https://doi.org/10.1016/j.neuron.2017.08.030' newline ...
                       '   noise_overlap -- (0 - 1) 0 indicating no overlap with generated noise cluster, 1 indicating unit indistinguishable from the generated noise cluster' newline ...
                       '   peak_amp -- average value at peak/trough of unit waveform' newline ...
                       '   peak_noise -- sd of values at peak/trough of unit waveform' newline ...
                       '   peak_snr -- peak_amp/peak_noise' newline ...
                       'ISIcontamination   -- vector of fraction of ISIs that were less than absolute refractory period of 2ms' newline ...
                       'spkRate_eachMinute -- cell array with time series of spike rate (spikes/sec) per 1min interval for duration of recording. This is a first order stability measure.' ];
    
    metrics = struct;
    metrics.readme = metrics_readme;

    metrics.fracISIbelowThr = [];
    metrics.spkRate_eachMinute = {};

    metrics.mountainsort = struct;
%     metrics.mountainsort.isolation = [];
    metrics.mountainsort.noise_overlap = [];
    metrics.mountainsort.peak_amp = [];
    metrics.mountainsort.peak_noise = [];
    metrics.mountainsort.peak_snr = [];
%     metrics.mountainsort.pair_overlap = [];
    
    aux = struct;
    aux.multi_unit_channel_present = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    num_channels = length(loaded_chan);
    num_total_units_added = 0;
    
    for iChan = 1:num_channels
        
        % how many units found on this channel
        channel_unique_units = unique(loaded_chan(iChan).firings(3,:));        
        num_units = length(channel_unique_units);
        
        selected_units_mountainsort_ID = [];
        selected_units_combinedNum = [];
        selected_units_total_count_inx = [];
        
        chan_SampFreq = split_jacksheet{iChan,'SampFreq'};
        sampFreq2MS_factor = chan_SampFreq/1000;
        
        % isol_metrics labels
        channel_isol_metrics_labels = [ loaded_chan(iChan).isol_metrics.clusters.label ];
        
        % incrementers
        num_channel_units_added = 0;
        num_channel_unit_spikes = 0;
        
        first_channel_unit = 1;
        
        
        for iUnit = 1:num_units
            
            unit_firings_name = channel_unique_units(iUnit);
            
            which_isol_metrics_idx = (channel_isol_metrics_labels == unit_firings_name);
            
            %does this unit pass the filter
            
            unit_amp = abs(loaded_chan(iChan).isol_metrics.clusters(iUnit).metrics.peak_amp);
            unit_snr = loaded_chan(iChan).isol_metrics.clusters(iUnit).metrics.peak_snr;
            unit_noise_overlap = loaded_chan(iChan).isol_metrics.clusters(iUnit).metrics.noise_overlap;
            
            if (unit_snr >= snr_min) && (unit_noise_overlap <= noise_overlap_max) && (unit_amp < amp_thresh)
                
                % increment unit counter for this channel
                num_channel_units_added = num_channel_units_added + 1;

                % increment unit counter for entire session
                num_total_units_added = num_total_units_added + 1;
                
                
                % add unit to sessUniqueUnitID
                
                PhysChanNum = split_jacksheet{iChan, 'PhysicalChan'};
                UnitNum = num_channel_units_added;
                DeviceNum = split_jacksheet{iChan, 'MicroDevNum'};
                CombinedNum = (DeviceNum * 1e6) + (PhysChanNum * 1e3) + (UnitNum); 
                ChanUnitName = [ split_jacksheet{iChan, 'ChanNameNew'}{1} '_' sprintf('%03d', UnitNum)];
                NSxChanName = split_jacksheet{iChan, 'ChanName'}{1};
                NSxFileName = split_jacksheet{iChan, 'FileName'}{1};
                NSPsuffix = split_jacksheet{iChan, 'NSPsuffix'}{1};
                ChanNameNew = split_jacksheet{iChan, 'ChanNameNew'}{1};
                
                sessUniqueUnitID = [ sessUniqueUnitID;  table(PhysChanNum, UnitNum, DeviceNum, CombinedNum, {ChanUnitName}, {NSxChanName}, {NSxFileName}, {NSPsuffix}, {ChanNameNew}, 'VariableNames', sessUniqueUnitID_varNames)];
                
                
                % check if this channel needs to be added to jackTableUsed
                
                if first_channel_unit
                   jackTableUsed = [ jackTableUsed ; split_jacksheet(iChan,:) ];
                   first_channel_unit = 0;
                end
                
                
                % store unit timestamps, converted to milliseconds
                
                firings_filt = (loaded_chan(iChan).firings(3,:) == unit_firings_name);
                
                        % increment unit spikes counter for this channel
                        num_channel_unit_spikes = num_channel_unit_spikes + sum(firings_filt);
                      
                unit_timeStamp_ms = loaded_chan(iChan).firings(2,firings_filt) ./ sampFreq2MS_factor;
                timeStamp{num_total_units_added, 1} = unit_timeStamp_ms;
                
                
                
                % store waveForms
                
                clip_size = size( loaded_chan(iChan).clips_spike , 2);
                clip_midpoint = ceil(clip_size/2);
                
                clip_timeMS = ((1:clip_size) - clip_midpoint) ./ sampFreq2MS_factor;
                clip_timeMS_boundsFilt = (clip_timeMS >= clip_timeMS_leftBound) & (clip_timeMS <= clip_timeMS_rightBound);
                
                
                % first "spike" waveForms
                
                waveForm_all{num_total_units_added, 1} = loaded_chan(iChan).clips_spike(firings_filt, clip_timeMS_boundsFilt);
                
                waveForm.timeMS = clip_timeMS(clip_timeMS_boundsFilt);
                waveForm.mean(num_total_units_added, 1:length(waveForm.timeMS)) = mean(waveForm_all{num_total_units_added, 1}, 1);
                waveForm.sd(num_total_units_added, 1:length(waveForm.timeMS)) = std(waveForm_all{num_total_units_added, 1}, 0, 1);
                
                
                % raw waveForms
                
                waveForm_raw_all{num_total_units_added, 1} = loaded_chan(iChan).clips_raw(firings_filt, clip_timeMS_boundsFilt);
                
                waveForm_raw.timeMS = clip_timeMS(clip_timeMS_boundsFilt);
                waveForm_raw.mean(num_total_units_added, 1:length(waveForm_raw.timeMS)) = mean(waveForm_raw_all{num_total_units_added, 1}, 1);
                waveForm_raw.sd(num_total_units_added, 1:length(waveForm_raw.timeMS)) = std(waveForm_raw_all{num_total_units_added, 1}, 0, 1);
                
                
                % sorted waveForms
                
                waveForm_sort_all{num_total_units_added, 1} = loaded_chan(iChan).clips_whiten(firings_filt, clip_timeMS_boundsFilt);
                
                waveForm_sort.timeMS = clip_timeMS(clip_timeMS_boundsFilt);
                waveForm_sort.mean(num_total_units_added, 1:length(waveForm_sort.timeMS)) = mean(waveForm_sort_all{num_total_units_added, 1}, 1);
                waveForm_sort.sd(num_total_units_added, 1:length(waveForm_sort.timeMS)) = std(waveForm_sort_all{num_total_units_added, 1}, 0, 1);
                
                
                % add ISI and spikerate metrics to struct
                
                unit_ISI = diff(unit_timeStamp_ms);
                metrics.fracISIbelowThr(num_total_units_added, 1) = sum(unit_ISI < 2) / length(unit_timeStamp_ms) ;
                
                
                ms_per_min = 1000 * 60;
                unit_timeStamp_min = unit_timeStamp_ms ./ ms_per_min ;
                max_whole_min = max(floor(unit_timeStamp_min));
                
                spkRate_eachMinute = [];
                
                for iMin = 0:max_whole_min
                    spkRate_eachMinute(iMin + 1) = sum((unit_timeStamp_min >= iMin) & (unit_timeStamp_min < (iMin+1))) / 60; % convert to spikes/sec
                end
                
                metrics.spkRate_eachMinute{num_total_units_added, 1} = spkRate_eachMinute;
                
                

                % add pre-calculated unit metrics to struct

                metrics.mountainsort.noise_overlap(num_total_units_added) = loaded_chan(iChan).isol_metrics.clusters(which_isol_metrics_idx).metrics.noise_overlap;
                metrics.mountainsort.peak_amp(num_total_units_added) = loaded_chan(iChan).isol_metrics.clusters(which_isol_metrics_idx).metrics.peak_amp;
                metrics.mountainsort.peak_noise(num_total_units_added) = loaded_chan(iChan).isol_metrics.clusters(which_isol_metrics_idx).metrics.peak_noise;
                metrics.mountainsort.peak_snr(num_total_units_added) = loaded_chan(iChan).isol_metrics.clusters(which_isol_metrics_idx).metrics.peak_snr;

                
                sessUnitSummary = [sessUnitSummary; table({NSxChanName},PhysChanNum,UnitNum,{ChanUnitName},DeviceNum, metrics.mountainsort.peak_snr(num_total_units_added),metrics.mountainsort.noise_overlap(num_total_units_added), sum(firings_filt) ,'VariableNames', sessUnitSummary_varNames)];


                % remember this unit passed the filter
                
                selected_units_mountainsort_ID(num_channel_units_added) = unit_firings_name;
                selected_units_combinedNum(num_channel_units_added) = CombinedNum;
                selected_units_total_count_inx(num_channel_units_added) = num_total_units_added;
                
                
            end
            
        end
        
        if num_channel_units_added > 0
        
            % store extract info string after all units on this channel have been filtered
            extractInfoStr{length(extractInfoStr) + 1, 1} = [ split_jacksheet{iChan, 'ChanName'}{1} ' --> ' num2str(num_channel_units_added) ' units, ' num2str(num_channel_unit_spikes) ' total spikes' ];

        end
        
        
        % add in pairwise isolation scores for units passing the filter
       
        if num_channel_units_added > 1
            aux.multi_unit_channel_present = 1;           
        end 
    end

end


function loaded_chan = load_chan_dirs(channel_dirs, split_path, parser)

    split_fname_suffix = parser.Results.split_fname_suffix;
    firings_fname = parser.Results.firings_fname;
    metrics_fname = parser.Results.metrics_fname;
    isol_metrics_fname = parser.Results.isol_metrics_fname;
    isol_pair_metrics_fname = parser.Results.isol_pair_metrics_fname;
    clip_features_fname = parser.Results.clip_features_fname;
    clips_raw_fname = parser.Results.clips_raw_fname;
    clips_spike_fname = parser.Results.clips_spike_fname;
    clips_whiten_fname = parser.Results.clips_whiten_fname;

    num_channels = length(channel_dirs);

    loaded_chan = struct;
    
    for iChan = 1:num_channels

        fprintf('loading %s\n', channel_dirs{iChan});

        channel_dir = [ split_path '/' channel_dirs{iChan} ];

        hp_reref_whiten_fpath = [ channel_dir '/' channel_dirs{iChan} '.' split_fname_suffix];
        hp_reref_fpath = [ split_path '/../splits_spike/' channel_dirs{iChan} '/' channel_dirs{iChan} '.' split_fname_suffix];
        lp_reref_fpath = [ split_path '/../splits_raw/' channel_dirs{iChan} '/' channel_dirs{iChan} '.' split_fname_suffix];

        firings_fpath = [channel_dir '/' firings_fname];
        firings = readmda(firings_fpath);

        metrics_fpath = [channel_dir '/' metrics_fname];
        metrics = readjson(metrics_fpath);

        isol_metrics_fpath = [channel_dir '/' isol_metrics_fname];
        isol_metrics = readjson(isol_metrics_fpath);

        isol_pair_metrics_fpath = [channel_dir '/' isol_pair_metrics_fname];
        isol_pair_metrics = readjson(isol_pair_metrics_fpath);

        clips_raw_fpath = [channel_dir '/' clips_raw_fname];
        clips_raw = squeeze(readmda(clips_raw_fpath))';

        clips_spike_fpath = [channel_dir '/' clips_spike_fname];
        clips_spike = squeeze(readmda(clips_spike_fpath))';

        clips_whiten_fpath = [channel_dir '/' clips_whiten_fname];
        clips_whiten = squeeze(readmda(clips_whiten_fpath))';

        clip_features_fpath = [channel_dir '/' clip_features_fname];

        loaded_chan(iChan).chan_dir = channel_dir;
        loaded_chan(iChan).firings = firings;
        loaded_chan(iChan).clips_raw = clips_raw;
        loaded_chan(iChan).clips_spike = clips_spike;
        loaded_chan(iChan).clips_whiten = clips_whiten;
        loaded_chan(iChan).metrics = metrics;
        loaded_chan(iChan).isol_metrics = isol_metrics;
        loaded_chan(iChan).isol_pair_metrics = isol_pair_metrics;
        loaded_chan(iChan).clips_whiten_features_fpath = clip_features_fpath;
        loaded_chan(iChan).lp_reref_fpath = lp_reref_fpath;
        loaded_chan(iChan).hp_reref_fpath = hp_reref_fpath;
        loaded_chan(iChan).hp_reref_whiten_fpath = hp_reref_whiten_fpath;
        
    end

end


function jsonstring = readjson(fname)

fid = fopen(fname);
s = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
jsonstring = jsondecode(strjoin(s{1}, ' '));

end



