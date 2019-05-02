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
p.addParameter('full_jacksheet_fpath', '', @ischar);
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
full_jacksheet_fpath = p.Results.full_jacksheet_fpath;
used_jacksheet_fpath = p.Results.used_jacksheet_fpath;

removeLargeAmpUnits = eval(p.Results.removeLargeAmpUnits);
skip_plots = eval(p.Results.skip_plots);

noise_overlap_max = eval(p.Results.noise_overlap_max);
snr_min = eval(p.Results.snr_min);

if ~exist(saveRoot, 'dir')
    mkdir(saveRoot);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path contents

session_path_ls = dir(session_path);
session_path_ls_names = {session_path_ls.name};

split_path_ls = dir(split_path);
split_path_ls_names = {split_path_ls.name};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load full jacksheet + used jacksheets

used_jacksheet = readtable(used_jacksheet_fpath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% how many channels were split for processing

num_orig_split_chan = size(used_jacksheet, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get session start time + session name

session_name = used_jacksheet{1, 'RawDir'}{1};
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
% collect the number of spikes across all channels

fprintf('count units\n');

channelTotals = zeros(num_channels, 1);
num_spikes = 0;
num_unique_units = 0;

for iChan = 1:num_channels
    
    current_firings = loaded_chan(iChan).firings;
    
    num_spikes = num_spikes + size(current_firings, 2);
    num_unique_units = num_unique_units + length(unique(current_firings(3,:)));
    channelTotals(iChan) = length(unique(current_firings(3,:)));
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the channel fig

if ~skip_plots

    fprintf('plot figs\n');

    for iChan = 1:num_channels

        chan_name = used_jacksheet{iChan, 'ChanName'}{1};
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

[sessUniqueUnitID, timeStamp, jackTableUsed, extractInfoStr, waveForm, waveForm_all, waveForm_raw, waveForm_raw_all, metrics ] = filter_units(loaded_chan, noise_overlap_max, snr_min, removeLargeAmpUnits, used_jacksheet);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now we can allocate the spikeInfo struct
% 
% spikeInfo_readme = '';
% metrics_readme = '';
% 
% spikeInfo = struct;
% spikeInfo.readme = spikeInfo_readme;
% spikeInfo.createdDate = datestr(datetime('now'), 'mmm.dd,yyyy HH:MM:SS');
% spikeInfo.sessUniqueUnitID = [];
% spikeInfo.timeStamp = cell(num_unique_units, 1);
% spikeInfo.waveForm = struct;
% spikeInfo.waveForm_raw = struct;
% spikeInfo.metrics = struct;
% spikeInfo.sessStr = session_start_time_str;
% spikeInfo.sessDurMin = used_jacksheet{1, 'DurationMin'};
% spikeInfo.startTime_datenum = 0;
% spikeInfo.extractInfoStr = {};
% spikeInfo.pulses = {};
% spikeInfo.alignedTo = '';
% spikeInfo.alignmentChan = '';
% spikeInfo.sort_dirname = 'reref_mountainsort';
% spikeInfo.sort_filenames = {};
% spikeInfo.jackTableSplit = used_jacksheet;
% spikeInfo.jackTableUsed = [];
% spikeInfo.sortNoteTable = {};
% 
% spikeInfo.waveForm.infoStr = 'mean and sd in uV, time in ms;  spike-band filtered (600-6000Hz)';
% spikeInfo.waveForm.timeMS = [];
% spikeInfo.waveForm.mean = [];
% spikeInfo.waveForm.sd = [];
% 
% spikeInfo.waveForm_raw.infoStr = 'mean and sd in uV, time in ms;  unfiltered (1-5000Hz) after global reref';
% spikeInfo.waveForm_raw.timeMS = [];
% spikeInfo.waveForm_raw.mean = [];
% spikeInfo.waveForm_raw.sd = [];

spikeInfo.metrics.fracISIbelowThr = zeros(num_unique_units, 1);
spikeInfo.metrics.spkRate_eachMinute = zeros(num_unique_units, 1);

spikeInfo.metrics.readme = metrics_readme;
spikeInfo.metrics.mountainsort = struct;
spikeInfo.metrics.mountainsort.isolation = zeros(num_unique_units, 1);
spikeInfo.metrics.mountainsort.noise_overlap = zeros(num_unique_units, 1);
spikeInfo.metrics.mountainsort.peak_amp = zeros(num_unique_units, 1);
spikeInfo.metrics.mountainsort.peak_noise = zeros(num_unique_units, 1);
spikeInfo.metrics.mountainsort.peak_snr = zeros(num_unique_units, 1);
spikeInfo.metrics.mountainsort.pair_overlap = zeros(num_unique_units, 1);

%                readme: 'this spikeInfo file, generated Mar.21,2019 15:38:59, ??
%           createdDate: 'Mar.21,2019 15:38:59'
%      sessUniqueUnitID: [4�9 table]
%             timeStamp: {4�1 cell}
%              waveForm: [1�1 struct]
%          waveForm_raw: [1�1 struct]
%               metrics: [1�1 struct]
%               sessStr: '190117_1720'
%            sessDurMin: 4.9000
%     startTime_datenum: 7.3744e+05
%        extractInfoStr: {2�1 cell}
%                pulses: {'INST0'  'YfGGmZkGMfb-20190117-172012-INST0.ns5'  [1�1 struct]}
%             alignedTo: ''
%         alignmentChan: ''
%          sort_dirname: 'reref(stimMask5)_sortedByJW'
%        sort_filenames: {'[INST0]utah_amtg01.txt'  '[INST0]utah_pmtg01.txt'}
%        jackTableSplit: [128�16 table]
%         jackTableUsed: [2�16 table]
%         sortNoteTable: {'sortNotes(190117_1720)_sortedByJW.xlsx'  [128�7 table]}
% 


sessUniqueUnitID = [];
timeStamp = {};
jackTableUsed = [];
extractInfoStr = {};

waveForm = struct;
waveForm_raw = struct;

waveForm.infoStr = 'mean and sd in uV, time in ms;  spike-band filtered (600-6000Hz)';
waveForm.timeMS = [];
waveForm.mean = [];
waveForm.sd = [];

waveForm_raw.infoStr = 'mean and sd in uV, time in ms;  unfiltered (1-5000Hz) after global reref';
waveForm_raw.timeMS = [];
waveForm_raw.mean = [];
waveForm_raw.sd = [];

metrics = struct;
metrics.readme = metrics_readme;

metrics.fracISIbelowThr = [];
metrics.spkRate_eachMinute = [];

metrics.mountainsort = struct;
metrics.mountainsort.isolation = [];
metrics.mountainsort.noise_overlap = [];
metrics.mountainsort.peak_amp = [];
metrics.mountainsort.peak_noise = [];
metrics.mountainsort.peak_snr = [];
metrics.mountainsort.pair_overlap = [];

              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
processed_units = 1;

%for each channel file, add the data to the struct
for iChan = 1:size(data_cmat, 1)
    
    data_cmat(iChan,:);
    
    currentFirings = data_cmat{iChan, firings_col};
    currentClips = data_cmat{iChan, clips_col};    
    currentMetrics = data_cmat{iChan, metrics_col}; 
    currentIsolMetrics = data_cmat{iChan, isol_metrics_col};
    currentIsolPairMetrics = data_cmat{iChan, isol_pair_metrics_col};
    currentChannelData = readmda(data_cmat{iChan, bp_col});
        
    spikeInfo.sessDurSec = length(currentChannelData)/nsxSamplingFreq;
    
    current_channel_dname = channel_dirs{iChan};
    fprintf('adding to struct: %s\n', current_channel_dname);
    
    %get the electrode name (NIHXXX_date_time_elec -> elec)
    channel_dname_splits = strsplit(current_channel_dname, '_');
    current_channel_ename = channel_dname_splits{end};
    
    %how many units are there ?
    current_channel_num_unit = length(unique(currentFirings(3,:)));
    
    %create the unit names (NIHXXX_date_time_channel --> channel_3digitnum
    current_channel_unit_names = cellstr(cellfun(@(x) [current_channel_ename '_' sprintf( '%03d', x )] , num2cell(1:current_channel_num_unit), 'UniformOutput', 0));
    
    %how many total spikes?
    current_channel_num_spikes = size( currentFirings , 2);
    spikeInfo.extractInfoStr{iChan,1} = [current_channel_ename ' --> ' int2str(current_channel_num_unit) ' units, ' int2str(current_channel_num_spikes) ' total spikes'];
        
    %make sure we even have a single unit
    if ~isempty(currentFirings)
        
        for iChanUnit = 1:current_channel_num_unit
            
            current_unitID = current_channel_unit_names{iChanUnit};
            
            spikeInfo.uniqueUnitID{ processed_units } = current_unitID;
            
            current_unit_filter = ( iChanUnit == currentFirings(3,:) );
            
            %fill in timeStamp values
            spikeInfo.timeStamp{ processed_units } = currentFirings(2, current_unit_filter);
            
            %fill in the waveform values
            spikeInfo.waveForm{ processed_units } = currentClips(1:clip_size, current_unit_filter)';

            %fill in the metrics
            %combine the metrics
            
            dur_sec = currentMetrics(iChanUnit).metrics.dur_sec;
            if ~isempty(dur_sec)
                spikeInfo.metrics.dur_sec(processed_units , 1) = dur_sec;
            else
                spikeInfo.metrics.dur_sec(processed_units , 1) = 0;
            end
            
            firing_rate = currentMetrics(iChanUnit).metrics.firing_rate;
            if ~isempty(firing_rate)
                spikeInfo.metrics.firing_rate(processed_units , 1) = firing_rate;
            else
                spikeInfo.metrics.firing_rate(processed_units , 1) = 0;
            end
            
            num_events = currentMetrics(iChanUnit).metrics.num_events;
            if ~isempty(num_events)
                spikeInfo.metrics.num_events(processed_units , 1) = num_events;
            else
                spikeInfo.metrics.num_events(processed_units , 1) = 0;
            end
            
            t1_sec = currentMetrics(iChanUnit).metrics.t1_sec;
            if ~isempty(t1_sec)
                spikeInfo.metrics.t1_sec(processed_units , 1) = t1_sec;
            else
                spikeInfo.metrics.t1_sec(processed_units , 1) = 0;
            end
            
            t2_sec = currentMetrics(iChanUnit).metrics.t2_sec;
            if ~isempty(t2_sec)
                spikeInfo.metrics.t2_sec(processed_units , 1) = t2_sec;
            else
                spikeInfo.metrics.t2_sec(processed_units , 1) = 0;
            end
            
            burst_parent = currentIsolMetrics(iChanUnit).metrics.bursting_parent;
            if ~isempty(burst_parent)
                if burst_parent > 0
                    spikeInfo.metrics.bursting_parent{processed_units , 1} = current_channel_unit_names{burst_parent};
                else
                    spikeInfo.metrics.bursting_parent{processed_units , 1} = burst_parent;
                end
            else
                spikeInfo.metrics.bursting_parent{processed_units , 1} = 0;
            end
            
            isolation = currentIsolMetrics(iChanUnit).metrics.isolation;
            if ~isempty(isolation)
                spikeInfo.metrics.isolation(processed_units , 1) = isolation;
            else
                spikeInfo.metrics.isolation(processed_units , 1) = 0;
            end
            
            noise_overlap = currentIsolMetrics(iChanUnit).metrics.noise_overlap;
            if ~isempty(noise_overlap)
                spikeInfo.metrics.noise_overlap(processed_units , 1) = noise_overlap;
            else
                spikeInfo.metrics.noise_overlap(processed_units , 1) = 0;
            end
            
            peak_amp = currentIsolMetrics(iChanUnit).metrics.peak_amp;
            if ~isempty(peak_amp)
                spikeInfo.metrics.peak_amp(processed_units , 1) = peak_amp;
            else
                spikeInfo.metrics.peak_amp(processed_units , 1) = 0;
            end
            
            peak_noise = currentIsolMetrics(iChanUnit).metrics.peak_noise;
            if ~isempty(peak_noise)
                spikeInfo.metrics.peak_noise(processed_units , 1) = peak_noise;
            else
                spikeInfo.metrics.peak_noise(processed_units , 1) = 0;
            end
            
            peak_snr = currentIsolMetrics(iChanUnit).metrics.peak_snr;
            if ~isempty(peak_snr)
                spikeInfo.metrics.peak_snr(processed_units , 1) = peak_snr;
            else
                spikeInfo.metrics.peak_snr(processed_units , 1) = 0;
            end
            
            pair_struct = struct;
            pair_struct.pair = {};
            pair_struct.overlap = [];
            
            iChanUnit_count = 0;
            
            if ~isempty(currentIsolPairMetrics)
                
                for iPair = 1:length(currentIsolPairMetrics)
                    
                    current_pair = cellfun(@str2num , cellstr(strsplit(currentIsolPairMetrics(iPair).label , ',')));
                    current_pair_filter = (current_pair == iChanUnit);
                    
                    if sum(current_pair_filter) > 0
                        
                        iChanUnit_count = iChanUnit_count + 1;
                        
                        overlap_val = currentIsolPairMetrics(iPair).metrics.overlap;
                        
                        pair_struct.pair{ length(pair_struct.pair) + 1, 1} = [current_channel_unit_names{current_pair(1)} ',' current_channel_unit_names{current_pair(2)}];
                        pair_struct.overlap{ length(pair_struct.overlap)  + 1, 1} = overlap_val;
                        
                    end
                    
                end
                
            end
            
            spikeInfo.metrics.pair_overlap{processed_units, 1} = pair_struct;
            
            fprintf('length of iChanUnit_count: %d -- length of pair_struct.pair %d\n', iChanUnit_count, length(pair_struct.pair));
            
            if iChanUnit_count ~= length(pair_struct.pair)
               keyboard; 
            end
            
            processed_units = processed_units + 1;
        end
        
    end
    
end

clear data_cmat;

    
fprintf('filtering spikeInfo\n');



orig_unitID_col = 1;
filt_unitID_col = 2;
filter_type_col = 3;
unit_channel_col = 4;
num_spikes_col = 5;
snr_col = 6;
isolation_col = 7;
noise_overlap_col = 8;

%liberal filter, filtering out the blatanly obvious noise unit
noise_overlap_max = 0.1;
snr_min = 1;
%not relevant for filtering blantant noise
isolation_min = 0.95;

filter_string = sprintf('spike clusters output from mountainsort were 1) filtered out as "noise" if noise overlap > %0.2f or SNR < %0.2f or 2) labelled as "non-isolated" if between-cluster isolation < %0.2f', noise_overlap_max, snr_min, isolation_min );

[spikeInfo_filt, noiseInfo, filter_dmat] = filterSpikes(spikeInfo, noise_overlap_max, snr_min, isolation_min, removeLargeAmpUnits);



%get the pulses
fprintf('getting the sync pulses from ns3: %s and nev: %s\n', ns3_pulse_fpath, nev_fpath);
spikeInfo_filt.pulses = getBlackrockPulses_DC_AN('ns3_fpath', ns3_pulse_fpath, 'nev_fpath', nev_fpath);


if ~exist(saveRoot, 'dir')
    mkdir(saveRoot);
end

traceFigs_saveRoot = [saveRoot '/traceFigs'];
if ~exist(traceFigs_saveRoot, 'dir')
    mkdir(traceFigs_saveRoot);
end

%plot the quantile stat figs
for iRefset = 1:num_refsets
    
%    raw_quantile_mat_fpath = [ saveRoot sprintf('/../raw_quantiles.refset%d.mat', iRefset) ];
%    bandpass_quantile_mat_fpath = [ saveRoot sprintf('/../bandpass_quantiles.refset%d.mat', iRefset) ];
   reref_quantile_mat_fpath = [ saveRoot sprintf('/../reref_quantiles.refset%d.mat', iRefset) ];
%    whiten_quantile_mat_fpath = [ saveRoot sprintf('/../whiten_quantiles.refset%d.mat', iRefset) ];

   quantile_save_root = [traceFigs_saveRoot '/' sess];  
    
   plotWindowStats('', '', reref_quantile_mat_fpath, '', quantile_save_root);
end


spikeWaveform = spikeInfo_filt.waveForm;
spikeInfo_filt = rmfield(spikeInfo_filt, 'waveForm');

about_me_spikeWaveform = ['this file contains the spike waveforms extracted from the bandpass filtered, whitened, and common average re-referenced data' newline ...
                          'the indices of the cell array correspond to unit IDs found in the spikeInfo file'];

timeStamp           = spikeInfo_filt.timeStamp;
sessUniqueUnitID    = expandUnitIDInfo(spikeInfo_filt.uniqueUnitID, referencing_info, used_chan_by_refset);
metrics             = spikeInfo_filt.metrics;     
extractInfoStr      = spikeInfo_filt.extractInfoStr;
sessDurSec          = spikeInfo_filt.sessDurSec;
sessStr             = spikeInfo_filt.sessStr;
alignedTo           = spikeInfo_filt.alignedTo;
alignmentChan       = spikeInfo_filt.alignmentChan;
startTime_mstime    = spikeInfo_filt.startTime_mstime;
pulses              = spikeInfo_filt.pulses;


about_me_spikeInfo = ['this spikeInfo file contains the following fields:' newline newline ...
                      '     sessStr - the session name the spikes were extracted from' newline ...
                      '     sessDurSec - the duration of session in seconds' newline ...
                      '     startTime_mstime - the start time of session in mstime (note that these values are taken from the timestamp in the nsx filename, not the incorrectly offset time values in the original nsx file)'  newline ...
                      '     extractInfoStr - a string describing how many potential units were found by mountainsort on each channel, and the subset remaining after filtering on unit metrics' newline ...
                      '     timeStamp - a ( #unit x 1 ) cell array with each cell containing the 30kHz timestamps of the correspondingly indexed unit in sessUniqueUnitID' newline ...
                      '     sessUniqueUnitID - a ( #unit x 6 ) cell array with each row containing ID info for a unit ( channel_num, unit_num, array_location_num, combined_num, array_num, nsx_channel_name)' newline ...
                      '     referencing_info - a cell array containing: the referencing set name, subject name, channel string, channel num set, array number, array location number ( temporal = 1, parietal = 2, frontal = 9 ), and an array description' newline ...
                      '     metrics - a struct containing various quality scores for each unit ( potential units with SNR < 1 and noise_overlap > 0.1 were considered "non-units" and were not included in this file' newline ...
                      '     pulses - a struct containing 30kHz uptimes for all "ain" and "din" channels, as well as 1kHz downsampled "ain" timeseries' newline ...
                      '     alignedTo - a string indicating which file the spikes have been aligned to' newline ...
                      '     alignmentChan - a string indicating which channel in the pulses struct was used for alignment' newline];


%save these out without plotting
fprintf('writing out mat files with liberal filter\n');
saveDir_spikeInfo = [saveRoot '/' sess '_spikeInfo.mat'];
save(saveDir_spikeInfo, '-v7.3', 'timeStamp', 'sessUniqueUnitID', 'referencing_info', 'metrics', 'extractInfoStr', 'sessDurSec', 'sessStr', 'filter_string', 'alignedTo', 'alignmentChan', 'startTime_mstime', 'pulses', 'about_me_spikeInfo');

saveDir_spikeWaves = [saveRoot '/' sess '_spikeWaveform.mat'];
save(saveDir_spikeWaves, '-v7.3', 'spikeWaveform', 'about_me_spikeWaveform');




saveDir = [saveRoot '/../' sess '_fullSortSummary.csv'];
filt_sum_fid = fopen(saveDir,'w');
fprintf(filt_sum_fid, '"unfilter_unitID","filter_unitID","filtered_type","channel","num_spikes","snr","isolation","noise_overlap"\n');
for i = 1:size(filter_dmat,1)
    fprintf(filt_sum_fid, '"%s","%s","%s","%s","%s","%s","%s","%s"\n', ... 
                          filter_dmat{i,orig_unitID_col}, ...
                          filter_dmat{i,filt_unitID_col}, ...
                          filter_dmat{i,filter_type_col}, ...
                          filter_dmat{i,unit_channel_col}, ...
                          num2str(filter_dmat{i,num_spikes_col}), ...
                          num2str(filter_dmat{i,snr_col}), ...
                          num2str(filter_dmat{i,isolation_col}), ...
                          num2str(filter_dmat{i,noise_overlap_col}) ...
                          );
end
fclose(filt_sum_fid);




saveDir = [saveRoot '/' sess '_sortSummary.csv'];
filt_sum_fid = fopen(saveDir,'w');
fprintf(filt_sum_fid, '"unitID","filtered_type","channel","num_spikes","snr","isolation","noise_overlap"\n');
for i = 1:size(filter_dmat,1)
    
    if ~isequal(filter_dmat{i,filter_type_col}, 'noise')
        
        fprintf(filt_sum_fid, '"%s","%s","%s","%s","%s","%s","%s"\n', ... 
                          filter_dmat{i,filt_unitID_col}, ...
                          filter_dmat{i,filter_type_col}, ...
                          filter_dmat{i,unit_channel_col}, ...
                          num2str(filter_dmat{i,num_spikes_col}), ...
                          num2str(filter_dmat{i,snr_col}), ...
                          num2str(filter_dmat{i,isolation_col}), ...
                          num2str(filter_dmat{i,noise_overlap_col}) ...
                          );
    end
end
fclose(filt_sum_fid);


fprintf('construct_spikeInfoMS -- done\n');


end


function [sessUniqueUnitID, timeStamp, jackTableUsed, extractInfoStr, waveForm, waveForm_all, waveForm_raw, waveForm_raw_all, waveForm_sort, waveForm_sort_all, metrics ] = filter_units(loaded_chan, noise_overlap_max, snr_min, removeLargeAmpUnits, used_jacksheet)
    
    if removeLargeAmpUnits 
        amp_thresh = 5000;
    else
        amp_thresh = Inf;
    end

    clip_timeMS_leftBound = -2;
    clip_timeMS_rightBound = 5;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define outputs 
    
    sessUniqueUnitID_varNames = {'PhysChanNum' 'UnitNum' 'DeviceNum' 'CombinedNum' 'ChanUnitName' 'NSxChanName' 'NSxFileName' 'NSPsuffix' 'ChanNameNew'};
    sessUniqueUnitID = table([], [], [], [], [], [], [], [], [], 'VariableNames', sessUniqueUnitID_varNames);

    jackTableUsed = table([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[], 'VariableNames', used_jacksheet.Properties.VariableNames);

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
    metrics.mountainsort.isolation = [];
    metrics.mountainsort.noise_overlap = [];
    metrics.mountainsort.peak_amp = [];
    metrics.mountainsort.peak_noise = [];
    metrics.mountainsort.peak_snr = [];
    metrics.mountainsort.pair_overlap = [];
    
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
        
        chan_SampFreq = used_jacksheet{iChan,'SampFreq'};
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
                
                PhysChanNum = used_jacksheet{iChan, 'PhysicalChan'};
                UnitNum = num_channel_units_added;
                DeviceNum = used_jacksheet{iChan, 'MicroDevNum'};
                CombinedNum = (DeviceNum * 1e6) + (PhysChanNum * 1e3) + (UnitNum); 
                ChanUnitName = [ used_jacksheet{iChan, 'ChanNameNew'}{1} '_' sprintf('%03d', UnitNum)];
                NSxChanName = used_jacksheet{iChan, 'ChanName'}{1};
                NSxFileName = used_jacksheet{iChan, 'FileName'}{1};
                NSPsuffix = used_jacksheet{iChan, 'NSPsuffix'}{1};
                ChanNameNew = used_jacksheet{iChan, 'ChanNameNew'}{1};
                
                sessUniqueUnitID = [ sessUniqueUnitID;  table(PhysChanNum, UnitNum, DeviceNum, CombinedNum, {ChanUnitName}, {NSxChanName}, {NSxFileName}, {NSPsuffix}, {ChanNameNew}, 'VariableNames', sessUniqueUnitID_varNames)];
                
                
                % check if this channel needs to be added to jackTableUsed
                
                if first_channel_unit
                   jackTableUsed = [ jackTableUsed ; used_jacksheet(iChan,:) ];
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



                % remember this unit passed the filter
                
                selected_units_mountainsort_ID(num_channel_units_added) = unit_firings_name;
                selected_units_combinedNum(num_channel_units_added) = CombinedNum;
                selected_units_total_count_inx(num_channel_units_added) = num_total_units_added;
                
                
            end
            
        end
        
        if num_channel_units_added > 0
        
            % store extract info string after all units on this channel have been filtered
            extractInfoStr{length(extractInfoStr) + 1, 1} = [ used_jacksheet{iChan, 'ChanName'}{1} ' --> ' num2str(num_channel_units_added) ' units, ' num2str(num_channel_unit_spikes) ' total spikes' ];

        end
        
        
        % add in pairwise isolation scores for units passing the filter
        
        if num_channel_units_added > 1
           keyboard; 
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



