function construct_spikeInfoMS(varargin)

echo construct_spikeInfoMS on;
 %NIH034_151001_1317_utah2_u
% construct_spikeInfoMS_edits('saveRoot','/Volumes/CZ_FRNU/spike_test/NIH034/NIH034_151001_1317_utah2/spike/outputs', 'nev_fpath', '/Volumes/CZ_FRNU/spike_test/NIH034/NIH034_151001_1317_utah2/NIH034_151001_1317_utah2.nev', 'elementInfo_fpath', '/Volumes/CZ_FRNU/spike_test/NIH034/NIH034_151001_1317_utah2/NIH034_151001_1317_utah2_elementInfo.txt', 'nsx_physio_fpath', '/Volumes/CZ_FRNU/spike_test/NIH034/NIH034_151001_1317_utah2/NIH034_151001_1317_utah2.ns6', 'bp_fname_suffix', 'mda_chan', 'subj_str', 'NIH034', 'time_str', '151001_1317', 'nsp_str', 'utah2', 'sessRoot', '/Volumes/CZ_FRNU/spike_test/NIH034/NIH034_151001_1317_utah2/spike/splits');

fprintf('****************************************************\n');
fprintf('* %s\n', mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse inputs

p = inputParser;

p.addParameter('subj_str', '', @ischar);
p.addParameter('time_str', '', @ischar);
p.addParameter('nsp_str', '', @ischar);

p.addParameter('sessRoot', '', @ischar);
p.addParameter('bp_fname_suffix', '', @ischar);
p.addParameter('nsx_physio_fpath', '', @ischar);
p.addParameter('ns3_pulse_fpath', '', @ischar);
p.addParameter('elementInfo_fpath', '', @ischar);
p.addParameter('nev_fpath', '', @ischar);
p.addParameter('saveRoot', '', @ischar);
p.addParameter('removeLargeAmpUnits', '0', @ischar);

p.addParameter('firings_fname', 'firings.mda', @ischar);
p.addParameter('metrics_fname', 'metrics.json', @ischar);
p.addParameter('isol_metrics_fname', 'isol_metrics.json', @ischar);
p.addParameter('isol_pair_metrics_fname', 'isol_pair_metrics.json', @ischar);
p.addParameter('clips_fname', 'clips.mda', @ischar);

parse(p, varargin{:});

disp(p.Results);

subj_str = p.Results.subj_str;
time_str = p.Results.time_str;
nsp_str = p.Results.nsp_str;

sessRoot = p.Results.sessRoot;
bp_fname_suffix = p.Results.bp_fname_suffix;
nsx_physio_fpath = p.Results.nsx_physio_fpath;
ns3_pulse_fpath = p.Results.ns3_pulse_fpath;
elementInfo_fpath = p.Results.elementInfo_fpath;
nev_fpath = p.Results.nev_fpath;
saveRoot = p.Results.saveRoot;
removeLargeAmpUnits = eval(p.Results.removeLargeAmpUnits);

firings_fname = p.Results.firings_fname;
metrics_fname = p.Results.metrics_fname;
isol_metrics_fname = p.Results.isol_metrics_fname;
isol_pair_metrics_fname = p.Results.isol_pair_metrics_fname;
clips_fname = p.Results.clips_fname;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get session start time

nsx = openNSx( nsx_physio_fpath , 'noread' );
nsxSamplingFreq = nsx.MetaTags.SamplingFreq;


nsxDateTime_mstime = datenum(time_str, 'yymmdd_HHMM');

sess = [subj_str '_' time_str '_' nsp_str];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the elementInfo

referencing_info = readtextcsv(elementInfo_fpath);
num_refsets = size(referencing_info, 1);

% read in the used channels
used_chan_dir = dir([sessRoot '/..']);
used_chan_dir_fnames = { used_chan_dir.name };

used_chan = {};
used_chan_by_refset = cell(num_refsets, 1);

for iRefset = 1:num_refsets
    
    used_chans_fname_string = sprintf('refset%d_used_chans.txt', iRefset);
    fprintf('used_chans_fname_string: %s\n', used_chans_fname_string);
    
    for iFile = 1:length(used_chan_dir_fnames)
        
       fprintf('looking for used_chans: %s\n', used_chan_dir_fnames{iFile});
       if contains(used_chan_dir_fnames{iFile}, used_chans_fname_string)
          
           fprintf('choosing: %s\n', used_chan_dir_fnames{iFile});
           used_chan_fpath = [sessRoot '/../' used_chan_dir_fnames{iFile}];
           
           used_chan_fid = fopen(used_chan_fpath);
           current_used_chan = textscan(used_chan_fid,'%s\n');
           current_used_chan = current_used_chan{1};
           fclose(used_chan_fid);
           
           used_chan = [ used_chan ; current_used_chan ];
           used_chan_by_refset{iRefset} = current_used_chan;
       end
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 1 - get list of session channels in this dir

channel_ls = dir(sessRoot);

%how many session channels are there?
% is it a directory -- AND -- does it have a firings file?
filtered_channel_ls = cellfun( @(f) ...
                                length(f) >= 3 ... %make sure its not a hidden file
                                && exist([sessRoot '/' f], 'dir') ... %make sure its a directory
                                && exist([sessRoot '/' f '/' firings_fname], 'file') ... %make sure there is a ms output file present, although this does not gauruntee that ms ran successfully
                                , {channel_ls.name});
               
%sort channel directories                            
channel_dirs = sort({channel_ls(filtered_channel_ls).name});
num_channels = sum( filtered_channel_ls );


if length(used_chan) ~= num_channels
   error('length(used_chan) ~= num_channels --- the number of split files with results does not equal the number of channels in all used_chan.txt');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make cell array for holding each channel's sort outputs filenames

firings_col = 1;
clips_col = 2;
metrics_col = 3;
isol_metrics_col = 4;
isol_pair_metrics_col = 5;
bp_col = 6;

fields_per_channel = 6;

data_cmat = cell( num_channels , fields_per_channel);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%collect the number of spikes across all channels as we load in data

channelTotals = zeros(num_channels, 1);
num_spikes = 0;
num_unique_units = 0;

if num_channels == 0
    error('0 channels detected, make sure input directory is correct');
end

for i = 1:num_channels
    
    channel_dname = channel_dirs{i};
    channel_path = [sessRoot '/' channel_dname];
    
    fprintf('loading %s\n', channel_dname);
    
    bp_pathf_dir = dir([channel_path '/*' bp_fname_suffix]);
    bp_pathf = [ bp_pathf_dir.folder '/' bp_pathf_dir.name ];
    
    fprintf('%s\n', bp_pathf);
    
    firings_pathf = [channel_path '/' firings_fname];
    metrics_pathf = [channel_path '/' metrics_fname];
    isol_metrics_pathf = [channel_path '/' isol_metrics_fname];
    isol_pair_metrics_pathf = [channel_path '/' isol_pair_metrics_fname];
    clips_pathf = [channel_path '/' clips_fname];
    
    data_cmat{i, bp_col} = bp_pathf;
    data_cmat{i, firings_col} = readmda(firings_pathf);
    data_cmat{i, clips_col} = squeeze(readmda(clips_pathf));
    
    metrics_struct = readjson(metrics_pathf);
    data_cmat{i, metrics_col} = metrics_struct.clusters;
    
    isol_metrics_struct = readjson(isol_metrics_pathf);
    data_cmat{i, isol_metrics_col} = isol_metrics_struct.clusters;
    
    isol_pair_metrics_struct = readjson(isol_pair_metrics_pathf);
    data_cmat{i, isol_pair_metrics_col} = isol_pair_metrics_struct.cluster_pairs;
    
    %record some details
    num_spikes = num_spikes + size(data_cmat{i, firings_col}, 2);
    num_unique_units = num_unique_units + length(unique(data_cmat{i, firings_col}(3,:)));
    channelTotals(i) = length(unique(data_cmat{i, firings_col}(3,:)));
    
    clip_size = size( data_cmat{i, clips_col} , 1);
    
end

fprintf('%d channels in session %s\n', size(data_cmat, 1), sess);
fprintf('total spikes across all channels for session %s is %d\n', sess, num_spikes);
fprintf('number of unique units across all channels for session %s is %d\n', sess, num_unique_units);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now we can allocate the spikeInfo struct

spikeInfo                   = struct;
spikeInfo.timeStamp         = cell(num_unique_units, 1);
spikeInfo.waveForm          = cell(num_unique_units, 1);
spikeInfo.uniqueUnitID      = cell(num_unique_units, 1);

spikeInfo.metrics                   = struct;
spikeInfo.metrics.dur_sec           = zeros(num_unique_units, 1);
spikeInfo.metrics.firing_rate       = zeros(num_unique_units, 1);
spikeInfo.metrics.num_events        = zeros(num_unique_units, 1);
spikeInfo.metrics.t1_sec            = zeros(num_unique_units, 1);
spikeInfo.metrics.t2_sec            = zeros(num_unique_units, 1);
spikeInfo.metrics.bursting_parent   = cell(num_unique_units, 1);
spikeInfo.metrics.isolation         = zeros(num_unique_units, 1);
spikeInfo.metrics.noise_overlap     = zeros(num_unique_units, 1);
spikeInfo.metrics.peak_amp          = zeros(num_unique_units, 1);
spikeInfo.metrics.peak_noise        = zeros(num_unique_units, 1);
spikeInfo.metrics.peak_snr          = zeros(num_unique_units, 1);
spikeInfo.metrics.pair_overlap      = cell(num_unique_units, 1);

%- cell array saved for each channel
spikeInfo.extractInfoStr    = {};

spikeInfo.sessDurSec        = 0;
spikeInfo.sessStr           = sess;
spikeInfo.alignedTo         = '';
spikeInfo.alignmentChan     = '';
spikeInfo.startTime_mstime = nsxDateTime_mstime;


              
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

function sessUniqueUnitID = expandUnitIDInfo(uniqueUnitID, referencing_info, used_chan_by_refset)


    sessUniqueUnitID = cell(length(uniqueUnitID), 6);
    
    for iUnit = 1:length(uniqueUnitID)
       
        current_unit = uniqueUnitID{iUnit};
        
        current_unit_splits = strsplit(current_unit, '_');
        
        nsx_chan_name = current_unit_splits{1};
        
        %get the array location and array num
        %which reference set is it a part of?
        
        refset_array_num = 0;
        refset_array_location = 0;
        
        for iRefset = 1:length(used_chan_by_refset)
            
            if any(cellfun(@(x) isequal(x, nsx_chan_name), used_chan_by_refset{iRefset}))
                % match!
            
                refset_array_num = str2num(referencing_info{iRefset, 5});
                refset_array_location = str2num(referencing_info{iRefset, 6});
            end
        end
        
        
        fprintf('refset_array_num: %d\n', refset_array_num);
        fprintf('refset_array_location: %d\n', refset_array_location);
    
        
        %get the unit num
        current_unit_num = str2num(current_unit_splits{2});
        
        %get the channelnum
        [startIndex,endIndex] = regexp(nsx_chan_name,'\d+');
        current_chan_num = str2num(nsx_chan_name(startIndex:endIndex));
        
        combined_unit_num = (refset_array_location * 1e6) + (current_chan_num * 1e3) + current_unit_num;
        
        sessUniqueUnitID{iUnit, 1} = current_chan_num;
        sessUniqueUnitID{iUnit, 2} = current_unit_num;
        sessUniqueUnitID{iUnit, 3} = refset_array_location;
        sessUniqueUnitID{iUnit, 4} = combined_unit_num;
        sessUniqueUnitID{iUnit, 5} = refset_array_num;
        sessUniqueUnitID{iUnit, 6} = current_unit;
      
    end

end

function [spikeInfo_filt, noiseInfo, filter_dmat] = filterSpikes(spikeInfo, noise_overlap_max, snr_min, isolation_min, removeLargeAmpUnits)

    orig_unitID_col = 1;
    filt_unitID_col = 2;
    filter_type_col = 3;
    unit_channel_col = 4;
    num_spikes_col = 5;
    snr_col = 6;
    isolation_col = 7;
    noise_overlap_col = 8;
    
    
    filter_dmat = cell(length(spikeInfo.uniqueUnitID), 8);
    %first column is original unit name
    % then filtered unit name
    % filter type
    % channel of unit
    % spikes for unit
    
    if removeLargeAmpUnits == 1

        largeAmp = 5000;
    else

        largeAmp = Inf;
    end
    
    current_renaming_chan = '';
    noise_spikes = 0;
    noise_units = 0;
    
    num_unique_units = size(spikeInfo.uniqueUnitID, 1);
    
    %label units as noise, non-isolated, and unit
    for iUnit = 1:length(spikeInfo.uniqueUnitID)
        
        current_unitID = spikeInfo.uniqueUnitID{iUnit};  
                
        current_unitID_splits = strsplit(current_unitID, '_');
        current_chan = current_unitID_splits{1};
        %current_ID_num = current_unitID_splits{2};
        
        filter_dmat{iUnit, orig_unitID_col} = current_unitID;
        filter_dmat{iUnit, unit_channel_col} = current_chan;
        filter_dmat{iUnit, num_spikes_col} = length(spikeInfo.timeStamp{iUnit});
        
        filter_dmat{iUnit, snr_col} = spikeInfo.metrics.peak_snr(iUnit);
        filter_dmat{iUnit, isolation_col} = spikeInfo.metrics.isolation(iUnit);
        filter_dmat{iUnit, noise_overlap_col} = spikeInfo.metrics.noise_overlap(iUnit);
                
        %determine unit type
        if ~(spikeInfo.metrics.noise_overlap(iUnit) < noise_overlap_max ...
                && spikeInfo.metrics.peak_snr(iUnit) > snr_min ...
                && spikeInfo.metrics.peak_amp(iUnit) < largeAmp)
            %this is noise unit
            
            filter_dmat{iUnit, filter_type_col} = 'noise';
            
            noise_spikes = noise_spikes + length(spikeInfo.timeStamp{iUnit});
            noise_units = noise_units + 1;
            
        elseif ~(spikeInfo.metrics.isolation(iUnit) > isolation_min)
            
            filter_dmat{iUnit, filter_type_col} = 'non-isolated';
        else
            filter_dmat{iUnit, filter_type_col} = 'unit';
        end
        
        if ~isequal(current_chan, current_renaming_chan)
            
            noise_unit_count = 1;
            preserved_unit_count = 1;
            current_renaming_chan = current_chan;
        end
        
        if isequal(filter_dmat{iUnit, filter_type_col}, 'noise')
            
            %noise_unit_ID = [current_chan '_nz0' num2str(noise_unit_count)];
            noise_unit_ID = [current_chan '_nz' sprintf('%03d', noise_unit_count) ];
            filter_dmat{iUnit, filt_unitID_col} = noise_unit_ID;
            
            noise_unit_count = noise_unit_count + 1;
        else
            
            %filt_unit_ID = [current_chan '_0' num2str(preserved_unit_count)];
            filt_unit_ID = [current_chan '_' sprintf('%03d',preserved_unit_count) ];
            filter_dmat{iUnit, filt_unitID_col} = filt_unit_ID;
            
            preserved_unit_count = preserved_unit_count + 1;
        end
        
    end
    
    
    %update extractInfoStrs
    
    infoStrs = spikeInfo.extractInfoStr;
    
    for iStr = 1:length(infoStrs)
        
        infoStrs_splits = strsplit(infoStrs{iStr}, ' --> ');
        full_channel_name = infoStrs_splits{1};
        full_channel_name_splits = strsplit(full_channel_name, '_');
        channel_name = full_channel_name_splits{end};
        
        per_channel_noise_units = 0;
        per_channel_noise_spikes = 0;
        
        per_channel_non_noise_units = 0;
        per_channel_non_noise_spikes = 0;
        
        for jUnit = 1:size(filter_dmat,1)
            
            if isequal(filter_dmat{jUnit, unit_channel_col} , channel_name)
                
                if isequal(filter_dmat{jUnit, filter_type_col} , 'noise')
                    
                    per_channel_noise_units = per_channel_noise_units + 1;
                    per_channel_noise_spikes = per_channel_noise_spikes + filter_dmat{jUnit, num_spikes_col};
                else
                    
                    per_channel_non_noise_units = per_channel_non_noise_units + 1;
                    per_channel_non_noise_spikes = per_channel_non_noise_spikes + filter_dmat{jUnit, num_spikes_col};
                end
            end
            
        end
        
        new_infoStr = [infoStrs{iStr} ' --> filter --> ' num2str(per_channel_noise_units) ' noise units, ' num2str(per_channel_noise_spikes) ' spikes; ' num2str(per_channel_non_noise_units) ' units, ' num2str(per_channel_non_noise_spikes) ' spikes'];
        infoStrs{iStr} = new_infoStr;
    end
    
    disp(filter_dmat);
    
    non_noise_units = num_unique_units - noise_units;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %create noise info struct
    
    noiseInfo                   = struct;
    noiseInfo.timeStamp         = cell(noise_units, 1);
    noiseInfo.waveForm          = cell(noise_units, 1);
    noiseInfo.uniqueUnitID      = cell(noise_units, 1);
    
    noiseInfo.metrics                   = struct;
    noiseInfo.metrics.dur_sec           = zeros(noise_units, 1);
    noiseInfo.metrics.firing_rate       = zeros(noise_units, 1);
    noiseInfo.metrics.num_events        = zeros(noise_units, 1);
    noiseInfo.metrics.t1_sec            = zeros(noise_units, 1);
    noiseInfo.metrics.t2_sec            = zeros(noise_units, 1);
    noiseInfo.metrics.bursting_parent   = cell(noise_units, 1);
    noiseInfo.metrics.isolation         = zeros(noise_units, 1);
    noiseInfo.metrics.noise_overlap     = zeros(noise_units, 1);
    noiseInfo.metrics.peak_amp          = zeros(noise_units, 1);
    noiseInfo.metrics.peak_noise        = zeros(noise_units, 1);
    noiseInfo.metrics.peak_snr          = zeros(noise_units, 1);
    noiseInfo.metrics.pair_overlap      = cell(noise_units, 1);
    
    noiseInfo.extractInfoStr    = infoStrs;
    noiseInfo.sessDurSec        = spikeInfo.sessDurSec;
    noiseInfo.sessStr           = spikeInfo.sessStr;
    noiseInfo.alignedTo         = '';
    noiseInfo.alignmentChan     = '';    
    noiseInfo.startTime_mstime = spikeInfo.startTime_mstime;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %create filtered spike info struct
    
    spikeInfo_filt                   = struct;
    spikeInfo_filt.timeStamp         = cell(non_noise_units, 1);
    spikeInfo_filt.waveForm          = cell(non_noise_units, 1);
    spikeInfo_filt.uniqueUnitID      = cell(non_noise_units, 1);
    
    spikeInfo_filt.metrics                   = struct;
    spikeInfo_filt.metrics.dur_sec           = zeros(non_noise_units, 1);
    spikeInfo_filt.metrics.firing_rate       = zeros(non_noise_units, 1);
    spikeInfo_filt.metrics.num_events        = zeros(non_noise_units, 1);
    spikeInfo_filt.metrics.t1_sec            = zeros(non_noise_units, 1);
    spikeInfo_filt.metrics.t2_sec            = zeros(non_noise_units, 1);
    spikeInfo_filt.metrics.bursting_parent   = cell(non_noise_units, 1);
    spikeInfo_filt.metrics.isolation         = zeros(non_noise_units, 1);
    spikeInfo_filt.metrics.noise_overlap     = zeros(non_noise_units, 1);
    spikeInfo_filt.metrics.peak_amp          = zeros(non_noise_units, 1);
    spikeInfo_filt.metrics.peak_noise        = zeros(non_noise_units, 1);
    spikeInfo_filt.metrics.peak_snr          = zeros(non_noise_units, 1);
    spikeInfo_filt.metrics.pair_overlap      = cell(non_noise_units, 1);
    
    spikeInfo_filt.extractInfoStr    = infoStrs;
    spikeInfo_filt.sessDurSec        = spikeInfo.sessDurSec;    
    spikeInfo_filt.sessStr           = spikeInfo.sessStr;
    spikeInfo_filt.alignedTo         = '';        
    spikeInfo_filt.alignmentChan     = '';    
    spikeInfo_filt.startTime_mstime  = spikeInfo.startTime_mstime;
    
    %%%%%%%%
    % seperate the spikeInfo
    
    noise_units_added = 1;
    filt_units_added = 1;
    
    for iUnit = 1:length(spikeInfo.uniqueUnitID)
        
        current_unit_filt_type = filter_dmat{iUnit, filter_type_col};
        fprintf('filtering unit %s\n', spikeInfo.uniqueUnitID{iUnit});
        
        if isequal(current_unit_filt_type, 'noise')
            
            %get the new unit name
            noise_unit_ID = filter_dmat{iUnit, filt_unitID_col};
            
            noiseInfo.timeStamp{noise_units_added} = spikeInfo.timeStamp{iUnit};
            noiseInfo.waveForm{noise_units_added}  = spikeInfo.waveForm{noise_units_added};
            
            noiseInfo.uniqueUnitID{noise_units_added}              = noise_unit_ID;
            
            noiseInfo.metrics.dur_sec(noise_units_added)           = spikeInfo.metrics.dur_sec(iUnit);
            noiseInfo.metrics.firing_rate(noise_units_added)       = spikeInfo.metrics.firing_rate(iUnit);
            noiseInfo.metrics.num_events(noise_units_added)        = spikeInfo.metrics.num_events(iUnit);
            noiseInfo.metrics.t1_sec(noise_units_added)            = spikeInfo.metrics.t1_sec(iUnit);
            noiseInfo.metrics.t2_sec(noise_units_added)            = spikeInfo.metrics.t2_sec(iUnit);
            
            %find the new name for the bursting parent
            %should return one index
            if ~isequal(spikeInfo.metrics.bursting_parent{iUnit}, 0)
                burst_parent_idx = cellfun(@(x) isequal(x, spikeInfo.metrics.bursting_parent{iUnit}), spikeInfo.uniqueUnitID);
                noiseInfo.metrics.bursting_parent{noise_units_added}   = filter_dmat{burst_parent_idx, filt_unitID_col};
            else
                noiseInfo.metrics.bursting_parent{noise_units_added} = 0;
            end
            
            noiseInfo.metrics.isolation(noise_units_added)         = spikeInfo.metrics.isolation(iUnit);
            noiseInfo.metrics.noise_overlap(noise_units_added)     = spikeInfo.metrics.noise_overlap(iUnit);
            noiseInfo.metrics.peak_amp(noise_units_added)          = spikeInfo.metrics.peak_amp(iUnit);
            noiseInfo.metrics.peak_noise(noise_units_added)        = spikeInfo.metrics.peak_noise(iUnit);
            noiseInfo.metrics.peak_snr(noise_units_added)          = spikeInfo.metrics.peak_snr(iUnit);

            current_unit_pair_overlap = spikeInfo.metrics.pair_overlap{iUnit};
            
            if isempty(current_unit_pair_overlap.pair)
                
                noiseInfo.metrics.pair_overlap{noise_units_added}      = current_unit_pair_overlap;
            
            else    
                
                overlap_struct = struct;
                overlap_struct.pair = cell(length(current_unit_pair_overlap.pair),1);
                overlap_struct.overlap = current_unit_pair_overlap.overlap;
                               
                for ov = 1:length(current_unit_pair_overlap.pair)
                    
                    ov_split = strsplit(current_unit_pair_overlap.pair{ov}, ',');
                    
                    el1 = ov_split{1};
                    el2 = ov_split{2};
                    
                    
                    el1_idx = cellfun(@(x) isequal(x, el1), spikeInfo.uniqueUnitID);
                    el2_idx = cellfun(@(x) isequal(x, el2), spikeInfo.uniqueUnitID);
                    
                    renamed_el1 = filter_dmat{el1_idx, filt_unitID_col};
                    renamed_el2 = filter_dmat{el2_idx, filt_unitID_col};
                    
                    ov_unsplit = [renamed_el1 ',' renamed_el2];
                                        
                    overlap_struct.pair{ov} = ov_unsplit;
                    
                end
                
                noiseInfo.metrics.pair_overlap{noise_units_added} = overlap_struct;
                
            end
            
            noise_units_added = noise_units_added + 1;
        
        else %it is a unit then
            
            filt_unit_ID = filter_dmat{iUnit, filt_unitID_col};
                        
            spikeInfo_filt.timeStamp{filt_units_added}         = spikeInfo.timeStamp{iUnit};
            spikeInfo_filt.waveForm{filt_units_added}          = spikeInfo.waveForm{iUnit};
            
            spikeInfo_filt.uniqueUnitID{filt_units_added}             = filt_unit_ID;
            
            spikeInfo_filt.metrics.dur_sec(filt_units_added)           = spikeInfo.metrics.dur_sec(iUnit);
            spikeInfo_filt.metrics.firing_rate(filt_units_added)       = spikeInfo.metrics.firing_rate(iUnit);
            spikeInfo_filt.metrics.num_events(filt_units_added)        = spikeInfo.metrics.num_events(iUnit);
            spikeInfo_filt.metrics.t1_sec(filt_units_added)            = spikeInfo.metrics.t1_sec(iUnit);
            spikeInfo_filt.metrics.t2_sec(filt_units_added)            = spikeInfo.metrics.t2_sec(iUnit);
            
            %find the new name for the bursting parent
            %should return one index
            if ~isequal(spikeInfo.metrics.bursting_parent{iUnit}, 0)
                burst_parent_idx = cellfun(@(x) isequal(x, spikeInfo.metrics.bursting_parent{iUnit}), spikeInfo.uniqueUnitID);
                spikeInfo_filt.metrics.bursting_parent{filt_units_added}   = filter_dmat{burst_parent_idx, filt_unitID_col};
            else
                spikeInfo_filt.metrics.bursting_parent{filt_units_added} = 0;
            end
            
            spikeInfo_filt.metrics.isolation(filt_units_added)         = spikeInfo.metrics.isolation(iUnit);
            spikeInfo_filt.metrics.noise_overlap(filt_units_added)     = spikeInfo.metrics.noise_overlap(iUnit);
            spikeInfo_filt.metrics.peak_amp(filt_units_added)          = spikeInfo.metrics.peak_amp(iUnit);
            spikeInfo_filt.metrics.peak_noise(filt_units_added)        = spikeInfo.metrics.peak_noise(iUnit);
            spikeInfo_filt.metrics.peak_snr(filt_units_added)          = spikeInfo.metrics.peak_snr(iUnit);

            current_unit_pair_overlap = spikeInfo.metrics.pair_overlap{iUnit};
            
            if isempty(current_unit_pair_overlap.pair)
                
                spikeInfo_filt.metrics.pair_overlap{filt_units_added}      = current_unit_pair_overlap;
            else    
                
                overlap_struct = struct;
                overlap_struct.pair = cell(length(current_unit_pair_overlap.pair),1);
                overlap_struct.overlap = current_unit_pair_overlap.overlap;
                               
                for ov = 1:length(current_unit_pair_overlap.pair)
                    
                    ov_split = strsplit(current_unit_pair_overlap.pair{ov}, ',');
                    
                    el1 = ov_split{1};
                    el2 = ov_split{2};
                                        
                    el1_idx = cellfun(@(x) isequal(x, el1), spikeInfo.uniqueUnitID);
                    el2_idx = cellfun(@(x) isequal(x, el2), spikeInfo.uniqueUnitID);
                    
                    renamed_el1 = filter_dmat{el1_idx, filt_unitID_col};
                    renamed_el2 = filter_dmat{el2_idx, filt_unitID_col};
                    
                    ov_unsplit = [renamed_el1 ',' renamed_el2];
                                        
                    overlap_struct.pair{ov} = ov_unsplit;
                    
                end
                
                spikeInfo_filt.metrics.pair_overlap{filt_units_added} = overlap_struct;
                
            end
            
            
            filt_units_added = filt_units_added + 1;
        end
        
    end

end

function jsonstring = readjson(fname)

fid = fopen(fname);
s = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
jsonstring = jsondecode(strjoin(s{1}, ' '));

end



