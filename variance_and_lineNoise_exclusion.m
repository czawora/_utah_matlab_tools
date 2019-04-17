function variance_and_lineNoise_exclusion(varargin)

    fprintf('****************************************************\n');
    fprintf('* %s\n', mfilename);

    p = inputParser;
    
    p.addParameter('session_path', '', @ischar);
    p.addParameter('jacksheet_fpath', '', @ischar);
    p.addParameter('analog_pulse_fpath', '', @ischar);
    p.addParameter('nev_fpath', '', @ischar);
    
    p.addParameter('completed_channel_string', 'filt', @ischar);
    p.addParameter('var_amp_sigma', '3', @ischar);

    parse(p, varargin{:});

    disp(p.Results);
 
    session_path = p.Results.session_path;
    jacksheet_fpath = p.Results.jacksheet_fpath;
    analog_pulse_fpath = p.Results.analog_pulse_fpath;
    nev_fpath = p.Results.nev_fpath;
    
    completed_channel_string = p.Results.completed_channel_string;
    var_amp_sigma = eval(p.Results.var_amp_sigma);
        
    if ~exist(session_path, 'dir')
       error('%s is not a valid dir\n', session_path);
    end   

    if ~exist(jacksheet_fpath, 'file')
       error('%s is not a valid file\n', jacksheet_fpath);
    end
            
    split_path = [session_path '/lfp_splits'];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find all the filtered channels in the split_dir
    
    channel_fnames = {};
    
    split_files = dir(split_path);
    
    for iFile = 1:length(split_files)
        if contains(split_files(iFile).name, completed_channel_string)        
            channel_fnames{ length(channel_fnames) + 1 } = split_files(iFile).name;
        end
    end
    
    if isempty(channel_fnames)
       error('0 filtered channel mat files are present in split dir\n'); 
    end
    
    if length(channel_fnames) ~= sum(cellfun(@(x) length(x) > 3 && contains(x, 'mat'), {split_files.name}))
        error('length(channel_fnames) ~= sum(~[split_files.isdir]) : there seem to unprocessed files left in the split_dir');
    end
    
    channel_fnames = sort(channel_fnames);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read jacksheet
    
    full_jacksheet = readtable(jacksheet_fpath);
    used_jacksheet = readtable([session_path '/used_jacksheet.csv']);
   
    original_split_channels = used_jacksheet{:, 'ChanName'};
   
    if length(original_split_channels) ~= length(channel_fnames)
        error('number of channels initially split and the number just read back in does not match');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %read in the channel data
    
    temp_struct = load([split_path '/' channel_fnames{1}]);
    
    num_timepoints = length(temp_struct.noreref);
    
    noreref = zeros( length(channel_fnames), num_timepoints);
    processed = zeros( length(channel_fnames), num_timepoints);
        
    loaded_samplingFreqs = zeros(1, length(channel_fnames));
    
    for iFile = 1:length(channel_fnames)
        
        temp_struct = load([split_path '/' channel_fnames{iFile}]);
        
        temp_channel_name = temp_struct.channel_name;
        
        %which row is this channel name in the jacksheet
        match2jacksheet = cellfun( @(x) isequal(x, temp_channel_name) , original_split_channels);
        
        if ~any(match2jacksheet)
           error('current loaded channel does not match a channel in the jacksheet: %s', channel_fnames{iFile}); 
        end
        
        if sum(match2jacksheet) > 1
            error('current loaded channel matches more than one channel in the jacksheet: %s, %s', channel_fnames{iFile}, temp_channel_name); 
        end
        
        
        fprintf('%s\n', [split_path '/' channel_fnames{iFile}]);

        noreref(match2jacksheet, :) = temp_struct.noreref;
        processed(match2jacksheet, :) = temp_struct.detrended_noLineNoise;
        loaded_samplingFreqs(match2jacksheet) = temp_struct.samplingFreq;
        
        clear('temp_struct');
    
    end
        
    save_dir = [session_path '/outputs/'];
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %create channel name sets to filter channels into rereferencing groups
    
    %how many referencing groups on this nsp
    reference_groups = unique(used_jacksheet{:, 'MicroDevNum'});
    num_reference_groups = length(reference_groups);

 
    % make sure all the loaded channels were downsampled to the same value
    if length(unique(loaded_samplingFreqs)) > 1
        error('more than one downsample rate was used in the loaded channels');
    else
        samplingFreq = loaded_samplingFreqs(1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %iteratively calculate variance and amplitude of signal to exclude
    %outliers from global average     
    
    glob_sig_all = cell(1, max(reference_groups));
    glob_sig_good = cell(1, max(reference_groups));
    
    variance_chanName = used_jacksheet{:, 'ChanName'};
    variance_is_good = ones(1, length(used_jacksheet{:, 'ChanName'}));
    variance_variance = nan(1, length(used_jacksheet{:, 'ChanName'}));
    variance_amplitude = nan(1, length(used_jacksheet{:, 'ChanName'}));
    variance_line_noise = nan(1, length(used_jacksheet{:, 'ChanName'}));
   
    for iRef = 1:num_reference_groups
        
        current_reference_num = reference_groups(iRef);
        
        current_reference_match = (used_jacksheet{:, 'MicroDevNum'} == current_reference_num);
        
        current_processed = processed( current_reference_match , :);
        current_chan_inds = 1:size(current_processed, 1);
        current_chan_names = used_jacksheet{current_reference_match, 'ChanName'};
        
        refset_savedir = [save_dir sprintf('microDev%d', current_reference_num) '/'];
        
        if ~exist(refset_savedir, 'dir')
           mkdir(refset_savedir); 
        end
        
        plotNaN(current_chan_names, current_processed, samplingFreq, refset_savedir);
        
        % check if all values are NaN, meaning the entire array is flat the whole session
        if sum(isnan(current_processed(:))) == (size(current_processed, 1) * size(current_processed, 2))
        
            current_bad_chan_names = current_chan_names';
            
            % mark these channels as bad in variance struct
            variance_is_good(cellfun(@(x) any(cellfun(@(y) isequal(x, y), current_bad_chan_names)) , variance_chanName)) = 0;    
            
        else
            
            % remove any channels that are totally NaN
            
            current_processed_nanChan = ( sum(isnan(current_processed), 2) == size(current_processed, 2) );
            excluded_chan_names = current_chan_names( current_processed_nanChan );
            
            if ~isempty(excluded_chan_names)
                
                % mark these channels as bad in variance struct
                variance_is_good(cellfun(@(x) any(cellfun(@(y) isequal(x, y), excluded_chan_names)) , variance_chanName)) = 0;    

            end
            
            % quality check only non NaN channels
            
            current_processed = current_processed(~current_processed_nanChan, :);
            current_chan_inds = current_chan_inds(~current_processed_nanChan);
            current_chan_names = current_chan_names(~current_processed_nanChan);
            
            [current_bad_chans_new, ...
             current_bad_chans_old, ...
             ~, ... %eeg_raw_proc
             current_glob_sig_all, ...
             current_glob_sig_good, ...
             current_mvar_ref_Z_last, ...
             current_mamp_ref_Z_last, ...
             current_line_rel_amp, ...
             ~] = eeg_noise_metrics_CZ([], [], [], [], [], [] , ...
                                                        'biowulf', '1',...
                                                        'biowulf_eeg_raw_proc', current_processed',...
                                                        'biowulf_Fs', samplingFreq ,...
                                                        'biowulf_chan_inds', current_chan_inds',...
                                                        'biowulf_chan_names', current_chan_names',...
                                                        'outputdir', refset_savedir ,...
                                                        'z_thresh', var_amp_sigma);                                       

            % add amp and var stats to variance struct
            variance_chanName_filt = cellfun(@(x) any(cellfun(@(y) isequal(x, y), current_chan_names)) , variance_chanName);
            
            variance_variance(variance_chanName_filt) = current_mvar_ref_Z_last(:);
            variance_amplitude(variance_chanName_filt) = current_mamp_ref_Z_last(:);
            variance_line_noise(variance_chanName_filt) = current_line_rel_amp(:);
   
            % mark these channels as bad in variance struct                                        
            current_bad_chan_names = union(current_bad_chans_new, current_bad_chans_old);
            if isempty(current_bad_chan_names), current_bad_chan_names = {}; end

            variance_is_good(cellfun(@(x) any(cellfun(@(y) isequal(x, y), current_bad_chan_names)) , variance_chanName)) = 0;    

            
            glob_sig_all{current_reference_num} = int16(current_glob_sig_all);
            glob_sig_good{current_reference_num} = int16(current_glob_sig_good);
            
        end                     
    
    end
                                            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save out the lfpInfo
    
    % write variance.csv file
    ndec = 4;
    variance = struct();
    variance.chanName = variance_chanName;
    variance.is_good = variance_is_good';
    variance.variance = round( variance_variance' , ndec, 'decimals');
    variance.amplitude = round( variance_amplitude' , ndec, 'decimals');
    variance.line_noise = round( variance_line_noise' , ndec, 'decimals');
    
    var_table = struct2table(variance);
    filename = [save_dir 'variance.csv'];
    writetable(var_table,filename);

    %get the pulses
    analog_splits = strsplit(analog_pulse_fpath, '/');
    fprintf('getting the sync pulses from ns3: %s and nev: %s\n', analog_pulse_fpath, nev_fpath);
    pulses_struct = getBlackrockPulses_DC_AN('ns3_fpath', analog_pulse_fpath, 'nev_fpath', nev_fpath);
    pulses = {used_jacksheet{1, 'NSPsuffix'} analog_splits{end} pulses_struct};
    
    gain_bin2uV = 0.25;
    
    load([session_path '/nsx_postProc.mat']);
    physio_nsx_postProc = nsx_postProc;
    physio_nsx_postProc = {used_jacksheet{1, 'NSPsuffix'} used_jacksheet{1, 'FileName'} physio_nsx_postProc};
    
    %extract start time from jacksheet
    dateInfo_column = 'FileName';
    old_fmt_regex_match = regexp( used_jacksheet{1, dateInfo_column} , '\d\d\d\d\d\d_\d\d\d\d', 'match');
    old_fmt_match = old_fmt_regex_match{1};
    
    new_fmt_regex_match = regexp( used_jacksheet{1, dateInfo_column} , '\d\d\d\d\d\d\d\d-\d\d\d\d\d\d', 'match');
    new_fmt_match = new_fmt_regex_match{1};
    
    if isempty(old_fmt_match) && isempty(new_fmt_match)
        
        error('no old or new format datestring found in used_jacksheet{1, "FileName"}');
    elseif ~isempty(old_fmt_match)
        
        startTime_datenum = datenum(old_fmt_match, 'yymmdd_HHMM');
    elseif ~isempty(new_fmt_match)
        
        startTime_datenum = datenum(new_fmt_match, 'yyyymmdd-HHMMSS');
    end
        
    
    createdDate = datestr(datetime('now'), 'mmm.dd,yyyy HH:MM:SS');
    chanIDperNSP = {used_jacksheet};
    sessDurMin = used_jacksheet{1, 'DurationMin'};
    sessStr = used_jacksheet{1, 'RawDir'}{1};
        
    
    readme = ['this lfp.mat file, generated ' sprintf('%s', createdDate) ', contains the following fields:' newline newline ...
              '     createdDate          - a string indicating when this mat file was created' newline ...
              '     chanIDperNSP         - a ( #NSPs used x 1 ) cell array. Each enty is a table corresponding to the channel data in lfp' newline ...
              '     lfp                  - a ( #NSPs used x 1 ) cell array. Each cell entry is #channels x #time points matrix (int16). Time point dimension might be slightly different prior to NSP alignment' newline ...
              '     rerefType            - a string indicating either "lfp_noreref" or "lfp_processed"' newline ...
              '     gain_bin2uV          - multiplicative factor to convert lfp_noreref to uV' newline ...
              '     samplingFreq         - samlple frequency of the data. typically 1000 Hz' newline ...
              '     filterSettings       - a string describing the processing steps applied to data in lfp field' newline ...
              '     sessStr              - the session name the lfps were extracted from (e.g., 190117_1336)' newline ...
              '     sessDurMin           - session duration in minutes' newline ...
              '     pulses               - (a #nsp x 3) cell with "nspSuffix" "nsxFilename" and "pulse" struct containing 30kHz uptimes for all "ain" and "din" channels, as well as 1kHz downsampled "ain" timeseries for that NSP' newline ...
              '     alignedTo            - a string indicating which file the spikes have been aligned to' newline ...
              '     alignmentChan        - a string indicating which channel in the pulses struct was used for alignment' newline ...
              '     jackTableFull        - complete jacksheetBR table from this session (all channels) with device numbers and new channel names' newline ...
              '     jackTableUsed        - just the jacksheet for the channels incorporated into this lfpStruct (combined across NSPs). this table can be used to go back and forth between original and new channel names' newline ...
              '     startTime_datenum    - the start time of session as output from datenum function (note that these values are taken from the timestamp in the nsx filename, not the incorrectly offset time values in the original nsx file)'  newline ...
              '     physio_nsx_postProc  - (a #nsp x 3) cell with "nspSuffix" "nsxFilename" and struct with info on how the original nsx file was modified using concatOpenNSx' ];
          

    lfp = {int16(noreref)};  
    rerefType = 'lfp_noreref';
          
    
    lfpStruct = struct;
    lfpStruct.readme = readme;          
    lfpStruct.createdDate = createdDate;
    lfpStruct.chanIDperNSP = chanIDperNSP;
    lfpStruct.lfp = lfp;
    lfpStruct.rerefType = rerefType;
    lfpStruct.gain_bin2uV = gain_bin2uV;
    lfpStruct.samplingFreq = samplingFreq;
    lfpStruct.filterSettings = ['LFP Band: butterworth low pass at 500 Hz; filtfilt for zero phase shift (effective order 4)' newline ...
                                ' + downsampled to lfpStruct.samplingFreq '];
    lfpStruct.sessStr = sessStr;
    lfpStruct.sessDurMin = sessDurMin;
    lfpStruct.pulses = pulses;
    lfpStruct.alignedTo = '';
    lfpStruct.alignmentChan = '';
    lfpStruct.jackTableFull = full_jacksheet;
    lfpStruct.jackTableUsed = used_jacksheet;
    lfpStruct.startTime_datenum = startTime_datenum;
    lfpStruct.physio_nsx_postProc = physio_nsx_postProc;

    %write out data mat files
    save([save_dir sessStr '_noreref.mat'],  '-v7.3', 'lfpStruct');    
    

    lfp = {int16(processed)};
    rerefType = 'lfp_processed';

    readme = ['this lfp.mat file, generated ' sprintf('%s', createdDate) ', contains the following fields:' newline newline ...
          '     createdDate          - a string indicating when this mat file was created' newline ...
          '     chanIDperNSP         - a ( #NSPs used x 1 ) cell array. Each enty is a table corresponding to the channel data in lfp' newline ...
          '     lfp                  - a ( #NSPs used x 1 ) cell array. Each cell entry is #channels x #time points matrix (int16). Time point dimension might be slightly different prior to NSP alignment' newline ...
          '     glob_sig_all         - 1-d cell array containing a global mean using all channels from a device. Index this array using MicroDevNum from the jacksheet' newline ...
          '     glob_sig_good        - 1-d cell array containing a global mean using only channels from a device that pass an amplitude and variance quality filter (see variance.csv). Index this array using MicroDevNum from the jacksheet' newline ... 
          '     rerefType            - a string indicating either "lfp_noreref" or "lfp_processed"' newline ...
          '     gain_bin2uV          - multiplicative factor to convert lfp_noreref to uV' newline ...
          '     samplingFreq         - samlple frequency of the data. typically 1000 Hz' newline ...
          '     filterSettings       - a string describing the processing steps applied to data in lfp field' newline ...
          '     sessStr              - the session name the lfps were extracted from (e.g., 190117_1336)' newline ...
          '     sessDurMin           - session duration in minutes' newline ...
          '     pulses               - (a #nsp x 3) cell with "nspSuffix" "nsxFilename" and "pulse" struct containing 30kHz uptimes for all "ain" and "din" channels, as well as 1kHz downsampled "ain" timeseries for that NSP' newline ...
          '     alignedTo            - a string indicating which file the spikes have been aligned to' newline ...
          '     alignmentChan        - a string indicating which channel in the pulses struct was used for alignment' newline ...
          '     jackTableFull        - complete jacksheetBR table from this session (all channels) with device numbers and new channel names' newline ...
          '     jackTableUsed        - just the jacksheet for the channels incorporated into this lfpStruct (combined across NSPs). this table can be used to go back and forth between original and new channel names' newline ...
          '     startTime_datenum    - the start time of session as output from datenum function (note that these values are taken from the timestamp in the nsx filename, not the incorrectly offset time values in the original nsx file)'  newline ...
          '     physio_nsx_postProc  - (a #nsp x 3) cell with "nspSuffix" "nsxFilename" and struct with info on how the original nsx file was modified using concatOpenNSx' ];

    lfpStruct = struct;
    lfpStruct.readme = readme;          
    lfpStruct.createdDate = createdDate;
    lfpStruct.chanIDperNSP = chanIDperNSP;
    lfpStruct.lfp = lfp;
    lfpStruct.glob_sig_all = glob_sig_all;
    lfpStruct.glob_sig_good = glob_sig_good;
    lfpStruct.rerefType = rerefType;
    lfpStruct.gain_bin2uV = gain_bin2uV;
    lfpStruct.samplingFreq = samplingFreq;
    lfpStruct.filterSettings = ['LFP Band: butterworth low pass at 500 Hz; filtfilt for zero phase shift (effective order 4)' newline ...
                                ' + downsampled to lfpStruct.samplingFreq ' newline ...
                                ' + detrended with locdetrend ' newline ...
                                ' + line noise removed with rmlinesmovingwinc ' newline ...
                                ' + constant voltage > 10 ms set to NaN'];
    lfpStruct.sessStr = sessStr;
    lfpStruct.sessDurMin = sessDurMin;
    lfpStruct.pulses = pulses;
    lfpStruct.alignedTo = '';
    lfpStruct.alignmentChan = '';
    lfpStruct.jackTableFull = full_jacksheet;
    lfpStruct.jackTableUsed = used_jacksheet;
    lfpStruct.startTime_datenum = startTime_datenum;
    lfpStruct.physio_nsx_postProc = physio_nsx_postProc;          
              
    
    save([save_dir sessStr '_processed.mat'],  '-v7.3', 'lfpStruct');

    rmdir(split_path, 's');                        

    
end
    
