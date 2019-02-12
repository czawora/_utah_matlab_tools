function variance_and_lineNoise_exclusion(varargin)

    fprintf('****************************************************\n');
    fprintf('* %s\n', mfilename);

    p = inputParser;

    p.addParameter('subj_str', '', @ischar);
    p.addParameter('time_str', '', @ischar);
    p.addParameter('nsp_str', '', @ischar);
    
    p.addParameter('split_path', '', @ischar);
    p.addParameter('completed_channel_string', 'filt', @ischar);
    p.addParameter('elementInfo_fpath', '', @ischar);
    p.addParameter('var_amp_sigma', '3', @ischar);
    
    p.addParameter('ns3_pulse_fpath', '', @ischar);
    p.addParameter('nev_fpath', '', @ischar);
    p.addParameter('nsx_physio_fpath', '', @ischar);

    parse(p, varargin{:});

    disp(p.Results);
    
    subj_str = p.Results.subj_str;
    time_str = p.Results.time_str;
    nsp_str = p.Results.nsp_str;
 
    split_path = p.Results.split_path;
    completed_channel_string = p.Results.completed_channel_string;
    var_amp_sigma = eval(p.Results.var_amp_sigma);
    
    elementInfo_fpath = p.Results.elementInfo_fpath;
    
    ns3_pulse_fpath = p.Results.ns3_pulse_fpath;
    nev_fpath = p.Results.nev_fpath;
    nsx_physio_fpath = p.Results.nsx_physio_fpath;
    
    nsx_physio_fpath_splits = strsplit(nsx_physio_fpath, '/');
    nsx_fname = nsx_physio_fpath_splits{end};
    
    session_name = [subj_str '_' time_str '_' nsp_str];
    sessStr = [subj_str '_' time_str '_' nsp_str];

    startTime_mstime = datenum(time_str, 'yymmdd_HHMM');
    
    load([split_path '/../nsx_postProc.mat']);
    physio_nsx_postProc = nsx_postProc;

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
    
    channel_fnames = sort(channel_fnames);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %read in the channel data
    
    temp_struct = load([split_path '/' channel_fnames{1}]);
    
    num_timepoints = length(temp_struct.noreref);
    samplingFreq = temp_struct.samplingFreq;
    
    noreref = zeros( length(channel_fnames), num_timepoints);
    processed = zeros( length(channel_fnames), num_timepoints);
    
    channel_ranges = zeros(1, length(channel_fnames));
    chan_names = cell(1, length(channel_fnames));
    chan_inds = zeros(1, length(channel_fnames));
    
    for iFile = 1:length(channel_fnames)
        
        temp_struct = load([split_path '/' channel_fnames{iFile}]);
        
        fprintf('%s\n', [split_path '/' channel_fnames{iFile}]);
        noreref(iFile, :) = temp_struct.noreref;
        processed(iFile, :) = temp_struct.detrended_noLineNoise;
        channel_ranges(iFile) = temp_struct.channel_range;
        chan_names{iFile} = temp_struct.channel_name;
        chan_inds(iFile) = iFile;
        
        clear('temp_struct');
    
    end
        
    save_dir = [split_path '/../outputs/'];
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %create channel name sets to filter channels into rereferencing groups
    
    chan_name_length = length(chan_names{1});
    
    referencing_info = readtextcsv(elementInfo_fpath);
    
    chan_names_by_refset = cell(size(referencing_info, 1), 1);
    
    for iRef = 1:size(referencing_info, 1)
       
        chan_string_length = length(referencing_info{iRef, 3});
        chan_digit_length = chan_name_length - chan_string_length;
        
        digit_replace_string = ['%0' num2str(chan_digit_length) 'd'];
        
        chan_names_by_refset{iRef} = strcat(referencing_info{iRef, 3} , cellfun(@(n) sprintf(digit_replace_string, n), num2cell(eval(referencing_info{iRef, 4})), 'UniformOutput', 0));
    end
    
    
    % sort channels into reference sets for processing
    chan_inds_by_refset = cell(length(chan_names_by_refset), 1);
    
    for iChan = 1:length(chan_names)
        
       current_channel = chan_names{iChan};
       
       for iRef = 1:length(chan_names_by_refset)
          
           current_refset_chan_names = chan_names_by_refset{iRef};
           
           if any(cellfun(@(x) isequal(x, current_channel), current_refset_chan_names))
              
               chan_inds_by_refset{iRef} = [ chan_inds_by_refset{iRef} chan_inds(iChan) ]; 
           end
       end
       
    end

    % seperate the timeseries into reference channel groups
    ts_by_refset = cell(length(chan_names_by_refset), 1);

    for iRef = 1:length(chan_inds_by_refset)

        ts_by_refset{iRef} = processed(chan_inds_by_refset{iRef}, :);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %iteratively calculate variance and amplitude of signal to exclude
    %outliers from global average
    
    bad_chan_names = {};
    
    mvar_ref_Z_last = [];
    mamp_ref_Z_last = [];
    line_rel_amp = [];
    
    glob_sig_all = zeros(length(ts_by_refset), size(processed, 2));
    glob_sig_good = zeros(length(ts_by_refset), size(processed, 2));

    for iRef = 1:length(ts_by_refset)
        
        current_processed = ts_by_refset{iRef};
        current_chan_inds = 1:size(current_processed, 1);
        current_chan_names = chan_names_by_refset{iRef};
        
        refset_savedir = [save_dir sprintf('refset%d', iRef) '/'];
        
        if ~exist(refset_savedir, 'dir')
           mkdir(refset_savedir); 
        end
        
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

        current_bad_chan_names = union(current_bad_chans_new, current_bad_chans_old);
        if isempty(current_bad_chan_names), current_bad_chan_names = {}; end
       
        bad_chan_names = [ bad_chan_names ; current_bad_chan_names];
       
        glob_sig_all(iRef, :) = current_glob_sig_all;
        glob_sig_good(iRef, :) = current_glob_sig_good;
        
        mvar_ref_Z_last = [mvar_ref_Z_last current_mvar_ref_Z_last];
        mamp_ref_Z_last = [mamp_ref_Z_last current_mamp_ref_Z_last];
        line_rel_amp = [line_rel_amp current_line_rel_amp];
    
    end
                                            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save out the lfpInfo
    
    % write variance.csv file
    ndec = 4;
    variance = struct();
    variance.chanName = chan_names(:);
    variance.is_good  = double(~ismember(chan_names', bad_chan_names));
    variance.variance = round(mvar_ref_Z_last(:), ndec, 'decimals');
    variance.amplitude = round(mamp_ref_Z_last(:), ndec, 'decimals');
    variance.line_noise = round(line_rel_amp(:), ndec, 'decimals');
    var_table = struct2table(variance);
    filename = [save_dir 'variance.csv'];
    writetable(var_table,filename);

    %get the pulses
    fprintf('getting the sync pulses from ns3: %s and nev: %s\n', ns3_pulse_fpath, nev_fpath);
    pulses = getBlackrockPulses_DC_AN('ns3_fpath', ns3_pulse_fpath, 'nev_fpath', nev_fpath);
    
    alignmentChan = '';
    alignedTo = '';
    
    gain_bin2uV = 0.25;
    
    lfp = int16(noreref);
    
    about_me_noreref = ['this noreref file contains the following fields:' newline newline ...
                      '     sessStr - the session name the spikes were extracted from' newline ...
                      '     startTime_mstime - the start time of session in mstime (note that these values are taken from the timestamp in the nsx filename, not the incorrectly offset time values in the original nsx file)'  newline ...
                      '     lfp - the chan x timepoint matrix of recorded micro data (noreref data has only been low-pass filtered at 500 Hz and downsampled to 1 kHz' newline ...
                      '     chan_names - a cell array of channel names corresponding to rows in the lfp matrix' newline ...
                      '     chan_names_by_refset - a cell array containing groups of channel names from the same array' newline ...
                      '     samplingFreq - the sampling frequency of the downsamples lfp' newline ...
                      '     referencing_info - a cell array containing: the referencing set name, subject name, channel string, channel num set, array number, array location number ( temporal = 1, parietal = 2, frontal = 9 ), and an array description' newline ...
                      '     pulses - a struct containing 30kHz uptimes for all "ain" and "din" channels, as well as 1kHz downsampled "ain" timeseries' newline ...
                      '     alignedTo - a string indicating which file the spikes have been aligned to' newline ...
                      '     alignmentChan - a string indicating which channel in the pulses struct was used for alignment' newline ...
                      '     physio_nsx_postProc - contains info on how the original nsx file was edited' newline ...
                      '     nsx_fname - the original nsx filename'];


    %write out data mat files
    save([save_dir session_name '_noreref.mat'],  '-v7.3', 'about_me_noreref', 'sessStr', 'startTime_mstime', 'lfp', 'pulses', 'chan_names', 'chan_names_by_refset', 'channel_ranges', 'samplingFreq', 'referencing_info', 'alignmentChan', 'alignedTo', 'physio_nsx_postProc', 'nsx_fname', 'gain_bin2uV');

    lfp = int16(processed);
    glob_sig_all = int16(glob_sig_all);
    glob_sig_good = int16(glob_sig_good);
    
    about_me_processed = ['this processed file contains the following fields:' newline newline ...
                  '     sessStr - the session name the spikes were extracted from' newline ...
                  '     startTime_mstime - the start time of session in mstime (note that these values are taken from the timestamp in the nsx filename, not the incorrectly offset time values in the original nsx file)'  newline ...
                  '     lfp - the chan x timepoint matrix of recorded micro data (processed data has been low-pass filtered at 500 Hz, downsampled to 1 kHz, detrended, and has had line noise removed)' newline ...
                  '     chan_names - a cell array of channel names corresponding to rows in the lfp matrix' newline ...
                  '     chan_names_by_refset - a cell array containing groups of channel names from the same array that were referenced together. Cell array indices correspond to computed global averages in "glob_sig_all" and "glob_sig_good"' newline ...
                  '     samplingFreq - the sampling frequency of the downsamples lfp' newline ...
                  '     glob_sig_all - computed per referecing group (usually by array), referencing group x timepoint matrix of the mean value at every timepoint in the lfp data computed using all the channels' newline ...
                  '     glob_sig_good - computed per referecing group (usually by array), referencing group x timepoint matrix of the mean value at every timepoint in the lfp data computed using only the channels marked as good in variance.csv' newline ...
                  '     referencing_info - a cell array containing: the referencing set name, subject name, channel string, channel num set, array number, array location number ( temporal = 1, parietal = 2, frontal = 9 ), and an array description' newline ...
                  '     pulses - a struct containing 30kHz uptimes for all "ain" and "din" channels, as well as 1kHz downsampled "ain" timeseries' newline ...
                  '     alignedTo - a string indicating which file the spikes have been aligned to' newline ...
                  '     alignmentChan - a string indicating which channel in the pulses struct was used for alignment' newline ...
                  '     physio_nsx_postProc - contains info on how the original nsx file was edited' newline ...
                  '     nsx_fname - the original nsx filename'];

    
    save([save_dir session_name '_processed.mat'],  '-v7.3', 'about_me_processed', 'sessStr', 'startTime_mstime', 'lfp', 'glob_sig_all', 'glob_sig_good',  'pulses', 'chan_names', 'chan_names_by_refset', 'channel_ranges', 'samplingFreq', 'referencing_info', 'alignmentChan', 'alignedTo', 'physio_nsx_postProc', 'nsx_fname', 'gain_bin2uV');

    rmdir(split_path, 's');                        

    
end
    
