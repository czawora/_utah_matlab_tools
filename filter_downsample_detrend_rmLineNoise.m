function filter_downsample_detrend_rmLineNoise(varargin)
    
    fprintf('****************************************************\n');
    fprintf('* %s\n', mfilename);

    p = inputParser;

    p.addParameter('channel_fpath', '', @ischar);
    p.addParameter('max_freq', '500', @ischar);
    p.addParameter('downsample_freq', '1000', @ischar);
    
    parse(p, varargin{:});

    disp(p.Results);
    
    channel_fpath = p.Results.channel_fpath;
    max_freq = eval(p.Results.max_freq);
    downsample_freq = eval(p.Results.downsample_freq);
    
    if ~exist(channel_fpath, 'file')
       error('%s is not a valid file\n', channel_fpath);
    end
    
    channel_fpath_splits = strsplit(channel_fpath, '/');
    channel_path = strjoin(channel_fpath_splits(1:end -1 ), '/');
    channel_fname = channel_fpath_splits{end};
    
    error_save_fpath = [ channel_path '/filt_' channel_fname '.log' ];
    mat_save_fpath = [ channel_path '/filt_' channel_fname ];
    
    if exist(error_save_fpath, 'file')
        delete(error_save_fpath);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load channel data
    fprintf('load\n');
    
    channelInfo = load(channel_fpath);
    channelInfo.channel_data = double(channelInfo.channel_data);
    
    channel_name = channelInfo.channel_name;
    
    if isempty(channelInfo.channel_data)
       
        error_fid = fopen(error_save_fpath, 'w');
        fprintf(error_fid, 'channel_data was empty when opened\n');
        fclose(error_fid);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %find samples with constant values for > 300 samples and set to NaN
    
    [rem_starts, rem_ends] = removeSat(channelInfo.channel_data);
    
    fprintf('detect chuncks of constants\n');
    
    fprintf('length(channel_data): %d NaNs: %d\n', length(channelInfo.channel_data), sum(isnan(channelInfo.channel_data)));

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %low pass filter
    fprintf('low pass filter\n');    
    
    [filtered_data, ~] = buttfilt(channelInfo.channel_data, max_freq, channelInfo.samplingFreq, 'low', 2);
    
    nan_detect_filtered_data = filtered_data;
    for j = 1:length(rem_starts)
        nan_detect_filtered_data(rem_starts(j):rem_ends(j))=NaN;
    end  
    
    fprintf('length(filtered_data): %d NaNs: %d\n', length(filtered_data), sum(isnan(filtered_data)));

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %downsample
    fprintf('downsample\n');    
    
    downsample_increment = floor( channelInfo.samplingFreq / downsample_freq );
    num_timepoints = length( channelInfo.channel_data );
    subsample_timepoints = 1:downsample_increment:num_timepoints;
    
    noreref = filtered_data(subsample_timepoints);
    
    fprintf('length(noreref): %d NaNs: %d\n', length(noreref), sum(isnan(noreref)));
    
    nan_detect_noreref = nan_detect_filtered_data(subsample_timepoints);
    set2nan = isnan(nan_detect_noreref);
    
    filtered_data = [];
    channelInfo.channel_data = [];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %detrend
    fprintf('detrend\n');    

    detrended = locdetrend(noreref, downsample_freq, [1 0.5]);
    
    fprintf('length(detrended): %d NaNs: %d\n', length(detrended), sum(isnan(detrended)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %remove line noise 
    fprintf('remove line noise\n');    

    line_params = struct('tapers', [3 5], 'Fs', downsample_freq, 'pad', -1);
    line_freqs = 60.*[1 2];
    t_params = [4 1];

    [detrended_noLineNoise, ~, ~, ~] = rmlinesmovingwinc(detrended , t_params , 10, line_params, [], 'n', line_freqs);
    
    fprintf('length(detrended_noLineNoise): %d NaNs: %d\n', length(detrended_noLineNoise), sum(isnan(detrended_noLineNoise)));
%     
%     %set the NaNs
%     
%     detrended_noLineNoise(set2nan) = NaN;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %save out
    
    detrended_noLineNoise = detrended_noLineNoise';
    samplingFreq = downsample_freq;
    
    fprintf('save mat file\n');    
    save( mat_save_fpath , '-v7.3', 'noreref', 'detrended_noLineNoise', 'samplingFreq' , 'channel_name', 'set2nan');
    delete(channel_fpath);
end


function [rem_starts, rem_ends] = removeSat(ts)

    res_ts = ts;

    rem_sat = 300;
    
    run_starts = [0 find(diff(ts)~=0)] + 1;
    run_lengths = [diff(run_starts) (numel(ts) - run_starts(end) + 1)];
    rem_starts = run_starts(run_lengths>rem_sat);
    rem_ends = rem_starts + run_lengths(run_lengths>rem_sat);
    edge_ind= rem_ends > length(ts);
            
    rem_ends(edge_ind) = length(ts);
    
end