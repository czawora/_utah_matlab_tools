function filter_downsample_detrend_rmLineNoise(varargin)

    echo filter_downsample_detrend_rmLineNoise on;
    
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %load channel data
    fprintf('load\n');
    
    channelInfo = load(channel_fpath);
    channelInfo.channel_data = double(channelInfo.channel_data);
    
    channel_name = channelInfo.channel_name;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %low pass filter
    fprintf('low pass filter\n');    
    
    [filtered_data, ~] = buttfilt(channelInfo.channel_data, max_freq , channelInfo.samplingFreq, 'low', 2);
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %downsample
    fprintf('downsample\n');    
    
    downsample_increment = floor( channelInfo.samplingFreq / downsample_freq );
    num_timepoints = length( channelInfo.channel_data );
    subsample_timepoints = 1:downsample_increment:num_timepoints;
    
    noreref = filtered_data(subsample_timepoints);
    
    filtered_data = [];
    channelInfo.channel_data = [];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %detrend
    fprintf('detrend\n');    

    detrended = locdetrend(noreref, downsample_freq, [1 0.5]);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %remove line noise 
    fprintf('remove line noise\n');    

    line_params = struct('tapers', [3 5], 'Fs', downsample_freq, 'pad', -1);
    line_freqs = 60.*[1 2];
    t_params = [4 1];

    [detrended_noLineNoise, ~, ~, ~] = rmlinesmovingwinc(detrended , t_params , 10, line_params, [], 'n', line_freqs);
    
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %remove artifacts where V is constant for at least 10 samples, not
    %impletemented yet
    
    %in order to figure out the constant jitter, save out the data range
    channel_range = max(detrended_noLineNoise) - min(detrended_noLineNoise);
    
    
    detrended_noLineNoise = detrended_noLineNoise';
    samplingFreq = downsample_freq;
    
    channel_fpath_splits = strsplit(channel_fpath, '/');
    channel_path = strjoin(channel_fpath_splits(1:end -1 ), '/');
    channel_fname = channel_fpath_splits{end};
    
    fprintf('save mat file\n');    
    save( [ channel_path '/filt_' channel_fname ] , '-v7.3', 'noreref', 'detrended_noLineNoise', 'channel_range', 'samplingFreq' , 'channel_name');
    delete(channel_fpath);
end