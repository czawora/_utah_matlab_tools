function calcWindowStats(data_mat, sampling_freq, output_fpath, chan_names)

    if size(data_mat, 1) ~= length(chan_names)
        error('number of chan_names must match number of channels in data_mat (rows)');
    end

    % calculate the quantiles of the raw data  
    
    samples_per_ms = sampling_freq / 1000;
    num_data_samples = size(data_mat, 2);
    
    ms_window_size = 1000 * 60 * 5; % 5 minute windows
    sample_window_size = ms_window_size * samples_per_ms;

    window_slide_frac = 1;
    window_move_num = floor(sample_window_size/window_slide_frac);

    std_mats = [];
    rms_mats = [];

    start_window_idx = 1;
    stop_window_idx = start_window_idx + sample_window_size;

    window_center_sec = [];
    
    while stop_window_idx < num_data_samples

       window_center_sec = [ window_center_sec ((stop_window_idx - start_window_idx)/2 * (1/sampling_freq))];
        
       %window_data is a time x channel matrix 
       window_data = double(data_mat(:, start_window_idx:stop_window_idx)');

       std_res = std(window_data, 0, 1);
       rms_res = rms(window_data, 1);
       
       std_mats = [ std_mats std_res' ];
       rms_mats = [ rms_mats rms_res' ];

       start_window_idx = start_window_idx + window_move_num;
       stop_window_idx = stop_window_idx + window_move_num;

    end
    
    save(output_fpath, 'std_mats', 'rms_mats', 'chan_names', 'window_center_sec');

end