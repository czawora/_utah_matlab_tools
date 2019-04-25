function calcWindowVarStats(data_mat, sampling_freq, output_fpath, jacksheet)

    if size(data_mat, 1) ~= size(jacksheet, 1)
        error('number of jacksheet rows must match number of channels in data_mat (rows)');
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
    
    var_stats = struct;
    var_stats.std_mats = std_mats;
    var_stats.rms_mats = rms_mats;
    var_stats.jacksheet = jacksheet;
    var_stats.window_center_sec = window_center_sec;
    
    save(output_fpath, 'var_stats');

end