function calcWindowStats(data_mat, sampling_freq, output_fpath, chan_names)

    if size(data_mat, 1) ~= length(chan_names)
        error('number of chan_names must match number of channels in data_mat (rows)');
    end

    % calculate the quantiles of the raw data  

    samples_per_sec = sampling_freq * 1;
    num_data_samples = size(data_mat, 2);

    window_move_num = floor(samples_per_sec/2);

    quantile_cutoffs = [0.001 0.01 0.99 0.999];

    quantile_mats = cell(length(quantile_cutoffs), 1);
    mean_mats = [];
    std_mats = [];
    rms_mats = [];

    start_window_idx = 1;
    stop_window_idx = start_window_idx + samples_per_sec;

    start_windows = [];
    
    while stop_window_idx < num_data_samples

       start_windows = [ start_windows start_window_idx ];
        
       %fprintf('total samples: %d, start: %d - stop: %d -- %0.2f \n', num_data_samples, start_window_idx, stop_window_idx, (stop_window_idx/num_data_samples) * 100); 

       %window_data is a time x channel matrix 
       window_data = double(data_mat(:, start_window_idx:stop_window_idx)');

       %quantiles will be calculated on all columns of window_data
       %quantile_res will be num_quantiles x channel matrix
       quantile_res = quantile(window_data, quantile_cutoffs, 1);

       for iQuantile = 1:length(quantile_cutoffs)

           quantile_mats{iQuantile} = [ quantile_mats{iQuantile} quantile_res(iQuantile, :)' ];
       end
       
       mean_res = mean(window_data, 1);
       std_res = std(window_data, 0, 1);
       rms_res = rms(window_data, 1);
       
       mean_mats = [ mean_mats mean_res' ];
       std_mats = [ std_mats std_res' ];
       rms_mats = [ rms_mats rms_res' ];

       start_window_idx = start_window_idx + window_move_num;
       stop_window_idx = stop_window_idx + window_move_num;

    end
    
    
    %calculate histogram bins
    
    hist_edges = cell(size(data_mat, 1), 1);
    hist_n = cell(size(data_mat, 1), 1);
    
    for iChan = 1:size(data_mat, 1)
                
        fprintf('calculating histogram for channel indx %d\n', iChan);
        
        [n , edges ] = histcounts( data_mat(iChan,:) ); 
        hist_edges{iChan} = edges;
        hist_n{iChan} = n;
        
    end
    
    
    

    fprintf('saving quantiles mat\n');
    save(output_fpath, 'hist_edges', 'hist_n', 'quantile_mats', 'mean_mats', 'std_mats', 'rms_mats', 'chan_names', 'start_windows', 'sampling_freq', 'samples_per_sec', 'num_data_samples', 'quantile_cutoffs');

end