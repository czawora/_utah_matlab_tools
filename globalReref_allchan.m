function globalReref_allchan(input_fpath, output_fpath)

    fprintf('starting\n');

    if ~exist(input_fpath, 'file')
       error('not a valid path, %s', input_fpath); 
    end
          
    fprintf('readmda\n');
    input_data = readmda(input_fpath);
        
    fprintf('mean_calc\n');
    global_avg = mean(input_data, 1);

    for iRow = 1:size(input_data, 1)
       
        input_data(iRow, :) = input_data(iRow, :) - global_avg;
    end
    
    writemda(input_data, output_fpath, 'float32');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate some variance stats
    
    output_fpath_splits = strsplit(output_fpath, '/');
    output_path = strjoin(output_fpath_splits(1:(end-1)), '/');
    output_fname = output_fpath_splits{end};
    
    output_var_stats_fpath = [ output_path '/' output_fname '_var_stats.mat'];
    
    chan_var = var(input_data, 0, 2);
    chan_rms = rms(input_data, 2);
    
    chan_var_stats = struct;
    chan_var_stats.chan_var = chan_var;
    chan_var_stats.chan_rms = chan_rms;
    
    save(output_var_stats_fpath, 'chan_var_stats');
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fprintf('globalReref_allchan -- done');

end