function globalReref_allchan(input_fpath, output_fpath, refset, used_chan_fpath)

    echo globalReref_allchan on;
    fprintf('starting\n');

    if ~exist(input_fpath, 'file')
       error('not a valid path, %s', input_fpath); 
    end
    
    if ~exist(used_chan_fpath, 'file')
       error('not a valid path, %s', used_chan_fpath); 
    end
    
    [input_path, ~, ~] = fileparts(input_fpath);
    
    bandpass_quantile_output_fpath = [input_path '/bandpass_'  sprintf('quantiles.refset%s.mat', refset)];
    reref_quantile_output_fpath = [input_path '/reref_'  sprintf('quantiles.refset%s.mat', refset)];
    
    fprintf('readmda\n');
    input_data = readmda(input_fpath);
    size(input_data)
    
    load([input_path '/samplingFreq.mat']);
    
    used_chans_fid = fopen(used_chan_fpath);
    used_chans = textscan(used_chans_fid, '%s\n');
    used_chans = used_chans{1};
    fclose(used_chans_fid);
    
    %calcWindowStats(input_data, sampling_freq, bandpass_quantile_output_fpath, used_chans);
    
    fprintf('mean_calc\n');
    global_avg = mean(input_data, 1);

    for iRow = 1:size(input_data, 1)
       
        input_data(iRow, :) = input_data(iRow, :) - global_avg;
    end
    
    size(input_data)
    writemda(input_data, output_fpath, 'float32');
    
    calcWindowStats(input_data, sampling_freq, reref_quantile_output_fpath, used_chans);
    
    fprintf('globalReref_allchan -- done');

end