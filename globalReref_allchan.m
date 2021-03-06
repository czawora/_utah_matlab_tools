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
    
    fprintf('globalReref_allchan -- done');

end