function csv_dat = readtextcsv(fpath)

    fid = fopen(fpath);
    
    top_line = fgets(fid);
    
    if contains(top_line, '","')
        
       splitter_string = '","';
    else
        
       splitter_string = ',';
    end
    
    top_line = strip(top_line, 'both', '"');
    splits = strsplit(top_line, splitter_string);
    
    csv_dat = splits;
    
    next_line = fgets(fid);
    
    while next_line ~= -1

        next_line = strip(next_line, 'both', '"');
        splits = strsplit(next_line, splitter_string);
        
        csv_dat = [csv_dat ; splits];
        next_line = fgets(fid);
    end
    
    fclose(fid);

end