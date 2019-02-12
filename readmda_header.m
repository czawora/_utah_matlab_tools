function readmda_header(fname)

    if ~exist(fname, 'file')
        error('%s - file does not exist', fname);
    end
    
    F=fopen(fname,'rb');

    try
    code=fread(F,1,'int32');
    catch
        error('Problem reading file: %s',fname);
    end
    
    if (code>0) 
        
        num_dims=code;
        code=-1;
    else
        
        fread(F,1,'int32');
        num_dims=fread(F,1,'int32');    
    end

end