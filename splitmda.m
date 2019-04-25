function splitmda(varargin)

    % parse inputs and validate

    fprintf('****************************************************\n');
    fprintf('* %s\n', mfilename);

    p = inputParser;
    p.addParameter('input_filename', '', @ischar);
    p.addParameter('refset', '', @ischar);
    p.addParameter('output_dir', '', @ischar);
        
    parse(p, varargin{:});
    
    disp(p.Results);
    
    input_filename = p.Results.input_filename;
    output_dir = p.Results.output_dir;
    refset = p.Results.refset;
    
    if isequal(refset, '') 
        error('refset must be an string integer');
    else
        refset = eval(refset);
    end
    
     if ~exist(input_filename, 'file')
       fprintf('%s is not a valid file\n', input_fname); 
     end
    
    if ~exist(output_dir, 'dir')
       fprintf('%s is not a valid dir\n', output_dir); 
    end
        
    input_fname_splits = strsplit(input_filename, '/');
    input_fname_no_path = input_fname_splits{end};
    
    input_fname_no_path_split = strsplit(input_fname_no_path, '.');
    input_fname_no_path_no_ext = input_fname_no_path_split{1};
    
    fprintf('reading mda file\n');
    mda = readmda(input_filename);
    
    for iChan = 1:size(mda, 1)
       
        fprintf('writing channel %d \n', iChan);
        
        channel_output_dir = [output_dir '/' input_fname_no_path_no_ext sprintf('_%d', iChan)];
        mkdir(channel_output_dir);
        
        channel_output_fname = [channel_output_dir '/' input_fname_no_path_no_ext sprintf('_%03d', iChan) '.mda_chan'];
        writemda(mda(iChan,:), channel_output_fname);
    end
    
    fprintf('done\n');

end