function splitmda(varargin)

    echo splitmda on;
    % parse inputs and validate

    fprintf('****************************************************\n');
    fprintf('* %s\n', mfilename);

    p = inputParser;
    p.addParameter('input_filename', '', @ischar);
    p.addParameter('used_chans_fpath', '', @ischar);
    p.addParameter('refset', '', @ischar);
    p.addParameter('output_dir', '', @ischar);
        
    parse(p, varargin{:});
    
    disp(p.Results);
    
    input_filename = p.Results.input_filename;
    used_chans_fpath = p.Results.used_chans_fpath;
    output_dir = p.Results.output_dir;
    refset = p.Results.refset;
    
    if isequal(refset, '') 
        error('refset must be an string integer corresponding to a row in elementInfo');
    else
        refset = eval(refset);
    end
    
     if ~exist(input_filename, 'file')
       fprintf('%s is not a valid file\n', input_fname); 
     end
    
     if ~exist(used_chans_fpath, 'file')
       fprintf('%s is not a valid file\n', used_chans_fpath); 
     end
    
    if ~exist(output_dir, 'dir')
       fprintf('%s is not a valid dir\n', output_dir); 
    end
    
    used_chans_fid = fopen(used_chans_fpath);
    used_chans = textscan(used_chans_fid, '%s\n');
    used_chans = used_chans{1};
    fclose(used_chans_fid);
    
    [input_path, ~, ~] = fileparts(input_filename);
    
    input_fname_splits = strsplit(input_filename, '/');
    input_fname_no_path = input_fname_splits{end};
    
    input_fname_no_path_split = strsplit(input_fname_no_path, '.');
    input_fname_no_path_no_ext = input_fname_no_path_split{1};
    
    fprintf('reading mda file\n');
    mda = readmda(input_filename);
    
    load([input_path '/samplingFreq.mat']);
    
    whiten_quantile_output_fpath = [input_path '/whiten_'  sprintf('quantiles.refset%d.mat', refset)];
    %calcWindowStats(mda, sampling_freq, whiten_quantile_output_fpath, used_chans);

    
    for iChan = 1:size(mda, 1)
       
        fprintf('writing channel %d %s\n', iChan, used_chans{iChan});
        
        channel_output_dir = [output_dir '/' input_fname_no_path_no_ext sprintf('_%s', used_chans{iChan})];
        mkdir(channel_output_dir);
        
        channel_output_fname = [channel_output_dir '/' input_fname_no_path_no_ext sprintf('_%s', used_chans{iChan}) sprintf('.refset%d_mda_chan', refset)];
        writemda(mda(iChan,:), channel_output_fname);
    end
    
    fprintf('done\n');

end