function nsx2mda(varargin)

    % parse inputs and validate

    fprintf('****************************************************\n');
    fprintf('* %s\n', mfilename);

    p = inputParser;
    p.addParameter('nsx_fpath', '/Volumes/72A/UTAH_A/NIH036/data_raw/151211_1041/151211_1041_utah1.ns6', @ischar);
    p.addParameter('jacksheet_fpath', '/Volumes/72A/UTAH_A/NIH036/data_raw/151211_1041/jacksheetBR_complete.csv', @ischar);
    p.addParameter('refset', '1', @ischar);
    p.addParameter('session_dir', '/Users/zaworaca/Desktop/Desktop/debug/ms_debug/NIH036_151211_1041', @ischar);
   
    parse(p, varargin{:});
    
    disp(p.Results);
    
    nsx_fpath = p.Results.nsx_fpath;  
    jacksheet_fpath = p.Results.jacksheet_fpath;
    refset = p.Results.refset;
    session_dir = p.Results.session_dir;
    
    if isequal(refset, '') 
        error('refset must be an string integer matching a microDevNum in the jacksheet');
    else
        refset = eval(refset);
    end
        
    if ~exist(nsx_fpath, 'file')
       fprintf('%s is not a valid file\n', nsx_fpath); 
    end
    
    if ~exist(jacksheet_fpath, 'file')
       fprintf('%s is not a valid file\n', jacksheet_fpath); 
    end
    
    if ~exist(session_dir, 'dir')
       fprintf('%s is not a valid dir\n', session_dir); 
    end
     
    % read the jacksheet    
    used_jacksheet = readtable(jacksheet_fpath);

    % read the data
    nsx = concatOpenNSx(nsx_fpath);
    
    % save this for use in the final step of the pipeline
    nsx_postProc = nsx.postProc;
    save([session_dir '/nsx_postProc.mat'], 'nsx_postProc');

    
    % fiddle with electrodeLabels
    electrodeLabels = { nsx.ElectrodesInfo.Label };
        
    % remove null spaces from electrodeLabels
    for iElec = 1:length(electrodeLabels)

        currentElec = electrodeLabels{iElec};
        endIdx = 1;

        while (endIdx + 1) <= length(currentElec) && double(currentElec(endIdx + 1)) ~= 0
            endIdx = endIdx + 1;     
        end
        
        electrodeLabels{iElec} = electrodeLabels{iElec}(1:endIdx);
    end
       
    % select channels from data based on filtered jacksheet
    % make sure data channels are ordered as they appear in the jacksheet
    reorder_idx = [];
    
    for iRow_jack = 1:size(used_jacksheet, 1)
        
        % current chan name in jacksheet
        chan_name_jack = used_jacksheet{iRow_jack, 'ChanName'}{1};
        
        % find this channel in electrodeLabels
        for iDataChan = 1:length(electrodeLabels)
            
            % does this index indicate the channel I am looking for
            if isequal(chan_name_jack, electrodeLabels{iDataChan})
                
                reorder_idx = [ reorder_idx iDataChan ];
            end
        end
    end
   
    output_mda_fpath = [session_dir '/' sprintf('refset%d.mda', refset)];
    writemda(nsx.Data(reorder_idx, :), output_mda_fpath, 'int16');
    
end