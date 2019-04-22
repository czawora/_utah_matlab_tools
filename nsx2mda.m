function nsx2mda(varargin)

    % parse inputs and validate

    fprintf('****************************************************\n');
    fprintf('* %s\n', mfilename);

    p = inputParser;
    p.addParameter('nsx_fpath', '', @ischar);
    p.addParameter('jacksheet_fpath', '', @ischar);
    p.addParameter('refset', '', @ischar);
    p.addParameter('session_dir', '', @ischar);
    p.addParameter('min_range_cutoff_microvolt', '10', @ischar);
    p.addParameter('min_duration_data_min', '5', @ischar);
    
    p.addParameter('ignoreShortFile', '1', @ischar);
    
    parse(p, varargin{:});
    
    disp(p.Results);
    
    nsx_fpath = p.Results.nsx_fpath;  
    jacksheet_fpath = p.Results.jacksheet_fpath;
    refset = p.Results.refset;
    session_dir = p.Results.session_dir;
    min_range_cutoff_ungained = eval(p.Results.min_range_cutoff_ungained);
    min_duration_data_min = eval(p.Results.min_duration_data_min);
    
    if isequal(refset, '') 
        error('refset must be an string integer matching a microDevNum in the jacksheet');
    else
        refset = eval(refset);
    end
    
    ignoreShortFile = logical(str2num(p.Results.ignoreShortFile));
    
    if ~exist(nsx_fpath, 'file')
       fprintf('%s is not a valid file\n', nsx_fpath); 
    end
    
    if ~exist(jacksheet_fpath, 'file')
       fprintf('%s is not a valid file\n', jacksheet_fpath); 
    end
    
    if ~exist(session_dir, 'dir')
       fprintf('%s is not a valid dir\n', session_dir); 
    end
     
    ignore_me_fname = [session_dir sprintf('/_ignore_me%d.txt', refset)];

    %convert min_range_cutoff_microvolt to millivolt
    min_range_cutoff_millivolt = min_range_cutoff_ungained * 0.25 * (1/1000); 
    
    % read the jacksheet    
    jacksheet = readtable(jacksheet_fpath);
    
    % jacksheet for this refset
    used_jacksheet = jacksheet( jacksheet{:,'MicroDevNum'} == refset && jacksheet{:,'RangeMilliV'} > min_range_cutoff_millivolt , : );
    
    % check if all the channels failed to pass the range filter
    if isempty(used_jacksheet)
       
        ignore_me_fid = fopen(ignore_me_fname, 'w');
        fprintf(ignore_me_fid, 'all channels have data range < %0.4f milliV', min_range_cutoff_millivolt);
        fclose(ignore_me_fid);
    end
    
    % check if the file is too short to care about
    if ignoreShortFile
        
        if used_jacksheet{1, 'DurationMin'} < min_duration_data_min
            
            ignore_me_fid = fopen(ignore_me_fname, 'w');
            fprintf(ignore_me_fid, 'data length less than 5 min ( %0.2f )', nsx.MetaTags.DataDurationSec/60);
            fclose(ignore_me_fid);
        end
        
    end
        

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
    
    % which electrodeLabels are in the refset
    refset_nsx_chan_filt = cellfun( @(x) any(cellfun( @(y) isequal(x, y) , used_jacksheet{:, 'ChanName'})), electrodeLabels);
   
    output_mda_fpath = [session_dir '/' sprintf('refset%d.mda', refset)];
    output_used_jacksheet_fpath = [session_dir '/' sprintf('jacksheet_refset%d.csv', refset)];
    
    writemda(nsx.Data(refset_nsx_chan_filt, :), output_mda_fpath, 'int16');
    writetable(used_jacksheet, output_used_jacksheet_fpath);

end