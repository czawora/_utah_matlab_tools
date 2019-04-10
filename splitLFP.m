function splitLFP(varargin)
    
    fprintf('****************************************************\n');
    fprintf('* %s\n', mfilename);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % params and input checks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p = inputParser;

    p.addParameter('nsx_fpath', '', @ischar);
    p.addParameter('jacksheet_fpath', '', @ischar);
    p.addParameter('nsp_suffix', '', @ischar);
    p.addParameter('save_dir', '', @ischar);
    
    p.addParameter('debug', '0', @ischar);
    
    parse(p, varargin{:});

    disp(p.Results);
    
    nsx_fpath = p.Results.nsx_fpath;
    jacksheet_fpath = p.Results.jacksheet_fpath;
    nsp_suffix = p.Results.nsp_suffix;
    save_dir = p.Results.save_dir;
    
    debug = eval(p.Results.debug);

    if ~exist(nsx_fpath, 'file')
        error('%s is not a valid filepath\n', nsx_fpath);
    end
    
    if ~exist(jacksheet_fpath, 'file')
        error('%s is not a valid filepath\n', jacksheet_fpath);
    end
    
    if ~exist(save_dir, 'dir')
        error('%s is not a valid directory\n', save_dir);
    end
    
    applyGain = 0;
    nsx = concatOpenNSx(nsx_fpath, applyGain);
    
    nsx_postProc = nsx.postProc;
    save([save_dir '/nsx_postProc.mat'], 'nsx_postProc');
    
    % if there is a < 5 min session, don't bother
    if ~debug && sum(nsx.MetaTags.DataDurationSec)/60 < 5
       
        ignore_fid = fopen([save_dir '/_ignore_me.txt'], 'w');   
        fprintf(ignore_fid, '%s', 'session duration less than 5 min');
        fclose(ignore_fid);
    
        error('ignore condition for this session');
    end
    
    % make split dir
    split_dir = [save_dir '/lfp_splits'];
    mkdir(split_dir);
    
    
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
        

    % read jacksheet and find micro channels to split
    jacksheet = readtable(jacksheet_fpath);
    
    nsp_suffix_filter = cellfun(@(x) isequal(x, nsp_suffix), jacksheet{:, 'NSPsuffix'});
    micro_chan_filter = jacksheet{:, 'MicroDevNum'} > 0;
    
    current_channel_jacksheet = jacksheet(nsp_suffix_filter & micro_chan_filter, :);
    
    % write out used jacksheet
    writetable(current_channel_jacksheet, [save_dir '/used_jacksheet.csv']);
    
    
    samplingFreq = current_channel_jacksheet{1, 'SampFreq'};
    
    
    for iChan = 1:size(current_channel_jacksheet, 1)
       
        fprintf('splitting channel %d\n', iChan);
        
        channel_name = current_channel_jacksheet{ iChan, 'ChanName' }{1};
        
        nsp_chan_index = cellfun(@(x) isequal(channel_name, x), electrodeLabels);
        
        channel_data = nsx.Data(nsp_chan_index, :);

        save( [ split_dir '/' sprintf('ch_%03d', iChan) '.mat'] , '-v7.3', 'channel_data', 'channel_name', 'samplingFreq');
        pause(2);
        
    end
    

end