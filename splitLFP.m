function splitLFP(varargin)

    echo splitLFP on;
    
    fprintf('****************************************************\n');
    fprintf('* %s\n', mfilename);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %params and input checks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p = inputParser;

    p.addParameter('nsx_fpath', '', @ischar);
    p.addParameter('save_dir', '', @ischar);
    
    parse(p, varargin{:});

    disp(p.Results);
    
    nsx_fpath = p.Results.nsx_fpath;
    save_dir = p.Results.save_dir;

    if ~exist(nsx_fpath, 'file')
        error('%s is not a valid filepath\n', nsx_fpath);
    end
    
    if ~exist(save_dir, 'dir')
        error('%s is not a valid directory\n', save_dir);
    end
    
    applyGain = 0;
    nsx = concatOpenNSx(nsx_fpath, applyGain);
    
    nsx_postProc = nsx.postProc;
    save([save_dir '/nsx_postProc.mat'], 'nsx_postProc');
    
    %if there is a < 5 min session, don't bother
    if sum(nsx.MetaTags.DataDurationSec)/60 < 5
       
        ignore_fid = fopen([save_dir '/_ignore_me.txt'], 'w');   
        fprintf(ignore_fid, '%s', 'session duration less than 5 min');
        fclose(ignore_fid);
    
        error('ignore condition for this session');
    end
    
    split_dir = [save_dir '/lfp_splits'];
    mkdir(split_dir);
    
    %fiddle with electrodeLabels
    electrodeLabels = { nsx.ElectrodesInfo.Label };
        
    %remove null spaces from electrodeLabels
    for iElec = 1:length(electrodeLabels)

        currentElec = electrodeLabels{iElec};
        endIdx = 1;

        while (endIdx + 1) <= length(currentElec) && double(currentElec(endIdx + 1)) ~= 0
            endIdx = endIdx + 1;     
        end
        
        electrodeLabels{iElec} = electrodeLabels{iElec}(1:endIdx);
    end
        
    
    samplingFreq = nsx.MetaTags.SamplingFreq;
    
    channels_split = 1;
    
    for iChan = 1:size(nsx.Data, 1)
       
        fprintf('splitting channel %d\n', iChan);
        
        channel_name = electrodeLabels{iChan};
        
        if ~contains(channel_name, 'ain') && ~contains(lower(channel_name), 'stim') 
        
            channel_data = nsx.Data(iChan, :);

            save( [ split_dir '/' sprintf('ch_%03d', channels_split) '.mat'] , '-v7.3', 'channel_data', 'channel_name', 'samplingFreq');
            
            channels_split = channels_split + 1;

        end
    end
    

end