function nsx2mda(varargin)

    % this function reads an nsx file as well as an elementInfo (containing a list of channel name inclusion strings)
    % channels from the nsx file are included in the output mda file if they pass two criteria: 
    % 1) the channel name must include one of the channel name inclusion strings in elementInfo as a substring
    % 2) the data range of that channel must be greater than min_range_cutoff (currently hardcoded as 10)

    % Inputs: input_filename - the file path to the input nsx file to be converted
    %         elementInfo_filename - the file path of the text file containing a list of channel inclusion strings
    %         output_dir - the path of the directory for outputs
    %         ignoreShortFile - (default 1) set to 0 to ignore the "min_duration_data" criteria
    
    % Outputs: mda file - the converted nsx file data
    %          chans.csv - two column file containing all channels present in nsx file and a string indicating whether they were selected in elementInfo and if that channels data range was adequate
    %          used_chans.txt - a text file containing a list of channels that were included in the converted mda file
    %          (_ignore_me.txt) - a file created when the nsx file includes no desirable channels
    
    
    %some constants
    min_range_cutoff = 10; % range of data must be > 10 to be considered, good. Very conservative filter
    min_duration_data = 5; % minimum number of minutes the recording must be care about

    % parse inputs and validate

    fprintf('****************************************************\n');
    fprintf('* %s\n', mfilename);

    p = inputParser;
    p.addParameter('nsx_fpath', '', @ischar);
    p.addParameter('jacksheet_fpath', '', @ischar);
    p.addParameter('refset', '', @ischar);
    p.addParameter('session_dir', '', @ischar);
    
    p.addParameter('ignoreShortFile', '1', @ischar);
    
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
        
    %read the jacksheet    
    jacksheet = readtable(jacksheet_fpath);
 
    %create output fname
    
    nsx_fpath_splits = strsplit(nsx_fpath, '/');
    nsx_fpath_no_path = nsx_fpath_splits{end};
    
    nsx_fpath_no_path_split = strsplit(nsx_fpath_no_path, '.');
    nsx_fpath_no_path_ext = nsx_fpath_no_path_split{1};
     
    ignore_me_fname = [session_dir sprintf('/_ignore_me%d.txt', refset)];
    %output_chans_csv_fname = [output_dir '/' input_fname_no_path_no_ext '_chans.csv'];
    
    
    %read the data
    nsx = concatOpenNSx(nsx_fpath);
    
    %save this for use in the final step of the pipeline
    nsx_postProc = nsx.postProc;
    save([session_dir '/nsx_postProc.mat'], 'nsx_postProc');

    %check if the file is too short to care about
    if ignoreShortFile
       
        if nsx.MetaTags.DataDurationSec/60 < min_duration_data
           
            ignore_me_fid = fopen(ignore_me_fname, 'w');
            fprintf(ignore_me_fid, 'data length less than 5 min ( %0.2f )', nsx.MetaTags.DataDurationSec/60);
            fclose(ignore_me_fid);
        end
    end

    
    %calculate data ranges
    ranges = range(nsx.Data, 2);
    
    if refset == 1

        % file records all channel names we start with
        output_chans_csv_fid = fopen(output_chans_csv_fname, 'w');

        %print header to csv
        fprintf(output_chans_csv_fid, '"channel","info"\n');  

    end
    
    %filter for channels in elementInfo_micro and channel range
    nsx_channelNames = {nsx.ElectrodesInfo.Label};
    nsx_channelNames_trim = {};
    nsx_channelNames_match = {};
    nsx_channelNames_matchIdx = [];
    
    for iChan = 1:length(nsx_channelNames)
        
        current_channel = nsx_channelNames{iChan};
        
        % trim null characters at end of channel name
        current_channel_trim = ''; 
        for iChar = 1:length(current_channel)
            
            if ~isequal(sprintf('%d', current_channel(iChar)), '0')
                current_channel_trim = [current_channel_trim current_channel(iChar)];
            end
        end

        nsx_channelNames_trim{length(nsx_channelNames_trim) + 1} = current_channel_trim;
        
        fprintf('checking channel %s\n', current_channel_trim);

        %check whether channel was selected and if it has adequate range
        selected_status = 'unselected';
        quality_status = 'good';
        
        if any(cellfun(@(x) isequal(x, current_channel_trim), channel_set))          
            selected_status = 'selected';
        end
        
        if ranges(iChan) < min_range_cutoff
            quality_status = 'bad';
        end
        
        if isequal(selected_status, 'selected') && isequal(quality_status, 'good')
            
            nsx_channelNames_match{length(nsx_channelNames_match) + 1} = current_channel_trim;
            nsx_channelNames_matchIdx = [nsx_channelNames_matchIdx iChan];
            
        end
        
        fprintf(sprintf('"%s","%s"\n', current_channel_trim, [selected_status '_' quality_status]));        
        
        if refset == 1
            fprintf(output_chans_csv_fid, sprintf('"%s","%s"\n', current_channel_trim, [selected_status '_' quality_status]));
        end
    end
    
    if refset == 1
        fclose(output_chans_csv_fid);
    end
    
     % if no channel pass, write and _ignore_me, else write the data out, reference set
    if isempty(nsx_channelNames_match)
       
        ignore_me_fid = fopen(ignore_me_fname, 'w');
        fprintf(ignore_me_fid, 'all channels unselected or have data range < 10');
        fclose(ignore_me_fid);
    else
        
        for iRef = 1:size(elInfo, 1)

            output_chans_tmp_fname = [output_dir '/' input_fname_no_path_no_ext sprintf('_refset%d_used_chans.tmp%d', iRef, refset)];
            output_chans_txt_fname = [output_dir '/' input_fname_no_path_no_ext sprintf('_refset%d_used_chans.txt', iRef)];
            
            %find names in nsx_channelNames_match that are in chan_names_by_refset{iRef}
            current_refset_chan_names = chan_names_by_refset{iRef};
            current_refset_idx = [];

            output_chans_tmp_fid = fopen(output_chans_tmp_fname, 'w');        

            for iRefChan = 1:length(current_refset_chan_names)

                current_refset_chan = current_refset_chan_names{iRefChan};
                current_refset_chan_binar = cellfun(@(x) isequal(x, current_refset_chan), nsx_channelNames_match);

                if any(current_refset_chan_binar)

                    fprintf(output_chans_tmp_fid, sprintf('%s\n', current_refset_chan));
                    current_refset_idx(length(current_refset_idx) + 1) = find(current_refset_chan_binar);
                end

            end

            fclose(output_chans_tmp_fid);
            movefile(output_chans_tmp_fname, output_chans_txt_fname);
            

            if iRef == refset
            
                output_mda_fname = [output_dir '/' input_fname_no_path_no_ext sprintf('.refset%d.mda', iRef)];

                writemda(nsx.Data(current_refset_idx, :), output_mda_fname, 'int16');
                
                %calcWindowStats(nsx.Data(current_refset_idx, :), nsx.MetaTags.SamplingFreq, [output_dir '/raw_'  sprintf('quantiles.refset%d.mat', iRef)], current_refset_chan_names);
                
                sampling_freq = nsx.MetaTags.SamplingFreq;
                save([output_dir '/samplingFreq.mat'], 'sampling_freq');
            
            end  
            
        end
    end
    
   
end