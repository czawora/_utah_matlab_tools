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
    p.addParameter('input_filename', '', @ischar);
    p.addParameter('elementInfo_filename', '', @ischar);
    p.addParameter('refset', '', @ischar);
    p.addParameter('output_dir', '', @ischar);
    
    p.addParameter('ignoreShortFile', '1', @ischar);
    
    parse(p, varargin{:});
    
    disp(p.Results);
    
    input_fname = p.Results.input_filename;  
    elementInfo_filename = p.Results.elementInfo_filename;
    output_dir = p.Results.output_dir;
    refset = p.Results.refset;
    
    ignoreShortFile = logical(str2num(p.Results.ignoreShortFile));
    
    if ~exist(input_fname, 'file')
       fprintf('%s is not a valid file\n', input_fname); 
    end
    
    if ~exist(elementInfo_filename, 'file')
       fprintf('%s is not a valid file\n', elementInfo_filename); 
    end
    
    if ~exist(output_dir, 'dir')
       fprintf('%s is not a valid dir\n', output_dir); 
    end
        
    %read the elementInfo_micro
    channel_set = [];
    
    elInfo = readtextcsv(elementInfo_filename);
    
    for iRow = 1:size(elInfo, 1)
       
        channel_set = [ channel_set eval(elInfo{iRow, 4})];
    end
    
    %create output fname
    
    input_fname_splits = strsplit(input_fname, '/');
    input_fname_no_path = input_fname_splits{end};
    
    input_fname_no_path_split = strsplit(input_fname_no_path, '.');
    input_fname_no_path_no_ext = input_fname_no_path_split{1};
     
    ignore_me_fname = [output_dir '/_ignore_me.txt'];
    output_mda_fname = [output_dir '/' input_fname_no_path_no_ext '.mda'];
    output_chans_csv_fname = [output_dir '/' input_fname_no_path_no_ext '_chans.csv'];
    output_chans_txt_fname = [output_dir '/' input_fname_no_path_no_ext '_used_chans.txt'];
    
    
    %read the data
    nsx = concatOpenNSx(input_fname);
    
    %check if the file is too short to care about
    if ignoreShortFile
       
        if nsx.MetaTags.DataDurationSec/60 < min_duration_data
           
            ignore_me_fid = fopen(ignore_me_fname, 'w');
            fprintf(ignore_me_fid, 'data length less than 5 min ( %0.2f )', nsx.MetaTags.DataDurationSec/60);
            fclose(ignore_me_fid);
        end
    end
    
    %open the channel name recording files
    output_chans_txt_fid = fopen(output_chans_txt_fname, 'w');        
    output_chans_csv_fid = fopen(output_chans_csv_fname, 'w');
   
    %print header to csv
    fprintf(output_chans_csv_fid, '"channel","info"\n');  
   
    
    %calculate data ranges
    ranges = range(nsx.Data, 2);
    
    %filter for channels in elementInfo_micro and channel range
    nsx_channelNames = {nsx.ElectrodesInfo.Label};
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

        fprintf('checking channel %s\n', current_channel_trim);

        %check whether channel was selected and if it has adequate range
        selected_status = 'unselected';
        quality_status = 'good';
        
        if any(iChan == channel_set)          
            selected_status = 'selected';
        end
        
        if ranges(iChan) < min_range_cutoff
            quality_status = 'bad';
        end
        
        if isequal(selected_status, 'selected') && isequal(quality_status, 'good')
            
            nsx_channelNames_match{length(nsx_channelNames_match) + 1} = current_channel_trim;
            nsx_channelNames_matchIdx = [nsx_channelNames_matchIdx iChan];
            
            fprintf(output_chans_txt_fid, sprintf('%s\n', current_channel_trim));
        end
        
        fprintf(sprintf('"%s","%s"\n', current_channel_trim, [selected_status '_' quality_status]));        
        fprintf(output_chans_csv_fid, sprintf('"%s","%s"\n', current_channel_trim, [selected_status '_' quality_status]));
    end
    
    fclose(output_chans_csv_fid);
    fclose(output_chans_txt_fid);
    
    % if no channel pass, write and _ignore_me, else write the data
    if isempty(nsx_channelNames_match)
       
        ignore_me_fid = fopen(ignore_me_fname, 'w');
        fprintf(ignore_me_fid, 'all channels unselected or have data range < 10');
        fclose(ignore_me_fid);
    else
    
        writemda(nsx.Data(nsx_channelNames_matchIdx, :), output_mda_fname, 'int16');
    end
end