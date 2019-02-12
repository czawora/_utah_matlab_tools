function [returnData, returnChannels] = modifySessionData(varargin)

    fprintf('****************************************************\n');
    fprintf('* %s\n', mfilename);

    fprintf('new log\n');

    p = inputParser;
    p.addParameter('input_type', 'nsx', @ischar);
    p.addParameter('input_filename', '', @ischar);
    
    p.addParameter('output_dir', '', @ischar);
    p.addParameter('split_output', '0', @ischar);
    p.addParameter('chans', 'all', @ischar);
    p.addParameter('writeOut', '1', @ischar);
    
    p.addParameter('ignoreShortFile', '0', @ischar);
    
    parse(p, varargin{:});
    
    disp(p.Results);
    
    input_type = p.Results.input_type;
    input_fname = p.Results.input_filename;
        
    output_dir = p.Results.output_dir;
    split_output = eval(p.Results.split_output); % eval to int
    chans = p.Results.chans;
    writeOut = logical(str2num(p.Results.writeOut));
    
    ignoreShortFile = logical(str2num(p.Results.ignoreShortFile));
    
%     output_dir = '/Users/zaworaca/Desktop/test_output';
%     input_type = 'nsx'; 
%     input_fname = '/Users/zaworaca/dev/biowulf/modifySessionData/NIH037_160125_1705_utah_m_beh.ns6';
%     split_output = 1;
%     chans = 'all';
%     
    %####################################################################################################
    %####################################################################################################
    % checking inputs
    
    
    
    if ~isequal(input_type, 'nsx') && ~isequal(input_type, 'mda')
       error('input_type parameter should be either "nsx" or "mda"\n');
    end
   
    if ~exist(input_fname, 'file')   
        error('%s is not a valid file\n', input_fname);
    end
 
    if ~exist(output_dir, 'dir')   
        error('%s is not a valid directory\n', output_dir);
    end
    
%     if length(subj) < 3 || (length(subj) >=3 && ~isequal(subj(1:3), 'NIH'))
%        error('subject name should be of format NIHXXX ** this is mandatory right now but may not be in the future **'); 
%     end
    
    
    %####################################################################################################
    %####################################################################################################
    % read in input file
    
    subj = regexp(input_fname, 'NIH\d\d\d', 'once', 'match')
    
    if isempty(subj)
       
        error('no NIHXXX subject name in input filename');
    end
    
    if isequal(input_type, 'nsx') 
       
        applyGain = 1;
        nsx_data = concatOpenNSx(input_fname, applyGain);
        
        disp(nsx_data);
        
        electrodeLabels = { nsx_data.ElectrodesInfo.Label };

        session_min = sum(nsx_data.MetaTags.DataDurationSec)/60;
        
        %if file duration is less than 5 min, we should ignore this session
        if ~ignoreShortFile && session_min < 5
           
            ignore_file = [output_dir '/_ignore_me.txt'];
            write_ignore('file duration is less than 5 minutes, too short to care about\n', ignore_file);
    
        end
        
        %get filename of nsx file
        fpath_split = strsplit(input_fname, '/');
        fname = fpath_split{end};
        
        fname_split = strsplit(fname, '.');
        fname_noext = fname_split{1};
        
        %get the start time from 
        input_datetime_split = strsplit(input_fname, '/');
        input_datetime = input_datetime_split{end};
        
        fprintf('input_datetime %s\n', input_datetime);
        
        input_datetime_match = regexp(input_datetime, '\d\d\d\d\d\d_\d\d\d\d', 'match');
        input_starttime = datetime(input_datetime_match, 'InputFormat', 'yyMMdd_HHmm');
        
        input_endtime = input_starttime + minutes(session_min);
        
        data = nsx_data.Data;
        input_duration_min = session_min;
        
        input_chans_internal = 1:(nsx_data.MetaTags.ChannelCount);
        input_chans_external = input_chans_internal;
        
        fprintf('inspecting channels 1\n');
        %check out the channels to write out
        % chans refers to external chan names
        if ~isequal(chans, 'all')
            out_chans = eval(chans);
        else
            out_chans = input_chans_external;
        end


        %determine which channels are bad ones

        chan_info_cell = {};

        output_starttime = datestr(input_starttime, 'yymmdd-HHMM');
        output_endtime = datestr(input_endtime, 'yymmdd-HHMM'); 

        bad_chans_internal = [];

        fprintf('inspecting channels 2\n');
        for iRow = 1:size(data, 1)

            info_string = '';

            iRow_external = input_chans_external(iRow);

            if sum( iRow_external == out_chans ) > 0

                info_string = [ info_string 'selected' ];
            else

                info_string = [ info_string 'unselected' ];
            end

            %identify bad channels
            data_range = range(data(iRow,:));
            data_range_min = 10;

            if data_range < data_range_min

               info_string = [ info_string '_bad' ]; 
               bad_chans_internal = [bad_chans_internal iRow]; 

            else

               info_string = [ info_string '_good' ]; 

            end

            chan_info_cell = [ chan_info_cell ; {electrodeLabels{iRow} info_string} ];
        end
        
        disp(chan_info_cell);

        fprintf('inspecting channels 3\n');
        %if they are all bad, write an ignore me
        if length(bad_chans_internal) == size(data,1)

            ignore_file = [output_dir '/_ignore_me.txt'];
            write_ignore( sprintf('all the channels are bad ( range < data_range_min %d )\n', data_range_min), ignore_file);

        end


        good_selected_chans_internal_logic = [];

        for iRow = 1:size(data, 1)

            iRow_external = input_chans_external(iRow);
            good_selected_chans_internal_logic = [ good_selected_chans_internal_logic (( sum( iRow_external == out_chans ) > 0 ) && ( sum(bad_chans_internal == iRow) == 0)) ];
        end

        if sum( good_selected_chans_internal_logic ) == 0

            ignore_file = [output_dir '/_ignore_me.txt'];
            write_ignore('no non-bad channels were selected\n', ignore_file);

        end

        good_selected_chans_internal = input_chans_internal(logical(good_selected_chans_internal_logic));
        good_selected_chans_external = input_chans_external(logical(good_selected_chans_internal_logic));

        
        fprintf('writing out\n');
        
        if split_output == 1

            if writeOut

                for iRow = 1:size(data, 1)

                    iRow_external = input_chans_external(iRow);

                    if sum( iRow == good_selected_chans_internal ) > 0

                        split_name = [ output_dir '/' subj '_' output_starttime '_' output_endtime '_' sprintf('%03d', iRow_external) ];
                        mkdir(split_name);

                        split_fname = [split_name '/' subj '_' output_starttime '_' output_endtime '_' sprintf('%03d', iRow_external) '.mda'];
                        fprintf('%s\n', split_fname);
                        writemda(data(iRow,:), split_fname, 'float32');

                    end

                end

            end

            returnData = 'split';

        else

            data = data(good_selected_chans_internal, :);

            chan_str = shrink_num_array(good_selected_chans_external);
            chan_str = strrep(chan_str, ' ' , '--');
            chan_str = strrep(chan_str, ':' , '-');

            unsplit_fname = [output_dir '/' subj '_' output_starttime '_' output_endtime '_' chan_str '.mda'];

            fprintf('%s\n', unsplit_fname);

            if writeOut
                writemda(data, unsplit_fname, 'float32');
            end

            returnData = data;

        end

        returnChannels = good_selected_chans_external;

        chan_info = fopen([output_dir '/' fname_noext '_chans.csv'], 'w');
        fprintf(chan_info, '"channel","info"\n'); 

        for iChanInfo = 1:size(chan_info_cell, 1)
            fprintf(chan_info, '"%s","%s"\n', chan_info_cell{iChanInfo, 1}, chan_info_cell{iChanInfo, 2} );
        end

        fclose(chan_info);
        
        save([output_dir '/' fname_noext '_chans.mat'], '-v7.3', 'chan_info_cell');
        
        
    elseif isequal(input_type, 'mda')
        
        fprintf('reading mda\n');
        data = readmda_single(input_fname);
        
        fname = strsplit(input_fname, '/');
        fname_splits = strsplit(fname{end}, '.mda');
        fname_info_splits = strsplit(fname_splits{1}, '_')
                
        input_starttime = datetime(fname_info_splits{2}, 'InputFormat', 'yyMMdd-HHmm');
        input_endtime = datetime(fname_info_splits{3}, 'InputFormat', 'yyMMdd-HHmm');
        
        %input_duration_min = minutes(input_endtime - input_starttime);
        
        %input_chans_internal = 1:size(data, 1);
        input_chans_external_string = [ '[' strrep(strrep(fname_info_splits{4}, '--', ' '), '-', ':') ']' ]
        input_chans_external = eval(input_chans_external_string);
      
        if ~isequal(chans, 'all')
            out_chans = eval(chans);
        else
            out_chans = input_chans_external;
        end
             
        
        if split_output == 1

            if writeOut
                
                for iRow = 1:size(data, 1)

                    iRow_external = input_chans_external(iRow);

                    if sum( iRow_external == out_chans ) > 0

                        split_name = [ output_dir '/' subj '_' fname_info_splits{2} '_' fname_info_splits{3} '_' sprintf('%02d', iRow_external) ];
                        mkdir(split_name);

                        split_fname = [split_name '/' subj '_' fname_info_splits{2} '_' fname_info_splits{3} '_' sprintf('%02d', iRow_external) '.mda_chan'];
                        fprintf('%s\n', split_fname);
                        writemda(data(iRow,:), split_fname, 'float32');

                    end

                end
                
            end

            returnData = 'split';

        else

            data = data(out_chans, :);

            chan_str = shrink_num_array(chans);
            chan_str = strrep(chan_str, ' ' , '--');
            chan_str = strrep(chan_str, ':' , '-');

            unsplit_fname = [output_dir '/' subj '_' input_starttime '_' input_endtime '_' chan_str '.mda_chan'];

            fprintf('%s\n', unsplit_fname);

            if writeOut
                writemda(data, unsplit_fname, 'float32');
            end
            
            returnData = data;

        end
        
        returnChannels = out_chans;
        
    end
         
    fprintf('modifySesionData -- done');
        
end


function write_ignore(message, ignore_fname)
    
    ignore_fid = fopen(ignore_fname, 'w');
            
    fprintf(ignore_fid, '%s', message);
            
    fclose(ignore_fid);
    
    error('ignore condition for this session');
end
 
