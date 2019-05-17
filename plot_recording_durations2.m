
rerun = 0;
save_fpath = '/Users/zaworaca/Dropbox/biowulf_weeds/_sync_folders/_utah_matlab_tools/subj_min_coverage.mat';
% paths = { '/Volumes/72A/UTAH_A/NIH029' };

paths = { ...
          '/Volumes/72A/UTAH_A/NIH029' ...
          '/Volumes/72A/UTAH_A/NIH030' ...
          '/Volumes/72A/UTAH_A/NIH034' ...
          '/Volumes/72A/UTAH_A/NIH036' ...
          '/Volumes/72A/UTAH_A/NIH037' ...
          ...
          '/Volumes/72B/UTAH_B/NIH039' ...
          '/Volumes/72B/UTAH_B/NIH042' ...
          '/Volumes/72B/UTAH_B/NIH046' ...
          '/Volumes/72B/UTAH_B/NIH047' ...
          '/Volumes/72B/UTAH_B/NIH050' ...
          ...
          '/Volumes/72C/UTAH_C/NIH054' ...
          '/Volumes/72C/UTAH_C/NIH056' ...
          '/Volumes/72C/UTAH_C/NIH057' ...
          '/Volumes/72C/UTAH_C/NIH058' ...
          '/Volumes/72C/UTAH_C/NIH059' ...
          '/Volumes/72C/UTAH_C/NIH060' ...
          '/Volumes/72C/UTAH_C/NIH061' ...
          '/Volumes/72C/UTAH_C/NIH062' ...
          '/Volumes/72C/UTAH_C/NIH063' ...
          '/Volumes/72C/UTAH_C/NIH064' ...
          ...
          '/Volumes/72D/UTAH_D/NIH065' ...
          '/Volumes/72D/UTAH_D/NIH066' ...
          '/Volumes/72D/UTAH_D/NIH067' ...
          '/Volumes/72D/UTAH_D/NIH068' ...
          '/Volumes/72D/UTAH_D/NIH069' ...
          '/Volumes/72D/UTAH_D/NIH071' ...
          '/Volumes/72D/UTAH_D/NIH072' ...
          };
    
      
 min_range_cutoff_ungained = 10;
 min_range_cutoff_millivolt = min_range_cutoff_ungained * 0.25 * (1/1000);

 min_duration_minutes = 5;


 min_per_day = 24 * 60;
 num_days_plot = 25;

 total_min = min_per_day * num_days_plot;

 num_subj = length(paths);

 
 if rerun

     subj_dur_min = cell(num_subj, 1);
     subj_min_coverage = cell(num_subj, 1);
    %  subj_min_coverage_skip = cell(num_subj, 1);
     subj_min_coverage_microDevNum = cell(num_subj, 1);

     for iSubj = 1:num_subj

        current_subj_path = [paths{iSubj} '/data_raw'];
        fprintf('collecting %s jacksheets\n', current_subj_path);

        subj_ls = dir(current_subj_path);
        subj_ls_names = {subj_ls.name}; 
        num_subj_ls_names = length(subj_ls_names);

        smallest_datenum = inf;
        smallest_datetime = [];

        % how many devices are there?

        good_sess_jacksheets = {};
        microDevNums = [];

        for iSess = 1:num_subj_ls_names

           sess_path = [current_subj_path '/' subj_ls_names{iSess}];
           sess_jacksheet_fpath = [sess_path '/jacksheetBR_complete.csv'];

           % is there a jacksheet
           if exist(sess_jacksheet_fpath, 'file')

               jacktable = readtable(sess_jacksheet_fpath);
               jacktable_micro = jacktable(jacktable{:,'MicroDevNum'} > 0, :);

               session_rawdir = jacktable{1,'RawDir'};
               session_rawdir_splits = strsplit(session_rawdir{1}, '_');

               session_start_datestr = [session_rawdir_splits{1} '_' session_rawdir_splits{2}];
               session_start_datenum = datenum(session_start_datestr, 'yymmdd_HHMM');
               session_start_datetime = datetime(session_start_datestr, 'InputFormat', 'yyMMdd_HHmm');

               if session_start_datenum < smallest_datenum
                   smallest_datenum = session_start_datenum;
                   smallest_datetime = session_start_datetime;
               end

               % are there micro devices recorded in this session?
               if ~isempty(jacktable_micro)

                   sess_unique_microDevNum = unique(jacktable_micro{:, 'MicroDevNum'});
                   microDevNums = [ microDevNums ; sess_unique_microDevNum];
                   good_sess_jacksheets{length(good_sess_jacksheets) + 1} =  jacktable;
               end    
           end
        end

        unique_microDevNums = unique(microDevNums);
        num_unique_microDevNums = length(unique_microDevNums);

        if num_unique_microDevNums < 1
            continue;
        end

        reference_datetime = datetime(year(smallest_datetime), month(smallest_datetime), day(smallest_datetime));

        subj_dur_min{iSubj} = zeros(num_unique_microDevNums, 1);
        subj_min_coverage{iSubj} = zeros(num_unique_microDevNums, total_min);
    %     subj_min_coverage_skip{subj_min_coverage_idx} = zeros(num_unique_microDevNums, total_min);
        subj_min_coverage_microDevNum{iSubj} = unique_microDevNums;


        % collect start and stop datenums for each session
        fprintf('checking %s jacksheets\n', current_subj_path);

        for iSess = 1:length(good_sess_jacksheets)

            jacktable = good_sess_jacksheets{iSess};

            session_rawdir = jacktable{1,'RawDir'};
            session_rawdir_splits = strsplit(session_rawdir{1}, '_');

            session_start_datestr = [session_rawdir_splits{1} '_' session_rawdir_splits{2}];
            session_start_datetime = datetime(session_start_datestr, 'InputFormat', 'yyMMdd_HHmm');


            % process each device

            for iDevNum = 1:num_unique_microDevNums

                current_dev = unique_microDevNums(iDevNum);
                jacktable_current_dev = jacktable(jacktable{:, 'MicroDevNum'} == current_dev,:);

                % is this dev recorded in this session?
                if isempty(jacktable_current_dev)
                    continue;
                end

                % does this device pass voltage threshold
                if all(jacktable_current_dev{:, 'RangeMilliV'} < min_range_cutoff_millivolt)
                    continue;
                end

                % does this session pass time criteria
                session_duration_min = unique(ceil(jacktable_current_dev{:,'DurationMin'}));

                if session_duration_min < min_duration_minutes
                    continue; 
                end

                % passed filters, this is a device in a session I care about
                start_min_abs = minutes(session_start_datetime - reference_datetime);
                stop_min_abs = minutes((minutes(session_duration_min) + session_start_datetime) - reference_datetime);

                subj_dur_min{iSubj}(iDevNum) = subj_dur_min{iSubj}(iDevNum) + (stop_min_abs - start_min_abs);        
                subj_min_coverage{iSubj}(iDevNum, start_min_abs:stop_min_abs) = 1;

            end

        end

     end

     save(save_fpath, 'paths', 'subj_dur_min', 'subj_min_coverage', 'subj_min_coverage_microDevNum');
 
 else
     
     load(save_fpath);
     
     num_colors = 10;
     colors = distinguishable_colors(num_colors, {'w'});
     
     figure();
     set(gcf,'color','w');
     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01, 0.01, 1, 0.96]);
     
     minute_vals = 1:size(subj_min_coverage{1}, 2);
     
     yticklabels_cell = {};
     yicks_mat = [];
     
     plotted_count = 0;
     
     for iSubj = 1:length(paths)
         
         paths_splits = strsplit(paths{iSubj}, '/');
         subj_name = paths_splits{end};
                  
         % plot the device durations
         
         current_subj_dur_min = subj_dur_min{iSubj};
         
         current_subj_min_coverage = subj_min_coverage{iSubj};
         current_subj_min_coverage_microDevNum = subj_min_coverage_microDevNum{iSubj};
         
         if isempty(current_subj_min_coverage)
            fprintf('skipping %s\n', subj_name);
            continue; 
         end
         
         plotted_count = plotted_count + 1;
         current_subj_min_coverage(current_subj_min_coverage(:) == 1) = plotted_count;
         current_subj_min_coverage(current_subj_min_coverage(:) == 0) = NaN;
         
         for iDev = 1:length(current_subj_min_coverage_microDevNum)
         
            if iDev == 1
                yticklabels_cell{length(yticklabels_cell) + 1} = [sprintf('%0.2f days', max(current_subj_dur_min)/min_per_day) '   ' subj_name];
            else
                yticklabels_cell{length(yticklabels_cell) + 1} = '';
            end
            
            yicks_mat(length(yicks_mat) + 1) = max(current_subj_min_coverage(iDev,:)) - ((iDev-1)*0.25);

            plot(minute_vals, current_subj_min_coverage(iDev,:) - ((iDev-1)*0.25), 'Color', colors(current_subj_min_coverage_microDevNum(iDev), :), 'LineWidth', 7);
            hold on;
         
         end
         
     end
     
     [sorted_yticks_mat, sorted_idx] = sort(yicks_mat);
     
     yticks(sorted_yticks_mat);
     yticklabels(yticklabels_cell(sorted_idx));
     
     ylim([0 plotted_count+1]);
     
     xticks(0:min_per_day:total_min);
     xticklabels(strsplit(num2str(1:(total_min/min_per_day)), ' '));
     
     xlabel('time (days)');
     
     set(gca, 'FontSize', 20);
     
     figure();
     set(gcf,'color','w');

     for iColor = 1:num_colors
        plot(1:10, repmat(iColor, 10, 1), 'Color', colors(iColor,:), 'LineWidth', 7);
        hold on;
     end
     
     title('microDevNum');
     xticklabels({});
     
     set(gca, 'FontSize', 14);
 end 
      
 
 