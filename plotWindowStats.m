function plotWindowStats(raw_stats_fpath, bandpass_stats_fpath, reref_stats_fpath, whiten_stats_fpath, output_fpath)

    raw_stats = load(raw_stats_fpath);
    bandpass_stats = load(bandpass_stats_fpath);
    reref_stats = load(reref_stats_fpath);
    whiten_stats = load(whiten_stats_fpath);
    
    num_channels = size(raw_stats.mean_mats, 1);
    
    subplot_rows = 11;
    subplot_cols = 12;
    
    trace_subplot_cols = [ 1 2 3 4 5 6 7 8 ];
    hist_subplot_cols = [ 9 10 11 12 ];
    
    minute_interval = 5;
    
    for iChan = 1:num_channels

        current_raw_channel = raw_stats.chan_names{iChan};
        
        save_fpath = [output_fpath '_' current_raw_channel '.png'];
        
        figure('Units', 'Normalized', 'OuterPosition', [0 0 1 0.95]);
        set(gcf,'color','w');

        % raw_stats
        
        current_stats = raw_stats;
        
        start_windows_minutes = floor(current_stats.start_windows/(current_stats.samples_per_sec * 60));
        [labels, label_intervals] = findMinuteIntervals(start_windows_minutes, minute_interval);
        
        string_quantiles = strjoin(cellfun(@num2str, num2cell(current_stats.quantile_cutoffs), 'UniformOutput', 0), ' , ');
        [chan_quantile_mat, chan_mean, chan_std] = calcChanMats(current_stats, iChan);   
        
        num_calc_points = length(chan_mean);
        
        xmax = num_calc_points + 10;
        xmin = - 10;
        
        ymax = max(chan_quantile_mat(:)) + 1;
        ymin = min(chan_quantile_mat(:)) - 1;
        
        
        %%%%%%%%%%%%
        % subplot 1 - raw plot quantiles
        %%%%%%%%%%%%
        subplot(subplot_rows, subplot_cols, trace_subplot_cols );
        hold on;
  
        set(gca, 'FontSize', 6);
        
        title(sprintf('%s raw recording - top: quantiles [ %s ] -- bottom: mean +/- std', current_raw_channel, string_quantiles));
        
        for iQuantile = 1:size(chan_quantile_mat, 1)
        
            plot(1:num_calc_points, chan_quantile_mat(iQuantile,:), 'k-');
        end
        
        ylim([ymin ymax]);
        xlim([xmin xmax]);
        xticks(label_intervals);
        xticklabels(labels);
        
        
        %%%%%%%%%%%%
        % subplot 2 - raw plot mean +/- sd
        %%%%%%%%%%%%
        subplot(subplot_rows, subplot_cols, trace_subplot_cols + ( subplot_cols  * 1 ) );
        hold on;
     
        set(gca, 'FontSize', 6);
        
        plot(1:num_calc_points, chan_mean, 'k-');        
        
        upper_std = chan_mean + chan_std;
        lower_std = chan_mean - chan_std;
        
        fill( [ (1:num_calc_points)' ; flipud((1:num_calc_points)')], [ lower_std'; flipud((upper_std)')], 'r', 'linestyle','none');
        alpha(0.3);
        
        ylim([ymin ymax]);
        xlim([xmin xmax]);
        xticks(label_intervals);
        xticklabels(labels);
        
        
        %%%%%%%%%%%%
        % subplot 3 - histogram of voltage
        %%%%%%%%%%%%
        ax = subplot(subplot_rows, subplot_cols, [ hist_subplot_cols(:) ; ( hist_subplot_cols(:) + ( subplot_cols  * 1 )) ]' );
        hold on;
   
        set(gca, 'FontSize', 6);
        
        current_hist_edges = current_stats.hist_edges{iChan};
        current_hist_n = log10(current_stats.hist_n{iChan} + 1);

        x_axis_points = [];

        for iPoint = 1:(length(current_hist_edges)-1)  
            x_axis_points(iPoint) = mean([ current_hist_edges(iPoint) current_hist_edges(iPoint + 1) ]);
        end

        plot(ax, x_axis_points, current_hist_n, 'r-', 'LineWidth', 3);

        xlabel('voltage value');
        ylabel('log count');

        
        %bandpassed
        
        current_stats = bandpass_stats;
        
        current_stats_name_match_idx = find(cellfun(@(x) isequal(x, current_raw_channel), current_stats.chan_names));
        
        if ~isempty(current_stats_name_match_idx)
            
            current_chan_idx = current_stats_name_match_idx(1);

            start_windows_minutes = floor(current_stats.start_windows/(current_stats.samples_per_sec * 60));
            [labels, label_intervals] = findMinuteIntervals(start_windows_minutes, minute_interval);

            string_quantiles = strjoin(cellfun(@num2str, num2cell(current_stats.quantile_cutoffs), 'UniformOutput', 0), ' , ');
            [chan_quantile_mat, chan_mean, chan_std] = calcChanMats(current_stats, current_chan_idx);   

            num_calc_points = length(chan_mean);

            xmax = num_calc_points + 10;
            xmin = - 10;

            ymax = max(chan_quantile_mat(:)) + 1;
            ymin = min(chan_quantile_mat(:)) - 1;


            %%%%%%%%%%%%
            % subplot 4 - bandpassed plot quantiles
            %%%%%%%%%%%%
            subplot(subplot_rows, subplot_cols, trace_subplot_cols + ( subplot_cols  * 3 ));
            hold on;

            set(gca, 'FontSize', 6);
            
            title(sprintf('+bandpass recording - top: quantiles [ %s ] -- bottom: mean +/- std', string_quantiles));

            for iQuantile = 1:size(chan_quantile_mat, 1)

                plot(1:num_calc_points, chan_quantile_mat(iQuantile,:), 'k-');
            end

            ylim([ymin ymax]);
            xlim([xmin xmax]);
            xticks(label_intervals);
            xticklabels(labels);

            %%%%%%%%%%%%
            % subplot 5 - bandpassed plot mean +/- sd
            %%%%%%%%%%%%
            subplot(subplot_rows, subplot_cols, trace_subplot_cols + ( subplot_cols  * 4 ));
            hold on;

            set(gca, 'FontSize', 6);
            
            plot(1:num_calc_points, chan_mean, 'k-');        

            upper_std = chan_mean + chan_std;
            lower_std = chan_mean - chan_std;

            fill( [ (1:num_calc_points)' ; flipud((1:num_calc_points)')], [ lower_std'; flipud((upper_std)')], 'r', 'linestyle','none');
            alpha(0.3);

            ylim([ymin ymax]);
            xlim([xmin xmax]);
            xticks(label_intervals);
            xticklabels(labels);

            %%%%%%%%%%%%
            % subplot 6 - histogram of voltage
            %%%%%%%%%%%%
            ax = subplot(subplot_rows, subplot_cols, [ ( hist_subplot_cols(:) + ( subplot_cols  * 3 )) ; ( hist_subplot_cols(:) + ( subplot_cols  * 4 )) ]');
            hold on;

            set(gca, 'FontSize', 6);
            
            current_hist_edges = current_stats.hist_edges{current_chan_idx};
            current_hist_n = log10(current_stats.hist_n{current_chan_idx} + 1);

            x_axis_points = [];

            for iPoint = 1:(length(current_hist_edges)-1)  
                x_axis_points(iPoint) = mean([ current_hist_edges(iPoint) current_hist_edges(iPoint + 1) ]);
            end

            plot(ax, x_axis_points, current_hist_n, 'r-', 'LineWidth', 3);

            xlabel('voltage value');
            ylabel('log count');
        
        end
        
        % reref
                
        current_stats = reref_stats;
        
        current_stats_name_match_idx = find(cellfun(@(x) isequal(x, current_raw_channel), current_stats.chan_names));
        
        if ~isempty(current_stats_name_match_idx)
            
            current_chan_idx = current_stats_name_match_idx(1);
        
            start_windows_minutes = floor(current_stats.start_windows/(current_stats.samples_per_sec * 60));
            [labels, label_intervals] = findMinuteIntervals(start_windows_minutes, minute_interval);

            string_quantiles = strjoin(cellfun(@num2str, num2cell(current_stats.quantile_cutoffs), 'UniformOutput', 0), ' , ');
            [chan_quantile_mat, chan_mean, chan_std] = calcChanMats(current_stats, current_chan_idx);   

            num_calc_points = length(chan_mean);

            xmax = num_calc_points + 10;
            xmin = - 10;

            ymax = max(chan_quantile_mat(:)) + 1;
            ymin = min(chan_quantile_mat(:)) - 1;


            %%%%%%%%%%%%
            % subplot 7 - reref plot quantiles
            %%%%%%%%%%%%        
            subplot(subplot_rows, subplot_cols, trace_subplot_cols + ( subplot_cols  * 6 ));
            hold on;

            set(gca, 'FontSize', 6);
            
            title(sprintf('+reref recording - top: quantiles [ %s ] -- bottom: mean +/- std', string_quantiles));

            for iQuantile = 1:size(chan_quantile_mat, 1)

                plot(1:num_calc_points, chan_quantile_mat(iQuantile,:), 'k-');
            end

            ylim([ymin ymax]);
            xlim([xmin xmax]);
            xticks(label_intervals);
            xticklabels(labels);


            %%%%%%%%%%%%
            % subplot 8 - refer plot mean +/- sd
            %%%%%%%%%%%%
            subplot(subplot_rows, subplot_cols, trace_subplot_cols + ( subplot_cols  * 7 ));
            hold on;

            set(gca, 'FontSize', 6);
            
            plot(1:num_calc_points, chan_mean, 'k-');        

            upper_std = chan_mean + chan_std;
            lower_std = chan_mean - chan_std;

            fill( [ (1:num_calc_points)' ; flipud((1:num_calc_points)')], [ lower_std'; flipud((upper_std)')], 'r', 'linestyle','none');
            alpha(0.3);

            ylim([ymin ymax]);
            xlim([xmin xmax]);
            xticks(label_intervals);
            xticklabels(labels);


            %%%%%%%%%%%%
            % subplot 9 - histogram of voltage
            %%%%%%%%%%%%
            ax = subplot(subplot_rows, subplot_cols, [ ( hist_subplot_cols(:) + ( subplot_cols  * 6 )) ; ( hist_subplot_cols(:) + ( subplot_cols  * 7 )) ]');
            hold on;

            set(gca, 'FontSize', 6);

            current_hist_edges = current_stats.hist_edges{current_chan_idx};
            current_hist_n = log10(current_stats.hist_n{current_chan_idx} + 1);

            x_axis_points = [];

            for iPoint = 1:(length(current_hist_edges)-1)  
                x_axis_points(iPoint) = mean([ current_hist_edges(iPoint) current_hist_edges(iPoint + 1) ]);
            end

            plot(ax, x_axis_points, current_hist_n, 'r-', 'LineWidth', 3);

            xlabel('voltage value');
            ylabel('log count');

        end        
        
        %whiten
        
        current_stats = whiten_stats;
        
        current_stats_name_match_idx = find(cellfun(@(x) isequal(x, current_raw_channel), current_stats.chan_names));
        
        if ~isempty(current_stats_name_match_idx)
            
            current_chan_idx = current_stats_name_match_idx(1);
        
            start_windows_minutes = floor(current_stats.start_windows/(current_stats.samples_per_sec * 60));
            [labels, label_intervals] = findMinuteIntervals(start_windows_minutes, minute_interval);

            string_quantiles = strjoin(cellfun(@num2str, num2cell(current_stats.quantile_cutoffs), 'UniformOutput', 0), ' , ');
            [chan_quantile_mat, chan_mean, chan_std] = calcChanMats(current_stats, current_chan_idx);   

            num_calc_points = length(chan_mean);

            xmax = num_calc_points + 10;
            xmin = - 10;

            ymax = max(chan_quantile_mat(:)) + 1;
            ymin = min(chan_quantile_mat(:)) - 1;

            %%%%%%%%%%%%
            % subplot 10 - whiten plot quantiles
            %%%%%%%%%%%%     
            subplot(subplot_rows, subplot_cols, trace_subplot_cols + ( subplot_cols  * 9 ));
            hold on;
            set(gca, 'FontSize', 6);

            title(sprintf('+whiten recording - top: quantiles [ %s ] -- bottom: mean +/- std', string_quantiles));

            for iQuantile = 1:size(chan_quantile_mat, 1)

                plot(1:num_calc_points, chan_quantile_mat(iQuantile,:), 'k-');
            end

            ylim([ymin ymax]);
            xlim([xmin xmax]);
            xticks(label_intervals);
            xticklabels(labels);

            %%%%%%%%%%%%
            % subplot 11 - whiten plot mean +/- sd
            %%%%%%%%%%%%            
            subplot(subplot_rows, subplot_cols, trace_subplot_cols + ( subplot_cols  * 10 ));
            hold on;
            set(gca, 'FontSize', 6);

            plot(1:num_calc_points, chan_mean, 'k-');        

            upper_std = chan_mean + chan_std;
            lower_std = chan_mean - chan_std;

            fill( [ (1:num_calc_points)' ; flipud((1:num_calc_points)')], [ lower_std'; flipud((upper_std)')], 'r', 'linestyle','none');
            alpha(0.3);

            ylim([ymin ymax]);
            xlim([xmin xmax]);
            xticks(label_intervals);
            xticklabels(labels);

            ylabel('ungained voltage');
            xlabel('time (min)');


            %%%%%%%%%%%%
            % subplot 12 - voltage histogram
            %%%%%%%%%%%%   
            ax = subplot(subplot_rows, subplot_cols, [ ( hist_subplot_cols(:) + ( subplot_cols  * 9 )) ; ( hist_subplot_cols(:) + ( subplot_cols  * 10 )) ]');
            hold on;

            set(gca, 'FontSize', 6);
    
            current_hist_edges = current_stats.hist_edges{current_chan_idx};
            current_hist_n = log10(current_stats.hist_n{current_chan_idx} + 1);

            x_axis_points = [];

            for iPoint = 1:(length(current_hist_edges)-1)  
                x_axis_points(iPoint) = mean([ current_hist_edges(iPoint) current_hist_edges(iPoint + 1) ]);
            end

            plot(ax, x_axis_points, current_hist_n, 'r-', 'LineWidth', 3);

            xlabel('voltage value');
            ylabel('log count');

        end
        
        tightfig();
        set(gcf, 'InvertHardcopy', 'off','PaperUnits','inches','PaperPosition', [0,0,25,19], 'PaperPositionMode','auto');

        print(gcf, save_fpath, '-dpng','-r200');
        
        fprintf('\t-- saving\n');
        close(gcf);
        
    end
    
end

function [chan_quantile_mat, chan_mean, chan_std] = calcChanMats(current_stats, current_chan)
    
    chan_quantile_mat = zeros(length(current_stats.quantile_cutoffs), size(current_stats.quantile_mats{1}, 2));

    for iQuantile = 1:length(current_stats.quantile_mats)
        chan_quantile_mat(iQuantile, :) = current_stats.quantile_mats{iQuantile}(current_chan, :);
    end

    chan_mean = current_stats.mean_mats(current_chan, :);
    chan_std = current_stats.std_mats(current_chan, :);

end


function [label, label_intervals] = findMinuteIntervals(start_windows_minutes, minute_interval)

    label = {};
    label_intervals = [];
    
    current_minute = 0;
    label{length(label) + 1} = '0';
    label_intervals = [ label_intervals 1 ];
    
    for idx = 1:length(start_windows_minutes)
    
        if start_windows_minutes(idx) - current_minute >= minute_interval
        
            label_intervals = [ label_intervals  idx ];
            label{length(label) + 1} = num2str(start_windows_minutes(idx));
            current_minute = start_windows_minutes(idx);
        end
    
    end

end


