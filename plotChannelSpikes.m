function plotChannelSpikes(varargin)

      fprintf('starting plotChannelSpikes %s\n', datetime('now'));

        %INPUTS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        samplingRateHz = 30000;
        
        p = inputParser;
        
        p.addParameter('session_name', '', @ischar);
        p.addParameter('channel_name', '', @ischar);
        
        p.addParameter('clip_features_fpath', '', @ischar);
        p.addParameter('clips_fpath', '', @ischar);
        p.addParameter('firings_fpath', '', @ischar);
        p.addParameter('isol_metrics_fpath', '', @ischar);
        p.addParameter('isol_pair_metrics_fpath', '', @ischar);
        p.addParameter('metrics_fpath', '', @ischar);
        p.addParameter('mda_fpath', '', @ischar);
        
        p.addParameter('clip_features', []);
        p.addParameter('clips', []);
        p.addParameter('firings', []);
        p.addParameter('isol_metrics', []);
        p.addParameter('isol_pair_metrics', []);
        p.addParameter('metrics', []);
        p.addParameter('mda', []);
        
        p.addParameter('unit_names_fpath', '');
        p.addParameter('good_units_fpath', '')
        
        p.addParameter('saveDir', '', @ischar);
        
        p.addParameter('removeLargeAmpUnits', '1', @ischar);
        
        p.addParameter('snr_min', '1', @ischar);
        p.addParameter('iso_min', '0.95', @ischar);
        p.addParameter('noise_overlap_max', '0.1', @ischar);
        
        parse(p, varargin{:});
        
        disp(p.Results);
        
        session_name = p.Results.session_name;
        current_chan_underscore = p.Results.channel_name;
        
        snr_filt = str2num(p.Results.snr_min);
        iso_filt = str2num(p.Results.iso_min);
        over_filt = str2num(p.Results.noise_overlap_max);
        
        removeLargeAmpUnits = str2num(p.Results.removeLargeAmpUnits);
        
        clip_features_fpath = p.Results.clip_features_fpath ;
        clips_fpath = p.Results.clips_fpath ;
        firings_fpath = p.Results.firings_fpath ;
        isol_metrics_fpath = p.Results.isol_metrics_fpath ;
        isol_pair_metrics_fpath = p.Results.isol_pair_metrics_fpath ;
        metrics_fpath = p.Results.metrics_fpath ;
        mda_fpath = p.Results.mda_fpath ;

        clip_features = p.Results.clip_features ;
        clips = p.Results.clips ;
        firings = p.Results.firings ;
        isol_metrics = p.Results.isol_metrics ;
        isol_pair_metrics = p.Results.isol_pair_metrics ;
        metrics = p.Results.metrics ;
        mda = p.Results.mda ;
        
        unit_names_fpath = p.Results.unit_names_fpath;
        good_units_fpath = p.Results.good_units_fpath;
        
        saveDir = p.Results.saveDir;
                
        
        session_name_underscore = session_name;
        session_name = strrep(session_name, '_', ' ');
        
        current_chan = strrep(current_chan_underscore, '_', ' ');
        
        if isempty(clip_features)
        
            if ~exist(clip_features_fpath, 'file')
               error('%s is not a valid file', clip_features_fpath);
            end

            clip_features = readmda(clip_features_fpath)';
        end
        
        if isempty(clips)

            if ~exist(clips_fpath, 'file')
               error('%s is not a valid file', clips_fpath);
            end

            clips = squeeze(readmda(clips_fpath))';

        end
        
        if isempty(firings)
        
            if ~exist(firings_fpath, 'file')
               error('%s is not a valid file', firings_fpath);
            end

            firings = readmda(firings_fpath);
        
        end
        
        if isempty(isol_metrics)
        
            if ~exist(isol_metrics_fpath, 'file')
               error('%s is not a valid file', isol_metrics_fpath);
            end

            isol_metrics = readjson(isol_metrics_fpath);

        end
        
        if length(isol_metrics) > 1
            if isempty(isol_pair_metrics)

                if ~exist(isol_pair_metrics_fpath, 'file')
                   error('%s is not a valid file', isol_pair_metrics_fpath);
                end

                isol_pair_metrics = readjson(isol_pair_metrics_fpath);

            end
        end
        
        if isempty(metrics)
        
            if ~exist(metrics_fpath, 'file')
               error('%s is not a valid file', metrics_fpath);
            end

            metrics = readjson(metrics_fpath);
        
        end
        
        if isempty(mda)
        
            if ~exist(mda_fpath, 'file')
               error('%s is not a valid file', mda_fpath);
            end

            mda = readmda(mda_fpath);
        
        end
        
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
        end
   
        
        if ~exist(good_units_fpath, 'file')
           error('%s is not a valid file', good_units_fpath);
        end
    
        load(good_units_fpath);
        
        if ~exist(unit_names_fpath, 'file')
           error('%s is not a valid file', unit_names_fpath);
        end
    
        load(unit_names_fpath);
        
        if length(good_units_filt) ~= length(unique(firings(3,:)))
           error('length(good_units_filt) ~= length(unique(firings(3,:)))'); 
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        example_clip_size = samplingRateHz * 3; %3 seconds of data
        num_example_clips = 2;
        
        example_clips = zeros(num_example_clips, example_clip_size);
        
        length_mda = length(mda);
    
        fprintf('getting example clips\n');
        %find the example  clips
        for b = 1:num_example_clips
        
            %fprintf('++++++ %d\n', b);

            min_point = floor((length_mda/num_example_clips)) * (b-1) + 1;
            max_point = floor((length_mda/num_example_clips)) * b;

            %fprintf('min point %d and max point %d\n', min_point, max_point);

            candidate_sample = zeros(1,example_clip_size);
            candidate_spike_count = 0;

            clip_count = min_point;

            while (clip_count + example_clip_size -1) < max_point

%                 fprintf('from %d to %d\n', clip_count, clip_count + example_clip_size);

                spikes_in_range = sum(...
                    (firings(2,:) > clip_count)...
                    &...
                    (firings(2,:) < (clip_count + example_clip_size))...
                    );

                if spikes_in_range > candidate_spike_count

                    candidate_spike_count = spikes_in_range;
                    candidate_sample = mda(clip_count:(clip_count + example_clip_size -1));

                end

                clip_count = clip_count + example_clip_size;

            end

            example_clips(b, :) = candidate_sample;

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        unique_units = unique(firings(3,:));
        
        good_units = [];
        good_units_names = {};

        noise_units = [];
        noise_units_names = {};
        
        for iFilt = 1:length(unique_units)
           
            if good_units_filt(iFilt) == 1
               
                good_units(length(good_units) + 1) = unique_units(iFilt);
                good_units_names{length(good_units_names) + 1} = unit_names{length(good_units)};
            else
                
                noise_units(length(noise_units) + 1) = unique_units(iFilt);
                noise_units_names{length(noise_units_names) + 1} = ['noise' num2str(length(noise_units_names) + 1)];
            end
            
        end
        
        total_unit_num = length(unique_units);
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %plotting settings and layout
        num_unique_colors = total_unit_num;
        unit_colors = distinguishable_colors( num_unique_colors , {'w', 'k'});
        nz_markers_set = ['.' '+' 'o' '*' 'x' 's' 'd' '^' '>' '<' 'h' 'p'];

        on_panels = [1 ... % avg waveforms
            1 ... % uncolored PCA
            1 ... % uncolored time-voltage
            1 ... % ISI metrics
            1 ... % colored PCA
            1 ... % colored time-voltage
            1 ... % metrics
            1 ... % example traces
            1 ... % spike rates;
            ];

%         full_subplot_width_scalar = 0.9;
%         full_subplot_height_scalar = 0.9;
% 
%         unit_hist_width_scalar = 0.9;
%         unit_hist_height_scalar = 0.9;
        full_subplot_width_scalar = 1;
        full_subplot_height_scalar = 1;

        unit_hist_width_scalar = 1;
        unit_hist_height_scalar = 1;

        %for each channel, produce the figures

        wf_shading_opacity = 0.3;

        font_size = 6;
        pca_msize = 2;
        time_volt_msize = 2;
        
        pca_sd_mult = 10;
        time_volt_sd_mult = 20;

        display_rows = 4;
        display_cols = 3;
        
        %unit ISI histograms need to plotted as seperate image
        if total_unit_num > 3

            subplot_rows_per_display_row = 2;
            subplot_cols_per_display_col = 1;
            
            figure(2); % declare the extra ISI plot
            set(gcf,'visible','off');

        %otherwise they can all fit    
        else
            
            subplot_rows_per_display_row = 2;
            subplot_cols_per_display_col = total_unit_num;

        end
        
        subplot_rows = display_rows * subplot_rows_per_display_row;
        subplot_cols = display_cols * subplot_cols_per_display_col;
        
        plot_layout = reshape(1:subplot_rows * subplot_cols, subplot_cols, subplot_rows)'
        subplot_sections = cell(display_rows, display_cols);

        for i = 1:size(plot_layout,1)
            for j = 1:size(plot_layout,2)

                current_display_row = ceil( i / subplot_rows_per_display_row);
                current_display_col = ceil( j / subplot_cols_per_display_col);

                subplot_sections{current_display_row, current_display_col} = [subplot_sections{current_display_row, current_display_col}  plot_layout(i,j)];
            end
        end
        
        fprintf('chan has %d units\n', total_unit_num);

        figure(1);%set figure 1 as gcf
        set(gcf,'color','w');
        set(gcf,'visible','off');
       
        %create unit_colors list to fraw from
        nz_markers = repmat(nz_markers_set, 1, ceil(length(noise_units)/length(nz_markers_set)));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%FIGURE PLOTTING%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %waveform figure - top left ( fine to run if unit contains only single spike)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('waveform fig %s\n', datetime('now'));
        current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{1,1}, 'replace');
        hold on;

        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01, 0.01, 1, 0.96]);

        s = get(gca, 'Position');
        set(gca, 'Position', [s(1), s(2), s(3)* full_subplot_width_scalar, s(4)* full_subplot_height_scalar ]);

        if on_panels(1) == 1

            %fprintf('\t-- plotting real panel %d for chan %s\n', 1, current_chan);

            ylabel('whitened units');
            %set bg color
            set(gca,'Color','k');
            set(gca, 'FontSize', font_size);

            min_y = 0;
            max_y = 0;

            legend_cell = {};
            legend_handles = [];

            %plot noise as white

            for n = 1:length(noise_units)

                current_nz = noise_units(n);
                %fprintf('\t\t-- plotting real panel %d for noise %d\n', 1, n);

                nzIDs_current_nz_condi = (firings(3,:) == current_nz);

                avg_wf = mean(clips(nzIDs_current_nz_condi,:), 1);
                sd_wf = std(clips(nzIDs_current_nz_condi,:),0, 1);

                if size(avg_wf, 1) < 2
                    avg_wf = avg_wf';
                end

                if size(sd_wf, 1) < 2
                    sd_wf = sd_wf';
                end

                if min(avg_wf) < min_y
                    min_y = min(avg_wf);
                end

                if max(avg_wf) > max_y
                    max_y = max(avg_wf);
                end

                fill( ...
                    [ (1:length(avg_wf))' ; flipud((1:length(avg_wf))')], ...
                    [ avg_wf - sd_wf ; flipud(avg_wf + sd_wf)],...
                    'w', 'linestyle','none');
                alpha(wf_shading_opacity);

                l = line( 1:1:length(avg_wf) , avg_wf);
                set(l, 'Color', 'k');

                legend_cell{length(legend_cell) + 1} = [noise_units_names{n} ' ( ' nz_markers( mod(n-1, length(nz_markers)) + 1 ) ' )'];
                legend_handles(length(legend_handles) + 1) = l;

            end
            %now plot the unit waveforms

            for u = 1:length(good_units)

                %fprintf('\t\t-- plotting real panel %d for unit %d\n', 1, u);


                current_unit = good_units(u);

                u_color = unit_colors( u, : );

                unitIDs_current_unit_condi = ( firings(3,:) == current_unit );

                avg_wf = mean(clips(unitIDs_current_unit_condi,:), 1);
                sd_wf = std(clips(unitIDs_current_unit_condi,:),0, 1);

                if size(avg_wf, 1) < 2
                    avg_wf = avg_wf';
                end

                if size(sd_wf, 1) < 2
                    sd_wf = sd_wf';
                end

                if min(avg_wf) < min_y
                    min_y = min(avg_wf);
                end

                if max(avg_wf) > max_y
                    max_y = max(avg_wf);
                end

                %plotting
                %plot std fill
                fill( ...
                    [ (1:length(avg_wf))' ; flipud((1:length(avg_wf))')], ...
                    [ avg_wf - sd_wf ; flipud(avg_wf + sd_wf)],...
                    u_color, 'linestyle','none');
                alpha(wf_shading_opacity);

                l = line( 1:size(clips,2) , avg_wf);
                set(l, 'Color', u_color);

                legend_cell{length(legend_cell) + 1} = good_units_names{u};
                legend_handles(length(legend_handles) + 1) = l;

                %fprintf('\t-- plotting\n');

            end

            %set axes
            y_axis_offset = 50;
            xmin = 1;
            xmax = size(clips,2);
            ymin = min_y - y_axis_offset;
            ymax = max_y + y_axis_offset;
            axis(current_plot, [ xmin xmax ymin ymax ]);

            title( [session_name ' ' current_chan] );
            
        else
            %fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 1, current_chan);
            plot(1:10, 1:10);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PCA figures ( fine to run if unit contains only single spike)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('plotting PCA %s\n', datetime('now'));

        %uncolored PCA
        hold off;
        current_plot = subplot(subplot_rows, subplot_cols, subplot_sections{1,2},'replace');
        hold on;

        s = get(gca, 'Position');
        set(gca, 'Position', [s(1), s(2), s(3)* full_subplot_width_scalar, s(4)* full_subplot_height_scalar ]);

        if on_panels(2) == 1

            %fprintf('\t-- plotting real panel %d for chan %s\n', 2, current_chan);

            xlabel('PC 1');
            ylabel('PC 2');
            %set bg color
            set(gca,'Color','k');
            set(gca, 'FontSize', font_size);

            for n = 1:length(noise_units)

                %fprintf('\t\t-- plotting real panel %d for noise %d\n', 2, n);

                current_nz = noise_units(n);
                
                nzIDs_current_nz_condi = ( firings(3,:) == current_nz );

                 %get waveform PC values
                pc1 = clip_features(nzIDs_current_nz_condi,1);
                pc2 = clip_features(nzIDs_current_nz_condi,2);

                plot(pc1, pc2, 'w.',  'MarkerSize', pca_msize);

            end

            for u = 1:length(good_units)

                %fprintf('\t\t-- plotting real panel %d for unit %d\n', 2, u);

                current_unit = good_units(u);
                unitIDs_current_unit_condi = ( firings(3,:) == current_unit );

                %keyboard;
                %get waveform PC values
                pc1 = clip_features(unitIDs_current_unit_condi,1);
                pc2 = clip_features(unitIDs_current_unit_condi,2);

                plot(pc1, pc2, 'w.',  'MarkerSize', pca_msize);

            end
            
            stdevs = std(clip_features);
            avgs = mean(clip_features);
            
            x_low_bound = stdevs(1) * -1 * pca_sd_mult + avgs(1);
            x_hi_bound = stdevs(1) * pca_sd_mult + avgs(1);
            
            y_low_bound = stdevs(2) * -1 * pca_sd_mult + avgs(2);
            y_hi_bound = stdevs(2) * pca_sd_mult + avgs(2);
            
            xlim([x_low_bound x_hi_bound]);
            ylim([y_low_bound y_hi_bound]);

        else
            %fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 2, current_chan);

            plot(1:10, 1:10);
        end

        %colored PCA ( fine to run if unit contains only single spike)
        hold off;
        current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{2,2} ,'replace');
        hold on;

        s = get(gca, 'Position');
        set(gca, 'Position', [s(1), s(2), s(3)* full_subplot_width_scalar, s(4)* full_subplot_height_scalar ]);


        if on_panels(5) == 1

            %fprintf('\t-- plotting real panel %d for chan %s\n', 5, current_chan);

            xlabel('PC 1');
            ylabel('PC 2');
            %set bg color
            set(gca,'Color','k');
            set(gca, 'FontSize', font_size);

            for n = 1:length(noise_units)

                %fprintf('\t\t-- plotting real panel %d for chan %s noise %d\n', 5, current_chan, n);

                n_marker = nz_markers( n );
                current_nz = noise_units(n);
                
                nzIDs_current_nz_condi = (firings(3,:) == current_nz);

                 %get waveform PC values
                pc1 = clip_features(nzIDs_current_nz_condi,1);
                pc2 = clip_features(nzIDs_current_nz_condi,2);

                plot(pc1, pc2, ['w' n_marker],  'MarkerSize', pca_msize);

            end

            for u = 1:length(good_units)

                %fprintf('\t\t-- plotting real panel %d for chan %s unit %d\n', 5, current_chan, u);

                u_color = unit_colors( u, : );
                current_unit = good_units(u);
                
                unitIDs_current_unit_condi = ( firings(3,:) == current_unit);

                %keyboard;
                %get waveform PC values
                pc1 = clip_features(unitIDs_current_unit_condi,1);
                pc2 = clip_features(unitIDs_current_unit_condi,2);

                plot(pc1, pc2, '.', 'MarkerSize', pca_msize, 'MarkerFaceColor', u_color, 'MarkerEdgeColor', u_color);

            end
            
            stdevs = std(clip_features);
            avgs = mean(clip_features);
            
            x_low_bound = stdevs(1) * -1 * pca_sd_mult + avgs(1);
            x_hi_bound = stdevs(1) * pca_sd_mult + avgs(1);
            
            y_low_bound = stdevs(2) * -1 * pca_sd_mult + avgs(2);
            y_hi_bound = stdevs(2) * pca_sd_mult + avgs(2);
            
            xlim([x_low_bound x_hi_bound]);
            ylim([y_low_bound y_hi_bound]);
            
        else
            %fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 5, current_chan);

            plot(1:10, 1:10);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %time-voltage figures ( fine to run if unit contains only single spike )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('plotting time-voltage %s\n', datetime('now'));

        hold off;
        current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{1,3} ,'replace');
        hold on;

        s = get(gca, 'Position');        
        set(gca, 'Position', [s(1), s(2), s(3)* full_subplot_width_scalar, s(4)* full_subplot_height_scalar ]);

        if on_panels(3) == 1
            
            %fprintf('\t-- plotting real panel %d for chan %s\n', 3, current_chan);

            xlabel('time (min)');
            ylabel('peak-to-trough');
            %set bg color
            set(gca,'Color','k');
            set(gca, 'FontSize', font_size);

            xmx = 0;
            
            all_pk2tr = [];
           
            for n = 1:length(noise_units)

                %fprintf('\t\t-- plotting real panel %d for chan %s noise %d\n', 3, current_chan, n);

                current_nz = noise_units(n);
                nzIDs_current_nz_condi = (firings(3,:) == current_nz);

                wv = clips(nzIDs_current_nz_condi,:);
                timestamps = firings(2, nzIDs_current_nz_condi)/30000/60;

                time_min = min(timestamps);
                time_max = max(timestamps);
                
                if time_max > xmx
                    xmx = time_max;
                end

                pk2tr = zeros(size(wv,1),1); 

                for iwv = 1:size(wv,1)

                    current_pk2tr = abs( max(wv(iwv,:)) - min(wv(iwv,:)) );

                    pk2tr(iwv) = current_pk2tr;

                end

                all_pk2tr = [all_pk2tr; pk2tr];
                plot(timestamps, pk2tr, 'w.', 'MarkerSize', time_volt_msize);

            end
            
            for u = 1:length(good_units)

                %fprintf('\t\t-- plotting real panel %d for chan %s unit %d\n', 3, current_chan, u);

                current_unit = good_units(u);
                unitIDs_current_unit_condi = ( firings(3,:) == current_unit);

                wv = clips(unitIDs_current_unit_condi,:);
                timestamps = firings(2, unitIDs_current_unit_condi)/30000/60;

                time_min = min(timestamps);
                time_max = max(timestamps);

                if time_max > xmx
                    xmx = time_max;
                end
                
                pk2tr = zeros(size(wv,1),1);

                for iwv = 1:size(wv,1)

                    current_pk2tr = abs( max(wv(iwv,:)) - min(wv(iwv,:)) );

                    pk2tr(iwv) = current_pk2tr;

                end

                all_pk2tr = [all_pk2tr; pk2tr];
                plot(timestamps, pk2tr, 'w.', 'MarkerSize', time_volt_msize);

            end

            axis(current_plot, [0 time_max + 1 min(all_pk2tr) max(all_pk2tr) + (max(all_pk2tr)/10)]);
            
            stdevs = std(all_pk2tr);
            avgs = mean(all_pk2tr);
            
            y_low_bound = stdevs * -1 * time_volt_sd_mult + avgs;
            y_hi_bound = stdevs * time_volt_sd_mult + avgs;
                        
            ylim([y_low_bound y_hi_bound]);

        else
            %fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 3, current_chan);

            plot(1:10, 1:10);
        end

        hold off;
        current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{2,3} ,'replace');
        hold on;
        
        s = get(gca, 'Position');
        set(gca, 'Position', [s(1), s(2), s(3)* full_subplot_width_scalar, s(4)* full_subplot_height_scalar ]);

        if on_panels(6) == 1
            %fprintf('\t-- plotting real panel %d for chan %s\n', 6, current_chan);

            xlabel('time (min)');
            ylabel('peak-to-trough');
            %set bg color
            set(gca,'Color','k');
            set(gca, 'FontSize', font_size);

            xmx = 0;

            all_pk2tr = [];
            
            for n = 1:length(noise_units)

                %fprintf('\t\t-- plotting real panel %d for chan %s noise %d\n', 6, current_chan, n);

                n_marker = nz_markers( mod(n-1, length(nz_markers)) + 1 );
                current_nz = noise_units(n);
                nzIDs_current_nz_condi = (firings(3,:) == current_nz);

                wv = clips(nzIDs_current_nz_condi,:);
                timestamps = firings(2, nzIDs_current_nz_condi)/30000/60;

                time_min = min(timestamps);
                time_max = max(timestamps);
                
                if time_max > xmx
                    xmx = time_max;
                end

                pk2tr = zeros(size(wv,1),1); 

                for iwv = 1:size(wv,1)

                    current_pk2tr = abs( max(wv(iwv,:)) - min(wv(iwv,:)) );

                    pk2tr(iwv) = current_pk2tr;

                end
                
                all_pk2tr = [all_pk2tr; pk2tr];
                plot(timestamps, pk2tr, ['w' n_marker], 'MarkerSize', time_volt_msize);

            end

            for u = 1:length(good_units)

                %fprintf('\t\t-- plotting real panel %d for chan %s unit %d\n', 6, current_chan, u);

                u_color = unit_colors( u , :);
                
                current_unit = good_units(u);
                unitIDs_current_unit_condi = ( firings(3,:) == current_unit);

                wv = clips(unitIDs_current_unit_condi,:);
                timestamps = firings(2, unitIDs_current_unit_condi)/30000/60;
                
                time_min = min(timestamps);
                time_max = max(timestamps);

                if time_max > xmx
                    xmx = time_max;
                end
                
                pk2tr = zeros(size(wv,1),1);

                for n = 1:size(wv,1)            

                    current_pk2tr = abs( max(wv(n,:)) - min(wv(n,:)) );

                    pk2tr(n) = current_pk2tr;            
                end
                
                all_pk2tr = [all_pk2tr; pk2tr];
                plot(timestamps, pk2tr, '.', 'MarkerSize', time_volt_msize, 'MarkerFaceColor', u_color, 'MarkerEdgeColor', u_color);

            end

            axis(current_plot, [0 time_max + 1 min(all_pk2tr) max(all_pk2tr) + (max(all_pk2tr)/10)]);
            
            stdevs = std(all_pk2tr);
            avgs = mean(all_pk2tr);
            
            y_low_bound = stdevs * -1 * time_volt_sd_mult + avgs;
            y_hi_bound = stdevs * time_volt_sd_mult + avgs;
                        
            ylim([y_low_bound y_hi_bound]);

        else
            %fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 6, current_chan);

            plot(1:10, 1:10);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %spike ISI distribution figures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('plot spike ISI %s\n', datetime('now'));

        if on_panels(4) == 1

            if length(noise_units) + length(good_units) < 10
                %a seperate figure will be made for these hist ISI plots
                if length(noise_units) + length(good_units) > 3

                    %fprintf('\t-- plotting SEPERATE panel %d for chan %s\n', 4, current_chan);

                    hold off;
                    current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{2,1} ,'replace');
                    hold on;

                    s = get(gca, 'Position');
                    set(gca, 'Position', [s(1), s(2), s(3)*unit_hist_width_scalar, s(4)*unit_hist_height_scalar ])

                    plot(1:10, zeros(1,10));
                    ylim([0 10]);

                    txt = '> 3 units, see extra figure';
                    text(5, 5, txt, 'FontSize', 14, 'HorizontalAlignment', 'center');

                    if total_unit_num < 10
                    %passing all the necessary items to make a plot
                        plotUnitHist(noise_units, good_units, unit_colors, firings, unit_hist_width_scalar, unit_hist_height_scalar, [saveDir '/' session_name_underscore '_' current_chan '_hist.png']);
                    end
                    figure(1);%set gcf back to main figure

                %no seperate figure needed, plot the unit ISI histograms as usual    
                else

                    %fprintf('\t-- plotting real panel %d for chan %s\n', 4, current_chan);
                    %fprintf('\t-- using plot points %d\n', subplot_sections{2,1});

                    mini_col_count = 1;

                     for n = 1:length(noise_units)

                        current_nz = noise_units(n);
                        nzIDs_current_nz_condi = (firings(3,:) == current_nz);     

                        %dont make image if there are barely any spikes
                        if sum(nzIDs_current_nz_condi) > 10

                            ms_diffs = diff(firings(2, nzIDs_current_nz_condi)/30);

                            %histogram of spike times ISI for this unit
                            hold off;
                            current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{2,1}(mini_col_count) ,'replace');
                            %fprintf('making hist for noise %d at plot point %d\n', n, subplot_sections{2,1}(mini_col_count));
                            hold on;

                             s = get(gca, 'Position');
                             set(gca, 'Position', [s(1), s(2), s(3)*unit_hist_width_scalar, s(4)*unit_hist_height_scalar ])

                             if n == 1
                                ylabel('count')
    %                             xlabel('ISI (ms)');
                             end

                             set(gca,'Color','k');
                             set(gca, 'FontSize', font_size);

                             histogram(current_plot, ms_diffs, 'FaceColor', 'w', 'EdgeColor', 'w', 'BinWidth', 1);
                             xlim([-5 inf]);

                             ax = gca;
                             expon = max([0 floor(min(log10(ax.YAxis.TickValues)))]);
                             ax.YAxis.Exponent = expon;

                             hold off;
                             current_plot = subplot(subplot_rows, subplot_cols, subplot_sections{2,1}(mini_col_count + subplot_cols_per_display_col) ,'replace');
                             hold on;

                             %fprintf('making hist for noise %d at plot point %d\n', n, subplot_sections{2,1}(mini_col_count + subplot_cols_per_display_col));

                             s = get(gca, 'Position');
                             set(gca, 'Position', [s(1), s(2), s(3)*unit_hist_width_scalar, s(4)*unit_hist_height_scalar ]);

                             if n == 1
                                ylabel('cumulative proportion');
                                xlabel('ISI (ms)');
                             end

                             set(gca,'Color','w');
                             set(gca, 'FontSize', font_size);

                             num_bins = ceil(max(ms_diffs));
                             hist_res = hist(ms_diffs, num_bins);

                             cum_dist = zeros(num_bins, 1);
                             cum_dist(1) = hist_res(1);            

                             for a = 2:num_bins
                                cum_dist(a) = cum_dist(a-1) + hist_res(a);
                             end
                             cum_dist = cum_dist./sum(hist_res);

                             txt = sprintf('3 ms (%0.2f)', sum(ms_diffs < 3)/length(ms_diffs));

                             plot(1:length(cum_dist), cum_dist, 'b');
                             ylim([0 1]);   
                             xlim([0 ceil(max(ms_diffs))]);

                             ax = gca;
                             expon = max([0 floor(min(log10(ax.YAxis.TickValues)))]);
                             ax.YAxis.Exponent = expon;

                             title(txt);

                        end
                        mini_col_count = mini_col_count + 1;

                     end

                     for u = 1:length(good_units)

                         u_color = unit_colors( u, : );
                         current_unit = good_units(u);
                         unitIDs_current_unit_condi = ( firings(3,:) == current_unit );

                         %dont make image if there are barely any spikes
                         if sum(unitIDs_current_unit_condi) > 10

                             ms_diffs = diff(firings(2,unitIDs_current_unit_condi)/30);

                             %histogram of spike times ISI for this unit
                             hold off;
                             current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{2,1}(mini_col_count) ,'replace');
                             %fprintf('making hist for unit %d at plot point %d\n', u, subplot_sections{2,1}(mini_col_count));
                             hold on;

                             s = get(gca, 'Position');
                             set(gca, 'Position', [s(1), s(2), s(3)*unit_hist_width_scalar, s(4)*unit_hist_height_scalar])

                             if n == 1
                                ylabel('count')
    %                             xlabel('ISI (ms)');
                             end

                             set(gca,'Color','k');
                             set(gca, 'FontSize', font_size);

                             histogram(current_plot, ms_diffs, 'FaceColor', u_color, 'EdgeColor', u_color, 'BinWidth', 1);
                             xlim([-5 inf]);

                             ax = gca;
                             expon = max([0 floor(min(log10(ax.YAxis.TickValues)))]);
                             ax.YAxis.Exponent = expon;

                             %then cdf of ISI for this unit

                             hold off;
                             current_plot = subplot(subplot_rows, subplot_cols, subplot_sections{2,1}(mini_col_count + subplot_cols_per_display_col) ,'replace');
                             hold on;

                             %fprintf('making hist for unit %d at plot point %d\n', u, subplot_sections{2,1}(mini_col_count + subplot_cols_per_display_col));

                             s = get(gca, 'Position');
                             set(gca, 'Position', [s(1), s(2),s(3)*unit_hist_width_scalar, s(4)*unit_hist_height_scalar ])

                             if n == 1
                                ylabel('cumulative proportion');
                                xlabel('ISI (ms)');
                             end

                             set(gca,'Color','w');
                             set(gca, 'FontSize', font_size);

                             num_bins = ceil(max(ms_diffs));
                             hist_res = hist(ms_diffs, num_bins);

                             cum_dist = zeros(num_bins, 1);
                             cum_dist(1) = hist_res(1);            

                             for a = 2:num_bins
                                cum_dist(a) = cum_dist(a-1) + hist_res(a);
                             end
                             cum_dist = cum_dist./sum(hist_res);

                             txt = sprintf('3 ms (%0.2f)', sum(ms_diffs < 3)/length(ms_diffs));

                             plot(1:length(cum_dist), cum_dist, 'b');
                             ylim([0 1]);   
                             xlim([0 ceil(max(ms_diffs))]);

                             ax = gca;
                             expon = max([0 floor(min(log10(ax.YAxis.TickValues)))]);
                             ax.YAxis.Exponent = expon;

                             title(txt);

                         end

                         mini_col_count = mini_col_count + 1;

                     end
                end
            end
        else

            %fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 4, current_chan);

            hold off;
            subplot(subplot_rows,subplot_cols, subplot_sections{2,1} ,'replace');
            hold on;

            plot(1:10, 1:10);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %metrics figures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('plot metrics %s\n', datetime('now'));

        hold off;
        current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{3,1},'replace');
        hold on;

        s = get(gca, 'Position');
        set(gca, 'Position', [s(1), s(2), s(3)* full_subplot_width_scalar, s(4)* full_subplot_height_scalar ]);
        %set(gca, 'FontSize', font_size);

        if on_panels(7) == 1
            %fprintf('\t-- plotting real panel %d for chan %s\n', 7, current_chan);
            
            %table_font = ;
            
            %will be known

            column_labels = {'' '' 'spikes' 'SNR' 'noise-overlap' 'isolation'};
            filter_labels = {'' '' '' num2str(snr_filt) num2str(over_filt) num2str(iso_filt)};
            unit_names = legend_cell; %get names from waveform plot
            lng_colors = zeros(length(unit_names) , 3);
            spikes = [];
            snr = [];
            iso = [];
            over = [];
            
            for n = 1:length(noise_units)

                current_nz = noise_units(n);
                
                lng_colors(n, :) = [1 1 1];
                spikes(length(spikes) + 1) = metrics.clusters(current_nz).metrics.num_events;
                snr(length(snr) + 1) = isol_metrics.clusters(current_nz).metrics.peak_snr;
                iso(length(iso) + 1) = isol_metrics.clusters(current_nz).metrics.isolation;
                over(length(over) + 1) = isol_metrics.clusters(current_nz).metrics.noise_overlap;
                    
            end
            
            for u = 1:length(good_units)

                u_color = unit_colors( u, : );
                current_unit = good_units(u);
               
                lng_colors(u + length(noise_units) , :) = u_color;
                spikes(length(spikes) + 1) = metrics.clusters(current_unit).metrics.num_events;
                snr(length(snr) + 1) = isol_metrics.clusters(current_unit).metrics.peak_snr;
                iso(length(iso) + 1) = isol_metrics.clusters(current_unit).metrics.isolation;
                over(length(over) + 1) = isol_metrics.clusters(current_unit).metrics.noise_overlap;
                     
            end
                        
            p = plot(0, 0);
            hold on;

            set(gca,'Color','k');
            set(gca,'xtick',[])
            set(gca,'xticklabel',[])
            set(gca,'ytick',[])
            set(gca,'yticklabel',[])

            xlim([0 1]);
            ylim([0 1]);

            col_coords = linspace(0.1, 0.9, length(column_labels));

            num_rows = length(unit_names) + 1;% + 1 for header
            row_coords = linspace(0.1,0.9,num_rows);

            %plot header
            header_row = row_coords(end);
            for i = 1:length(col_coords)
                text(col_coords(i), header_row, column_labels{i}, 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', font_size);
            end
            
            
            %plot header line
            header_line = (row_coords(end) - row_coords(end - 1))/4 + row_coords(end - 1);
            for i = linspace(0.1,0.9, 25)
                text(i, header_line, '-', 'Color', 'w', 'FontSize', font_size);
            end

            %plot filter cutoffs
            filter_line = (row_coords(end) - header_line)/2 + header_line;
            for i = 1:length(col_coords)
                if ~isempty(filter_labels{i})
                    text(col_coords(i), filter_line, ['(' filter_labels{i} ')'], 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', font_size);
                end
            end

            
            col_half_dist = (col_coords(end) - col_coords(end -1)) /3;
            row_quart_dist = (row_coords(end) - row_coords(end -1)) /6;
            
            for r = 1:length(unit_names)

                current_row = row_coords(length(row_coords)-r);

                fill_x = [col_coords(1) - col_half_dist; col_coords(1) - col_half_dist; col_coords(1) + col_half_dist; col_coords(1) + col_half_dist];
                fill_y =  [current_row  - row_quart_dist;current_row  + row_quart_dist; current_row  + row_quart_dist; current_row  - row_quart_dist];

                %fill
                fill(fill_x, fill_y, lng_colors(r, :));

                %name and metrics
                text(col_coords(2), current_row, strrep(unit_names{r}, '_', ' '), 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', font_size);

                text(col_coords(3), current_row, sprintf('%d', spikes(r)), 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', font_size);

                if snr(r) < snr_filt
                    text(col_coords(4), current_row, sprintf('%0.3f', snr(r)), 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', font_size, 'EdgeColor', 'r');
                else
                    text(col_coords(4), current_row, sprintf('%0.3f', snr(r)), 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', font_size);
                end

                if over(r) > over_filt
                    text(col_coords(5), current_row, sprintf('%0.3f', over(r)), 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', font_size, 'EdgeColor', 'r');
                else
                    text(col_coords(5), current_row, sprintf('%0.3f', over(r)), 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', font_size);
                end
                
                if iso(r) < iso_filt
                    text(col_coords(6), current_row, sprintf('%0.3f', iso(r)), 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', font_size, 'EdgeColor', 'g');
                else
                    text(col_coords(6), current_row, sprintf('%0.3f', iso(r)), 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', font_size);
                end

            end
            
        else
            %fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 7, current_chan);

            plot(1:10, 1:10);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %example trace figures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('plot example traces %s\n', datetime('now'));

        hold off;
        current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{3,2}(1:subplot_cols_per_display_col) ,'replace');
        hold on;

        s = get(gca, 'Position');
        set(gca, 'Position', [s(1), s(2), s(3)* full_subplot_width_scalar, s(4)* full_subplot_height_scalar ]);
        set(gca, 'FontSize', font_size);

        if on_panels(8) == 1
            %fprintf('\t-- plotting real panel %d for chan %s\n', 8, current_chan);

            top_clip = example_clips(1,:);
            plot(downsample(top_clip, 30), 'k-');
            
            xlabel(sprintf('time (%d sec)', ceil(length(top_clip)/30000)));
            ylim_pad = abs(max(top_clip) - min(top_clip))/2 + 1;
            ylim([min(top_clip) - ylim_pad max(top_clip) + ylim_pad]);
        
        else
            %fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 8, current_chan);

            plot(1:10, 1:10);
        end

        hold off;
        current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{3,2}(subplot_cols_per_display_col+1:end) ,'replace');
        hold on;

        s = get(gca, 'Position');
        set(gca, 'Position', [s(1), s(2), s(3)* full_subplot_width_scalar, s(4)* full_subplot_height_scalar ]);
        set(gca, 'FontSize', font_size);

        if on_panels(8) == 1
            %fprintf('\t-- plotting real panel %d for chan %s\n', 8, current_chan);

            bottom_clip = example_clips(2,:);
            plot(downsample(bottom_clip, 30), 'k-');
            
            xlabel(sprintf('time (%d sec)', ceil(length(bottom_clip)/30000)));
            ylim_pad = abs(max(bottom_clip) - min(bottom_clip))/2;
            
            ylim([ (min(bottom_clip) - ylim_pad - 1) (max(bottom_clip) + ylim_pad + 1) ]);

        else
            %fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 8, current_chan);

            plot(1:10, 1:10);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %spike rate figures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('plot spikerates %s\n', datetime('now'));

        hold off;
        current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{3,3} ,'replace');
        hold on;
        
        if ~isempty(good_units)

            if on_panels(9) == 1
                %fprintf('\t-- plotting real panel %d for chan %s\n', 9, current_chan);

                set(gca,'Color','k');
                set(gca, 'FontSize', font_size);  
                
                 s = get(gca, 'Position');
                set(gca, 'Position', [s(1), s(2), s(3)* full_subplot_width_scalar, s(4)* full_subplot_height_scalar ]);

                avg_secs = 180;

                xlabel('time (min)');
                ylabel(sprintf('spikes/sec (averaged over %d min)', avg_secs/60));

                ymn = 0;
                ymx = 0;

                xmn = 0;
                xmx = ceil(max(firings(2,:))/30000/60);


                for u = 1:length(good_units)

                    u_color = unit_colors( u, : );
                    current_unit = good_units(u);
                    unitIDs_current_unit_condi = (firings(3,:) == current_unit);

                    %dont include in plot if unit has few spikes
                    if sum(unitIDs_current_unit_condi) > 10
                    
                        unit_timestamps = firings(2,unitIDs_current_unit_condi);

                        first_timepoint = min(unit_timestamps);
                        last_timepoint = max(unit_timestamps);

                        first_sec = ceil(first_timepoint/30000);
                        last_sec = ceil(last_timepoint/30000);

                        spk_rate = zeros((last_sec - first_sec), 1);

                        %calculate spikes per second
                        for spr = first_sec:last_sec-1
                            spk_rate(spr-first_sec+1) = sum(unit_timestamps >= (spr*30000) & unit_timestamps < ((spr+1) * 30000));  
                        end


                        %average the rate over minutes
                        %if unit duration is less than averaging time, assign
                        %average spike/sec value to each second of that unit

                        if length(spk_rate) >= avg_secs

                            for spr = 1:length(spk_rate)    
                                if spr + avg_secs <= length(spk_rate)
                                    spk_rate(spr) = sum(spk_rate(spr:spr + avg_secs))/avg_secs;
                                else
                                    spk_rate(spr) = NaN;
                                end
                            end

                        else

                            spk_rate(:) = mean(spk_rate);

                        end

                        if max(spk_rate) > ymx
                            ymx = max(spk_rate);
                        end

                        plot(current_plot, (first_sec:last_sec-1)/60, spk_rate,  '-', 'LineWidth', 2, 'MarkerFaceColor', u_color, 'MarkerEdgeColor', u_color); 
                        hold on;

                    end
                end

                ylim([ymn ymx + 1]);
                xlim([xmn xmx]);

            else
                %fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 9, current_chan);

                plot(1:10, 1:10);
            end

        else
            %fprintf('\t-- plotting EMPTY panel %d for chan %s\n', 9, current_chan);

            plot(1:10, 1:10);
            hold on;
            plot(1:10, fliplr(1:10));
            hold off;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % data variance plot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('plot timeseries %s\n', datetime('now'));

        
        hold off;
        current_plot = subplot(subplot_rows, subplot_cols, [subplot_sections{4, :}] ,'replace');
        hold on;

        set(gca, 'FontSize', font_size);
        
        mda_sd = std(mda);
        sd_line = -3 * mda_sd;
        
        plot(1:length(mda), mda, 'k-');
        plot(1:length(mda), repmat(sd_line, 1, length(mda)), 'r-', 'LineWidth', 2);
        
        samples_per_five_minute = 30000 * 60 * 5;
        
        tickspots_index = 1:samples_per_five_minute:length(mda);
        ticklabels = cellstr(strsplit(num2str( (0:(length(tickspots_index)-1)) * 5 ), ' '));
        
        xticks(current_plot, tickspots_index);
        xticklabels(current_plot, ticklabels);
        
        ylim(quantile(mda,[0.00001 0.99999]));
        
        xlabel(current_plot, 'time (min)');
        ylabel(current_plot, 'ungained, whitened voltage');
        
        hold off;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % done plotting!
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

             
        
        %keyboard;
        
        %print([ms_figDir '/' session_name_no_spaces '_' current_chan '.png'], '-painters', '-dpng', '-r100');
        %saveas(gcf, [ms_figDir '/' session_name_no_spaces '_' current_chan '.pdf'], 'pdf');
        
        tightfig();
        set(gcf, 'InvertHardcopy', 'off','PaperUnits','inches','PaperPosition',[0,0,25,19],'PaperPositionMode','auto');
        fprintf('save fig %s\n', datetime('now'));

        print(gcf,[saveDir '/' session_name_underscore '_' current_chan_underscore '.png'],'-dpng','-r300');
        
        %export_fig([ms_figDir '/' session_name_no_spaces '_' current_chan '.bmp'], '-bmp', '-painters');
        fprintf('\t-- saving\n');
        close(gcf);

    
    
    fprintf('plotChannelSpikes -- done\n');

end




function plotUnitHist(noise_units, good_units, unit_colors, firings, unit_hist_width_scalar, unit_hist_height_scalar, saveFile)
 
     figure(2);
     
     font_size = 5;

     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01, 0.01, 1, 0.96]);
     
     display_cols = length(noise_units) + length(good_units);
     display_rows = 2;
     
     mini_col_count = 1;

     for n = 1:length(noise_units)

        current_nz = noise_units(n);
        nzIDs_current_nz_condi = ( firings(3,:) == current_nz );

        if sum(nzIDs_current_nz_condi) > 10
        
            ms_diffs = diff(firings(2, nzIDs_current_nz_condi)/30);

            %histogram of spike times ISI for this unit
            hold off;
            current_plot = subplot(display_rows, display_cols, mini_col_count ,'replace');
            %fprintf('making hist for noise %d at plot point %d\n', n, mini_col_count + display_cols * floor(mini_col_count/display_cols));
            hold on;

             s = get(gca, 'Position');
             set(gca, 'Position', [s(1), s(2), s(3)*unit_hist_width_scalar, s(4)*unit_hist_height_scalar ])

             if n == 1
                ylabel('count')
                xlabel('ISI (ms)');
             end

             set(gca,'Color','k');
             set(gca, 'FontSize', font_size);

             histogram(current_plot, ms_diffs, 'FaceColor', 'w', 'EdgeColor', 'w', 'BinWidth', 1);
             xlim([-5 inf]);

             ax = gca;
             expon = max([0 floor(min(log10(ax.YAxis.TickValues)))]);
             ax.YAxis.Exponent = expon;

             hold off;
             current_plot = subplot(display_rows, display_cols, mini_col_count + display_cols ,'replace');
             hold on;

             %fprintf('making hist for noise %d at plot point %d\n', n, (mini_col_count + display_cols * floor(mini_col_count/display_cols)) + display_cols);

             s = get(gca, 'Position');
             set(gca, 'Position', [s(1), s(2), s(3)*unit_hist_width_scalar, s(4)*unit_hist_height_scalar ])

             if n == 1
                ylabel('cumulative proportion');
                xlabel('ISI (ms)');
             end

             set(gca,'Color','w');
             set(gca, 'FontSize', font_size);

             num_bins = ceil(max(ms_diffs));
             hist_res = hist(ms_diffs, num_bins);

             cum_dist = zeros(num_bins, 1);
             cum_dist(1) = hist_res(1);            

             for a = 2:num_bins
                cum_dist(a) = cum_dist(a-1) + hist_res(a);
             end
             cum_dist = cum_dist./sum(hist_res);

             txt = sprintf('3 ms (%0.2f)', sum(ms_diffs < 3)/length(ms_diffs));

             plot(1:length(cum_dist), cum_dist, 'b');
             ylim([0 1]);   
             xlim([0 ceil(max(ms_diffs))]);

             ax = gca;
             expon = max([0 floor(min(log10(ax.YAxis.TickValues)))]);
             ax.YAxis.Exponent = expon;

             text(double(ceil(max(ms_diffs)))/2 , 0.5, txt, 'FontSize', 9);

        end
        
        mini_col_count = mini_col_count + 1;

     end

     
     for u = 1:length(good_units)
         
         u_color = unit_colors( u, :);
         current_unit = good_units(u);
         unitIDs_current_unit_condi = ( firings(3,:) == current_unit);

         if sum(unitIDs_current_unit_condi) > 10
         
             ms_diffs = diff( firings(2, unitIDs_current_unit_condi)/30 );

             %histogram of spike times ISI for this unit
             hold off;
             current_plot = subplot(display_rows, display_cols,  mini_col_count ,'replace');
             %fprintf('making hist for unit %d at plot point %d\n', u,  mini_col_count + display_cols * floor((mini_col_count-1)/display_cols));
             hold on;

             s = get(gca, 'Position');
             set(gca, 'Position', [s(1), s(2), s(3)*unit_hist_width_scalar, s(4)*unit_hist_height_scalar])

             if n == 1
                ylabel('count')
                xlabel('ISI (ms)');
             end

             set(gca,'Color','k');
             set(gca, 'FontSize', font_size);

             histogram(current_plot, ms_diffs, 'FaceColor', u_color, 'EdgeColor', u_color, 'BinWidth', 1);
             xlim([-5 inf]);

             ax = gca;
             expon = max([0 floor(min(log10(ax.YAxis.TickValues)))]);
             ax.YAxis.Exponent = expon;

             %then cdf of ISI for this unit

             hold off;
             current_plot = subplot(display_rows, display_cols, mini_col_count + display_cols ,'replace');
             hold on;

             %fprintf('making hist for unit %d at plot point %d\n', u, (mini_col_count + display_cols * floor((mini_col_count-1)/display_cols)) + display_cols);

             s = get(gca, 'Position');
             set(gca, 'Position', [s(1), s(2),s(3)*unit_hist_width_scalar, s(4)*unit_hist_height_scalar ])

             if n == 1
                ylabel('cumulative proportion');
                xlabel('ISI (ms)');
             end

             set(gca,'Color','w');
             set(gca, 'FontSize', font_size);

             num_bins = ceil(max(ms_diffs));
             hist_res = hist(ms_diffs, num_bins);

             cum_dist = zeros(num_bins, 1);
             cum_dist(1) = hist_res(1);            

             for a = 2:num_bins
                cum_dist(a) = cum_dist(a-1) + hist_res(a);
             end
             cum_dist = cum_dist./sum(hist_res);

             txt = sprintf('3 ms (%0.2f)', sum(ms_diffs < 3)/length(ms_diffs));

             plot(1:length(cum_dist), cum_dist, 'b');
             ylim([0 1]);   
             xlim([0 ceil(max(ms_diffs))]);

             ax = gca;
             expon = max([0 floor(min(log10(ax.YAxis.TickValues)))]);
             ax.YAxis.Exponent = expon;

             text(double(ceil(max(ms_diffs)))/2 , 0.5, txt, 'FontSize', 9);
         
         end
         
         mini_col_count = mini_col_count + 1;

     end
     
     fprintf('saving seperate ISI plot\n');
     
     tightfig();
     set(gcf, 'InvertHardcopy', 'off','PaperUnits','inches','PaperPosition', [0,0,25,19], 'PaperPositionMode','auto');
     print(gcf,saveFile,'-dpng','-r300');

     close(gcf);
    
end

function hfig = tightfig(hfig)
% tightfig: Alters a figure so that it has the minimum size necessary to
% enclose all axes in the figure without excess space around them.
%
% Note that tightfig will expand the figure to completely encompass all
% axes if necessary. If any 3D axes are present which have been zoomed,
% tightfig will produce an error, as these cannot easily be dealt with.
%
% hfig - handle to figure, if not supplied, the current figure will be used
% instead.

if nargin == 0
    hfig = gcf;
end

% There can be an issue with tightfig when the user has been modifying
% the contnts manually, the code below is an attempt to resolve this,
% but it has not yet been satisfactorily fixed
%     origwindowstyle = get(hfig, 'WindowStyle');
set(hfig, 'WindowStyle', 'normal');

% 1 point is 0.3528 mm for future use

% get all the axes handles note this will also fetch legends and colorbars as well
hax = findall(hfig, 'type', 'axes');

% get the original axes units, so we can change and reset these again later
origaxunits = get(hax, 'Units');

% change the axes units to cm
set(hax, 'Units', 'centimeters');

% get various position parameters of the axes
if numel(hax) > 1
    %         fsize = cell2mat(get(hax, 'FontSize'));
    ti = cell2mat(get(hax,'TightInset'));
    pos = cell2mat(get(hax, 'Position'));
else
    %         fsize = get(hax, 'FontSize');
    ti = get(hax,'TightInset');
    pos = get(hax, 'Position');
end

% ensure very tiny border so outer box always appears
ti(ti < 0.1) = 0.15;

% we will check if any 3d axes are zoomed, to do this we will check if they are not being viewed in any of the 2d directions
views2d = [0,90; 0,0; 90,0];

for i = 1:numel(hax)
    
    set(hax(i), 'LooseInset', ti(i,:));
    %         set(hax(i), 'LooseInset', [0,0,0,0]);
    
    % get the current viewing angle of the axes
    [az,el] = view(hax(i));
    
    % determine if the axes are zoomed
    iszoomed = strcmp(get(hax(i), 'CameraViewAngleMode'), 'manual');
    
    % test if we are viewing in 2d mode or a 3d view
    is2d = all(bsxfun(@eq, [az,el], views2d), 2);
    
    if iszoomed && ~any(is2d)
        error('TIGHTFIG:haszoomed3d', 'Cannot make figures containing zoomed 3D axes tight.')
    end
    
end

% we will move all the axes down and to the left by the amount necessary to just show the bottom and leftmost axes and labels etc.
moveleft = min(pos(:,1) - ti(:,1));
movedown = min(pos(:,2) - ti(:,2));

% we will also alter the height and width of the figure to just encompass the topmost and rightmost axes and lables
figwidth = max(pos(:,1) + pos(:,3) + ti(:,3) - moveleft);
figheight = max(pos(:,2) + pos(:,4) + ti(:,4) - movedown);

% move all the axes
for i = 1:numel(hax)
    set(hax(i), 'Position', [pos(i,1:2) - [moveleft,movedown], pos(i,3:4)]);
end

origfigunits = get(hfig, 'Units');
set(hfig, 'Units', 'centimeters');

% change the size of the figure
figpos = get(hfig, 'Position');
set(hfig, 'Position', [figpos(1), figpos(2), figwidth, figheight]);

% change the size of the paper
set(hfig, 'PaperUnits','centimeters');
set(hfig, 'PaperSize', [figwidth, figheight]);
set(hfig, 'PaperPositionMode', 'manual');
set(hfig, 'PaperPosition',[0 0 figwidth figheight]);

% reset to original units for axes and figure
if ~iscell(origaxunits)
    origaxunits = {origaxunits};
end

for i = 1:numel(hax)
    set(hax(i), 'Units', origaxunits{i});
end

set(hfig, 'Units', origfigunits);

end

function jsonstring = readjson(fname)

fid = fopen(fname);
s = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
jsonstring = jsondecode(strjoin(s{1}, ' '));

end