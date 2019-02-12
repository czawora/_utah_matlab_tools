
function plotSessionSpikes(varargin)

    %spikeInfo_path = '/Users/zaworaca/dev/biowulf/spikeInfos/NIH029_150608_1113_spikeInfo.mat';   

    %defining some randome helper functions
    split_wrap = @(s) strsplit(s, '_');
    get_first_el = @(s) s{1};
    get_ms_unit_channel = @(s) get_first_el(split_wrap(s));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % parsing varargin
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    MSfilt_snr_min = 1.5;
    MSfilt_iso_min = 0.95;
    MSfilt_noise_overlap_max = 0.03;
    
    p = inputParser;
    p.addParameter('spikeInfo_path', '', @ischar);
    p.addParameter('ms_spikeInfo', '');
    p.addParameter('ms_noiseInfo', '');
    p.addParameter('unit_filters', '');
    p.addParameter('noise_filters', '');
    p.addParameter('unit_filters_path', '', @ischar);
    p.addParameter('noise_filters_path', '', @ischar);
    p.addParameter('unit_filters_savePath', '', @ischar);
    p.addParameter('noise_filters_savePath', '', @ischar);
    p.addParameter('snr_min', MSfilt_snr_min, @isnumeric);
    p.addParameter('isolation_min', MSfilt_iso_min, @isnumeric);
    p.addParameter('noise_overlap_max', MSfilt_noise_overlap_max, @isnumeric);
    p.addParameter('debug', 0, @isint);
    parse(p, varargin{:});
    
    snr_filt = p.Results.snr_min
    iso_filt = p.Results.isolation_min
    over_filt = p.Results.noise_overlap_max
    
    fprintf('%s is %s\n', 'spikeInfo_path', p.Results.spikeInfo_path);
    
    if ~strcmp(p.Results.ms_spikeInfo, '')
        fprintf('%s was passed in\n', 'ms_spikeInfo');
    else
        fprintf('%s was NOT passed in\n', 'ms_spikeInfo');
    end
    
    if ~strcmp(p.Results.ms_noiseInfo, '')
        fprintf('%s was passed in\n', 'ms_noiseInfo');
    else
        fprintf('%s was NOT passed in\n', 'ms_noiseInfo');
    end
    
    if ~strcmp(p.Results.unit_filters, '')
        fprintf('%s was passed in\n', 'unit_filters');
    else
        fprintf('%s was NOT passed in\n', 'unit_filters');
    end
    
    if ~strcmp(p.Results.noise_filters, '')
        fprintf('%s was passed in\n', 'noise_filters');
    else
        fprintf('%s was NOT passed in\n', 'noise_filters');
    end
        
    fprintf('%s is %s\n', 'unit_filters_path', p.Results.unit_filters_path);
    fprintf('%s is %s\n', 'noise_filters_path', p.Results.noise_filters_path);
    fprintf('%s is %s\n', 'unit_filters_savePath', p.Results.unit_filters_savePath);
    fprintf('%s is %s\n', 'noise_filters_savePath', p.Results.noise_filters_savePath);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if debugging on desktop, use this default spikeInfo file
    
    if p.Results.debug == 1
    
        spikeInfo_path = '/Users/zaworaca/dev/biowulf/spikeInfos/NIH029_150608_1113_spikeInfo.mat';   
    
    else
        
        spikeInfo_path = p.Results.spikeInfo_path;
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % load the spikeInfo
    
    if strcmp(p.Results.ms_spikeInfo, '')
            
        % check for existence
        if ~exist(spikeInfo_path, 'file')
        
            error('%s does not exist\n', spikeInfo_path);
        end
        
        load(spikeInfo_path);
        ms_sim = spikeInfo_filt;

        clear spikeInfo;
        
    else
            
        ms_sim = p.Results.ms_spikeInfo;
            
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %make the session sortFig directory if it doesnt exist
    spikeInfo_path_splits = strsplit(spikeInfo_path, '.m');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the noiseInfo
    
    if strcmp(p.Results.ms_noiseInfo, '')
        
        %get the noise from the filtered spikeInfo
        nzInfo_path = [spikeInfo_path_splits{1} '_noise.mat'];

        use_nzInfo = 0;
        
        if ~exist(nzInfo_path, 'file')
           
            fprintf('no spikeInfo_noise.mat present in dir. Plotting without\n');
        else
            
            use_nzInfo = 1;
            load(nzInfo_path);
            ms_nim = noiseInfo;
            clear noiseInfo;
        end

    else
        
        use_nzInfo = 1;
        ms_nim = p.Results.ms_noiseInfo;
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the spike filters for units and spikes
    
    if ~strcmp(p.Results.unit_filters, '')
        
        unit_filters = p.Results.unit_filters;
        
    elseif ~strcmp(p.Results.unit_filters_path, '')
    
        if ~exist(p.Results.unit_filters_path, 'file')
            error('%s does not exist\n', p.Results.unit_filters_path);
        end
        
        load(p.Results.unit_filters_path);
        
    else
        
        %unit filters
        unit_filters = cell(length(ms_sim.uniqueUnitID),1);

        for uu = 1:length(ms_sim.uniqueUnitID)

            fprintf('making filter for unit %s\n', char(ms_sim.uniqueUnitID{ uu }));

            unit_filters{uu} = cellfun( @(x) isequal(x, char(ms_sim.uniqueUnitID{uu})), ms_sim.unitID);
            fprintf('\t-- had %d spikes\n', sum(unit_filters{uu}));

        end
        
        % if path exists, save them to prevent recomputation
        if ~strcmp(p.Results.unit_filters_savePath, '')
            save(p.Results.unit_filters_savePath, '-v7.3', 'unit_filters');
        end
        
    end
    
    if ~strcmp(p.Results.noise_filters, '')
        
        nz_filters = p.Results.noise_filters;
    
    elseif ~strcmp(p.Results.noise_filters_path, '')
       
        if ~exist(p.Results.noise_filters_path, 'file')
            error('%s does not exist\n', p.Results.noise_filters_path);
        end
        
        load(p.Results.noise_filters_path);
        
    else
        
        %noise filters
        nz_filters = cell(length(ms_nim.uniqueUnitID),1);

        for nu = 1:length(ms_nim.uniqueUnitID)

            fprintf('making filter for noise %s\n', char(ms_nim.uniqueUnitID{ nu }));

            nz_filters{ nu} = cellfun( @(x) isequal(x, char(ms_nim.uniqueUnitID{nu})), ms_nim.unitID);
            fprintf('\t-- had %d spikes\n', sum(nz_filters{nu}));
        end
        
        % if path exists, save them to prevent recomputation
        if ~strcmp(p.Results.noise_filters_savePath, '')
            save(p.Results.noise_filters_savePath, '-v7.3', 'nz_filters');
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % done parsing varargin
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    unit_colors = distinguishable_colors(15, {'w', 'k'});
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
    
    full_subplot_width_scalar = 0.9;
    full_subplot_height_scalar = 0.9;
    
    unit_hist_width_scalar = 0.9;
    unit_hist_height_scalar = 0.9;
    
    %make the session sortFig directory if it doesnt exist
    ms_figDir = [spikeInfo_path_splits{1} '_sortFigs'];
    
    if ~exist(ms_figDir, 'dir')
        mkdir(ms_figDir);
    end

    %get the session name
    slash_splits = strsplit(spikeInfo_path, '/');
    session_fname_splits = strsplit(slash_splits{end}, '_spikeInfo.mat');
    session_name_no_spaces = session_fname_splits{1};
    session_name = strrep(session_name_no_spaces, '_', ' ');

    unique_channels = unique(cellfun(get_ms_unit_channel, ms_sim.uniqueUnitID, 'UniformOutput', 0));
    
    if use_nzInfo == 1
        noise_channels = unique(cellfun(get_ms_unit_channel, ms_nim.uniqueUnitID, 'UniformOutput', 0));
        unique_channels = unique([unique_channels;noise_channels]);
    end
    %for each channel, produce the figures

    wf_shading_opacity = 0.3;

    font_size = 5;
    pca_msize = 2;
    time_volt_msize = 2;

    display_rows = 3;
    display_cols = 3;

    unique_channels
        
    for c = 1:length(unique_channels)

        current_chan = unique_channels{c};
        fprintf('top of the loop for chan %s\n', current_chan);

        uniq_unitIDs_current_channel_condi = contains( cellfun(@char, ms_sim.uniqueUnitID, 'UniformOutput',0), current_chan);
        uniq_unitIDs_current_channel = ms_sim.uniqueUnitID(uniq_unitIDs_current_channel_condi);

        %check if noise will be plotted for this channel
        uniq_nzID_current_channel = {};

        if use_nzInfo == 1

            nz_uniqueUnitID_match = contains(ms_nim.uniqueUnitID, current_chan);
            if ~isempty(nz_uniqueUnitID_match)

                uniq_nzID_current_channel = ms_nim.uniqueUnitID(nz_uniqueUnitID_match);
            end
        end

        %unit ISI histograms need to plotted as seperate image
        if length(uniq_unitIDs_current_channel) + length(uniq_nzID_current_channel) > 3

            subplot_rows_per_display_row = 2;
            subplot_cols_per_display_col = 1;
            
            figure(2); % declare the extra ISI plot
            
        %otherwise they can all fit    
        else
            
            subplot_rows_per_display_row = 2;
            subplot_cols_per_display_col = length(uniq_unitIDs_current_channel) + length(uniq_nzID_current_channel);

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
        
        fprintf('chan has %d units\n', sum(uniq_unitIDs_current_channel_condi));

        figure(1);%set figure 1 as gcf
        
        %create unit_colors list to fraw from
        nz_markers = repmat(nz_markers_set, 1, ceil(length(uniq_nzID_current_channel)/length(nz_markers_set)));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %waveform figure - top left ( fine to run if unit contains only single spike)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{1,1}, 'replace');
        hold on;

        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01, 0.01, 1, 0.96]);

        s = get(gca, 'Position');
        set(gca, 'Position', [s(1), s(2), s(3)* full_subplot_width_scalar, s(4)* full_subplot_height_scalar ]);

        if on_panels(1) == 1

            fprintf('\t-- plotting real panel %d for chan %s\n', 1, current_chan);

            ylabel('uV');
            %set bg color
            set(gca,'Color','k');
            set(gca, 'FontSize', font_size);

            min_y = 0;
            max_y = 0;

            legend_cell = {};
            legend_handles = [];

            %plot noise as white

            for n = 1:length(uniq_nzID_current_channel)

                current_nz = uniq_nzID_current_channel{n};
                fprintf('\t\t-- plotting real panel %d for chan %s noise %d\n', 1, current_chan, n);

                nzIDs_current_nz_condi = nz_filters{ cellfun( @(x) isequal(x, current_nz), ms_nim.uniqueUnitID) };

                avg_wf = mean(ms_nim.waveForm(nzIDs_current_nz_condi,:), 1);
                %fprintf('\t-- got avg\n');
                sd_wf = std(ms_nim.waveForm(nzIDs_current_nz_condi,:),0, 1);
                %fprintf('\t-- got std\n');

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

                legend_cell{length(legend_cell) + 1} = ['noise ' num2str(n) ' ( ' nz_markers( mod(n, length(nz_markers)) ) ' )'];
                legend_handles(length(legend_handles) + 1) = l;

            end
            %now plot the unit waveforms

            for u = 1:length(uniq_unitIDs_current_channel)

                fprintf('\t\t-- plotting real panel %d for chan %s unit %d\n', 1, current_chan, u);


                current_unit = uniq_unitIDs_current_channel{u};
                %fprintf('\t-- %s\n', current_unit{1});

                u_color = unit_colors( u, : );

                unitIDs_current_unit_condi = unit_filters{ cellfun( @(x) isequal(x, current_unit), ms_sim.uniqueUnitID) };
                %fprintf('\t-- %d spikes\n', sum(unitIDs_current_unit_condi));

                avg_wf = mean(ms_sim.waveForm(unitIDs_current_unit_condi,:), 1);
                %fprintf('\t-- got avg\n');
                sd_wf = std(ms_sim.waveForm(unitIDs_current_unit_condi,:),0, 1);
                %fprintf('\t-- got std\n');

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
                    [ (1:size(ms_sim.waveForm,2))' ; flipud((1:size(ms_sim.waveForm,2))')], ...
                    [ avg_wf - sd_wf ; flipud(avg_wf + sd_wf)],...
                    u_color, 'linestyle','none');
                alpha(wf_shading_opacity);

                l = line( 1:size(ms_sim.waveForm,2) , avg_wf);
                set(l, 'Color', u_color);

                legend_cell{length(legend_cell) + 1} = ['unit ' num2str(u)];
                legend_handles(length(legend_handles) + 1) = l;

                %fprintf('\t-- plotting\n');

            end

            %set axes
            y_axis_offset = 50;
            xmin = 1;
            xmax = size(ms_sim.waveForm,2);
            ymin = min_y - y_axis_offset;
            ymax = max_y + y_axis_offset;
            axis(current_plot, [ xmin xmax ymin ymax ]);

            title( [session_name ' ' current_chan] );
            
        else
            fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 1, current_chan);
            plot(1:10, 1:10);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %PCA figures ( fine to run if unit contains only single spike)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %uncolored PCA
        hold off;
        current_plot = subplot(subplot_rows, subplot_cols, subplot_sections{1,2},'replace');
        hold on;

        s = get(gca, 'Position');
        set(gca, 'Position', [s(1), s(2), s(3)* full_subplot_width_scalar, s(4)* full_subplot_height_scalar ]);

        if on_panels(2) == 1

            fprintf('\t-- plotting real panel %d for chan %s\n', 2, current_chan);

            xlabel('PC 1');
            ylabel('PC 2');
            %set bg color
            set(gca,'Color','k');
            set(gca, 'FontSize', font_size);

            for n = 1:length(uniq_nzID_current_channel)

                fprintf('\t\t-- plotting real panel %d for chan %s noise %d\n', 2, current_chan, n);

                current_nz = uniq_nzID_current_channel{n};
                nzIDs_current_nz_condi = nz_filters{ cellfun( @(x) isequal(x, current_nz), ms_nim.uniqueUnitID) };

                 %get waveform PC values
                pc1 = ms_nim.waveForm_PC(nzIDs_current_nz_condi,1);
                pc2 = ms_nim.waveForm_PC(nzIDs_current_nz_condi,2);

                plot(pc1, pc2, 'w.',  'MarkerSize', pca_msize);

            end

            for u = 1:length(uniq_unitIDs_current_channel)

                fprintf('\t\t-- plotting real panel %d for chan %s unit %d\n', 2, current_chan, u);

                current_unit = uniq_unitIDs_current_channel{u};
                unitIDs_current_unit_condi = unit_filters{ cellfun( @(x) isequal(x, current_unit), ms_sim.uniqueUnitID) };

                %keyboard;
                %get waveform PC values
                pc1 = ms_sim.waveForm_PC(unitIDs_current_unit_condi,1);
                pc2 = ms_sim.waveForm_PC(unitIDs_current_unit_condi,2);

                plot(pc1, pc2, 'w.',  'MarkerSize', pca_msize);

            end
            
            xlim([-2000 2000]);
            ylim([-2000 2000]);

        else
            fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 2, current_chan);

            plot(1:10, 1:10);
        end

        %colored PCA ( fine to run if unit contains only single spike)
        hold off;
        current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{2,2} ,'replace');
        hold on;

        s = get(gca, 'Position');
        set(gca, 'Position', [s(1), s(2), s(3)* full_subplot_width_scalar, s(4)* full_subplot_height_scalar ]);


        if on_panels(5) == 1

            fprintf('\t-- plotting real panel %d for chan %s\n', 5, current_chan);

            xlabel('PC 1');
            ylabel('PC 2');
            %set bg color
            set(gca,'Color','k');
            set(gca, 'FontSize', font_size);

            for n = 1:length(uniq_nzID_current_channel)

                fprintf('\t\t-- plotting real panel %d for chan %s noise %d\n', 5, current_chan, n);

                n_marker = nz_markers( n );
                current_nz = uniq_nzID_current_channel{n};
                nzIDs_current_nz_condi = nz_filters{ cellfun( @(x) isequal(x, current_nz), ms_nim.uniqueUnitID) };

                 %get waveform PC values
                pc1 = ms_nim.waveForm_PC(nzIDs_current_nz_condi,1);
                pc2 = ms_nim.waveForm_PC(nzIDs_current_nz_condi,2);

                plot(pc1, pc2, ['w' n_marker],  'MarkerSize', pca_msize);

            end

            for u = 1:length(uniq_unitIDs_current_channel)

                fprintf('\t\t-- plotting real panel %d for chan %s unit %d\n', 5, current_chan, u);

                u_color = unit_colors( u, : );
                current_unit = uniq_unitIDs_current_channel{u};
                unitIDs_current_unit_condi = unit_filters{ cellfun( @(x) isequal(x, current_unit), ms_sim.uniqueUnitID) };

                %keyboard;
                %get waveform PC values
                pc1 = ms_sim.waveForm_PC(unitIDs_current_unit_condi,1);
                pc2 = ms_sim.waveForm_PC(unitIDs_current_unit_condi,2);

                plot(pc1, pc2, '.', 'MarkerSize', pca_msize, 'MarkerFaceColor', u_color, 'MarkerEdgeColor', u_color);

            end
            
            xlim([-2000 2000]);
            ylim([-2000 2000]);

        else
            fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 5, current_chan);

            plot(1:10, 1:10);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %time-voltage figures ( fine to run if unit contains only single spike )
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        hold off;
        current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{1,3} ,'replace');
        hold on;

        s = get(gca, 'Position');        
        set(gca, 'Position', [s(1), s(2), s(3)* full_subplot_width_scalar, s(4)* full_subplot_height_scalar ]);

        if on_panels(3) == 1
            
            fprintf('\t-- plotting real panel %d for chan %s\n', 3, current_chan);

            xlabel('time (min)');
            ylabel('peak-to-trough');
            %set bg color
            set(gca,'Color','k');
            set(gca, 'FontSize', font_size);

            ymn = 0;
            ymx = 0;
            
            xmx = 0;

            for n = 1:length(uniq_nzID_current_channel)

                fprintf('\t\t-- plotting real panel %d for chan %s noise %d\n', 3, current_chan, n);

                current_nz = uniq_nzID_current_channel{n};
                nzIDs_current_nz_condi = nz_filters{ cellfun( @(x) isequal(x, current_nz), ms_nim.uniqueUnitID) };

                wv = ms_nim.waveForm(nzIDs_current_nz_condi,:);
                timestamps = ms_nim.timeStamp(nzIDs_current_nz_condi,:)/30000/60;

                time_min = min(timestamps);
                time_max = max(timestamps);
                
                if time_max > xmx
                    xmx = time_max;
                end

                pk2tr = zeros(size(wv,1),1); 

                for iwv = 1:size(wv,1)

                    current_pk2tr = abs( max(wv(iwv,:)) - min(wv(iwv,:)) );

                    if current_pk2tr < ymn
                        ymn = current_pk2tr;
                    elseif current_pk2tr > ymx 
                        ymx = current_pk2tr;
                    end

                    pk2tr(iwv) = current_pk2tr;

                end

                plot(timestamps, pk2tr, 'w.', 'MarkerSize', time_volt_msize);

            end

            for u = 1:length(uniq_unitIDs_current_channel)

                fprintf('\t\t-- plotting real panel %d for chan %s unit %d\n', 3, current_chan, u);

                current_unit = uniq_unitIDs_current_channel{u};
                unitIDs_current_unit_condi = unit_filters{ cellfun( @(x) isequal(x, current_unit), ms_sim.uniqueUnitID) };

                wv = ms_sim.waveForm(unitIDs_current_unit_condi,:);
                timestamps = ms_sim.timeStamp(unitIDs_current_unit_condi,:)/30000/60;

                time_min = min(timestamps);
                time_max = max(timestamps);

                if time_max > xmx
                    xmx = time_max;
                end
                
                pk2tr = zeros(size(wv,1),1);

                for iwv = 1:size(wv,1)

                    current_pk2tr = abs( max(wv(iwv,:)) - min(wv(iwv,:)) );

                    if current_pk2tr < ymn
                        ymn = current_pk2tr;
                    elseif current_pk2tr > ymx 
                        ymx = current_pk2tr;
                    end

                    pk2tr(iwv) = current_pk2tr;

                end

                plot(timestamps, pk2tr, 'w.', 'MarkerSize', time_volt_msize);

            end

            axis(current_plot, [0 time_max + 5 ymn ymx + (ymx/10)]);
            xlim([0 xmx]);
            ylim([0 2000]);

        else
            fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 3, current_chan);

            plot(1:10, 1:10);
        end

        hold off;
        current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{2,3} ,'replace');
        hold on;
        
        s = get(gca, 'Position');
        set(gca, 'Position', [s(1), s(2), s(3)* full_subplot_width_scalar, s(4)* full_subplot_height_scalar ]);

        if on_panels(6) == 1
            fprintf('\t-- plotting real panel %d for chan %s\n', 6, current_chan);

            xlabel('time (min)');
            ylabel('peak-to-trough');
            %set bg color
            set(gca,'Color','k');
            set(gca, 'FontSize', font_size);

            ymn = 0;
            ymx = 0;
            
            xmx = 0;

            for n = 1:length(uniq_nzID_current_channel)

                fprintf('\t\t-- plotting real panel %d for chan %s noise %d\n', 6, current_chan, n);

                n_marker = nz_markers( mod(n, length(nz_markers)) );
                current_nz = uniq_nzID_current_channel{n};
                nzIDs_current_nz_condi = nz_filters{ cellfun( @(x) isequal(x, current_nz), ms_nim.uniqueUnitID) };

                wv = ms_nim.waveForm(nzIDs_current_nz_condi,:);
                timestamps = ms_nim.timeStamp(nzIDs_current_nz_condi,:)/30000/60;

                time_min = min(timestamps);
                time_max = max(timestamps);
                
                if time_max > xmx
                    xmx = time_max;
                end

                pk2tr = zeros(size(wv,1),1); 

                for iwv = 1:size(wv,1)

                    current_pk2tr = abs( max(wv(iwv,:)) - min(wv(iwv,:)) );

                    if current_pk2tr < ymn
                        ymn = current_pk2tr;
                    elseif current_pk2tr > ymx 
                        ymx = current_pk2tr;
                    end

                    pk2tr(iwv) = current_pk2tr;

                end

                plot(timestamps, pk2tr, ['w' n_marker], 'MarkerSize', time_volt_msize);

            end

            for u = 1:length(uniq_unitIDs_current_channel)

                fprintf('\t\t-- plotting real panel %d for chan %s unit %d\n', 6, current_chan, u);

                u_color = unit_colors( u , :);
                current_unit = uniq_unitIDs_current_channel{u};
                unitIDs_current_unit_condi = unit_filters{ cellfun( @(x) isequal(x, current_unit), ms_sim.uniqueUnitID) };

                wv = ms_sim.waveForm(unitIDs_current_unit_condi,:);
                timestamps = ms_sim.timeStamp(unitIDs_current_unit_condi,:)/30000/60;

                time_min = min(timestamps);
                time_max = max(timestamps);

                if time_max > xmx
                    xmx = time_max;
                end
                
                pk2tr = zeros(size(wv,1),1);

                for n = 1:size(wv,1)            

                    current_pk2tr = abs( max(wv(n,:)) - min(wv(n,:)) );

                    if current_pk2tr < ymn
                        ymn = current_pk2tr;
                    elseif current_pk2tr > ymx 
                        ymx = current_pk2tr;
                    end

                    pk2tr(n) = current_pk2tr;            
                end

                plot(timestamps, pk2tr, '.', 'MarkerSize', time_volt_msize, 'MarkerFaceColor', u_color, 'MarkerEdgeColor', u_color);

            end

            axis(current_plot, [0  time_max + 5 ymn ymx + (ymx/10)]);
            xlim([0 xmx]);
            ylim([0 2000]);
            

        else
            fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 6, current_chan);

            plot(1:10, 1:10);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %spike ISI distribution figures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if on_panels(4) == 1

            %a seperate figure will be made for these hist ISI plots
            if length(uniq_nzID_current_channel) + length(uniq_unitIDs_current_channel) > 3
                
                fprintf('\t-- plotting SEPERATE panel %d for chan %s\n', 4, current_chan);
                
                hold off;
                current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{2,1} ,'replace');
                hold on;
                
                s = get(gca, 'Position');
                set(gca, 'Position', [s(1), s(2), s(3)*unit_hist_width_scalar, s(4)*unit_hist_height_scalar ])
                
                plot(1:10, zeros(1,10));
                ylim([0 10]);
                
                txt = '> 3 units, see extra figure';
                text(5, 5, txt, 'FontSize', 14, 'HorizontalAlignment', 'center');
             
                %passing all the necessary items to make a plot
                plotUnitHist(uniq_nzID_current_channel, uniq_unitIDs_current_channel, nz_filters, unit_filters, ms_sim, ms_nim, unit_colors, unit_hist_width_scalar, unit_hist_height_scalar, [ms_figDir '/' session_name_no_spaces '_' current_chan '_hist.png']);
                
                figure(1);%set gcf back to main figure
                
            %no seperate figure needed, plot the unit ISI histograms as usual    
            else

                fprintf('\t-- plotting real panel %d for chan %s\n', 4, current_chan);
                fprintf('\t-- using plot points %d\n', subplot_sections{2,1});

                mini_col_count = 1;

                 for n = 1:length(uniq_nzID_current_channel)

                    current_nz = uniq_nzID_current_channel{n};
                    nzIDs_current_nz_condi = nz_filters{ cellfun( @(x) isequal(x, current_nz), ms_nim.uniqueUnitID) };
                    
                    %dont make image if there are barely any spikes
                    if sum(nzIDs_current_nz_condi) > 10
                        
                        ms_diffs = diff(ms_nim.timeStamp(nzIDs_current_nz_condi,:)/30);

                        %histogram of spike times ISI for this unit
                        hold off;
                        current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{2,1}(mini_col_count) ,'replace');
                        fprintf('making hist for noise %d at plot point %d\n', n, subplot_sections{2,1}(mini_col_count));
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
                         current_plot = subplot(subplot_rows, subplot_cols, subplot_sections{2,1}(mini_col_count + subplot_cols_per_display_col) ,'replace');
                         hold on;

                         fprintf('making hist for noise %d at plot point %d\n', n, subplot_sections{2,1}(mini_col_count + subplot_cols_per_display_col));

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

                         text(mean([0 ceil(max(ms_diffs))]) , mean([0 1]), txt, 'FontSize', 6);

                    end
                    mini_col_count = mini_col_count + 1;

                 end

                 for u = 1:length(uniq_unitIDs_current_channel)

                     u_color = unit_colors( u, : );
                     current_unit = uniq_unitIDs_current_channel{u};
                     unitIDs_current_unit_condi = unit_filters{ cellfun( @(x) isequal(x, current_unit), ms_sim.uniqueUnitID) };
                    
                     %dont make image if there are barely any spikes
                     if sum(unitIDs_current_unit_condi) > 10
                     
                         ms_diffs = diff(ms_sim.timeStamp(unitIDs_current_unit_condi,:)/30);

                         %histogram of spike times ISI for this unit
                         hold off;
                         current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{2,1}(mini_col_count) ,'replace');
                         fprintf('making hist for unit %d at plot point %d\n', u, subplot_sections{2,1}(mini_col_count));
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
                         current_plot = subplot(subplot_rows, subplot_cols, subplot_sections{2,1}(mini_col_count + subplot_cols_per_display_col) ,'replace');
                         hold on;

                         fprintf('making hist for unit %d at plot point %d\n', u, subplot_sections{2,1}(mini_col_count + subplot_cols_per_display_col));

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

                         text(mean([0 ceil(max(ms_diffs))]) , mean([0 1]), txt, 'FontSize', 6);
                     
                     end
                     
                     mini_col_count = mini_col_count + 1;

                 end
            end
        else

            fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 4, current_chan);

            hold off;
            subplot(subplot_rows,subplot_cols, subplot_sections{2,1} ,'replace');
            hold on;

            plot(1:10, 1:10);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %metrics figures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        hold off;
        current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{3,1},'replace');
        hold on;

        s = get(gca, 'Position');
        set(gca, 'Position', [s(1), s(2), s(3)* full_subplot_width_scalar, s(4)* full_subplot_height_scalar ]);
        %set(gca, 'FontSize', font_size);

        if on_panels(7) == 1
            fprintf('\t-- plotting real panel %d for chan %s\n', 7, current_chan);
            
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
            
            for n = 1:length(uniq_nzID_current_channel)

                current_nz = uniq_nzID_current_channel{n};
                nzIDs_current_nz_condi = cellfun( @(x) isequal(x, current_nz), ms_nim.uniqueUnitID);

                
                lng_colors(n, :) = [1 1 1];
                spikes(length(spikes) + 1) = ms_nim.metrics.num_events(nzIDs_current_nz_condi);
                snr(length(snr) + 1) = ms_nim.metrics.peak_snr(nzIDs_current_nz_condi);
                iso(length(iso) + 1) = ms_nim.metrics.isolation(nzIDs_current_nz_condi);
                over(length(over) + 1) = ms_nim.metrics.noise_overlap(nzIDs_current_nz_condi);
                    
            end
            
            for u = 1:length(uniq_unitIDs_current_channel)

                u_color = unit_colors( u, : );
                current_unit = uniq_unitIDs_current_channel{u};
                unitIDs_current_unit_condi = cellfun( @(x) isequal(x, current_unit), ms_sim.uniqueUnitID);
               
                lng_colors(u + length(uniq_nzID_current_channel) , :) = u_color;
                spikes(length(spikes) + 1) = ms_sim.metrics.num_events(unitIDs_current_unit_condi);
                snr(length(snr) + 1) = ms_sim.metrics.peak_snr(unitIDs_current_unit_condi);
                iso(length(iso) + 1) = ms_sim.metrics.isolation(unitIDs_current_unit_condi);
                over(length(over) + 1) = ms_sim.metrics.noise_overlap(unitIDs_current_unit_condi);
                    
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
                text(col_coords(2), current_row, unit_names{r}, 'Color', 'w', 'HorizontalAlignment', 'center', 'FontSize', font_size);

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
            fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 7, current_chan);

            plot(1:10, 1:10);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %example trace figures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        hold off;
        current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{3,2}(1:subplot_cols_per_display_col) ,'replace');
        hold on;

        s = get(gca, 'Position');
        set(gca, 'Position', [s(1), s(2), s(3)* full_subplot_width_scalar, s(4)* full_subplot_height_scalar ]);
        set(gca, 'FontSize', font_size);

        if on_panels(8) == 1
            fprintf('\t-- plotting real panel %d for chan %s\n', 8, current_chan);

            top_bp_clip = squeeze(ms_sim.example_BPclip(c, 1, :));
            plot(downsample(top_bp_clip,30), 'k-');
            
            xlabel(sprintf('time (%d sec)', ceil(length(top_bp_clip)/30000)));
            ylim_pad = abs(max(top_bp_clip) - min(top_bp_clip))/2;
            ylim([min(top_bp_clip) - ylim_pad max(top_bp_clip) + ylim_pad]);
        
        else
            fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 8, current_chan);

            plot(1:10, 1:10);
        end

        hold off;
        current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{3,2}(subplot_cols_per_display_col+1:end) ,'replace');
        hold on;

        s = get(gca, 'Position');
        set(gca, 'Position', [s(1), s(2), s(3)* full_subplot_width_scalar, s(4)* full_subplot_height_scalar ]);
        set(gca, 'FontSize', font_size);

        if on_panels(8) == 1
            fprintf('\t-- plotting real panel %d for chan %s\n', 8, current_chan);

            bottom_bp_clip = squeeze(ms_sim.example_BPclip(c, 2, :));
            plot(downsample(bottom_bp_clip,30), 'k-');
            
            xlabel(sprintf('time (%d sec)', ceil(length(bottom_bp_clip)/30000)));
            ylim_pad = abs(max(bottom_bp_clip) - min(bottom_bp_clip))/2;
            ylim([min(bottom_bp_clip) - ylim_pad max(bottom_bp_clip) + ylim_pad]);

        else
            fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 8, current_chan);

            plot(1:10, 1:10);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %spike rate figures
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        hold off;
        current_plot = subplot(subplot_rows,subplot_cols, subplot_sections{3,3} ,'replace');
        hold on;
        
        if ~isempty(uniq_unitIDs_current_channel)

            if on_panels(9) == 1
                fprintf('\t-- plotting real panel %d for chan %s\n', 9, current_chan);

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
                xmx = ceil(max(ms_sim.timeStamp)/30000/60);


                for u = 1:length(uniq_unitIDs_current_channel)

                    u_color = unit_colors( u, : );
                    current_unit = uniq_unitIDs_current_channel{u};
                    unitIDs_current_unit_condi = unit_filters{ cellfun( @(x) isequal(x, current_unit), ms_sim.uniqueUnitID) };

                    %dont include in plot if unit has few spikes
                    if sum(unitIDs_current_unit_condi) > 10
                    
                        unit_timestamps = ms_sim.timeStamp(unitIDs_current_unit_condi);

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
                                    lower_idx = spr - avg_secs + (length(spk_rate) - spr);

                                    %fprintf('unit %d: spr - avg_secs + (length(spk_rate) - spr), %d - %d + (%d - %d) = %d\n',u, spr, avg_secs, length(spk_rate), spr, lower_idx); 

                                    spk_rate(spr) = sum(spk_rate(lower_idx:length(spk_rate)))/avg_secs;
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
                fprintf('\t-- plotting DUMMY panel %d for chan %s\n', 9, current_chan);

                plot(1:10, 1:10);
            end

        else
            fprintf('\t-- plotting EMPTY panel %d for chan %s\n', 9, current_chan);

            plot(1:10, 1:10);
            hold on;
            plot(1:10, fliplr(1:10));
            hold off;
        end

        hold off;

        %keyboard;
        
        %print([ms_figDir '/' session_name_no_spaces '_' current_chan '.png'], '-painters', '-dpng', '-r100');
        %saveas(gcf, [ms_figDir '/' session_name_no_spaces '_' current_chan '.pdf'], 'pdf');
        
        tightfig();
        set(gcf, 'InvertHardcopy', 'off','PaperUnits','inches','PaperPosition',[0,0,25,19],'PaperPositionMode','auto');
        print(gcf,[ms_figDir '/' session_name_no_spaces '_' current_chan '.png'],'-dpng','-r300');
        
        %export_fig([ms_figDir '/' session_name_no_spaces '_' current_chan '.bmp'], '-bmp', '-painters');
        fprintf('\t-- saving\n');
        close(gcf);

    end
    
    fprintf('done\n');
    %keyboard;
end

function plotUnitHist(uniq_nzID_current_channel, uniq_unitIDs_current_channel, nz_filters, unit_filters, ms_sim, ms_nim, unit_colors, unit_hist_width_scalar, unit_hist_height_scalar, saveFile)
 
     figure(2);
     
     font_size = 8;

     set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01, 0.01, 1, 0.96]);
     
     display_cols = 5;
     display_rows = 2 * ceil((length(uniq_nzID_current_channel) + length(uniq_unitIDs_current_channel))/display_cols);
     
     mini_col_count = 1;

     for n = 1:length(uniq_nzID_current_channel)

        current_nz = uniq_nzID_current_channel{n};
        nzIDs_current_nz_condi = nz_filters{ cellfun( @(x) isequal(x, current_nz), ms_nim.uniqueUnitID) };

        if sum(nzIDs_current_nz_condi) > 10
        
            ms_diffs = diff(ms_nim.timeStamp(nzIDs_current_nz_condi,:)/30);

            %histogram of spike times ISI for this unit
            hold off;
            current_plot = subplot(display_rows, display_cols, mini_col_count + display_cols * floor(mini_col_count/display_cols) ,'replace');
            fprintf('making hist for noise %d at plot point %d\n', n, mini_col_count + display_cols * floor(mini_col_count/display_cols));
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
             current_plot = subplot(display_rows, display_cols, (mini_col_count + display_cols * floor(mini_col_count/display_cols)) + display_cols ,'replace');
             hold on;

             fprintf('making hist for noise %d at plot point %d\n', n, (mini_col_count + display_cols * floor(mini_col_count/display_cols)) + display_cols);

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

             text(mean([0 ceil(max(ms_diffs))]) , mean([0 1]), txt, 'FontSize', 6);

        end
        
        mini_col_count = mini_col_count + 1;

     end

     
     for u = 1:length(uniq_unitIDs_current_channel)
         
         u_color = unit_colors( u, :);
         current_unit = uniq_unitIDs_current_channel{u};
         unitIDs_current_unit_condi = unit_filters{ cellfun( @(x) isequal(x, current_unit), ms_sim.uniqueUnitID) };

         if sum(unitIDs_current_unit_condi) > 10
         
             ms_diffs = diff(ms_sim.timeStamp(unitIDs_current_unit_condi,:)/30);

             %histogram of spike times ISI for this unit
             hold off;
             current_plot = subplot(display_rows, display_cols,  mini_col_count + display_cols * floor((mini_col_count-1)/display_cols) ,'replace');
             fprintf('making hist for unit %d at plot point %d\n', u,  mini_col_count + display_cols * floor((mini_col_count-1)/display_cols));
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
             current_plot = subplot(display_rows, display_cols, (mini_col_count + display_cols * floor((mini_col_count-1)/display_cols)) + display_cols ,'replace');
             hold on;

             fprintf('making hist for unit %d at plot point %d\n', u, (mini_col_count + display_cols * floor((mini_col_count-1)/display_cols)) + display_cols);

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

             text(mean([0 ceil(max(ms_diffs))]) , mean([0 1]), txt, 'FontSize', 6);
         
         end
         
         mini_col_count = mini_col_count + 1;

     end
     
     fprintf('saving seperate ISI plot\n');
     
     tightfig();
     set(gcf, 'InvertHardcopy', 'off','PaperUnits','inches','PaperPosition',[0,0,25,19],'PaperPositionMode','auto');
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

%%%%% old barplot code

%             xlabels = {'snr' 'isolation' 'noise overlap'};
%             ymx = 1;
% 
%             %get isolation, noise overlap, peak_snr
%             barplot_mat = zeros( 3, length(uniq_nzID_current_channel) + length(uniq_unitIDs_current_channel));
% 
%             cluster_count = 1;
% 
%             for n = 1:length(uniq_nzID_current_channel)
% 
%                 current_nz = uniq_nzID_current_channel{n};
%                 nzIDs_current_nz_condi = cellfun( @(x) isequal(x, current_nz), ms_nim.uniqueUnitID);
% 
%                 isolation = ms_nim.metrics.isolation(nzIDs_current_nz_condi);
%                 noise_overlap = ms_nim.metrics.noise_overlap(nzIDs_current_nz_condi);
%                 snr = ms_nim.metrics.peak_snr(nzIDs_current_nz_condi);
% 
%                 barplot_mat(1,cluster_count) = snr;
%                 barplot_mat(2,cluster_count) = isolation;
%                 barplot_mat(3,cluster_count) = noise_overlap;
% 
%                 if max(barplot_mat(:)) > ymx
%                     ymx = max(barplot_mat(:));
%                 end
% 
%                 cluster_count = cluster_count + 1;
%             end
% 
%              for u = 1:length(uniq_unitIDs_current_channel)
% 
%                  current_unit = uniq_unitIDs_current_channel{u};
%                  unitIDs_current_unit_condi = cellfun( @(x) isequal(x, current_unit), ms_sim.uniqueUnitID);
% 
%                  isolation = ms_sim.metrics.isolation(unitIDs_current_unit_condi);
%                  noise_overlap = ms_sim.metrics.noise_overlap(unitIDs_current_unit_condi);
%                  snr = ms_sim.metrics.peak_snr(unitIDs_current_unit_condi);
% 
%                  barplot_mat(1,cluster_count) = snr;
%                  barplot_mat(2,cluster_count) = isolation;
%                  barplot_mat(3,cluster_count) = noise_overlap;
% 
%                  if max(barplot_mat(:)) > ymx
%                      ymx = max(barplot_mat(:));
%                  end
% 
%                  cluster_count = cluster_count + 1;
% 
%              end
% 
%              b = bar(1:3,barplot_mat);
%              pause(0.5);
%              bar_width = get(b, 'BarWidth');
% 
%              set(gca, 'XTickLabel',xlabels, 'XTick',1:numel(xlabels));
% 
%              for g = 1:size(barplot_mat,1)
%                  for v = 1:size(barplot_mat,2)
% 
%                      current_val = barplot_mat(g,v);
% 
%                      txt_color = 'k';
% 
%                      if g == 1 %snr > 1.5
% 
%                          if current_val < 1.5
%                              txt_color = 'r';
%                          end
% 
%                      elseif g == 2 % isolation > 0.95
% 
%                          if current_val < 0.95
%                              txt_color = 'r';
%                          end
% 
%                      elseif g == 3 % noise overlap < 0.03
% 
%                          if current_val > 0.03
%                              txt_color = 'r';
%                          end
% 
%                      end
% 
%                      txt = sprintf('%0.3f', barplot_mat(g,v));
%                      tx = text(b(v).XData(g) + b(v).XOffset, b(v).YData(g) + 0.5, txt, 'Color', txt_color, 'HorizontalAlignment', 'center', 'FontSize', 6); 
%                      set(tx,'Rotation',45);
%                      
%                  end
%              end
% 
%              %set bar colors
%              for n = 1:length(uniq_nzID_current_channel)
%                  set(b(n), 'FaceColor', 'w');
%              end
% 
%              for u = (length(uniq_nzID_current_channel)+1):(length(uniq_nzID_current_channel) + length(uniq_unitIDs_current_channel))
%                  u_color = unit_colors( u - length(uniq_nzID_current_channel) );
%                  set(b(u), 'FaceColor', u_color);
%              end
% 
%              hold on;
%              plot([2 3], [1 1], 'k--', 'LineWidth', 1);
% 
%              fprintf('%d and %d\n', 0, ymx + 0.5);
%              ylim([0 ymx + 0.5]);


