function extract_psd_var(varargin)

fprintf('****************************************************\n');
fprintf('* %s\n', mfilename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse inputs

p = inputParser;

p.addParameter('sess_path', '', @ischar);
p.addParameter('min_per_window', '5', @ischar);

parse(p, varargin{:});

disp(p.Results);

sess_path = p.Results.sess_path;
min_per_window = eval(p.Results.min_per_window);

if ~exist(sess_path, 'dir')
    error('%s is not a valid path', sess_path);
end

psd_path = [ sess_path '/psd' ];

if ~exist(psd_path, 'dir')
   mkdir(psd_path); 
end


sess_dir_name = strsplit(sess_path, '/');
sess_path_splits = strsplit(sess_dir_name{end}, '_');
sess_time_str = strjoin(sess_path_splits(1:(end-1)), '_');
% sess_time_str = [sess_path_splits{1} '_' sess_path_splits{2}];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load noreref

noreref_fpath = [sess_path '/raw/' sess_time_str '_noreref.mat'];
lfpStruct = load(noreref_fpath);
lfpStruct = lfpStruct.lfpStruct;

% load processed

processed_fpath = [sess_path '/raw/' sess_time_str '_processed.mat'];
processed = load(processed_fpath);
processed = processed.lfpStruct;

% count NaNs on each device

frac_nan = [];
unique_microDevNums = unique(processed.chanIDperNSP{1}{:, 'MicroDevNum'});

for iDev = 1:length(unique_microDevNums)

    current_dev = unique_microDevNums(iDev);
    current_dev_filt = (processed.chanIDperNSP{1}{:, 'MicroDevNum'} == current_dev);
    current_dev_data = processed.lfp{1}(current_dev_filt, :);
    current_dev_data_frac_nan = sum(current_dev_data(:) == 0)/(size(current_dev_data, 1)*size(current_dev_data, 2));
    
    frac_nan(iDev) = current_dev_data_frac_nan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find intervals

% how many windows to use? 6 min windows = 10 per hour

if lfpStruct.sessDurMin < min_per_window
    error('lfpStruct.sessDurMin < min_per_window');
end


lfp = double(lfpStruct.lfp{1});
lfp_sample_num = size(lfp, 2);

samples_per_min = lfpStruct.samplingFreq * 60;
samples_per_window = samples_per_min * min_per_window;

samples_interval_starts = 1:samples_per_window:lfp_sample_num;
samples_interval_stops = samples_interval_starts + samples_per_window - 1;

if samples_interval_stops(end) > lfp_sample_num
   
    samples_interval_starts = samples_interval_starts(1:(end-1));
    samples_interval_stops = samples_interval_stops(1:(end-1));
end

samples_intervals = [samples_interval_starts' samples_interval_stops'];
samples_intervals_midpoints = samples_interval_starts + samples_per_window/2;

fprintf('sessDurMin: %0.2f\n', lfpStruct.sessDurMin);
for iInterval = 1:length(samples_interval_starts)
    fprintf('interval %d: %0.2f - %0.2f min\n', iInterval, samples_interval_starts(iInterval)/samples_per_min, samples_interval_stops(iInterval)/samples_per_min);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate power for intervals

[noreref_psd, f] = extract_psd(lfp, samples_intervals, lfpStruct.samplingFreq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate variance for intervals

[noreref_var, noreref_rms] = extract_var(lfp, samples_intervals, lfpStruct.samplingFreq);

createdDate = datestr(datetime('now'));

noreref_readme = ['this psd.mat file, generated ' sprintf('%s', createdDate) ', contains the following fields:' newline ...
          '     createdDate          - a string indicating when this mat file was created' newline ...
          '     chanID               - a table corresponding to the channel data in psd, var, and rms' newline ...
          '     frac_nan             - vector with value for each microDevNum (indexed in sorted order) indicating the fraction of NaNs found in processed.mat' newline ...
          '     psd                  - a (chan x timepoint x freq) matrix with power values in decibels at each time window (described by timepoints min and window_min), calculated with pwelch' newline ...
          '     var                  - a (chan x timepoint x subwindow) var(chan, timepoint, :) describes the distribution of variance values calculated using non-overlapping subwindows for the indexed timepoint and channel' newline ...
          '     rms                  - a (chan x timepoint x subwindow) rms(chan, timepoint, :) describes the distribution of rms values calculated using non-overlapping subwindows for the indexed timepoint and channel' newline ...
          '     freqs                - a vector containing frequency labels for values in the 3rd dimension of psd' newline ...
          '     timepoints_samples   - a vector with timepoints (in unit samples) corresponding the center of each time window used in psd, var, and rms' newline ...
          '     timepoints_min       - a vector with timepoints (in unit minutes) corresponding the center of each time window used in psd, var, and rms' newline ...          
          '     window_min           - the number of minutes of data used for each time window. Windows are non-overlapping' newline ...
          '     sessDurMin           - session duration in minutes' newline ...
          '     sessStr              - the session name the lfps were extracted from (e.g., 190117_1336)' newline ...
          '     samplingFreq         - samlple frequency of the data. typically 1000 Hz' newline ...
          '     rerefType            - a string indicating which lfp_type this psd.mat is derived from' newline ...
          '     jackTableUsed        - just the jacksheet for the channels incorporated into this psdStruct. this table can be used to go back and forth between original and new channel names' newline ...
          '     startTime_datenum    - the start time of session as output from datenum function (note that these values are taken from the timestamp in the nsx filename, not the incorrectly offset time values in the original nsx file)' ];

psdStruct = struct;
psdStruct.readme = noreref_readme;
psdStruct.createdDate = createdDate;
psdStruct.chanID = lfpStruct.chanIDperNSP{1};
psdStruct.frac_nan = frac_nan;
psdStruct.psd = noreref_psd;
psdStruct.var = noreref_var;
psdStruct.rms = noreref_rms;
psdStruct.freqs = f;
psdStruct.timepoints_samples = samples_intervals_midpoints;
psdStruct.timepoints_min = round(samples_intervals_midpoints ./ samples_per_min, 2);
psdStruct.window_min = min_per_window;
psdStruct.sessDurMin = lfpStruct.sessDurMin;
psdStruct.sessStr = lfpStruct.sessStr;
psdStruct.samplingFreq = lfpStruct.samplingFreq;
psdStruct.rerefType = 'lfp_noreref';
psdStruct.jackTableUsed = lfpStruct.jackTableUsed;
psdStruct.startTime_datenum = lfpStruct.startTime_datenum;

noreref_psd_fpath = [ psd_path '/noreref.mat' ];
save(noreref_psd_fpath, 'psdStruct');

%make plots

plot_psd(psdStruct, psd_path);


end

function plot_psd(psdStruct, psd_path)

num_timepoints = length(psdStruct.timepoints_min);

subplot_cols = 3;
subplot_rows = ceil(num_timepoints/subplot_cols);

num_unique_microDevNum = unique(psdStruct.chanID{:,'MicroDevNum'});

for iDevNum = 1:length(num_unique_microDevNum)
    
    current_microDevNum = num_unique_microDevNum(iDevNum);
    current_microDevNum_filt = (psdStruct.chanID{:,'MicroDevNum'} == current_microDevNum);

    figure();
    set(gcf,'color','w');
    set(gcf,'visible','off');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01, 0.01, 1, 0.96]);


    for iTime = 1:num_timepoints

        subplot(subplot_rows, subplot_cols, iTime);

        current_psd = squeeze(psdStruct.psd(current_microDevNum_filt, iTime, :));

        plot(psdStruct.freqs, current_psd);

        xlim([0 200]);

        if iTime == 1
            xlabel('freq (Hz)');
            ylabel('power (decibels)');
        end

        title(sprintf('minute %0.2f', psdStruct.timepoints_min(iTime) ));

    end

    tightfig();
    set(gcf, 'InvertHardcopy', 'off','PaperUnits','inches','PaperPosition',[0,0,25,19],'PaperPositionMode','auto');
    print(gcf,[psd_path '/' sprintf('microDev%d_psd.png', current_microDevNum)],'-dpng','-r300');

end

end



function [psd, f] = extract_psd(lfp, samples_intervals, samplingFreq)

lfp_channel_num = size(lfp, 1);
num_intervals = size(samples_intervals, 1);

pwelch_window = samplingFreq * 3; % 3 second window size
pwelch_overlap = fix(pwelch_window/2); % 50 percent overlap


% get frequency vector

current_interval = samples_intervals(1, :);
current_interval_data = lfp(:, current_interval(1):current_interval(2) );
[~,f] = pwelch(current_interval_data', pwelch_window, pwelch_overlap, [], samplingFreq);



psd = zeros(lfp_channel_num, num_intervals, length(f));

for iInterval = 1:num_intervals
   
    current_interval = samples_intervals(iInterval, :);
    
    current_interval_data = lfp(:, current_interval(1):current_interval(2) );
   
    fprintf('computing psd: interval %d of %d\n', iInterval, num_intervals);
    [pxx,~] = pwelch(current_interval_data', pwelch_window, pwelch_overlap, [], samplingFreq);
    
    power_decibels = 10*log10(pxx);
    
    for iChan = 1:lfp_channel_num
        psd(iChan, iInterval, :) = power_decibels(:, iChan);
    end
end

end


function [noreref_var, noreref_rms] = extract_var(lfp, samples_intervals, samplingFreq)

window_size = samplingFreq * 3; % 3 second window size

% how many subwindows per interval

subwindow_count = 0;
window_start = samples_intervals(1,1);
while (window_start + window_size) <= samples_intervals(1,2)
    
    window_start = window_start + window_size;
    subwindow_count = subwindow_count + 1;
end

fprintf('sample interval length: %d samples (%0.2f min), window_size %d samples, subwindow count %d\n', samples_intervals(1,2) - samples_intervals(1,1) ,(samples_intervals(1,2) - samples_intervals(1,1))/(samplingFreq*60), window_size, subwindow_count);

lfp_channel_num = size(lfp, 1);
num_intervals = size(samples_intervals, 1);

noreref_var = zeros(lfp_channel_num, num_intervals, subwindow_count);
noreref_rms = zeros(lfp_channel_num, num_intervals, subwindow_count);

for iInterval = 1:num_intervals

    fprintf('computing variance and rms: interval %d of %d\n', iInterval, num_intervals);
    
    current_interval = samples_intervals(iInterval, :);
    
    window_start = current_interval(1);
    window_end = window_start + window_size;
    
    for iSubWin = 1:subwindow_count
    
        current_interval_data = lfp(:, window_start:window_end );

        noreref_var(:, iInterval, iSubWin) = var(current_interval_data, 0, 2);

        noreref_rms(:, iInterval, iSubWin) = rms(current_interval_data, 2);

    end
end

end