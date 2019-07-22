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


% sess_dir_name = strsplit(sess_path, '/');
% sess_path_splits = strsplit(sess_dir_name{end}, '_');
% sess_time_str = strjoin(sess_path_splits(1:2), '_');
% sess_time_str = [sess_path_splits{1} '_' sess_path_splits{2}];

noreref_ls = dir([sess_path '/raw/*_noreref.mat']);
processed_ls = dir([sess_path '/raw/*_processed.mat']);
variance_ls = dir([sess_path '/cleaning/variance.csv']);

if isempty(variance_ls)
   error('command dir([sess_path "/cleaning/variance.csv"]) found no results'); 
end

if isempty(noreref_ls)
   error('command dir([sess_path "/raw/*_noreref.mat"]) found no results'); 
end

if isempty(processed_ls)
   error('command dir([sess_path "/raw/*_processed.mat"]) found no results'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load variance.csv

cleaning_info_fpath = [variance_ls.folder '/' variance_ls.name];
cleaning_info = readtable(cleaning_info_fpath);

if all(~cleaning_info.is_good)
    ignore_fid = fopen([psd_path '/_ignore_me.txt'], 'w');
    fprintf(ignore_fid, 'all channels marked as not not_good in variance.csv\n');
    fclose(ignore_fid);
    return;
end

% load noreref

noreref_fpath = [noreref_ls.folder '/' noreref_ls.name];
noreref = load(noreref_fpath);
noreref = noreref.lfpStruct;

% load processed

processed_fpath = [processed_ls.folder '/' processed_ls.name];
processed = load(processed_fpath);
processed = processed.lfpStruct;


% count NaNs on each device

old_NaN_value = 0;

frac_nan = [];
unique_microDevNums = unique(processed.chanIDperNSP{1}{:, 'MicroDevNum'});
microDevFilts = [];

for iDev = 1:length(unique_microDevNums)

    current_dev = unique_microDevNums(iDev);
    current_dev_filt = (processed.chanIDperNSP{1}{:, 'MicroDevNum'} == current_dev);
    
    microDevFilts = [ microDevFilts current_dev_filt ];
    
    if ~isfield(processed, 'nan_mask')
    
        current_dev_data = processed.lfp{1}(current_dev_filt, :);
        current_dev_frac_nan = sum(current_dev_data(:) == old_NaN_value)/(size(current_dev_data, 1)*size(current_dev_data, 2));
        
    else
       
        current_dev_nan_mask = processed.nan_mask{1}(current_dev_filt, :);
        current_dev_frac_nan = sum(current_dev_nan_mask(:))/(size(current_dev_nan_mask, 1)*size(current_dev_nan_mask, 2));
        
    end
    
    frac_nan(iDev) = current_dev_frac_nan;
end

microDevFilts = logical(microDevFilts);

if ~isfield(processed, 'nan_mask')
   nan_mask = (processed.lfp{1} == old_NaN_value);
else
   nan_mask = processed.nan_mask{1}; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find intervals

% how many windows to use? 6 min windows = 10 per hour

if noreref.sessDurMin < min_per_window
    error('lfpStruct.sessDurMin < min_per_window');
end

noreref_lfp = double(noreref.lfp{1});
processed_lfp = double(processed.lfp{1});

lfp_sample_num = size(noreref_lfp, 2);

samples_per_min = noreref.samplingFreq * 60;

samples_per_window = samples_per_min * min_per_window;

minutes_per_day = (60*24);
% windows_per_day = minutes_per_day/min_per_window;

% these are the windows i want to use across all sessions
window_start_times_minute = 0:min_per_window:minutes_per_day;
window_start_times_minute(end) = [];

window_start_times_samples = window_start_times_minute .* samples_per_min;
window_stop_times_samples = window_start_times_samples + samples_per_window;

% what samples numbers does our data correspond to in the daily count
sessStart_datetime = datetime(noreref.startTime_datenum,'ConvertFrom','datenum');
sessStart_reference_datetime = datetime(year(sessStart_datetime), month(sessStart_datetime), day(sessStart_datetime));

sessStart_minutes_since_midnight = minutes(sessStart_datetime - sessStart_reference_datetime);
sessStart_samples_since_midnight = sessStart_minutes_since_midnight * samples_per_min;

sess_samples_since_midnight = (0:lfp_sample_num-1) + sessStart_samples_since_midnight;



timepoints_daily_window = [];
samples_interval_starts = [];
samples_interval_stops = [];

samples_captured = 0;

% which windows are present in this session?

minWindowCoverage = 0.95; % a window must have 95% coverage from the data in this session to be considered as filled

for iWin = 1:length(window_start_times_samples)
    
    current_window_start = window_start_times_samples(iWin);
    current_window_stop = window_stop_times_samples(iWin);
    
    % does any data fall into this window?
    current_window_start_filt = (sess_samples_since_midnight >= current_window_start);
    current_window_stop_filt = (sess_samples_since_midnight < current_window_stop);
    
    current_window_filt = current_window_start_filt & current_window_stop_filt;
    
    if sum(current_window_filt)/samples_per_window >= minWindowCoverage
       %found a window!
              
       samples_captured = samples_captured + sum(current_window_filt);
       
       current_window_start_hours = current_window_start/(samples_per_min * 60);
       current_window_start_hour_mark = floor(current_window_start_hours);
       current_window_start_min_mark = round(60 * (current_window_start_hours - current_window_start_hour_mark));
       
       current_window_stop_hours = current_window_stop/(samples_per_min * 60);
       current_window_stop_hour_mark = floor(current_window_stop_hours);
       current_window_stop_min_mark = round(60 * (current_window_stop_hours - current_window_stop_hour_mark));
       
       fprintf('found %d samples in window %02d:%02d - %02d:%02d\n', sum(current_window_filt), current_window_start_hour_mark, current_window_start_min_mark, current_window_stop_hour_mark, current_window_stop_min_mark);
       
       timepoints_daily_window(length(timepoints_daily_window) + 1) = window_start_times_minute(iWin);
       samples_interval_starts(length(samples_interval_starts) + 1) = find(current_window_start_filt, 1, 'first');
       samples_interval_stops(length(samples_interval_stops) + 1) = find(current_window_stop_filt, 1, 'last');
    end
    
end

frac_samples_covered = samples_captured/lfp_sample_num;

if samples_captured == 0
    ignore_fid = fopen([psd_path '/_ignore_me.txt'], 'w');
    fprintf(ignore_fid, 'no daily time windows reach %0.2f percent coverage in this session\n', minWindowCoverage);
    fclose(ignore_fid);
    return;
end

samples_intervals = [samples_interval_starts' samples_interval_stops'];
samples_intervals_midpoints = samples_interval_starts + samples_per_window/2;

fprintf('sessDurMin: %0.2f\n', noreref.sessDurMin);
for iInterval = 1:length(samples_interval_starts)
    fprintf('interval %d: %0.2f - %0.2f min\n', iInterval, samples_interval_starts(iInterval)/samples_per_min, samples_interval_stops(iInterval)/samples_per_min);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply reref to processed

for iDev = 1:length(unique_microDevNums)

    current_dev = unique_microDevNums(iDev);
    current_dev_filt = (processed.chanIDperNSP{1}{:, 'MicroDevNum'} == current_dev);

    current_global_mean = double(processed.glob_sig_good{current_dev});
    
    current_global_mean(isnan(current_global_mean)) = 0;
    
    for iFilt = 1:length(current_dev_filt)
        if current_dev_filt(iFilt)
            processed_lfp(iFilt, :) = processed_lfp(iFilt, :) - current_global_mean';
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate power for intervals

[noreref_psd, f] = extract_psd(noreref_lfp, samples_intervals, noreref.samplingFreq);
[processed_psd, ~] = extract_psd(processed_lfp, samples_intervals, noreref.samplingFreq, cleaning_info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate variance for intervals

[noreref_var, noreref_rms, interval_frac_nan] = extract_var(noreref_lfp, samples_intervals, noreref.samplingFreq, microDevFilts, nan_mask);
[processed_var, processed_rms, ~] = extract_var(processed_lfp, samples_intervals, processed.samplingFreq, microDevFilts, nan_mask);

createdDate = datestr(datetime('now'));

noreref_readme = ['this psd.mat file, generated ' sprintf('%s', createdDate) ', contains the following fields:' newline ...
          '     createdDate              - a string indicating when this mat file was created' newline ...
          '     chanID                   - a table corresponding to the channel data in psd, var, and rms' newline ...
          '     frac_nan                 - vector with value for each microDevNum (indexed in sorted order) indicating the fraction of NaNs found in processed.mat' newline ...
          '     interval_nan_frac        - a (num_microDev x timepoint) matrix with fraction of NaNs in each timepoint per microDev' newline ...
          '     psd (noreref, processed) - a (chan x timepoint x freq) matrix with power values in decibels at each time window (described by timepoints min and window_min), calculated with pwelch. Processed lfp had glob_sig_good subtracted before psd calc' newline ...
          '     var (noreref, processed) - a (chan x timepoint x subwindow) var(chan, timepoint, :) describes the distribution of variance values calculated using non-overlapping subwindows for the indexed timepoint and channel. Processed lfp had glob_sig_good subtracted before psd calc' newline ...
          '     rms (noreref, processed) - a (chan x timepoint x subwindow) rms(chan, timepoint, :) describes the distribution of rms values calculated using non-overlapping subwindows for the indexed timepoint and channel. Processed lfp had glob_sig_good subtracted before psd calc' newline ...
          '     cleaning_info            - copy of the variance.csv that is produced from the lfp pipeline' newline ...
          '     freqs                    - a vector containing frequency labels for values in the 3rd dimension of psd' newline ...
          '     timepoints_samples       - a vector with timepoints (in unit samples) corresponding the center of each time window used in psd, var, and rms' newline ...
          '     timepoints_min           - a vector with timepoints (in unit minutes) corresponding the center of each time window used in psd, var, and rms' newline ...
          '     timepoints_daily_window  - a vector with timepoints (in minutes) corresponding to the beginning of each window calculated since the beginning of the day (daily windows calculated 0:window_min:minutes_per_day)' newline ...
          '     frac_samples_covered     - fraction of samples in the lfp matrix that were used with the specfied window_min' newline ...
          '     minWindowCoverage        - minimum fraction of the daily window that must be present in the data for the data to count towards that window' newline ...
          '     window_min               - the number of minutes of data used for each time window. Windows are non-overlapping' newline ...
          '     sessDurMin               - session duration in minutes' newline ...
          '     sessStr                  - the session name the lfps were extracted from (e.g., 190117_1336)' newline ...
          '     samplingFreq             - samlple frequency of the data. typically 1000 Hz' newline ...
          '     rerefType                - a string indicating which lfp_type this psd.mat is derived from' newline ...
          '     jackTableUsed            - just the jacksheet for the channels incorporated into this psdStruct. this table can be used to go back and forth between original and new channel names' newline ...
          '     startTime_datenum        - the start time of session as output from datenum function (note that these values are taken from the timestamp in the nsx filename, not the incorrectly offset time values in the original nsx file)' ];

psdStruct = struct;
psdStruct.readme = noreref_readme;
psdStruct.createdDate = createdDate;
psdStruct.chanID = noreref.chanIDperNSP{1};
psdStruct.frac_nan = frac_nan;
psdStruct.interval_frac_nan = interval_frac_nan;
psdStruct.noreref_psd = noreref_psd;
psdStruct.noreref_var = noreref_var;
psdStruct.noreref_rms = noreref_rms;
psdStruct.processed_psd = processed_psd;
psdStruct.processed_var = processed_var;
psdStruct.processed_rms = processed_rms;
psdStruct.cleaning_info = cleaning_info;
psdStruct.freqs = f;
psdStruct.timepoints_samples = samples_intervals_midpoints;
psdStruct.timepoints_min = round(samples_intervals_midpoints ./ samples_per_min, 2);
psdStruct.timepoints_daily_window = timepoints_daily_window;
psdStruct.window_min = min_per_window;
psdStruct.frac_samples_covered = frac_samples_covered;
psdStruct.minWindowCoverage = minWindowCoverage;
psdStruct.sessDurMin = noreref.sessDurMin;
psdStruct.sessStr = noreref.sessStr;
psdStruct.samplingFreq = noreref.samplingFreq;
psdStruct.rerefType = 'lfp_noreref';
psdStruct.jackTableUsed = noreref.jackTableUsed;
psdStruct.startTime_datenum = noreref.startTime_datenum;

noreref_psd_fpath = [ psd_path '/psd.mat' ];
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % noreref    
    figure();
    set(gcf,'color','w');
    set(gcf,'visible','off');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01, 0.01, 1, 0.96]);


    for iTime = 1:num_timepoints

        subplot(subplot_rows, subplot_cols, iTime);

        current_psd = squeeze(psdStruct.noreref_psd(current_microDevNum_filt, iTime, :));
        current_timepoint_frac_nan = psdStruct.interval_frac_nan(iDevNum, iTime);

        plot(psdStruct.freqs, current_psd);

        xlim([0 200]);

        if iTime == 1
            xlabel('freq (Hz)');
            ylabel('power (decibels)');
        end

        title(sprintf('minute %0.2f (frac nan %0.2f)', psdStruct.timepoints_min(iTime), current_timepoint_frac_nan));

    end

    tightfig();
    set(gcf, 'InvertHardcopy', 'off','PaperUnits','inches','PaperPosition',[0,0,25,19],'PaperPositionMode','auto');
    print(gcf,[psd_path '/' sprintf('noreref_microDev%d_psd.png', current_microDevNum)],'-dpng','-r300');
    close(gcf);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % processed
    figure();
    set(gcf,'color','w');
    set(gcf,'visible','off');
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01, 0.01, 1, 0.96]);


    for iTime = 1:num_timepoints

        subplot(subplot_rows, subplot_cols, iTime);

        current_psd = squeeze(psdStruct.processed_psd(current_microDevNum_filt, iTime, :));
        current_timepoint_frac_nan = psdStruct.interval_frac_nan(iDevNum, iTime);

        plot(psdStruct.freqs, current_psd);

        xlim([0 200]);

        if iTime == 1
            xlabel('freq (Hz)');
            ylabel('power (decibels)');
        end

        title(sprintf('minute %0.2f (frac nan %0.2f)', psdStruct.timepoints_min(iTime), current_timepoint_frac_nan));

    end

    tightfig();
    set(gcf, 'InvertHardcopy', 'off','PaperUnits','inches','PaperPosition',[0,0,25,19],'PaperPositionMode','auto');
    print(gcf,[psd_path '/' sprintf('processed_microDev%d_psd.png', current_microDevNum)],'-dpng','-r300');
    close(gcf);
end

end



function [psd, f] = extract_psd(lfp, samples_intervals, samplingFreq, cleaning_info)

if nargin < 4
   cleaning_info = {}; 
end


lfp_channel_num = size(lfp, 1);
num_intervals = size(samples_intervals, 1);


pwelch_window = samplingFreq * 3; % 3 second window size
pwelch_overlap = fix(pwelch_window/2); % 50 percent overlap


% get frequency vector

current_interval = samples_intervals(1, :);
current_interval_data = lfp(:, current_interval(1):current_interval(2) );
[~,f] = pwelch(current_interval_data', pwelch_window, pwelch_overlap, [], samplingFreq);



psd = nan(lfp_channel_num, num_intervals, length(f));

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

% nan bad channels if passed in
if ~isempty(cleaning_info)
   psd(~cleaning_info.is_good, :, :) = NaN; 
end

end


function [noreref_var, noreref_rms, interval_nan_frac] = extract_var(lfp, samples_intervals, samplingFreq, microDevFilts, nan_mask)

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
num_dev = size(microDevFilts, 2);

noreref_var = nan(lfp_channel_num, num_intervals, subwindow_count);
noreref_rms = nan(lfp_channel_num, num_intervals, subwindow_count);
interval_nan_frac = zeros(num_dev, num_intervals);

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
    
    for iDev = 1:num_dev
    
        current_microDevFilt = microDevFilts(:, iDev);
        current_interval_nans = nan_mask(current_microDevFilt, samples_intervals(iInterval,1):samples_intervals(iInterval,2) );
    
        current_interval_frac_nan = sum(current_interval_nans(:))/length(current_interval_nans(:));
        interval_nan_frac(iDev, iInterval) = current_interval_frac_nan;
    end
    
end

end