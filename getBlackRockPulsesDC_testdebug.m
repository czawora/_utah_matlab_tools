
%test files
nsx = concatOpenNSx()
pulses = getBlackrockPulses_DC_AN('postProc', , 'ns3_fpath', '/Volumes/56D/UTAH_D/NIH068/data_raw/181209_0135/PiGlkGDWmEEiGjDfG-20181209-013554-INST1.ns3', 'nev_fpath', '/Volumes/56D/UTAH_D/NIH068/data_raw/181209_0135/PiGlkGDWmEEiGjDfG-20181209-013554-INST1.nev');
pulses = getBlackrockPulses_DC_AN('ns3_fpath', '/Volumes/56C/UTAH_C/NIH062/data_raw/180923_0027/YmCiLtiHHid-20180923-002702-INST1.ns3', 'nev_fpath', '/Volumes/56C/UTAH_C/NIH062/data_raw/180923_0027/YmCiLtiHHid-20180923-002702-INST1.nev');


max_din = max(pulses.din4_uptimes);
max_ain = max(pulses.ain20_uptimes);

ain_past_din = sum(pulses.ain20_uptimes > max_din);
din_past_ain = sum(pulses.din4_uptimes > max_ain);

fprintf('ain_past_din %d\n', ain_past_din);
fprintf('din_past_ain %d\n', din_past_ain);


%xcorr
[r,lags] = xcorr(pulses.ain20_ts, pulses.din4_ts, 100);

figure();
stem(lags,r);
title('xcorr');


% uptimes plot
npulse_ain = length(pulses.ain20_uptimes);
npulse_din = length(pulses.din4_uptimes);
figure();
plot(1:npulse_ain, pulses.ain20_uptimes, 'k.', 'MarkerSize', 13);
hold on;
plot(1:npulse_din, pulses.din4_uptimes , 'r*');
title('uptimes');


% shifted uptimes
npulse_ain = length(pulses.ain20_uptimes(2:end));
npulse_din = length(pulses.din4_uptimes);
figure();
plot(1:npulse_ain, pulses.ain20_uptimes(2:end), 'k.', 'MarkerSize', 13);
hold on;
plot(1:npulse_din, pulses.din4_uptimes , 'r*');
title('shifted uptimes');


% range plot
sample_plot_range = 1:2000;
%sample_plot_range = (5150744-1000):(5150744);

figure();
plot(sample_plot_range, pulses.ain20_ts(sample_plot_range) ./ max(pulses.ain20_ts(sample_plot_range)), 'k-');
hold on;
plot(sample_plot_range, pulses.din4_ts(sample_plot_range) ./ max( pulses.din4_ts(sample_plot_range)), 'r--', 'LineWidth',3);


%diff hist
%how close is an ain uptime to its nearest din uptime
min_diffs = [];

for iAin = 1:length(pulses.ain20_uptimes)
   
    pulse_diff = pulses.ain20_uptimes(iAin) - pulses.din4_uptimes;
    abs_pulse_diff = abs(pulse_diff);
    [~, min_diff_idx] = min(abs_pulse_diff);
    min_diff_val = pulse_diff(min_diff_idx);
    min_diffs(iAin) = min_diff_val;
    
end

figure();
histogram(min_diffs, 'BinWidth', 1);
title('diff hist');



%diff hist, skip 1 ain
%how close is an ain uptime to its nearest din uptime
min_diffs = [];

for iAin = 10:length(pulses.ain20_uptimes)
   
    pulse_diff = pulses.ain20_uptimes(iAin) - pulses.din4_uptimes;
    abs_pulse_diff = abs(pulse_diff);
    [~, min_diff_idx] = min(abs_pulse_diff);
    min_diff_val = pulse_diff(min_diff_idx);
    min_diffs(iAin) = min_diff_val;
    
end

figure();
histogram(min_diffs, 'BinWidth', 1);
title('diff hist, skip 1 ain');


