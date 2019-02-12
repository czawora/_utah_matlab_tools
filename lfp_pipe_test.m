splitLFP('nsx_fpath','/Volumes/56E/UTAH_E/NIH066/data_raw/181024_1539/QaHkESkEkA-20181024-153949-INST0.ns5', 'save_dir', '/Volumes/CZ_FRNU/lfp_test/NIH066');

for i = 1:128
   
    chan_fpath = sprintf('/Volumes/CZ_FRNU/lfp_test/NIH066/lfp_splits/ch_%03d.mat', i);
    filter_downsample_detrend_rmLineNoise('channel_fpath', chan_fpath);    
end


variance_and_lineNoise_exclusion('nsx_physio_fpath', '/Volumes/56E/UTAH_E/NIH066/data_raw/181024_1539/QaHkESkEkA-20181024-153949-INST0.ns5', 'nev_fpath', '/Volumes/56E/UTAH_E/NIH066/data_raw/181024_1539/QaHkESkEkA-20181024-153949-INST0.nev', 'ns3_pulse_fpath', '/Volumes/56E/UTAH_E/NIH066/data_raw/181024_1539/QaHkESkEkA-20181024-153949-INST0.ns3', 'subj_str', 'NIH066', 'time_str', '123456_7890', 'nsp_str', 'utah1', 'split_path', '/Volumes/CZ_FRNU/lfp_test/NIH066/lfp_splits', 'elementInfo_fpath', '/Volumes/CZ_FRNU/lfp_test/NIH066/NIH066_elementInfo.txt');


splitLFP('nsx_fpath','/Volumes/56B/UTAH_B/NIH047/data_raw/170311_1112/170311_1112_utah.ns5', 'save_dir', '/Volumes/CZ_FRNU/lfp_test/NIH047');

for i = 1:96
   
    chan_fpath = sprintf('/Volumes/CZ_FRNU/lfp_test/NIH047/lfp_splits/ch_%03d.mat', i);
    filter_downsample_detrend_rmLineNoise('channel_fpath', chan_fpath);    
end


variance_and_lineNoise_exclusion('nsx_physio_fpath', '/Volumes/56B/UTAH_B/NIH047/data_raw/170311_1112/170311_1112_utah.ns5', 'nev_fpath', '/Volumes/56B/UTAH_B/NIH047/data_raw/170311_1112/170311_1112_utah.nev', 'ns3_pulse_fpath', '/Volumes/56B/UTAH_B/NIH047/data_raw/170311_1112/170311_1112_utah.ns4', 'subj_str', 'NIH047', 'time_str', '123456_7890', 'nsp_str', 'utah', 'split_path', '/Volumes/CZ_FRNU/lfp_test/NIH047/lfp_splits', 'elementInfo_fpath', '/Volumes/CZ_FRNU/lfp_test/NIH047/NIH047_elementInfo.txt');




construct_spikeInfoMS_edits('saveRoot','/Volumes/CZ_FRNU/spike_test/NIH034/NIH034_151001_1317_utah2/spike/outputs', 'nev_fpath', '/Volumes/CZ_FRNU/spike_test/NIH034/NIH034_151001_1317_utah2/NIH034_151001_1317_utah2.nev', 'elementInfo_fpath', '/Volumes/CZ_FRNU/spike_test/NIH034/NIH034_151001_1317_utah2/NIH034_151001_1317_utah2_elementInfo.txt', 'nsx_physio_fpath', '/Volumes/CZ_FRNU/spike_test/NIH034/NIH034_151001_1317_utah2/NIH034_151001_1317_utah2.ns6', 'bp_fname_suffix', 'mda_chan', 'subj_str', 'NIH034', 'time_str', '151001_1317', 'nsp_str', 'utah2', 'sessRoot', '/Volumes/CZ_FRNU/spike_test/NIH034/NIH034_151001_1317_utah2/spike/splits');


p.addParameter('subj_str', '', @ischar);
p.addParameter('time_str', '', @ischar);
p.addParameter('nsp_str', '', @ischar);

p.addParameter('sessRoot', '', @ischar);
p.addParameter('bp_fname_suffix', '', @ischar);
p.addParameter('nsx_physio_fpath', '', @ischar);
p.addParameter('ns3_pulse_fpath', '', @ischar);
p.addParameter('elementInfo_fpath', '', @ischar);
p.addParameter('nev_fpath', '', @ischar);
p.addParameter('num_refset', '1', @ischar);
p.addParameter('saveRoot', '', @ischar);
p.addParameter('removeLargeAmpUnits', '0', @ischar);

p.addParameter('firings_fname', 'firings.mda', @ischar);
p.addParameter('metrics_fname', 'metrics.json', @ischar);
p.addParameter('isol_metrics_fname', 'isol_metrics.json', @ischar);
p.addParameter('isol_pair_metrics_fname', 'isol_pair_metrics.json', @ischar);
p.addParameter('clips_fname', 'clips.mda', @ischar);