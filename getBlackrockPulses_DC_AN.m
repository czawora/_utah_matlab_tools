function pulses = getBlackrockPulses_DC_AN(varargin)

    %output struct
    pulses = struct;
        
    p = inputParser;
        
    p.addParameter('ns3_fpath', '', @ischar);
    p.addParameter('nev_fpath', '', @ischar);
    p.addParameter('postProc', []);

    parse(p, varargin{:});
    
    %record filenames passed in
    ns3_fpath = p.Results.ns3_fpath;
    nev_fpath = p.Results.nev_fpath;
    postProc = p.Results.postProc;
    
    ns3_fpath_splits = strsplit(ns3_fpath, '/');
    nev_fpath_splits = strsplit(nev_fpath, '/');
    
    pulses.ain_filename = ns3_fpath_splits{end};
    pulses.din_filename = nev_fpath_splits{end};
        
    %first read ns3
    
    if ~isequal(ns3_fpath, '') && exist(ns3_fpath, 'file')
            
        micro_freq = 30000; % sampling rate for micros 
        output_timeseries_freq = 1000;
        
        pulses.timeseries_freq = output_timeseries_freq;
        pulses.uptimes_freq = micro_freq;
            
        rawData_NS3 = concatOpenNSx(ns3_fpath);
        pulses.pulse_nsx_postProc = rawData_NS3.postProc;
               
        ns3_freq = rawData_NS3.MetaTags.SamplingFreq;
        
        %create vector of indices for downsampling analog pulse timeseries 
        downsampled_step_factor = ns3_freq/output_timeseries_freq;
        
        
        
        electrodeLabels = { rawData_NS3.ElectrodesInfo.Label };
        
        %remove null spaces from electrodeLabels
        for iElec = 1:length(electrodeLabels)
           
            currentElec = electrodeLabels{iElec};
            endIdx = 1;
            
            while (endIdx + 1) <= length(currentElec) && double(currentElec(endIdx + 1)) ~= 0
           
                endIdx = endIdx + 1;     
            end
            
            
            electrodeLabels{iElec} = electrodeLabels{iElec}(1:endIdx);
        end
        
        fprintf('electrode labels\n');
        disp(electrodeLabels);
        
        
        
        % look for ain channels
        ains = electrodeLabels( find(cellfun( @(x) contains(x, 'ain'), electrodeLabels)) );
            
        fprintf('ain electrodes\n');
        disp(ains);
        
        %loop looking for ain 
        for iAin = 1:length(ains)
           
            current_ain = ains{iAin};
            channelIdx = find(cellfun( @(x) isequal(x, current_ain), electrodeLabels));
            
            if length(channelIdx) > 1
                error('more than one channel in %s with name %s. Is this possible?\n', ns3_fpath, current_ain);
            end
            
            fprintf('channelIdx %d\n', channelIdx);
            size(rawData_NS3.Data(channelIdx,:)')
            
            % call get_triggers for current channel, set pulse times to physio sample rate and store in pulses struct
            thisAinPulses = get_triggers(double(rawData_NS3.Data(channelIdx,:)'),ns3_freq);
            eval( sprintf('pulses.%s_uptimes = thisAinPulses{1}*micro_freq/ns3_freq;', current_ain ) );
            
            % downsample current channel and store in the pulses struct
            downsampled_channel = rawData_NS3.Data(channelIdx,1:downsampled_step_factor:size(rawData_NS3.Data(channelIdx,:),2));
            eval( sprintf('pulses.%s_ts = downsampled_channel;', current_ain ) );
            
            
        end
                
        
    end
    
    pulses.din1_uptimes = [];
    pulses.din2_uptimes = [];
    pulses.din3_uptimes = [];
    pulses.din4_uptimes = [];

    pulses.din1_ts = [];
    pulses.din2_ts = [];
    pulses.din3_ts = [];
    pulses.din4_ts = [];
    
    pulses.din1_flags = [];
    pulses.din2_flags = [];
    pulses.din3_flags = []; 
    pulses.din4_flags = [];
    pulses.what_is_IPIViolationFlag = 'tldr: when correctSplitNEV is called on din uptimes AND a clock reset is detected an estimated timestamp is calculated marking the end of a segment. If the time between that timestamp and the last preceeding pulse is greater than the channel-specifc IPI, this flag is set to 1. read correctSplitNEV';
    pulses.what_is_moreClockResetsThanNegDiff = 'based on the postProc from the nsx file more clock resets occurred than could be detected in digital channel pulse times using negative diffs. This condition prohibits clock reset correction for a channels pulses. Returned pulses are uncorrected';
    pulses.what_is_moreNegDiffThanClockReset = 'inverse of moreClockResetsThanNegDiff';
   
    
    if ~isequal(nev_fpath, '') && exist(nev_fpath, 'file')
     
        NEVdata = openNEV( nev_fpath , 'nosave' , 'nomat');
        
        [pulses.din1_uptimes, pulses.din1_ts, pulses.din1_flags] = getBlackRockPulsesDC(NEVdata, 9, postProc);
        [pulses.din2_uptimes, pulses.din2_ts, pulses.din2_flags] = getBlackRockPulsesDC(NEVdata, 10, postProc);
        [pulses.din3_uptimes, pulses.din3_ts, pulses.din3_flags] = getBlackRockPulsesDC(NEVdata, 11, postProc);
        [pulses.din4_uptimes, pulses.din4_ts, pulses.din4_flags] = getBlackRockPulsesDC(NEVdata, 12, postProc);

    end 
    
end