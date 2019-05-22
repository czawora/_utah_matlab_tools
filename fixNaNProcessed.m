function fixNaNProcessed(processed_mat_fpath)

    lfpStruct = load(processed_mat_fpath);
    lfpStruct = lfpStruct.lfpStruct;
    
    if isfield(lfpStruct, 'nan_mask')
        return;
    end
    
    fprintf('working on %s\n', processed_mat_fpath);
    
    lfp = lfpStruct.lfp{1};
    
    nan_mask = false(size(lfp));
    
    for iRow = 1:size(lfp, 1)
        
       fprintf('%d ...', iRow);
       
       row_vec = false(1, size(lfp, 2));
       
       for iCol = 1:(size(lfp, 2)-1)
           if lfp(iRow,iCol) == 0 && lfp(iRow,iCol+1) == 0
               row_vec(1, iCol) = true(1);
           end
       end
       
       if lfp(iRow,size(lfp, 2)-1) == 0 && lfp(iRow,size(lfp, 2)) == 0
           row_vec(end) = true(1);
       end
       
       nan_mask(iRow, :) = row_vec;
       
    end
    
    
    createdDate = lfpStruct.createdDate;
    
    new_readme = ['this lfp.mat file, generated ' sprintf('%s', createdDate) ', contains the following fields:' newline newline ...
          '     createdDate          - a string indicating when this mat file was created' newline ...
          '     chanIDperNSP         - a ( #NSPs used x 1 ) cell array. Each enty is a table corresponding to the channel data in lfp' newline ...
          '     lfp                  - a ( #NSPs used x 1 ) cell array. Each cell entry is #channels x #time points matrix (int16). Time point dimension might be slightly different prior to NSP alignment' newline ...
          '     nan_mask             - a ( #NSPs used x 1 ) cell array. Each cell entry is #channels x #time points matrix (logical) indicating where "saturation" NaNs should be applied to corresponding lfp matrix' newline ...
          '     glob_sig_all         - 1-d cell array containing a global mean using all channels from a device. Index this array using MicroDevNum from the jacksheet' newline ...
          '     glob_sig_good        - 1-d cell array containing a global mean using only channels from a device that pass an amplitude and variance quality filter (see variance.csv). Index this array using MicroDevNum from the jacksheet' newline ... 
          '     rerefType            - a string indicating either "lfp_noreref" or "lfp_processed"' newline ...
          '     gain_bin2uV          - multiplicative factor to convert lfp_noreref to uV' newline ...
          '     samplingFreq         - samlple frequency of the data. typically 1000 Hz' newline ...
          '     filterSettings       - a string describing the processing steps applied to data in lfp field' newline ...
          '     sessStr              - the session name the lfps were extracted from (e.g., 190117_1336)' newline ...
          '     sessDurMin           - session duration in minutes' newline ...
          '     pulses               - (a #nsp x 3) cell with "nspSuffix" "nsxFilename" and "pulse" struct containing 30kHz uptimes for all "ain" and "din" channels, as well as 1kHz downsampled "ain" timeseries for that NSP' newline ...
          '     alignedTo            - a string indicating which file the spikes have been aligned to' newline ...
          '     alignmentChan        - a string indicating which channel in the pulses struct was used for alignment' newline ...
          '     jackTableFull        - complete jacksheetBR table from this session (all channels) with device numbers and new channel names' newline ...
          '     jackTableUsed        - just the jacksheet for the channels incorporated into this lfpStruct (combined across NSPs). this table can be used to go back and forth between original and new channel names' newline ...
          '     startTime_datenum    - the start time of session as output from datenum function (note that these values are taken from the timestamp in the nsx filename, not the incorrectly offset time values in the original nsx file)'  newline ...
          '     physio_nsx_postProc  - (a #nsp x 3) cell with "nspSuffix" "nsxFilename" and struct with info on how the original nsx file was modified using concatOpenNSx' ];

    lfpStruct.readme = new_readme;
    lfpStruct.nan_mask = {nan_mask};
    
    save(processed_mat_fpath,  '-v7.3', 'lfpStruct');
      
    fprintf('\n');
end