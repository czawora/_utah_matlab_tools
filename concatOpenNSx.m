function nsx_data = concatOpenNSx(nsx_filepath, applyGain)
%this function will read an NSx file using openNSx
%it will remove the first data segment if it's a junk segment
%and concatenate the remaining data segments and other segmented struct fields

    if nargin < 2
        applyGain = 0;
    end

    if ~exist(nsx_filepath, 'file')
       error('%s is not a valid filepath', nsx_filepath); 
    end
        
    nsx_data = openNSx(nsx_filepath); 

    nsx_data.postProc = struct;
    nsx_data.postProc.junkSamplesRemoved = 0;
    nsx_data.postProc.junkSecRemoved = 0;
    nsx_data.postProc.appliedGain = 0;
    
    disp(nsx_data);

    if iscell(nsx_data.Data) 
           
        fprintf('data contains %d segments\n', length(nsx_data.Data));
        nsx_data.postProc.numRawSegments = length(nsx_data.Data);
        nsx_data.postProc.numNonJunkSegments = length(nsx_data.Data);
        
        %multi-segment data
        if length(nsx_data.Data) > 1
                       
            %if first segment is less than 5 sec, 
            %treat it as a pre-sync junk segment
            if nsx_data.MetaTags.DataDurationSec(1) < 5

                removedTimeSec = nsx_data.MetaTags.DataDurationSec(1);
                
                nsx_data.postProc.junkSamplesRemoved = size(nsx_data.Data{1}, 2);
                nsx_data.postProc.junkSecRemoved = removedTimeSec;
                
                fprintf('first data segment only %d seconds long, removing it\n', removedTimeSec);
                
                nsx_data.MetaTags.Timestamp(1) = [];
                nsx_data.MetaTags.DataPoints(1) = [];
                nsx_data.MetaTags.DataDurationSec(1) = [];
                nsx_data.Data(1) = [];

                
                
                nsx_data.postProc.removedJunk = 1;
                nsx_data.postProc.numNonJunkSegments = nsx_data.postProc.numRawSegments - 1;
            
            end
            
            %if there are still multple segments after the junk segment is
            %possible removed, we have run into the blackrock clock reset issue
            if length(nsx_data.Data) > 1    

                fprintf('concatenating remaining %d data segments', length(nsx_data.Data));
                
                nsx_data.Data = horzcat(nsx_data.Data{:});

                nsx_data.MetaTags.Timestamp = -1;
                nsx_data.MetaTags.DataPoints = sum(nsx_data.MetaTags.DataPoints);
                nsx_data.MetaTags.DataDurationSec = sum(nsx_data.MetaTags.DataDurationSec);

            else

                nsx_data.Data = nsx_data.Data{1};

            end

        
        else %data array is a cell array, but is not multi-segment, pull the data out of the cell array
                
            nsx_data.Data = nsx_data.Data{1};
        end
        
    end
    
    if applyGain == 1
       
        fprintf('apply 1/4 gain to data\n');
        nsx_data.postProc.appliedGain = 1;
        
        nsx_data.Data = nsx_data.Data/4;
    end
                            
end