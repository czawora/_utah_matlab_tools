function nsx_data = concatOpenNSx(nsx_filepath, applyGain, skipfactor, chanToOpen, nonClockResetTimeStampThresh)
%this function will read an NSx file using openNSx
%it will remove the first data segment if it's a junk segment
%and concatenate the remaining data segments and other segmented struct fields
%
%  skipfactor - (1) dont skip any sample, >1 skip samples
%  chanToOpen - [] for all channels, else specify channels


VERBOSE = 0; %- 0:  set to 1 to get more outputs about cells existing, etc

if ~exist(nsx_filepath, 'file')
    error('%s is not a valid filepath', nsx_filepath);
end



%- options for how openNSx is called
if nargin < 2,
    applyGain = 0;
end
if nargin < 3,
    skipfactor = 1;
end
if nargin < 4,
    chanToOpen = [];
end
if nargin<5;
    nonClockResetTimeStampThresh = 5000;
end

if isempty(chanToOpen)
    nsx_data = openNSx(nsx_filepath,'skipfactor',skipfactor);
else
    nsx_data = openNSx(nsx_filepath,'skipfactor',skipfactor,'channels',chanToOpen); %- JW tweaked so chanToOpen can actually be passed 12/2018
end

nsx_data.postProc = struct;
nsx_data.postProc.removedJunk = 0;
nsx_data.postProc.appliedGain = 0;

% disp(nsx_data);

if iscell(nsx_data.Data)
    
    if VERBOSE, fprintf('\n in concatOpenNSx, data contains %d segments ', length(nsx_data.Data)); end 
    nsx_data.postProc.numRawSegments     = length(nsx_data.Data);
    nsx_data.postProc.numNonJunkSegments = length(nsx_data.Data);
    
    %multi-segment data
    if length(nsx_data.Data) > 1
        
        %if first segment is less than 5 sec, treat it as a pre-sync junk segment
        if nsx_data.MetaTags.DataDurationSec(1) < 5
            
            removedTimeSec = nsx_data.MetaTags.DataDurationSec(1);
            
            if VERBOSE, fprintf('first data segment only %.1f seconds long, removing it\n', removedTimeSec); end
            
            nsx_data.MetaTags.Timestamp(1)  = [];
            nsx_data.MetaTags.DataPoints(1) = [];
            nsx_data.MetaTags.DataDurationSec(1) = [];
            nsx_data.Data(1) = [];
            
            nsx_data.postProc.removedJunk = 1;
            nsx_data.postProc.numNonJunkSegments = nsx_data.postProc.numRawSegments - 1;
            
        end
        
        %if there are still multiple segments after the junk segment is
        %possible removed, we have run into the blackrock clock reset issue
        if length(nsx_data.Data) > 1
            
            if VERBOSE, fprintf('concatenating remaining %d data segments\n\n', length(nsx_data.Data)); end
            
            
            original_timestamps = nsx_data.MetaTags.Timestamp;
            original_cell_lengths = cellfun(@(x) size(x,2), nsx_data.Data);
            
            nsx_data.postProc.nonJunkTimeStamps = original_timestamps;
            nsx_data.postProc.nonJunkCellLengths = original_cell_lengths;
            nsx_data.postProc.samplesAdded = [];
            
            postProcIter = 1;
            
            while length(nsx_data.Data)~=1
                
                
                timeStampDifference = original_timestamps(2) - original_timestamps(1);
                
                samplingFreq = nsx_data.MetaTags.SamplingFreq;
                samplingRateOfTimestamps = 30000;
                segmentLength = original_cell_lengths(1);
                samplesInGap = round(timeStampDifference * (samplingFreq/samplingRateOfTimestamps));
                
                numSamplesToInsert = floor((samplesInGap - segmentLength) / skipfactor);  
                    
                if (timeStampDifference > nonClockResetTimeStampThresh) && (numSamplesToInsert>=1) % if timestamps were close to one another,
                    % there was probably a clock reset
                    
                    if numSamplesToInsert==1
                        samplesToInsert = nsx_data.Data{1,1}(end);
                        
                    elseif mod(numSamplesToInsert,2)==1
                        takeFromBefore = nsx_data.Data{1,1}(:,end-round(numSamplesToInsert/2)+1:end);
                        takeFromAfter = nsx_data.Data{1,2}(:,1:round(numSamplesToInsert/2)-1);
                        samplesToInsert = cat(2,fliplr(takeFromBefore),fliplr(takeFromAfter));
                        
                    elseif mod(numSamplesToInsert,2)==0
                        takeFromBefore = nsx_data.Data{1,1}(:,end-round(numSamplesToInsert/2)+1:end);
                        takeFromAfter = nsx_data.Data{1,2}(:,1:round(numSamplesToInsert/2));
                        samplesToInsert = cat(2,fliplr(takeFromBefore),fliplr(takeFromAfter));
                        
                    end
                    
                    nsx_data.Data{1} = horzcat(nsx_data.Data{1},samplesToInsert,nsx_data.Data{2});
                    
                else
                    
                    nsx_data.Data{1} = horzcat(nsx_data.Data{1},nsx_data.Data{2});
                    
                end
                
                nsx_data.postProc.samplesAdded(postProcIter) = numSamplesToInsert;
                postProcIter = postProcIter + 1;
                
                nsx_data.Data(2) = [];
                original_timestamps(1) = [];
                original_cell_lengths(1) = [];
            end

            
        end   % real data segment end

    end   % junk segment end
    
    nsx_data.Data = nsx_data.Data{1};
    
    nsx_data.MetaTags.DataPoints = size(nsx_data.Data,2);
    nsx_data.MetaTags.DataDurationSec = size(nsx_data.Data,2) / nsx_data.MetaTags.SamplingFreq;
    
end


if applyGain == 1
    
    fprintf('apply 1/4 gain to data\n');
    nsx_data.postProc.appliedGain = 1;
    
    nsx_data.Data = nsx_data.Data/4;
end

end