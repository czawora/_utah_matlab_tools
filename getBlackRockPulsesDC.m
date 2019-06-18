function [output_br_timeStamp_up, output_ts, IPIViolationFlag] = getBlackRockPulsesDC(NEVdata,whichDC,postProc)
%Given the NEV data, outputs the pulse information from a certain DC
%channel
%
% whichDC = 9, 11, 12, etc.

% keyboard

output_br_timeStamp_up = [];
output_ts = [];
IPIViolationFlag = 0;

if nargin < 3,
    postProc = [];
end

% 2. Get raw pulse information
rawUnparsedData = NEVdata.Data.SerialDigitalIO.UnparsedData;
raw_br_timeStamp = double(NEVdata.Data.SerialDigitalIO.TimeStamp);

binaryChar = dec2bin(rawUnparsedData);
binaryMat = double(binaryChar);
binaryMat(binaryMat==49) = 1; binaryMat(binaryMat==48) = 0;

diffMat = nan(size(binaryMat));
for instance = 2:size(binaryMat,1)
    diffMat(instance,:) = binaryMat(instance,:) - binaryMat(instance-1,:);
end

timeStamps = cell(1,size(binaryMat,2));
for loc = 1:size(binaryMat,2)
    timeStamps{loc} = cat(2,diffMat(diffMat(:,loc)~=0,loc),raw_br_timeStamp(diffMat(:,loc)~=0)');
end

br_timeStamp_up = [];

if ~isempty(timeStamps)
    
    % DC09 is D3, which is first from the back (D0,D1,D2,D3)
    % DC11 is D1, which is third from the back (D0,D1,D2,D3)
    % DC12 is D0, which is fourth from the back (D0,D1,D2,D3)
    
    DCposition = length(timeStamps) - (whichDC - 9);
    
    br_timeStamp_up = timeStamps{DCposition}(timeStamps{DCposition}(:,1)==1,2);
        
end


if ~isempty(br_timeStamp_up) && ~any(isnan(br_timeStamp_up))
    
    if ~isempty(postProc) && isfield(postProc,'samplesAdded')
        
        [br_timeStamp_up, IPIViolationFlag] = correctSplitNEV(br_timeStamp_up, whichDC, postProc);
    end
    
end

output_br_timeStamp_up = br_timeStamp_up;


end

