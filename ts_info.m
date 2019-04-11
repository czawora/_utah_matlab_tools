function [m1, m2, m3, m4, amp, time, varargout] = ts_info(X,Fs,T,overlap,varargin)
%
%TS_INFO
%Written by: Julio I. Chapeton
%
%Calculates the first four moments, range, and entropy of a vector or matrix for the specified
%windowing. m3 and m4 are satandardized, i.e skewness and kurtosis.
%Segments can be overlapping. Matrix inputs should be of the form (time x channel)
%
%Usage:
%TS_INFO(X,Fs,T,overlap,varargin)
%
%Inputs:
%    Required:
%        X        = vector/matrix (contains time series, if matrix should be (samples x channels))
%        Fs       = scalar (sampling frequency in Hz)
%        T        = scalar (size of the time window in seconds, can use Inf as an input to use the whole time series with no overlap)
%        overlap  = scalar (size of overlap for adjacent windows, as a fraction of T)
%                   NOTE: To define T in samples as opposed to seconds set Fs=1;
%    Optional:
%        scale    = string (use 'ztrans_seg' to z-score the data at each segment, 'ztrans_full' to zscore the whole clip,...
%                   or 'bounded' to scale each segment to go from 0 to 1)
%        bins     = vector/scalar (vector of bin edges or scalar with number of bins to use. Given a scalar input data will be binned using that number of bins with limits...
%                   dependent on the scaling: 0 to 1 for bounded, -5 to 5 for both 'ztrans' options, and min(X(:)) to max(X(:)) if there is no scaling. The default is 30 bins...
%                   using the scaling dependent limits. If input as a vector it will define the bin edges explicitly)
%
%Outputs:
%        m1...m4  = vector/matrix (ith moment for each of the time segments)
%        amp      = vector/matrix (rms of the data)
%        times    = vector/matrix (provides start and end times for each time segment)
%
%    Optional:
%        H        = vector/matrix (entropy of the data)
%
%NOTES:
% 1. Add visualization options
% 2. May need to use nanmin and nanmax
%% parse inputs
inp_pars = inputParser;

def_scale = [];
expected_scale = {'ztrans_seg','ztrans_full','bounded'};%,'detrend'};
def_bins = [];

addParameter(inp_pars,'scale',def_scale,@(x) assert(any(strcmpi(x,expected_scale)) | isempty(x),...
    'Not a valid scaling option, the options are ''ztrans_seg'',''ztrans_full'' or ''bounded'''));
addParameter(inp_pars,'bins',def_bins,@(x) (isvector(x) || isscalar(x) || isempty(x)));

parse(inp_pars,varargin{:})
%inp_pars.Results

%% define limits and binning
if isempty(inp_pars.Results.bins)
    nbins=30;
    if strcmpi(inp_pars.Results.scale,'bounded')
        lb=0;ub=1;
    elseif strncmpi(inp_pars.Results.scale,'ztrans',6)
        lb=-5;ub=5;
    else
        lb=min(X(:)); ub=max(X(:));
    end
    
    bins=[lb+(ub-lb)/nbins/2:(ub-lb)/nbins:ub-(ub-lb)/nbins/2];
    bins=double(bins);
    d=diff(bins)/2;
    edges = [bins(1)-d(1), bins(1:end-1)+d, bins(end)+d(end)];
    edges(2:end) = edges(2:end)+eps(edges(2:end));
    edges(1)=-Inf;edges(end)=Inf;
elseif isscalar(inp_pars.Results.bins)
    nbins=inp_pars.Results.bins;
    
    if strcmpi(inp_pars.Results.scale,'bounded')
        lb=0;ub=1;
    elseif strncmpi(inp_pars.Results.scale,'ztrans',6)
        lb=-5;ub=5;
    else
        lb=min(X(:)); ub=max(X(:));
    end
    
    bins=[lb+(ub-lb)/nbins/2:(ub-lb)/nbins:ub-(ub-lb)/nbins/2];
    bins=double(bins);
    d=diff(bins)/2;
    edges = [bins(1)-d(1), bins(1:end-1)+d, bins(end)+d(end)];
    edges(2:end) = edges(2:end)+eps(edges(2:end));
    edges(1)=-Inf;edges(end)=Inf;
    
elseif isvector(inp_pars.Results.bins)
    edges=inp_pars.Results.bins;
    lb=min(edges(:)); ub=max(edges(:));
end

if isinf(T)
    Nsamples=size(X,1);
    ovrlp_samples=0;
else
    Nsamples=round(T*Fs);
    ovrlp_samples=round(overlap*Nsamples);
end

if isvector(X)
    X=X(:);
end

%% main calculation
nwins=floor((size(X,1)-Nsamples)/(Nsamples-ovrlp_samples))+1;
time=nan(nwins,2);
time(:) = [((1:nwins)-1)*(Nsamples-ovrlp_samples),((1:nwins)*(Nsamples-ovrlp_samples)+ovrlp_samples)-1]/Fs;

X_windowed=nan(Nsamples,nwins);
m1=nan(nwins,size(X,2));
m2=nan(nwins,size(X,2));
m3=nan(nwins,size(X,2));
m4=nan(nwins,size(X,2));
amp=nan(nwins,size(X,2));
H=nan(nwins,size(X,2));
for k=1:size(X,2)
    if strcmpi(inp_pars.Results.scale,'ztrans_seg')
        [X_windowed, Z] = buffer(X(:,k),Nsamples,ovrlp_samples,'nodelay');
        X_windowed=bsxfun(@rdivide,bsxfun(@minus,X_windowed,nanmean(X_windowed)),nanstd(X_windowed)); % may be better to check for nans and only use nanmean/nanstd when needed
    elseif  strcmpi(inp_pars.Results.scale,'ztrans_full')
        X(:,k)=(X(:,k)-nanmean(X(:,k)))./nanstd(X(:,k));
        [X_windowed, Z] = buffer(X(:,k),Nsamples,ovrlp_samples,'nodelay');
    elseif strcmpi(inp_pars.Results.scale,'bounded')
        [X_windowed, Z] = buffer(X(:,k),Nsamples,ovrlp_samples,'nodelay');
        X_windowed=bsxfun(@minus,X_windowed,min(X_windowed,[],1));
        X_windowed=bsxfun(@rdivide,X_windowed,max(X_windowed,[],1));
        %     elseif strcmpi(inp_pars.Results.scale,'detrend')
        %         [X_windowed, Z] = buffer(X(:,k),Nsamples,ovrlp_samples,'nodelay');
        %         X_windowed=detrend(X_windowed);
    elseif isempty(inp_pars.Results.scale)
        [X_windowed, Z] = buffer(X(:,k),Nsamples,ovrlp_samples,'nodelay');
    end
    
    assert(~all(isnan(X_windowed(:))), 'Bad error: all NaN');
    
    m1(:,k)=nanmean(X_windowed);
    m2(:,k)=nanvar(X_windowed);
    %var_check=(mean(m2+m1.^2)-mean(m1).^2-var(X))/var(X);
    m3(:,k)=skewness(X_windowed);
    m4(:,k)=kurtosis(X_windowed);
    %     amp(:,k)=range(X_windowed);
%     amp(:,k)=rms(X_windowed);
    amp(:,k)=sqrt(nanmean(X_windowed.^2));
    if nargout==7
        %entropy calculation
        for l=1:size(X_windowed,2)
            %p=histc(X_windowed(:,l),edges); % histc produces n bins for n edges, it adds a bin at the end which is only non-zero if x=edges(end) exactly, which almost never happens.
            % use histc only if you don't have ndHistc and you have a matlab version prior to R2014b
            
            p=ndHistc(X_windowed(:,l),edges); % ndHistc produces n bins for n+1 edges, behaviour is similar to histcounts, which has replaced histc in newer versions of matlab. Requires a mex function but is very fast.
            %p=histcounts(X_windowed(:,l),edges)'; % Use histcounts if you don't have ndHistc.
            p=p./sum(p);
            H(l,k)=-sum(p.*log2(p+eps))+log2((ub-lb)/(nbins)); %% the second term is a correction for discretization which won't matter if the binning is fixed, the Miller correction can be added if nbins is not << Nsamples ((nbins-1)/2N)
        end
    end
end
varargout{1}=H;
