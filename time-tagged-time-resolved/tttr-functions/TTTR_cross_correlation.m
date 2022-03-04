function T2data = TTTR_cross_correlation(parsedData, options)
% Inputs:
%   parsedData, structure returned from TTTR_extract_channel_times.m with
%   the following fields:
%       chan0Times:             column of arrival times for channel 0
%       chan1Times:             column of arrival times for channel 1
%       nChan0:                 total number of events for channel 0
%       nChan1:                 total number of events for channel 1
%       globalRes:              global resolution, should be 4e-12 s
%       acquisitionTime:        total acquisition time in s
%    options, structure with following fields and default values:
%       tRes = []:              resolution for count rates calculation. omit
%                               if you don't wish to have count rates calculated
%       tauRes = 4e-12:         resolution over which cross-correlation
%                               is calculated
%       tauLimits = 1e-9*[-1 1]:limits for the tau axis
%       tauAxis = []:           tau axis over which cross-correlation is
%                               calculated. if omitted, tau axis will be 
%                               calculated from tauRes and
%                               tauLimits
%       countRateRanges = []:   time axis range over which count
%                               rates is calculated. omit if there's no 
%                               partitioning of the time axis.
%       tFlag = []:             flags for paritioning the time axis for
%                               count rates calculation
%       verbose = false:        whether to enable command window
%                               updates
%       statusbar = false:      whether to enable pop-up progress bar
%                               for cross-correlation calculation
% Output:
%   T2data, a structure with the following fields:
%       nBins:                  number of bins in tau axis
%       nPartitions:            number of time partitions
%       tauAxis:                axis of delay times
%       xCorr:                  raw correlated data
%       g2:                     normalized g2
%       g2norm:                 normalization factor 
%       g2Lerr:                 lower Poisson error on g2
%       g2Uerr:                 upper Poisson error on g2
%       nSamples:               number of samples per each bin
%   if count rates are calculated:
%       tAxis:                  time axis for count rates
%       countRates:             average count rates over acquisition time
%       if there are partitions for calculating count rate:
%       tFlag:                  flags for partitioning the time axis


% Parse the inputs
if nargin == 1
    options = struct;
end
options = parse_options(options);
TreatIntegerTau = false; % this controls an option to treat bin widths that 
% are power-of-two multiples of the global resolution differently from the
% more general case.  In principle it should be more efficient but in
% practice this is much slower, possibly due to a bug in the code or my
% misunderstanding of how MATLAB handles integer arithmetic.  So disable
% for now.
globalRes = parsedData.globalRes; % should be 4e-12 s
acquisitionTime = parsedData.acquisitionTime; % given in ms, convert to s

% Unpack validated options
if isempty(options.tauAxis) % no axis given -use fixed resolution
    uniformTau = true; % uniform bin centers at (...,-2*tau_res,-tau_res,0,tau_res,2*tau_res,...)
    tauRes = options.tauRes; 
    tauLimits = options.tauLimits; % in this case tau_limits must also be supplied
    % check if tau_res/2 is a multiple of GlobalResolution*2^n, in which case can use
    % integer rather than floating point representations and speed up histogramming
    integerTau =  TreatIntegerTau & fix(log2(tauRes/(2*globalRes)))==log2(tauRes/(2*globalRes));
    if integerTau % make sure limits are also integer multiples of tau_res
        tauLimits = tauRes*fix(tauLimits/tauRes);
%         tau_bitshift = log2(tau_res/(2*GlobalResolution)); % can reduce resolution by 2^tau_bitshift to speed comparisons
    end
    tauAxis = (tauLimits(1):tauRes:tauLimits(2));
else % tauAxis given.  tauLimits ignored if provided
    tauAxis = options.tauAxis; 
    tauLimits = tauAxis([1,end]); 
    tauRes = tauAxis(2)-tauAxis(1);
    tolTau = 1e-12;
    uniformTau = all(abs(diff(tauAxis)./tauRes-1)<tolTau); % check if uniform vector
    if ~uniformTau % nonuniform vector of delay values
        tauWidths = zeros(size(tauAxis));
        dtau = diff(tauAxis);
        tauWidths(1) = dtau(1);
        tauWidths(2:end-1) = min(dtau(1:end-1),dtau(2:end));
        tauWidths(end) = dtau(end);
        integerTau = TreatIntegerTau & (all(fix(log2(tauAxis/globalRes))==log2(tauAxis/globalRes)) &&...
            all(fix(log2(tauWidths/(2*globalRes)))==log2(tauWidths/(2*globalRes)))); % all bin boundaries must fall on integer clock ticks
    else
        integerTau =  TreatIntegerTau & fix(log2(tauRes/(2*globalRes)))==log2(tauRes/(2*globalRes)); % check if can use integer values
        if integerTau
            if fix(tauAxis(1)/tauRes)~= tauAxis(1)/tauRes
                tauAxis = tauRes*floor(tauAxis/tauRes); % shift delay positions slightly to align with hardware resolution
            end
%             tau_bitshift = log2(tau_res/(2*GlobalResolution));
        end
    end 
end

if ~isempty(options.tRes) % Perform count-rate calculations
    if isscalar(options.tRes) % fixed resolution
        tRes = options.tRes;
        tAxis = (0.5*tRes):tRes:acquisitionTime; % centered bins
%         t_axis = 0:t_res:AcquisitionTime; % beginning of each bin
%         t_uniform = true;
    else % assume this is tAxis
        tAxis = options.tRes;
        tRes = tAxis(2)-tAxis(1);
        tolT = 1e-2*tRes;
%         t_uniform = all(abs(diff(t_axis)./t_res-1)<tol_t); % check if uniform vector
        if ~issorted(tAxis) || any(abs(diff(tAxis)./tRes-1)>tolT) % check if uniform and sorted
            error('List of time values provided in tRes must be uniform and sorted!');
        end
        if any(~inrange(tAxis,[0,acquisitionTime]))
            warning('Provided values of tAxis extend outside the acquisition time');
        end
    end
    calcCountRates = true;
else
    calcCountRates = false;
end

if calcCountRates % otherwise ignore these fields
    if ~isempty(options.tFlag) % this takes precedence if provided
        tFlag = options.tFlag;
        if length(tFlag)~=length(tAxis)
            error('Optional input "tFlag" must be the same length as tAxis!');
        end
        tPartition = 'flags';
        nPartitions = max(tFlag);
    elseif ~isempty(options.countRateRanges)
        countRateRanges = options.countRateRanges;
        tPartition = 'CountRate';
        nPartitions = size(countRateRanges,1);
    else
        tPartition = 'none';
        nPartitions = 1;
        tFlag = ones(size(tAxis)); % set all flags to same value
    end
else
    tPartition = 'none';
    nPartitions = 1;
%     tFlag = ones(size(tAxis)); % set all flags to same value
end

% check that tauLimits<tRes
if ~strcmp(tPartition,'none')
    if max(abs(tauLimits))>tRes
        warning('Maximum delay time exceeds resolution of time-partition vector.  Partitioning might not be valid.');
    end
end

verbose = options.verbose;
plotWaitbar = options.statusbar;

if verbose
    fprintf(1,'Unpacking T2 data...');
end
ChAtimes = parsedData.chan0Times;
ChBtimes = parsedData.chan1Times;
nChA = parsedData.nChan0;
nChB = parsedData.nChan1;

if integerTau
    tMax = acquisitionTime/globalRes;
else
    tMax = acquisitionTime;
end

% Calculate count rates in each channel as a function of time

if calcCountRates
    if verbose
        fprintf(1,'Calculating count rates...');
    end
    countRates = zeros(2,length(tAxis));
    counts = zeros(2,length(tAxis));
    validTBins = find(inrange(tAxis,[0,acquisitionTime]));
    nTBins = length(validTBins);
    if integerTau
        tValues = tAxis/globalRes; % convert to clock tick units
        tRes = tRes/globalRes;
    else
        tValues = tAxis; % leave in seconds
    end
    tIdx = validTBins(1); % first point in range
    tCenter = tValues(tIdx); 
    tRange = [max([tCenter-tRes/2,0]),tCenter+tRes/2]; % covers case where bin extends earlier than t=0
    tWidth = diff(tRange); 
    [lIdxA,uIdxA] = BinarySearch(ChAtimes,tRange);    
    [lIdxB,uIdxB] = BinarySearch(ChBtimes,tRange);
    counts(:,tIdx) = [(1+uIdxA-lIdxA); (1+uIdxB-lIdxB)];
    countRates(:,tIdx) = counts(:,tIdx)/tWidth; % will need to convert to real units later
%     CountRates(1,t_ix) = (1+uixA-lixA)/t_width; 
%     CountRates(2,t_ix) = (1+uixB-lixB)/t_width;
    lastIdxA = uIdxA;
    lastIdxB = uIdxB;
   
    for tBinIdx = 2:nTBins % remaining bins
        tIdx = validTBins(tBinIdx);
        tUpper = tValues(tIdx)+tRes/2; % upper bound of current bin
        uIdxA = BinaryFindUpperBound(ChAtimes,tUpper);
        uIdxB = BinaryFindUpperBound(ChBtimes,tUpper);
        if tBinIdx==nTBins % deal with width of last bin separately
            tWidth = tRes - max([tUpper-tMax,0]);
        else
            tWidth = tRes;
        end
        counts(:,tIdx) = [(uIdxA-lastIdxA); (uIdxB-lastIdxB)];
        countRates(:,tIdx) = counts(:,tIdx)/tWidth;
%         CountRates(1,t_ix) = (uixA-last_ixA)/t_width;
%         CountRates(2,t_ix) = (uixB-last_ixB)/t_width;
        lastIdxA = uIdxA;
        lastIdxB = uIdxB;
    end
    if integerTau
        countRates = countRates./globalRes; % convert back to Hz units
    end   

    if strcmp(tPartition,'CountRate') % need to determine flags, otherwise should already be set.
        tFlag=NaN(size(tAxis));
        for tAxIdx=1:length(tAxis)
            totalR = sum(countRates(:,tAxIdx));
            in_range = totalR>countRateRanges(:,1) & totalR<=countRateRanges(:,2);
            gaveWarning1 = false;
            gaveWarning2 = false;
            if any(in_range)
                if ~gaveWarning1 && sum(in_range)>1 % only give warning once
                    warning('countRateRanges are not mutually exclusive. First matching range is used.');
                    gaveWarning1 = true;
                end
                tFlag(tAxIdx) = find(in_range,1);
            else
                tFlag(tAxIdx) = -1;
                if ~gaveWarning2
                    warning('Some measurements fall outside CountRate_ranges and will be ignored.');
                    gaveWarning2 = true;
                end
            end
        end
    end
    avgRateB = NaN(1,nPartitions); % Normalization factor needed to calculate g2
    avgRateA = NaN(1,nPartitions);
    for eventIdx=1:nPartitions
        currFlag = (tFlag==eventIdx)&validTBins; % logical vector
        totalCtsB = sum(counts(2,currFlag));
        totalCtsA = sum(counts(1,currFlag));
        totalTime = sum(currFlag)*tRes;
        avgRateB(eventIdx) = totalCtsB/totalTime; % units are cts/time where time is either s or clock ticks
        avgRateA(eventIdx) = totalCtsA/totalTime;
    end
    if verbose
        fprintf(1,'Done!\n');
    end
else
    avgRateA = nChA/tMax;
    avgRateB = nChB/tMax;
end


% Calculate cross correlation

nBins = length(tauAxis);
xCorr = zeros(nBins,nPartitions);

% Determine limits of each time-lag bin
if uniformTau
    tauMin = tauAxis-tauRes/2;
    tauMax = tauAxis+tauRes/2;
else
    tauMin = tauAxis-tauWidths/2;
    tauMax = tauAxis+tauWidths/2;
end

if integerTau % convert to integer type (same as ChAtimes, ChBtimes)
    tauMin = uint64(tauMin/globalRes); 
    tauMax = uint64(tauMax/globalRes);
    tMax = uint64(tMax);
end

nSamples = zeros(nBins,nPartitions); % integer array that will give the number of samples for each bin, accounting for edge effects and partitioning

% Generate vectors l_ix,m_ix corresponding to time-lag limits of each bin
lVec = zeros(nBins,1);
mVec = zeros(nBins,1);

% Use channel with fewest counts as Ch1 for efficiency
if nChA>nChB
    Ch1times = ChBtimes; % I believe Matlab assignes this via pointer and does not copy the arrays...
    Ch2times = ChAtimes;
    nChan0 = nChB;
    nChan1 = nChA;
    avgRate = avgRateA; % average count rate of higher-rate channel, needed for g2 normalization
else
    Ch1times = ChAtimes;
    Ch2times = ChBtimes;
    nChan0 = nChA;
    nChan1 = nChB;
    avgRate = avgRateB;
end

if verbose
    fprintf(1,'Calculating cross-correlation...');
end

if plotWaitbar
    wbh = waitbar(0,'Calculating cross-correlation...'); % update message
%     wbinc = floor(N1/20);
    wThresh = 1e5;
    wCounter = 0;
end

for eventIdx=1:nChan0 % photon number in channel with lower photon counts
    t0 = Ch1times(eventIdx);
    for binIdx=1:nBins
        if binIdx>1 % check if new bin is adjacent to previous
            if integerTau
                adjacentBins = uniformTau || tauMin(binIdx)==tauMax(binIdx-1);
            else
                adjacentBins = uniformTau || abs(tauMin(binIdx)-tauMax(binIdx-1))<globalRes;
            end
        else
            adjacentBins = false;
        end
        if adjacentBins % if adjacent bin limits are equivalent
            lb = mVec(binIdx-1); % these will be the same
        else
            if eventIdx==1 % need to find limit for first time
                lb = BinaryFindLowerBound(Ch2times,t0+tauMin(binIdx));
            else % start with value from last iteration
                lb = LinearFindLowerBound(Ch2times,t0+tauMin(binIdx),lVec(binIdx)); 
            end
        end
        if isempty(lb) || lb>nChan1 % means no lower bound exists
            lVec(binIdx) = nChan1+1;
            mVec(binIdx) = nChan1+1; % this will necessarily also be out of range
        else
            lVec(binIdx) = lb;
            if eventIdx==1 % need to find limit for first time
                lb2 = BinaryFindLowerBound(Ch2times,t0+tauMax(binIdx));
            else % start with value from last iteration
                lb2 = LinearFindLowerBound(Ch2times,t0+tauMax(binIdx),mVec(binIdx)); 
            end
            if isempty(lb2) % out of range
                mVec(binIdx) = nChan1+1;
            else
                mVec(binIdx) = lb2;
            end
        end
    end
    
    validBins = (t0+tauMin>0) & (t0+tauMax<tMax);
    xCorrIncre = (mVec-lVec).*validBins(:); % new xcorr points to add, ignoring bins that are out of range
    
    switch tPartition
        case 'none'
            nSamples = nSamples+validBins(:);
            xCorr = xCorr+xCorrIncre;
        otherwise % using flags
            tIdx = LinearFindLowerBound(tValues+tRes/2,t0,1); % find index
            flag = tFlag(tIdx);
            if flag>0
                nSamples(:,flag) = nSamples(:,flag)+validBins(:);
                xCorr(:,flag) = xCorr(:,flag)+xCorrIncre;
            end
    end
    
    if plotWaitbar
        if wCounter>wThresh
            waitbar(eventIdx/nChan0,wbh);
            wCounter = 0;
        else
            wCounter = wCounter+1;
        end
    end
end

if plotWaitbar
    close(wbh);
end

if verbose
    fprintf(1,'Done!\n');
end

% Now calculate uncertainties and normalized g2

% "Usual" Poisson error calculation, but this is incorrect for bins with
% zero counts.
% xCorrErr = sqrt(xCorr); % Poissonian uncertainties on counts in each bin

% Calculate asymmetric errors for Poisson-distributed data following:
%   lower-error = -0.5+sqrt(0.25+N)
%   upper-error = +0.5+sqrt(0.25+N)
PoissErrs = sqrt(0.25+xCorr); 
xCorrErr = cat(3,-0.5+PoissErrs,0.5+PoissErrs);

if integerTau
    binwidths = globalRes.*repmat(tauMax(:)-tauMin(:),1,nPartitions);
else
    binwidths = repmat(tauMax(:)-tauMin(:),1,nPartitions);
end

g2norm = nSamples.*repmat(avgRate,nBins,1).*binwidths; % Normalization factor to calculate g2
g2 = xCorr./g2norm;
g2Lerr = xCorrErr(:,:,1)./g2norm;
g2Uerr = xCorrErr(:,:,2)./g2norm;

% Pack and return outputs

T2data.nBins = nBins;
T2data.nPartitions = nPartitions;
T2data.tauAxis = tauAxis(:);
T2data.xCorr = xCorr;
T2data.g2 = g2;
T2data.g2norm = g2norm;
T2data.g2Lerr = g2Lerr;
T2data.g2Uerr = g2Uerr;
T2data.Nsamples = nSamples;
if calcCountRates
    T2data.tAxis = tAxis;
    T2data.countRates = countRates;
    if ~strcmp(tPartition,'none')
        T2data.tFlag = tFlag;
    end
end
end

%% 
function validatedOptions = parse_options(inputOptions)
parseInputs = inputParser; %RF changed parseInputs to parse_inputs

% Default parameters
defaultOptions.tRes = []; % this will disable calculation of count rates if t_res is not provided
defaultOptions.tauRes = 4e-12;
defaultOptions.tauLimits = 1e-9*[-10 10];
defaultOptions.tauAxis = [];
defaultOptions.countRateRanges = [];
defaultOptions.tFlag = [];
defaultOptions.verbose = false;
defaultOptions.statusbar = false;

if (nargin == 0)
        inputOptions = defaultOptions;
end

addParameter(parseInputs,'tRes',defaultOptions.tRes,...
    @(x) validateattributes(x,{'numeric'},{'nonempty','increasing'}));

addParameter(parseInputs,'tauRes',defaultOptions.tauRes,...
    @(x) validateattributes(x,{'numeric'},{'scalar','nonempty'}));

addParameter(parseInputs,'tauLimits',defaultOptions.tauLimits,...
    @(x) validateattributes(x,{'numeric'},{'vector','increasing','numel',2}));

addParameter(parseInputs,'tauAxis',defaultOptions.tauAxis,...
    @(x) validateattributes(x,{'numeric'},{'vector','nonempty','increasing'}));

addParameter(parseInputs,'countRateRanges',defaultOptions.countRateRanges,...
    @(x) validateattributes(x,{'numeric'},{'nonnegative','ncols',2}));

addParameter(parseInputs,'tFlag',defaultOptions.tFlag,...
    @(x) validateattributes(x,{'numeric'},{'vector','nonnegative'}));

addParameter(parseInputs,'verbose',defaultOptions.verbose,...
    @(x) validateattributes(x,{'logical'},{'scalar'}));

addParameter(parseInputs,'statusbar',defaultOptions.statusbar,...
    @(x) validateattributes(x,{'logical'},{'scalar'}));
parse(parseInputs,inputOptions);
validatedOptions = parseInputs.Results;
end
function [b,c]=BinarySearch(x,range)
%Binary search replacement for ismember(A,B) for the
%special case where the first input argument is sorted.  Returns indices
% [b,c] such that x(b:c) gives the elements contained within
% [range(1),range(end)].
%
% Based on findInSorted() function by Daniel Roeske <danielroeske.de>
% 
% Includes additional check that range is contained within input vector

A=range(1);
B=range(end);

if A>x(end) || B<x(1) || B<A
    % no indices satify bounding conditions
    b = [];
    c = [];
    return;
end

a=1;
b=numel(x);
c=1;
d=numel(x);
if A<=x(1)
   b=a;
end
if B>=x(end)
    c=d;
end
while (a+1<b)
    lw=(floor((a+b)/2));
    if (x(lw)<A)
        a=lw;
    else
        b=lw;
    end
end
while (c+1<d)
    lw=(floor((c+d)/2));
    if (x(lw)<=B)
        c=lw;
    else
        d=lw;
    end
end
end

%% 

function imax = BinaryFindUpperBound(x,s)
% Use binary search to locate the largest index imax such that x(imax)<=s.
% Similar to BinarySearch(x,range) above, but only looks for one bound to
% save time.

N = numel(x);
b = N;
if s>x(end)
    imax = N;
    return;
elseif s<x(1)
    imax = [];
    return;
else
    imax = 1; % initialize
end

while (imax+1<b)
    lw=(floor((imax+b)/2));
    if (x(lw)<=s)
        imax=lw;
    else
        b=lw;
    end
end

end

%%

function imin = BinaryFindLowerBound(x,s)
% Use binary search to locate the smallest index imin such that x(imin)>=s.
% Similar to BinarySearch(x,range) above, but only looks for one bound to
% save time.

N = numel(x);
a = 1;
if s>x(end)
    imin = [];
    return;
elseif s<x(1)
    imin = 1;
    return;
else
    imin = N; % initialize
end

while (a+1<imin)
    lw=(floor((a+imin)/2));
    if (x(lw)<s)
        a=lw;
    else
        imin=lw;
    end
end

end

%% 

function imin = LinearFindLowerBound(x,s,istart)
% Linear search to locate the smallest index imin such that x(imin)>=s.
% Start with imin = istart and only consider values >istart
% Note that this function returns imin = length(x)+1 if no lower bound
% exists, i.e. if x(end)<s

N = numel(x);
if isempty(istart)
    if x(end)<s % no lower bound exists
        imin = N+1;
        return;
    elseif x(1)>s % lower bound is 1
        imin = 1; 
        return;
    else
        imin = 1; % initialize
    end
else
    imin = istart;
end


while imin<=N && x(imin)<s 
    imin = imin+1;
%     if imin>N
%         imin = [];
%         return
%     end
end
end

function [lvec coercedX] = inrange(X,xrange,option)
% returns logical vector of size(X) of the equivalent
% X>=min(xrange)&X<=max(xrange). Also optionally returns a vector of
% coerced values coercedX which lies within the given range.
% Additional option
%   inrange(X,xrange,'inclusive')
% or
%   inrange(X,xrange,'exclusive')
% controls whether equality is allowed at the bounds.  The default behavior
% is 'inclusive'.

xrange = sort(xrange); % be sure input range is sorted

if length(xrange)~=2
    error('xrange must be a vector with 2 elements.\n');
end

if nargin>2
    if strcmpi(option,'inclusive')
        inclusive = true;
    elseif strcmpi(option,'exclusive')
        inclusive = false;
    else
        error('Unrecognized option: %s',option);
    end
else
    inclusive = true;
end
        
if inclusive
    lvec = X>=min(xrange)&X<=max(xrange);
else
    lvec = X>min(xrange)&X<max(xrange);
end

if nargout>1
    maxX = xrange(2)*ones(size(X));
    minX = xrange(1)*ones(size(X));
    coercedX = min(max(X,minX),maxX);
end
end