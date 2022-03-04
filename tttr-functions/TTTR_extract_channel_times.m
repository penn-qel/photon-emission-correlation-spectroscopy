function parsedData = TTTR_extract_channel_times(TTTRData, options)
% Function to parse a TTTRData structure (from TTTR_ImportPTU) into the
% corresponding event times, with overflows, along with the number of
% syncs and total counts for averaging and other information.
%
% Inputs:
%   TTTRData - the matlab structure returned from TTTR_ImportPTU 
%   options - options structure with fields corresponding to
%       integerTime - (T/F) whether to keep the times as a 64-bit integer
%                     default is double precision
%       verbose     - (T/F) whether to print updates to the console
%                     default is to be verbose
%       vectorized  - (T/F) whether input is vectorized
%                     default is vectorized
%       
% Output Structure:
%   chan1Times = chan1Times;
%   chan0Times = chan0Times;
%   nChan0     = nChan0;
%   nChan1     = nChan1;
%   nOverflows = nOverflows;
%   globalResolution = globalResolution
%   acquisitionTime = acquisitionTime
%
%
% Last updated by DAH 12/18/2017
%--------------------------------------------------------------------------

% Parse the inputs
if nargin == 1
    options = struct;
end
options = parse_inputs(options);

% necessary constants
%OVERFLOW_HEADER = uint8(binaryVectorToDecimal([1 1 1 1]));  %changing this
%for use on mac -RF
%OVERFLOW_HEADER = uint8(binaryVectorToDecimal([1 1 1 1])); 

%WRAPAROUND = uint64(binaryVectorToDecimal([1 1 0 0 1 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0])); % time to add on each overflow event

OVERFLOW_HEADER = uint8(bin2dec(num2str([1 1 1 1]))); 

WRAPAROUND = uint64(bin2dec(num2str([1 1 0 0 1 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]))); % time to add on each overflow event


% Important values from Tags
globalResolution = TTTRData.Tags.MeasDesc_GlobalResolution; % should be 4e-12 s
acquisitionTime = 1e-3*double(TTTRData.Tags.MeasDesc_AcquisitionTime); % given in ms, convert to s
binningFactor = TTTRData.Tags.MeasDesc_BinningFactor; % this is provided, but does it matter for T2 measurements?

if options.verbose
    fprintf(1,'Extracting times...\n');
    tic
end

% Expand event records into photon time events flagged by channel
nEvents = length(TTTRData.Events);

% Only first 4 bits determine the channel flags. I think the uint8 is to
% reduce this thing in memory?
eventFlags = uint8(bitshift(TTTRData.Events,-28));

% Grab the indices for channel 0, channel1, and any overflows. == is faster
% than ismember
chan0Markers = eventFlags == 0;
chan1Markers = eventFlags == 1;
overflowMarkers = uint64(eventFlags == OVERFLOW_HEADER);

% find the number of syncs, counts and overflows
nChan0 = sum(chan0Markers);
nChan1 = sum(chan1Markers);
nOverflows = sum(overflowMarkers);

if options.verbose
   fprintf(1, '# Events: %d, # Channel 0: %d, # Channel 1: %d, # Overflows: %d\n', nEvents, nChan0, nChan1, nOverflows) 
end

if options.vectorized
    % Fast overflow calculation - convert overflowMarker indicators to a
    % cumuluative sum
    % example:
    %
    % overflowMarker: 0 0 0 0 0 1 0 0 0 1 0 0
    % cumsum ->
    % overflowMArker: 0 0 0 0 0 1 1 1 1 2 2 2
    
    % I am so very clever.
    % This is Richard's reinvention of cumsum -- keep it here as a reminder
    % to not be clever. - DAH
    %     Nevents = length(TTTRData.Events);
    %     Nmarkers = [find(overflowMarkers); Nevents+1];
    %     for loop_index = 1:(length(Nmarkers)-1)
    %             overflowMarkers(Nmarkers(loop_index):Nmarkers(loop_index+1)-1) = loop_index.*ones(size(overflowMarkers(Nmarkers(loop_index):Nmarkers(loop_index+1)-1)));
    %     end
    
    overflowMarkers = cumsum(overflowMarkers);
    
    % Find the global timeData vector as a 64-bit vector. The overflowMarkers
    % are multipled to take into account the WRAPAROUND efficiently
    %timeData = uint64(bitand(TTTRData.Events, binaryVectorToDecimal(ones(1, 28)))) ...
    %    + overflowMarkers.*WRAPAROUND;
    % changing for use on mac -RF
    timeData = uint64(bitand(TTTRData.Events, bin2dec(num2str(ones(1, 28))))) ...
        + overflowMarkers.*WRAPAROUND;
    
    
    % Pull out the two channel time vectors from the Event vector using logical indexing
    chan0Times = timeData(chan0Markers);
    chan1Times = timeData(chan1Markers);
    
    % Convert to double if specified
    if ~options.integerTime
        chan0Times = double(chan0Times);
        chan1Times = double(chan1Times);
    end
    
    % scale for time resolution
    chan0Times = chan0Times * globalResolution;
    chan1Times = chan1Times * globalResolution;
    
else % old for loop method for efficiency comparison
    % Expand event records into photon time events flagged by channel
    Nevents = length(TTTRData.Events);
    
    % Iterate through each event
    % % Parse event data
    event_flag = zeros(Nevents,1,'uint8'); % preallocate 1-byte array
    % T2_data = zeros(Nevents,1,'unit32'); % preallocate
    for ii=1:Nevents
        event_flag(ii) = uint8(bitand(bitshift(TTTRData.Events(ii),-28),15)); % highest 4 bits encode channel/overflow/marker information
        %     T2_data(ii) = bitand(TTTRData.Events(ii),268435455);
    end
    % event_flag = uint8(bitand(bitshift(TTTRData.Events,-28),15)); % highest 4 bits encode channel/overflow/marker information
    
    nChannel0 = sum(ismember(event_flag,0));
    nChannel1 = sum(ismember(event_flag,1));
    Nmarkers = sum(ismember(event_flag,15));
    
    if nChannel0+nChannel1+Nmarkers ~= Nevents
        error('Error interpreting T2 event data. Unrecognized marker flag in data stream.');
    end
    
    % Initialize arrays for counts
    if options.integerTime% use integer types
        chan0Times = zeros(nChannel0,1,'uint64');
        chan1Times = zeros(nChannel1,1,'uint64');
    else % double precision
        chan0Times = zeros(nChannel0,1); % initialize arrays (double precision)
        chan1Times = zeros(nChannel1,1);
    end
    
    count0 = 0;
    count1 = 0;
    nOverflows = 0;
    markerCounts = 0;
    
    ofltime = uint64(0); % (in GlobalResolution units)
    WRAPAROUND = uint64(210698240); % time to add on each overflow event
    
    for ii=1:Nevents
        time_data = uint64(bitand(TTTRData.Events(ii),268435455)); % lowest 28 bits of record, converted to 64bit integer
        switch event_flag(ii)
            case 0
                chan0Times(count0+1) = time_data+ofltime;
                count0 = count0+1;
            case 1
                chan1Times(count1+1) = time_data+ofltime;
                count1 = count1+1;
            case 15
                markers = bitand(time_data,15); % lowest 4 bits are marker bits
                if markers ==0 % overflow record
                    ofltime = ofltime + WRAPAROUND;
                    nOverflows = nOverflows+1;
                else % otherwise this is a true marker
                    markerCount = markerCount+1; % for now do not do anything with markers, but record their total number
                end
            otherwise
                error('Unrecognized channel flag in data stream.');
        end
    end
    
end

if options.verbose
   fprintf(1, 'Extraction completed in %.1f seconds.\n', toc) 
end

% Pack results for returning
parsedData.chan0Times = chan0Times;
parsedData.chan1Times = chan1Times;
parsedData.nChan0 = nChan0;
parsedData.nChan1 = nChan1;
parsedData.nOverflows = nOverflows;
parsedData.globalRes = globalResolution;
parsedData.acquisitionTime = acquisitionTime;

end

function validatedOptions = parse_inputs(options)
% parse the function inputs

defaultOptions.integerTime = false;
defaultOptions.verbose = true;
defaultOptions.vectorized = true;

% Parse inputs
parseInputs = inputParser;

addParameter(parseInputs,'integerTime',defaultOptions.integerTime,...
@(x) validateattributes(x,{'logical'},{'scalar'}));

addParameter(parseInputs,'verbose',defaultOptions.verbose,...
@(x) validateattributes(x,{'logical'},{'scalar'}));

addParameter(parseInputs, 'vectorized', defaultOptions.vectorized,...
    @(x)validateattributes(x, {'logical'},{'scalar'}));

parse(parseInputs, options);

validatedOptions = parseInputs.Results;



end