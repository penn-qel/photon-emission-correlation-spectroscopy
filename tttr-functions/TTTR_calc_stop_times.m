function outData = TTTR_calc_stop_times(parsedData, syncChannel, nVars)
% Function to calculate the start-stop times for lifetime/transient
% measuremetns with the picoHarp. Takes in parsed data from
% TTTR_extract_channel_time and has a few options that lets you choose what
% the start and stop channels are.
%
% Inputs
%   parsedData - matlab structure from TTTR_extract_channel_time. Can be an
%   array for batch processing
%   syncChannel - (0/1) which channel is the sync channel 0 default
%
% Output structure array fields:
%   stopTimes - a vector of all of the relative stop times for each sweep
%
% Last updated by DAH 12/18/2017
%--------------------------------------------------------------------------

% Handle defaults
if nargin == 1
    syncChannel = 0;
    nVars = 1;
elseif nargin == 2
    nVars = 1;
end

nSweeps = length(parsedData);
% outData(nSweeps) = struct();

for iSweep = 1:nSweeps
    
    % break into the sync and count times.
    switch syncChannel
        case 0
            syncTimes = parsedData(iSweep).chan0Times;
            countTimes = parsedData(iSweep).chan1Times;
        case 1
            syncTimes = parsedData(iSweep).chan1Times;
            countTimes = parsedData(iSweep).chan0Times;
        otherwise
            error('invalid syncChannel. 0 or 1 only.')
    end
    
    % remove counts that were recorded before the first sync
    firstSync = syncTimes(1);
    syncedIndices = countTimes > firstSync;
    countTimes = countTimes(syncedIndices);
    
    % Pre-allocate arrays for speed and outputting
    outData(iSweep).stopTimes = zeros(length(countTimes), 2);
    relativeArrival = zeros(length(countTimes), 2);
    
    syncIdx = 1;
    
    % create the variable idx
    varIdx = mod(1:length(syncTimes), nVars);
    varIdx(varIdx == 0) = nVars;
    

    for iEvent = 1:length(countTimes)
        
        % Check if a new sync is needed, don't look for a new sync if we
        % are at the last one
        if syncIdx >= length(syncTimes)
            newSync = 0;
        else
            newSync = countTimes(iEvent) > syncTimes(syncIdx + 1);
        end
        
        % Iterate the syncIdx until it either reaches the end, or the
        % current countTime is aligned to the syncIdx
        while newSync
            syncIdx = syncIdx + 1;
            if syncIdx == length(syncTimes)
                % break out of the loop if the sync reach the end
                newSync = 0;
                
            elseif countTimes(iEvent) > syncTimes(syncIdx+1)
                % check the next sync
                newSync = 1; 
                
            else
                % valid sync found, break out of loop
                newSync = 0;
               
            end
        end
        
         % subtract and add
        relativeArrival(iEvent, 1) = countTimes(iEvent) - syncTimes(syncIdx);
        relativeArrival(iEvent, 2) = varIdx(syncIdx);
        
    end
    
    % Store data for output after iterating through. Assigning to this
    % structure array is time consuming compared to just copying at the
    % end
    outData(iSweep).stopTimes = relativeArrival;
    
    for i = 1:nVars
        outData(iSweep).nSyncs(i) = sum(varIdx == i);
    end

end