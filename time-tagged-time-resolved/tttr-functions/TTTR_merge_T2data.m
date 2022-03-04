function T2data = TTTR_merge_T2data(T2datasets)
% Function to merge correlation results from multiple exposures.
%
% Input structure 'T2datasets' should be a structure array with an element
% for each dataset, produced by calling TTTR_ProcessT2data with the same
% options (i.e., the correlation bin definitions should be identical).  The
% function does not check that this condition is met.
% 
% Output structure is of the same form, but containing the combined
% count-rate and correlation information.

nSets = length(T2datasets);

T2data = T2datasets(1); % initialize output structure

% Concatenate count rates
if isfield(T2datasets,'countRates')
    partitions = [T2datasets.nPartitions]; % these could possibly be different for each set 
    T2data.nPartitions = max(partitions);
    nTimes = arrayfun(@(x)length(x.tAxis),T2datasets,'UniformOutput',1); % vector giving length of each trace
    T2data.tAxis = zeros(1,sum(nTimes));
    T2data.countRates = zeros(2,sum(nTimes));
    if isfield(T2data,'tFlag')
        T2data.tFlag = zeros(1,sum(nTimes));
    end
    lastTime = 0;
    lastIdx = 0;
    for sIdx = 1:nSets
        T2data.tAxis(lastIdx+(1:nTimes(sIdx))) = lastTime + T2datasets(sIdx).tAxis;
        T2data.countRates(:,lastIdx+(1:nTimes(sIdx))) = T2datasets(sIdx).countRates;
        if isfield(T2data,'tFlag')
            T2data.tFlag(lastIdx+(1:nTimes(sIdx))) = T2datasets(sIdx).tFlag;
        end
        lastTime = lastTime + T2datasets(sIdx).tAxis(end) + diff(T2datasets(sIdx).tAxis(end-1:end))/2; % need to add half of final bin width to get ending time
        lastIdx = lastIdx+nTimes(sIdx);
    end
end
   
% Merge correlation data - this is the current version on the gitlab 1/2022
% but can't handle partitions...
%xCorr = sum([T2datasets.xCorr],2);
%g2norm = sum([T2datasets.g2norm],2);

% Merge correlation data
xCorr = 0;
g2norm = 0;
for sIdx = 1:nSets
    xCorr = xCorr + T2datasets(sIdx).xCorr;
    g2norm = g2norm + T2datasets(sIdx).g2norm;
end

% Calculate asymmetric errors for Poisson-distributed data following:
%   lower-error = -0.5+sqrt(0.25+N)
%   upper-error = +0.5+sqrt(0.25+N)
PoissErrs = sqrt(0.25+xCorr); 
xCorrErr = cat(3,-0.5+PoissErrs,0.5+PoissErrs);

g2 = xCorr./g2norm;
g2Lerr = xCorrErr(:,:,1)./g2norm;
g2Uerr = xCorrErr(:,:,2)./g2norm;

% Store merged results
T2data.xCorr = xCorr;
T2data.g2 = g2;
T2data.g2norm = g2norm;
T2data.g2Lerr = g2Lerr;
T2data.g2Uerr = g2Uerr;