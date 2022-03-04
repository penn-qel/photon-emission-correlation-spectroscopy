function T2data_out = TTTR_fold_T2data(T2data,tauZero,exclude)
% Function to fold correlation data about zero-delay symmetry point and
% recalculate g2 along with uncertainties.

sym_tol = 1e-12; 

if nargin<3
    exclude = false(size(T2data.tauAxis));
end

% Check that tau_axis is symmetric
tauFull = T2data.tauAxis(:)-tauZero;

if floor(length(tauFull)/2)==length(tauFull)/2 % if even number of elements (no bin at zero delay)
    zeroix = length(tauFull)/2+1;
    is_symmetric = flipud(tauFull(1:(zeroix-1)))+tauFull((zeroix):end)<sym_tol;
    contains_zero = false;
else
    zeroix = ceil(length(tauFull)/2);
    if abs(tauFull(zeroix))>sym_tol
        error('Central bin of provided tau_axis must equal tau_zero.');
    end
    is_symmetric = flipud(tauFull(1:(zeroix-1)))+tauFull((zeroix+1):end)<sym_tol;
    contains_zero = true;
end
   
if any(~is_symmetric)
    error('Field "tau_axis" is not symmetric about the provided tau_zero value');
end

tauAbs = tauFull(zeroix:end);
nP = T2data.nPartitions;
nBins = length(tauAbs);
xCorr = zeros(nBins,nP);
g2norm = zeros(nBins,nP);

if contains_zero
    xCorr(1) = T2data.xCorr(zeroix,:);
    g2norm(1) = T2data.g2norm(zeroix,:);
    for tt = 2:nBins
        %     ixs = zeroix+(tt-1)*[-1,1];
        ixm = zeroix - tt + 1;
        ixp = zeroix + tt -1;
        xCorr(tt,:) = ~exclude(ixm)*T2data.xCorr(ixm,:)+~exclude(ixp)*T2data.xCorr(ixp,:);
        g2norm(tt,:) = ~exclude(ixm)*T2data.g2norm(ixm,:)+~exclude(ixp)*T2data.g2norm(ixp,:);
    end
else
    for tt = 1:nBins
        %     ixs = zeroix+(tt-1)*[-1,1];
        ixm = zeroix - tt;
        ixp = zeroix + tt -1;
        xCorr(tt,:) = ~exclude(ixm)*T2data.xCorr(ixm,:)+~exclude(ixp)*T2data.xCorr(ixp,:);
        g2norm(tt,:) = ~exclude(ixm)*T2data.g2norm(ixm,:)+~exclude(ixp)*T2data.g2norm(ixp,:);
    end
end


% Calculate asymmetric errors for Poisson-distributed data following:
%   lower-error = -0.5+sqrt(0.25+N)
%   upper-error = +0.5+sqrt(0.25+N)
Poiss_errs = sqrt(0.25+xCorr); 
Xcorr_err = cat(3,-0.5+Poiss_errs,0.5+Poiss_errs);

g2 = xCorr./g2norm;
g2Lerr = Xcorr_err(:,:,1)./g2norm;
g2Uerr = Xcorr_err(:,:,2)./g2norm;

% Store merged results
T2data_out = T2data; % initialize with values for count rates, etc.

T2data_out.nBins = nBins;
T2data_out.tauAxis = tauAbs;
T2data_out.xCorr = xCorr;
T2data_out.g2 = g2;
T2data_out.g2norm = g2norm;
T2data_out.g2Lerr = g2Lerr;
T2data_out.g2Uerr = g2Uerr;

