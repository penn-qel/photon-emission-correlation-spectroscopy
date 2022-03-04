function [g2Bksub,g2BksubErr,rhoMean,rhoErr,T2data] = TTTR_bksub(T2data,rho_fit)
%function [g2Bksub,g2BksubErr,rhoMean,rhoErr,T2data] = TTTR_bksub(AllT2data,AllT2bkgd,rho_fit,T2data)

% background subtracted g2 is calculated based on equations described in
% Brouri et al., Optics Letters 25, 17 (2000).
% 
% Calculation of rho came from performing algebraic transformations on the
% typical definition of rho = S/(S + B)

%if T2data.bkgdFlag
%    B = [AllT2bkgd(:).countRates];
%    T = [AllT2data(:).countRates];
%    rho = 1 - sum(B)./sum(T);
%    rhoMean = mean(rho);
%    rhoErr = std(rho)/length(rho);
%    T2data.rhoMean = rhoMean;
%    T2data.rhoErr = rhoErr;
%elseif isfield(T2data,'rhoMeanTracking')
if isfield(T2data,'rhoMeanTracking')
    rhoMean = T2data.rhoMeanTracking;
    rhoErr = T2data.rhoErrTracking;
else
    rhoMean = rho_fit;
    rhoErr = 'unassigned';
end



if isfield(T2data,'rhoMeanTracking')
%if T2data.bkgdFlag || isfield(T2data,'rhoMeanTracking')

    % Calculate background corrected g2 using error propogation (Lee'
    % ErrorPropogation function)
    
    g2Err = mean([T2data.g2Lerr,T2data.g2Uerr],2);
    
    for ii = 1:length(T2data.g2)
    n = ErrorPropagation(@(A,B)((A - (1 - B^2))/B^2),...
        [T2data.g2(ii),g2Err(ii)],[rhoMean,rhoErr]);
            g2Bksub(ii,1) = n(1); 
            g2BksubErr(ii,1) = diff(n(2:3))/2;
    end

    T2data.g2Bksub = g2Bksub;
    T2data.g2BksubErr = g2BksubErr;
    
else
    g2Bksub = (T2data.g2 - (1 - (rhoMean^2)))./(rhoMean^2);
    g2BksubErr = mean([(T2data.g2Lerr/rhoMean^2),(T2data.g2Uerr/rhoMean^2)],2);
end    

end