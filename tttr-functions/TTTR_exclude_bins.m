function exclude = TTTR_ExcludeBins(tau_axis,exclude_range)
% Function to determine which (if any) bins defined by the vector
% 'tau_axis' overlap with the range of values given in 'exclude_range'

% Assumes both tau_axis and exlude_range are sorted.  Returns logical
% vector of size(tau_axis).


tau_widths = zeros(size(tau_axis));
dtau = diff(tau_axis);
tau_widths(1) = dtau(1);
tau_widths(2:end-1) = min(dtau(1:end-1),dtau(2:end));
tau_widths(end) = dtau(end);

tau_min = tau_axis - tau_widths./2;
tau_max = tau_axis + tau_widths./2;

exclude = (exclude_range(2)>tau_min)&(exclude_range(1)<tau_max);