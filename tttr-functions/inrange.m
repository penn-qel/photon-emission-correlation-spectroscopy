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
