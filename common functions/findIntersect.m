function x0 = findIntersect(y,x, y0, t0)
% function to find the intersect of y(x) and y0 (constant)
% t0, earliest time for a possible rise time
yIntp = interp1(x,y,x(1):0.01:x(end));
tIntp = x(1):0.01:x(end);
t = yIntp-y0 > 0;
% find the intersection point that is closest to the cue and increase
i0 = find(diff(t(:))==1);
if isempty(i0)
    x0 = NaN;
else
    intersect = tIntp(i0);
    % only keep the t>0 part
    i0 = i0(intersect>t0);
    [~,ind] = min(intersect(intersect>t0));
    x0 = (tIntp(i0(ind))+tIntp(i0(ind)+1))/2;
end

