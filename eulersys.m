function [t,x,y]=eulersys(e, f, tspan, x0, y0, h, varargin)
if nargin<4, error("at least 4 inputs required"), end
ti = tspan(1);tf = tspan(2);
if ~(tf>ti), error('upper limit must be greater than lower'), end
t = (ti:h:tf)'; n = length(t);
if t(n)<tf
    t(n+1) = tf;
    n = n + 1;
end
y = y0*ones(n,1);
x = x0*ones(n,1);
for i = 1:n-1
    x(i+1) = x(i) + e(t(i), x(i), y(i))*h;
    y(i+1) = y(i) + f(t(i), x(i), y(i))*h;
end
figure(1);plot(t,x)
figure(2);plot(t,y)
end