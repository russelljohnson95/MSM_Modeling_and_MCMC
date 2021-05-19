function [b] = deriv(a,deltat)
% DERIV  Numerical first derivative.
%   B = DERIV(A,DELTAT) returns the first derivative of the vector A.
%   The interval between samples is specified by deltat (in sec), and
%   sampling rate is 1/deltat. Derivatives at first & last points are
%   estimated using second order forward & backwards differences, and
%   all other points use first order central differences.
%
%   Written by Brian R. Umberger, June 2000

% make sure we have been passed a row or column vector
if isvector(a) == 0
   error('Argument to function deriv must be a vector');
end

delta2 = deltat*2;   % twice the sampling interval
b = zeros(size(a));  % space to store the results

% second order forward difference approximation for sample 1
b(1) = (-3*a(1) + 4*a(2) - a(3))/delta2;

%first central difference approximation for samples >1 and <n
n = length(a);
nMinusOne = length(a) - 1;
for i = 2:nMinusOne
	b(i) = (a(i+1) - a(i-1))/delta2;
end

% second order backwards difference approximation for sample n
b(n) = (a(n-2) - 4*a(n-1) + 3*a(n))/delta2;
