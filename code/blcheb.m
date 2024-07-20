% This function transforms the chebyshev polynomials into the boundary 
% layer domain using an exponential map

function [eta, S, dS, d2S] = blcheb(z, T, dT, d2T, N, Y_max)

% The mapping has the form z = A * (e ^{a * eta} - 1) - 1
% To ensure that [-1, 1] -> [0, Y_max] A must have the form
% A = 2 / (e ^{a * eta_max} - 1)
% The value of 'a' essentially defines the exponential map

a = - 1 / 5;
A = 2 / (exp(a * Y_max) - 1);
eta = log((z + 1 + A) / A) / a;
 
dzdeta = a * A * exp(a * eta)';
d2zdeta2 = a ^ (2) * A * exp(a * eta)';

S = T;
dS = dT .* repmat(dzdeta, 1, N + 1);
d2S = d2T .* repmat(dzdeta, 1, N + 1) .^ (2) ... 
                                        + dT .* repmat(d2zdeta2, 1, N + 1);
