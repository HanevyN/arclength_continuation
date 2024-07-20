% This function computes the chebyshev polynomials and their derivatives

function [z, T, dT, d2T] = cheb(N)

% Use - cos(pi * (0 : N) / N) to go from [-1, 1]

z = - cos(pi * (0 : N) / N);

T = zeros(N + 1, N + 1);
dT = T; d2T = T;

T(1, :) = 1;
T(2, :) = z;

dT(2, :) = 1;
dT(3, :) = 4 * T(2, :);

d2T(3, :) = 4 * dT(2, :);

for i = 2 : N
    T(i + 1, :) = 2 * z .* T(i, :) - T(i - 1, :);
    if i > 2
        dT(i + 1, :) = 2 * T(i, :) + 2 * z .* dT(i, :) - dT(i - 1, :);
        d2T(i + 1, :) = 4 * dT(i, :) + 2 * z .* d2T(i, :) - d2T(i - 1, :);
    end
end    

T = T';
dT = dT';
d2T = d2T';
