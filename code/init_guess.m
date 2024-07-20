
N = 50;

% set up chebychev basis function
[z, T, dT, d2T] = cheb(N);
[eta, S, dS, d2S] = blcheb(z, T, dT, d2T, N, Y_max);

load('flow.mat')

% map base flow onto chebychev domain
U = spline(YY, F(2, :), eta);
dU = spline(YY, F(3, :), eta);

R = 2500;
omega = 0.01:0.01:0.1;
eig_vals = zeros(size(omega));

for i = 1:length(omega)

    [A0, A1, A2] = D(U,dU,Y_max, delta, N, S, dS, d2S, R, omega(i));
    eig_vals(i) = filter_eigs(polyeig(A0,A1, A2), omega(i));

end

% find sign changes
alpha_i  = imag(eig_vals);
sc = alpha_i(1:end-1).*alpha_i(2:end);

index= find(sc < 0);

% dont have this issue here but I did for Cranes flow 
% picked up spurious eigenvalues
if length(index) > 2
    index = index(2:3);
end

for i = 1:2
    x1 = omega(index(i));
    x2 = omega(index(i)+1);
    
    y1 = alpha_i(index(i) );
    y2 = alpha_i(index(i)+1 );
    
    init_omega(i) = find_zero(x1,x2,y1,y2);

    % check values
    [A0, A1, A2] = D(U,dU,Y_max, delta, N, S, dS, d2S, R, init_omega(i));
    init_alpha(i) = filter_eigs(polyeig(A0,A1, A2), omega(i));

end


