clearvars; clc
% supress matrix inversion warnings
warning('off','MATLAB:nearlySingularMatrix')
%%%% NB:: turn back on after !!!!!!!!!!
% warning('on','MATLAB:nearlySingularMatrix')
load("flow.mat")

% init_omega = 0.0962;
init_omega = 0.028;
iR = 2500;
N = 50;

% set up chebychev basis function
[z, T, dT, d2T] = cheb(N);
[eta, S, dS, d2S] = blcheb(z, T, dT, d2T, N, Y_max);

% map base flow onto chebychev domain
U = spline(YY, F(2, :), eta);
dU = spline(YY, F(3, :), eta);

% initial guesses for first iteration
R = iR;
omega = init_omega;

[A0, A1, A2] = D(U,dU,Y_max, delta, N, S, dS, d2S, R, omega);

% check initial input is ok
e1 = filter_eigs(polyeig(A0,A1,A2), omega);
imag(e1)

% want to loop back to initial point
alpha = e1;
count = 0;
dir = -1;

while R <= iR
    count = count + 1;
    % store answers after each iteration
    R_vec(count)  = R;
    alpha_vec(count) = alpha;
    omega_vec(count) = omega;
    
    plot(R_vec, omega_vec, 'k', 'LineWidth',3)
    xticks( 0 : iR / 5 : iR)
    yticks(0 : 0.04 : 0.2)
    axis([0 iR 0 0.2])
    ylabel("$\omega$", 'Interpreter', 'latex')
    xlabel("$R$", 'Interpreter', 'latex')
    pause(0.001)
    hold on

    % set up initial guess for next step
    [A0, A1, A2] = D(U,dU,Y_max, delta, N, S, dS, d2S, R, omega);

    
    DD = A0 + alpha*A1 + alpha^2*A2;
    [VR, VL] =eig_vecs(DD);
    
    AR = D_R(Y_max, delta, N, S, dS, d2S, R, alpha);
    Aw = D_omega(Y_max, N, S, dS, alpha);
    Aa = D_alpha(U,Y_max, N, S, dS, R, alpha);
       
    [~, ~, Rs, ws] = D_lambda(AR,Aw, Aa, VR, VL, dir);
    
    % previous guess
    x0 = [omega; R]; 
    % intial guess for next point

    ds =  abs(Rs/ws/2000);

    x1 = x0 + ds*[ ws; Rs]; 
    
    [A0, A1, A2] = D(U,dU,Y_max, delta, N, S, dS, d2S, x1(2), x1(1)); 
    alpha = filter_eigs(polyeig(A0, A1, A2 ), x1(1) );

    % second while loop to update eigenvalues
    % R and omega
    [R, alpha, omega, dir] =  Newt_iter(x0, x1, VR,VL, Y_max, delta, N,S, dS, d2S, U,dU, alpha, ws, Rs,ds, dir);
    imag_alpha  = abs(imag(alpha));

end

