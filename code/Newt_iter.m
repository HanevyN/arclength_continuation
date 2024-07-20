function [R, alpha, omega,dir] =  Newt_iter(x0, x1, VR,VL, Y_max, delta, N,S, dS,d2S, U,dU, alpha,w_s,R_s,ds, dir)
% VR , VL left and right eigenvecors -- won'tchange as evelauted at last
%         stable point
% x0 - given last stable point
% x1 - initial guess for new stable point
%    - updated until imag alpha gets small

% x1_init = x1;
count  = 0; 


H = [imag(alpha);
    ([w_s; R_s']'*(x1 - x0) - ds)];

 while abs(H(1)) > 1e-6
    if count > 2
        ds = ds/3;
        x1 = x0 + ds*[w_s; R_s];
        count = 0;
    end


    count = count + 1;
    AR = D_R(Y_max, delta, N, S, dS, d2S, x1(2), alpha);
    Aw = D_omega(Y_max, N, S, dS, alpha);
    Aa = D_alpha(U,Y_max, N, S, dS, x1(2), alpha);
% 
        % set up initial guess for next step
 
    [lambda_r, lambda_w, ~, ~] = D_lambda(AR,Aw, Aa, VR, VL, dir); 
 

    J = [ lambda_w, lambda_r; w_s, R_s ];

    x1 = x1 - J\H;
    % update guess H based on new guess
    [A0, A1, A2] = D(U,dU,Y_max, delta, N, S, dS, d2S, x1(2), x1(1));    
    alpha = filter_eigs(polyeig(A0, A1, A2 ), x1(1) );

    H = [imag(alpha);
         ([w_s; R_s']'*(x1 - x0) - ds)];

if count == 1
     if x1(2) < 100
        fprintf('%s %13s %10s %8s \n', 'R', 'omega', 'alpha_r', 'alpha_i')
    elseif x1(2) >= 100 && x1(2) < 1000
        fprintf('%s %14s %10s %8s \n', 'R', 'omega', 'alpha_r', 'alpha_i')
    elseif x1(2) >= 1000 && x1(2) < 10000
        fprintf('%s %15s %10s %8s \n', 'R', 'omega', 'alpha_r', 'alpha_i')
    else
        fprintf('%s %16s %10s %8s \n', 'R', 'omega', 'alpha_r', 'alpha_i')
     end
end
    fprintf('%.4f   %.4f   %.4f   %.4f \n', x1(2), x1(1), real(alpha), imag(alpha))

         if  abs(x1(1) - x0(1)) > 0.1
             dir = - dir;
             x1  = x0;
             x1(1) = x1(1) + 1e-3*sign(w_s);
             R = x1(2);
             omega = x1(1);
            [A0, A1, A2] = D(U,dU,Y_max, delta, N, S, dS, d2S, x1(2), x1(1));    
            alpha = filter_eigs(polyeig(A0, A1, A2 ), x1(1) );

            return
        end
end

% update direction vector one iteration has converged

R = x1(2);
omega = x1(1);
end

