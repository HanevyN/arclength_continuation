function [lambda_r, lambda_w, R_s, w_s] = D_lambda(AR,Aw, Aa, VR, VL, dir)
    % inital_guess to newton iteration         


    lambda_r =  -imag((VL'*AR*VR)/(VL'*Aa*VR));
    lambda_w = -imag((VL'*Aw*VR)/(VL'*Aa*VR));
    
    R_s = dir*sqrt(1 -  abs(lambda_r/lambda_w)^2);
    w_s = -dir*lambda_r/lambda_w;


end