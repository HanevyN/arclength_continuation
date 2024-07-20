function [VR,VL] = eig_vecs(D)
    % % test with  = ~ 1e-10
    % max(abs(D*VR))
    % max(abs(VL'*D'))    
    VL = ones(length(D),1 );
    VR = VL; 
    for i = 1:3
         newR = D\VR;
         newL = D'\VL;
         VL = newL/norm(newL);
         VR = newR/norm(newR);
    end    
    
%     N = (length(D)-3)/3;
%     [~, D0, ~, ~] = cheb(N);
%     VL =[D0'*VL(1:N+1);D0'*VL(N+2:2*N+2);D0'*VL(2*N+3:end) ];
%     VR = [D0*VR(1:N+1);D0*VR(N+2:2*N+2);D0*VR(2*N+3:end) ];

    end