% This function constructs the matrix A that defines the eigenvalue
% problem
function A = D_R(Y_max, delta, N, S, dS, d2S, R, alpha)

% interpolate f and f' onto shebychev domain
i = sqrt(- 1);

er = - Y_max * 100 * i;

A2 = zeros(3*N + 3, 3*N+3);
A1 = A2;
A0 = A1;

for n = 0 : N
    for k = 1 : N - 1
        J = S(k + 1, n + 1);
%         dJ = dS(k + 1, n + 1);
        d2J = d2S(k + 1, n + 1);
        
        % alpha ^ 2 terms
        % u
        A2(3 * k + 1, 3 * n + 1) = -J / R^2;  % x - momentum
%         A2(3 * k + 1, 3 * n + 2) = 0;      % y - momentum
%         A2(3 * k + 1, 3 * n + 3) = 0;      % continuity
        % v
%         A2(3 * k + 2, 3 * n + 1) = 0;     
        A2(3 * k + 2, 3 * n + 2) = -J / R^2; 
%         A2(3 * k + 2, 3 * n + 3) = 0; 
        % p
%         A2(3 * k + 3, 3 * n + 1) = 0;
%         A2(3 * k + 3, 3 * n + 2) = 0;
%         A2(3 * k + 3, 3 * n + 3) = 0;

        % alpha ^ 1 terms

%         A1(3 * k + 1, 3 * n + 1) = i * U(k + 1) * J; 
%         A1(3 * k + 1, 3 * n + 2) = 0;
%         A1(3 * k + 1, 3 * n + 3) = i * J;

%         A1(3 * k + 2, 3 * n + 1) = 0;
%         A1(3 * k + 2, 3 * n + 2) = i * U(k + 1) * J; 
%         A1(3 * k + 2, 3 * n + 3) = 0;

%         A1(3 * k + 3, 3 * n + 1) = i * J; 
%         A1(3 * k + 3, 3 * n + 2) = 0; 
%         A1(3 * k + 3, 3 * n + 3) = 0;

        % alpha ^ 0 terms

        A0(3 * k + 1, 3 * n + 1) =  delta ^ 2 * d2J / R^2;
%         A0(3 * k + 1, 3 * n + 2) = delta * dU(k + 1) * J; 
%         A0(3 * k + 1, 3 * n + 3) = 0;

%         A0(3 * k + 2, 3 * n + 1) = 0;
        A0(3 * k + 2, 3 * n + 2) =  delta ^ 2 * d2J / R^2;  
%         A0(3 * k + 2, 3 * n + 3) = delta * dJ;

%         A0(3 * k + 3, 3 * n + 1) = 0; 
%         A0(3 * k + 3, 3 * n + 2) = delta * dJ; 
%         A0(3 * k + 3, 3 * n + 3) = 0;         

    end  
    
%     % Boundary conditions at the wall
%     % alpha^2 u
%     A2(1, 3 * n + 1) = S(1,n+1); % x                                                             	
%     A2(1, 3 * n + 2) = 0; % y
%     A2(1, 3 * n + 3) = 0; %ce
%     % alpha^2 v
%     A2(2, 3 * n + 1) = 0; % x                                                             
%     A2(2, 3 * n + 2) = S(1,n+1); % y 
%     A2(2, 3 * n + 3) = 0; % ce
%     % alpha^2 p
%     A2(3, 3 * n + 1) = 0; % x
%     A2(3, 3 * n + 2) = dS(1,n+1); % y
%     A2(3, 3 * n + 3) = 0; % ce                                                                                                                  
%     % alpha u
%     A1(1, 3 * n + 1) = er * S(1, n + 1); % x                                         	
%     A1(1, 3 * n + 2) = 0; % y
%     A1(1, 3 * n + 3) = 0; % ce
%     % alpha v  where er = - Y_max * 100 * i;
%     A1(2, 3 * n + 1) = 0; % x
%     A1(2, 3 * n + 2) = er * S(1, n + 1);  % y 
%     A1(2, 3 * n + 3) = 0; % ce
%     % alpha p
%     A1(3, 3 * n + 1) = 0;   % x                                    
%     A1(3, 3 * n + 2) = er * dS(1, n + 1); % y 
%     A1(3, 3 * n + 3) = 0;  % ce                                                                                                                                                                                                   
%     % u
%     A0(1, 3 * n + 1) = S(1, n + 1);                                                                           
%     A0(1, 3 * n + 2) = 0;
%     A0(1, 3 * n + 3) = 0;
%     % v
%     A0(2, 3 * n + 1) = 0;                                                                    
%     A0(2, 3 * n + 2) = S(1, n + 1);
%     A0(2, 3 * n + 3) = 0; 
%     % p
%     A0(3, 3 * n + 1) = 0;                                                  
%     A0(3, 3 * n + 2) = dS(1, n + 1); 
%     A0(3, 3 * n + 3) = 0;
%                                                                                                                                                                    
%     % Boundary conditions at infinity
%     % alpha^2 u
%     A2(3 * N + 1, 3 * n + 1) = S(N+1, n+1);                           
%     A2(3 * N + 1, 3 * n + 2) = 0;
%     A2(3 * N + 1, 3 * n + 3) = 0;
%     % alpha^2 v
%     A2(3 * N + 2, 3 * n + 1) = 0;                                       
%     A2(3 * N + 2, 3 * n + 2) = S(N+1,n+1);
%     A2(3 * N + 2, 3 * n + 3) = 0;
%     % alpha^2 p
%     A2(3 * N + 3, 3 * n + 1) = 0;                                        
%     A2(3 * N + 3, 3 * n + 2) = 0;
%     A2(3 * N + 3, 3 * n + 3) = S(N+1,n+1);
%     % alpha u
%     A1(3 * N + 1, 3 * n + 1) = er * S(N + 1, n + 1);    
%     A1(3 * N + 1, 3 * n + 2) = 0;                       
%     A1(3 * N + 1, 3 * n + 3) = 0; 	
%     % alpha v
%     A1(3 * N + 2, 3 * n + 1) = 0;                       
%     A1(3 * N + 2, 3 * n + 2) = er * S(N + 1, n + 1);    
%     A1(3 * N + 2, 3 * n + 3) = 0;
%     % alpha p
%     A1(3 * N + 3, 3 * n + 1) = 0;                      
%     A1(3 * N + 3, 3 * n + 2) = 0;                       
%     A1(3 * N + 3, 3 * n + 3) = er * S(N + 1, n + 1);    
%     % u
%     A0(3 * N + 1, 3 * n + 1) = S(N + 1, n + 1);         
%     A0(3 * N + 1, 3 * n + 2) = 0;                       
%     A0(3 * N + 1, 3 * n + 3) = 0;                      	
%     % v
%     A0(3 * N + 2, 3 * n + 1) = 0;                       
%     A0(3 * N + 2, 3 * n + 2) = S(N + 1, n + 1);         
%     A0(3 * N + 2, 3 * n + 3) = 0;                       
%     % p
%     A0(3 * N + 3, 3 * n + 1) = 0;                       
%     A0(3 * N + 3, 3 * n + 2) = 0;                       
%     A0(3 * N + 3, 3 * n + 3) = S(N + 1, n + 1);         

end

A = A0 +alpha*A1 + alpha^2*A2;


end