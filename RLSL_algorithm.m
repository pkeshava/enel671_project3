function [e,gamma_s] = RLSL_algorithm(M,N,lamda,delta,d,u)
%% REQUIRES
    
    % u: tap-input vector
    % d: desired signal vector
    % mu: step size parameter
    % delta: delay of channel
    % M: filter order
% Returns   
    % W: Tap weight vector
    % e: estimation error 
%% Initialization
b = zeros(M+1,N+1);
f = zeros(M+1,N+1);
B = delta.*ones(M+1,N+1);
F = delta.*ones(M+1,N+1);
Delta = zeros(M+1,N+1);
gamma_s = ones(M+1,N+1);
rho = zeros(M+1,N+1);
e = zeros(M+1,N+1);
%u = u(:);

for n = 2:N %is that right?
    % should these be n-1? 
    %u_vec = u(n:-1:n-M+1);
    b(2,n) = u(n);
    f(2,n) = u(n);
    F(2,n) = lamda*F(2,n-1)+u(n)^2;
    e(1,n)= d(n);
    gamma_s(2,n) = 1;
    B(2,n) = F(2,n);
    for m = 2:M+1
        Delta(m,n) = Delta(m,n-1) + b(m,n-1)*f(m,n)/(gamma_s(m,n-1));
        gamma_f(m+1,n) = -Delta(m,n)/B(m,n-1);
        gamma_b(m+1,n) = -Delta(m,n)/F(m,n);
        rho(m-1,n) = lamda*rho(m-1,n-1)+ b(m,n)/gamma_s(m,n)*e(m-1,n); %might need work
        kap(m-1,n) = rho(m-1,n)/B(m,n);%might need work
        
        f(m+1,n) = f(m,n) + gamma_f(m+1,n)*b(m,n-1);
        b(m+1,n) = b(m,n-1)+ gamma_b(m+1,n)*f(m,n);
        F(m+1,n) = F(m,n) + gamma_f(m+1,n)*Delta(m,n);
        B(m+1,n)= B(m,n-1)+ gamma_b(m+1,n)*Delta(m,n);
        gamma_s(m+1,n) = gamma_s(m,n) - b(m,n)^2/B(m,n);
        % NOTE IN ALL THE ESTIMATION CALCS I HAVE LEFT THE PREDICTION
        % VARIABLE m AS IS FOR PREDICTION
        e(m,n) = e(m-1,n) - kap(m-1,n)*b(m,n);
    end
end

end
