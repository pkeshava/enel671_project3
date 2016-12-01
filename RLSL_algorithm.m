function [alpha,gamma_s] = RLSL_algorithm(M,N,lamda,delta,d,u)
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
b = zeros(M+1,N);
f = zeros(M+1,N);
B = delta.*ones(M+1,N);
F = delta.*ones(M+1,N);
Delta = zeros(M+1,N);
gamma_s = ones(M+1,N);
gamma_f = zeros(M+1,N);
gamma_b = zeros(M+1,N);
erroraverage = 0;
for n = 2:N 

    b(1,n) = u(n);
    f(1,n) = u(n);
    F(1,n) = lamda*F(1,n-1)+(u(n))^2;
    gamma_s(1,n) = 1;
    B(1,n) = F(1,n);
    for m = 2:M+1
        Delta(m-1,n) = Delta(m-1,n-1) + b(m-1,n-1)*f(m-1,n)/(gamma_s(m-1,n-1));
        gamma_f(m,n) = -Delta(m-1,n)/B(m-1,n-1);
        gamma_b(m,n) = -Delta(m-1,n)/F(m-1,n);
             
        f(m,n) = f(m-1,n) + gamma_f(m,n)*b(m-1,n-1);
        b(m,n) = b(m-1,n-1)+ gamma_b(m,n)*f(m-1,n);
        F(m,n) = F(m-1,n) + gamma_f(m,n)*Delta(m-1,n);
        B(m,n)= B(m-1,n-1)+ gamma_b(m,n)*Delta(m-1,n);
        gamma_s(m,n) = gamma_s(m-1,n) - (b(m-1,n))^2/B(m-1,n);
            
    end
end

rho = zeros(M+1,N);
e = zeros(M,N);
alpha = zeros(M,N);
Kappa = zeros(M+1,N);

for n = 2:N
    e(1,n)= d(n);
    rho(m,1) = 0;
    for m = 1:M+1
        rho(m,n) = lamda*rho(m,n-1)+ b(m,n)/gamma_s(m,n)*e(m,n);
        kap(m,n) = rho(m,n)/B(m,n);
        e(m+1,n) = e(m,n) - kap(m,n)*b(m,n);
    end

end
 alpha = e(11,:)./gamma_s(11,:);
 errorsquared = alpha.^2;
 %erroraverage = erroraverage + errorsquared; 
end
