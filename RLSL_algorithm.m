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
b = zeros(M+1,N);
B = delta.*ones(M+1,N);
F = delta.*ones(M+1,N);
Delta = zeros(M+1,N);
gamma_s = ones(M+1,N);
rho = zeros(M+1,N);
b(1,1) = u(1);
f(1,1) = u(1);
no = 1;
for n = 1:N 
    if n > 1
        no = n-1;
    end
    b(1,n) = u(n);
    f(1,n) = u(n);
    gamma_s(1,n) = 1;
    F(1,n) = lamda*F(1,no)+u(n)^2;
    B(1,n) = F(1,n);
    for m = 2:M+1
        
            Delta(m-1,n) = Delta(m-1,no) + b(m-1,no)*f(m-1,n)/(gamma_s(m-1,no));
            gamma_f(m,n) = -Delta(m-1,n)/B(m-1,no);  
            gamma_b(m,n) = -Delta(m-1,n)/F(m-1,n);

            f(m,n) = f(m-1,n) + gamma_f(m,n)*b(m-1,no);
            b(m,n) = b(m-1,no)+ gamma_b(m,n)*f(m-1,n);
            F(m,n) = F(m-1,n) + gamma_f(m,n)*Delta(m-1,n);
            B(m,n)= B(m-1,no)+ gamma_b(m,n)*Delta(m-1,n);
            gamma_s(m,n) = gamma_s(m-1,n) - b(m-1,n)^2/B(m-1,n);
            
    end
end

for n = 1:N
    e(1,n)= d(n);
    for m = 1:M
        %rho(m,1) = 0;
        rho(m,n) = lamda*rho(m,no)+ b(m,n)*e(m,n)/gamma_s(m,n);
        kap(m,n) = rho(m,n)/B(m,n);
        e(m+1,n) = e(m,n) - kap(m,n)*b(m,n);
    end

end
eM = e(11,:);
end
