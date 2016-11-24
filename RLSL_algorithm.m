function [e] = RLSL_algorithm(M,N,lamda,delta,d,u)
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
    b(:,n) = u(n);
    f(:,n) = u(n);
    F(:,n) = lamda*F(:,n-1)+u(n)^2;
    e(:,n)= d(n);
    for m = 2:M+1
        Delta(m-1,n) = Delta(m-1,n) + b(m-1,n)*f(m-1,n)/(gamma_s(m-1,n));
        gamma_f(m,n) = -Delta(m-1,n)/B(m-1,n-1);
        gamma_b(m,n) = -Delta(m-1,n)/F(m-1,n);
        rho(m,n) = lamda*rho(m,n-1)+ b(m,n)/gamma_s(m,n);
        kap(m,n) = rho(m,n)/B(m,n);
        
        f(m,n) = f(m-1,n) + gamma_f(m,n)*b(m-1,n-1);
        b(m,n) = b(m-1,n-1)+ gamma_b(m,n)*f(m-1,n);
        F(m,n) = F(m-1,n) + gamma_f(m,n)*Delta(m-1,n);
        B(m,n)= B(m-1,n-1)+ gamma_b(m,n)*Delta(m-1,n);
        gamma_s(m,n) = gamma_s(m-1,n) - b(m-1,n)^2/B(m-1,n);
        % don't forget to adjust for delay
        e(m,n) = e(m-1,n) - kap(m-1,n)*b(m-1,n);
    end
end

end
