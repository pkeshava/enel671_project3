%% REQUIRES
    
    % u: tap-input vector
    % d: desired signal vector
    % mu: step size parameter
    % delta: delay of channel
    % M: filter order
% Returns   
    % W: Tap weight vector
    % e: estimation error 
clear all 
close all
clc

M = 11;
delta = 0.01;
N = 600;
K = 350;
h = [0.2194 1.0 0.2194;0.2798 1.0 0.2798;0.3365 1.0 0.3365;0.3887 1.0 0.3887];
lamda = 1;
gamma_ave = zeros(N,1);
for k=1:K
    a = BPSK(N);
    u = filterinput(a,h); 
    u = u(:,1);
    %Initialization
    b = zeros(N,M);
    f = zeros(N,M);
    B = delta.*ones(N,M);
    F = delta.*ones(N,M);
    Delta = zeros(N,M);
    gamma_s = ones(N,M);
    gamma_f = zeros(N,M);
    gamma_b = zeros(N,M);
    
    for n = 2:N 
        b(n,1) = u(n);
        f(n,1) = u(n);
        F(n,1) = lamda*F(n-1,1)+(u(n))^2;
        gamma_s(n,1) = 1;
        B(n,1) = F(n,1);
        for m = 2:M
            Delta(n,m-1) = Delta(n-1,m-1) + b(n-1,m-1)*f(n,m-1)/(gamma_s(n-1,m-1));
            gamma_f(n,m) = -Delta(n,m-1)/B(n-1,m-1);
            gamma_b(n,m) = -Delta(n,m-1)/F(n,m-1);

            f(n,m) = f(n,m-1) + gamma_f(n,m)*b(n-1,m-1);
            b(n,m) = b(n-1,m-1)+ gamma_b(n,m)*f(n,m-1);
            F(n,m) = F(n,m-1) + gamma_f(n,m)*Delta(n,m-1);
            B(n,m)= B(n-1,m-1)+ gamma_b(n,m)*Delta(n,m-1);
            gamma_s(n,m) = gamma_s(n,m-1) - (b(n,m-1))^2/B(n,m-1);
        end
    end

rho = zeros(N,M);
e = zeros(N,M);
alpha = zeros(N,M);
kap = zeros(N,M);
for n = 2:N
    e(n,1)= a(n);
    for m = 1:M
        rho(1,m) = 0;
        rho(n,m) = lamda*rho(n-1,m)+ b(n,m)/gamma_s(n,m)*e(n,m);
        kap(n,m) = rho(n,m)/B(n,m);
        e(n,m+1) = e(n,m) - kap(n,m)*b(n,m);
    end

end
 alpha = e(:,11)./gamma_s(:,11);
 gamma_ave = (gamma_ave + gamma_s(:,11));
 alphasum(:,k) = alpha.^2;
 MSEE11 = sum(alphasum,2)/K;
 MSEE11n = MSEE11./gamma_s(:,11);
end
gamma_ave = gamma_ave/K;

figure(1)
plot(1:N,gamma_ave,'LineWidth',2)
legend('Channel 1','Channel 2','Channel 3','Channel 4')
grid on
xlabel('Time (s)');
ylabel('Gamma'); 
title('Likilihood');

figure(2)
semilogy(1:N,MSEE11n,'LineWidth',2)
title('Normalized MSEE');