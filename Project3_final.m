%% REQUIRES
clear all 
close all
clc

M = 11;
delta = 0.01;
N = 800;
K = 350;
h = [0.2194 1.0 0.2194;0.2798 1.0 0.2798;0.3365 1.0 0.3365;0.3887 1.0 0.3887];
lamda = 1;
gamma_ave = zeros(N,1);
gamma_ave_f = zeros(N,1);
gamma_ave_b = zeros(N,1);
kap_ave = zeros(M,1);
for channel = 1:4
for k=1:K
    a = BPSK(N);
    u = filterinput(a,h); 
    u = u(:,channel);
    u = u(:);

    b = zeros(N,M);
    B = delta.*ones(N,M);
    F = delta.*ones(N,M);
    Delta = zeros(N,M);
    gamma_s = ones(N,M);
    
    for n = 2:N 
        b(n,1) = u(n);
        f(n,1) = u(n);
        F(n,1) = lamda*F(n-1,1)+(u(n))^2;
        gamma_s(n,1) = 1;
        B(n,1) = F(n,1);
        B(1,1) = F(1,1);
        for m = 2:M
            Delta(n,m-1) = Delta(n-1,m-1) + b(n-1,m-1)*f(n,m-1)/(gamma_s(n-1,m-1));
            gamma_f(n,m) = Delta(n,m-1)/B(n-1,m-1);
            gamma_b(n,m) = Delta(n,m-1)/F(n,m-1);

            f(n,m) = f(n,m-1) - gamma_f(n,m)*b(n-1,m-1);
            b(n,m) = b(n-1,m-1)- gamma_b(n,m)*f(n,m-1);
            F(n,m) = F(n,m-1) - gamma_f(n,m)*Delta(n,m-1);
            B(n,m)= B(n-1,m-1)- gamma_b(n,m)*Delta(n,m-1);
            gamma_s(n,m) = gamma_s(n,m-1) - (b(n,m-1))^2/B(n,m-1);
        end
    end
    

    for n = 2:N
        d = a(:);
        d = [zeros(5, 1); d];
        e(n,1)= d(n-1);
        for m = 1:M
            rho(1,m) = 0;
            rho(n,m) = lamda*rho(n-1,m)+ b(n,m)/gamma_s(n,m)*e(n,m);
            kap(n,m) = rho(n,m)/B(n,m);
            e(n,m+1) = e(n,m) - kap(n,m)*b(n,m);
        end

    end
    if channel == 1
        gamma_ave = (gamma_ave + gamma_s(:,11));
        gamma_ave_f = (gamma_ave_f + gamma_f(:,11));
        gamma_ave_b = (gamma_ave_b + gamma_b(:,11));
    end
    alpha = e(:,12)./gamma_s(:,11);
    alphasum(:,k) = alpha.^2;
    MSEE11 = sum(alphasum,2)/K;
end

if channel == 1 %only calculate regression and ref coeffs for channel 1
    for i =1:M
    kap_ave(i) = mean(kap(:,i));
    end
    gamma_ave = gamma_ave/K;
    gamma_ave_f = gamma_ave_f/K;
    gamma_ave_b = gamma_ave_b/K;
end

figure(1)
semilogy(1:N,MSEE11,'LineWidth',1)
grid on
legend('Channel 1','Channel 2','Channel 3','Channel 4')
xlabel('Time (s)');
ylabel('MSSE'); 
title('A Posteriori MSEE');
hold on
end
hold off

figure(2)
plot(1:N,gamma_ave,'LineWidth',2)
legend('Likelihood Parameter')
grid on
xlabel('Time (s)');
ylabel('Averaged \gamma_M'); 
title('Likilihood Parameter vs Time');

figure(3)
plot(1:N,gamma_ave_f,'LineWidth',2.5)
title('Forward Reflection Coeff');
xlabel('Time (s)');
ylabel('Averaged \Gamma_f');
grid on

figure(4)
plot(1:N,gamma_ave_b,'LineWidth',2.5)
title('Backward Reflection Coeff');
xlabel('Time (s)');
ylabel('Averaged \Gamma_b');
grid on

figure(4); 
stem(kap_ave,'color','b','LineWidth',3)
grid on
xlabel('Filter order(M)')
title('Steady-State Values of Regression Coeff');
ylabel('Averaged \kappa_b')