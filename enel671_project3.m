%% ENEL 671 Project3
% Pouyan Keshavarzian
% FALL 2016
%% Effect of Eignevalue Spread
clear all 
close all
clc

M = 10;
delta = 0.01;
N = 600;
K = 100;
h = [0.2194 1.0 0.2194;0.2798 1.0 0.2798;0.3365 1.0 0.3365;0.3887 1.0 0.3887];
lamda = 1;
for k=1:K
      a = BPSK(N);
      u = filterinput(a,h); 
      [e, gamma_s] = RLSL_algorithm(M,N,lamda,delta,a,u(:,1));
      e9 = e(9,:);
      e9 = e9(:);
      ed9(:,k) = e9.^2;
      MSEE9 = sum(ed9,2)/K;
      MSEE9n = MSEE9./gamma_s(9,:)';
end
    figure(1)
    plot(1:N+1,gamma_s(10,:),'LineWidth',2)
    legend('Channel 1','Channel 2','Channel 3','Channel 4')
    grid on
    xlabel('Time (s)');
    ylabel('Mean Squared Error'); 
    title('Likilihood');

    figure(2)
    semilogy(1:N+1,MSEE9,'LineWidth',2)