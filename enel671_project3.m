%% ENEL 671 Project3
% Pouyan Keshavarzian
% FALL 2016
%% Effect of Eignevalue Spread
M = 10;
delta = 0.01;
N = 600;
K = 100;
h = [0.2194 1.0 0.2194;0.2798 1.0 0.2798;0.3365 1.0 0.3365;0.3887 1.0 0.3887];
lamda = 1;
for k=1:K
      a = BPSK(N);
      u = filterinput(a,h); 
      [e] = RLSL_algorithm(M,N,lamda,delta,a,u(:,1));
      e9 = e(9,:);
      e9 = e9(:);
      ed9(:,k) = e9.^2;
      MSEE9 = sum(ed9,2)/K;
end
    figure(1)
    semilogy(1:N+1,MSEE9,'LineWidth',2)
    legend('Channel 1')
    grid on
    xlabel('Time (s)');
    ylabel('Mean Squared Error'); 
    title('Effect of Eigenvalue Spread');
    hold on