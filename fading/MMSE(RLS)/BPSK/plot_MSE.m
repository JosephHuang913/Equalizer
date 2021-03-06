close all;
clear all;

load MSE_35.log;
P=1000;
i=1:1:P;
semilogy(MSE_35(i,1)-9, MSE_35(i,2), '-b');
hold on;
grid on;

load MSE_33.log;
P=1000;
i=1:1:P;
semilogy(MSE_33(i,1)-9, MSE_33(i,2), '-m');

load MSE_31.log;
P=1000;
i=1:1:P;
semilogy(MSE_31(i,1)-9, MSE_31(i,2), '-g');

load MSE_29.log;
P=1000;
i=1:1:P;
semilogy(MSE_29(i,1)-9, MSE_29(i,2), '-r');

load MSE.log;
P=1000;
i=1:1:P;
semilogy(MSE(i,1)-4, MSE(i,2), '-k');

axis([0 200 0.001 10.0]);
title('Learning Curve of Adaptive MMSE EQ (RLS Direct Form)');
xlabel('Number of Iterations ( i )');
ylabel('Mean Square Error  J( i )');

legend('W = 3.5', 'W = 3.3', 'W = 3.1', 'W = 2.9');
%print -djpeg100 MSE_BPSK_RLS.jpg;