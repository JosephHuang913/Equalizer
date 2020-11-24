close all;
clear all;

%load MSE120_95.log;
%P=1000;
%i=1:1:P;
%semilogy(MSE120_95(i,1)-4, MSE120_95(i,2), '-b');
%hold on;
%grid on;

%load MSE60_975.log;
%P=1000;
%i=1:1:P;
%semilogy(MSE60_975(i,1)-4, MSE60_975(i,2), '-m');

%load MSE30_99.log;
%P=1000;
%i=1:1:P;
%semilogy(MSE30_99(i,1)-4, MSE30_99(i,2), '-g');

%load MSE0_999.log;
%P=1000;
%i=1:1:P;
%semilogy(MSE0_999(i,1)-4, MSE0_999(i,2), '-r');

load MSE.log;
P=1000;
i=1:1:P;
semilogy(MSE(i,1)-4, MSE(i,2), '-k');

%axis([0 2000.0001 20.0]);
title('Learning Curve of Adaptive MMSE EQ (RLS Direct Form)');
xlabel('Number of Iterations ( i )');
ylabel('Mean Square Error  J( i )');

%legend('W = 3.5', 'W = 3.3', 'W = 3.1', 'W = 2.9');
%print -djpeg100 MSE_BPSK_RLS_fading1.jpg;