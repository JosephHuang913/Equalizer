close all;
clear all;

load MSE_0075.log;
P=1000;
i=1:1:P;
semilogy(MSE_0075(i,1)-9, MSE_0075(i,2), '-b');
hold on;
grid on;

load MSE_025.log;
P=1000;
i=1:1:P;
semilogy(MSE_025(i,1)-9, MSE_025(i,2), '-m');

load MSE_075.log;
P=1000;
i=1:1:P;
semilogy(MSE_075(i,1)-9, MSE_075(i,2), '-g');


title('Learning Curve of Adaptive MMSE EQ (LMS Algorithm)');
xlabel('Number of Iterations ( i )');
ylabel('Mean Square Error  J( i )');

legend('\mu = 0.0075', '\mu = 0.025', '\mu = 0.075');
%print -djpeg100 MSE_BPSK_LMS2.jpg;