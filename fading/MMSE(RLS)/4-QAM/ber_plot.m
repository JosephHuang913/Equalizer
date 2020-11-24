close; 
clear;
set(gca, 'fontsize', 12)

load ber_RLS_EQ.log;
semilogy(ber_RLS_EQ(:,1), ber_RLS_EQ(:,2), '-ro',  'LineWidth', 1.6, 'MarkerSIze', 8);
hold on;
grid on;

r=0:1:10;
Pb=0.5.*erfc(sqrt(10.^(r./10)));
semilogy(r, Pb, '-.kx', 'LineWidth', 1.6, 'MarkerSIze', 10);

gamma=0:1:15;
gamma_c=10.^(gamma./10)./3;
mju=sqrt(gamma_c./(1+gamma_c));
Pb=((1/2).*(1-mju)).^3.*(1+3.*((1/2).*(1+mju)).^1+6.*((1/2).*(1+mju)).^2);
semilogy(gamma, Pb, '-.k*',  'LineWidth', 1.6, 'MarkerSIze', 8);

title('BER Performance of Adaptive RLS EQ with 4-QAM in Multipath ISI Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
axis([0 30 1e-4 1]);

legend('4-QAM, 3-path ISI Channel', 'BPSK, AWGN', 'Diversity Order = 3', 1);

%print -djpeg100 MMSE_EQ_RLS_4QAM.jpg;
