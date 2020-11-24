close; 
clear;

load ber_MMSE120_95.log;
i=1:1:16;
semilogy(ber_MMSE120_95(i,1), ber_MMSE120_95(i,2), '-ro');
hold on;
grid on;

load ber_MMSE60_975.log;
i=1:1:16;
semilogy(ber_MMSE60_975(i,1), ber_MMSE60_975(i,2), '-g^');

load ber_MMSE30_99.log;
i=1:1:16;
semilogy(ber_MMSE30_99(i,1), ber_MMSE30_99(i,2), '-bx');

load ber_MMSE5_999.log;
i=1:1:16;
semilogy(ber_MMSE5_999(i,1), ber_MMSE5_999(i,2), '-ms');

load ber_MLSE_fading.log;
i=1:1:11;
semilogy(ber_MLSE_fading(i,1), ber_MLSE_fading(i,2), '-.k.');

gamma=0:1:30;
gamma_c=10.^(gamma./10)./3;
mju=sqrt(gamma_c./(1+gamma_c));
Pb=((1/2).*(1-mju)).^3.*(1+3.*((1/2).*(1+mju)).^1+6.*((1/2).*(1+mju)).^2);
semilogy(gamma, Pb, '-k.');

title('BER Performance of Adaptive MMSE Equalizer in Multipath Fading Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
axis([0 30 1e-6 1]);

legend('V = 120 km/h, \lambda = 0.95', 'V = 60 km/h, \lambda = 0.975', 'V = 30 km/h, \lambda = 0.99', 'V = 5 km/h, \lambda = 0.999', 'MLSE, \{0.577, 0.577, 0.577\}', 'Diversity Order: 3', 3);
%print -djpeg100 MMSE_EQ_RLS_fading.jpg;
