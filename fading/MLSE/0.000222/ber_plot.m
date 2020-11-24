close; 
clear;

load ber_MLSE.log;

i=1:1:11;
semilogy(ber_MLSE(i,1), ber_MLSE(i,2), '-bo');
hold on;
grid on;

gamma=0:1:10;
gamma_c=10.^(gamma./10)./3;
mju=sqrt(gamma_c./(1+gamma_c));
Pb=((1/2).*(1-mju)).^3.*(1+3.*((1/2).*(1+mju)).^1+6.*((1/2).*(1+mju)).^2);
semilogy(gamma, Pb, '-k.');

title('BER Performance of MLSE Equalizer in Multipath Fading Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
%axis([0 10 1e-6 1]);

legend('SOVA Algorithm', 'Diversity Order: 3 (Theoretical)', 'Path Weighting Factor: \{0.577,0.577,0.577\}', 'f_d x t = 0.000222', 3);
%print -djpeg100 MLSE_EQ_fading.jpg;
