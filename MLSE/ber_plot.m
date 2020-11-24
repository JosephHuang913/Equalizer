close; 
clear;

load ber.log;

i=1:1:11;
semilogy(ber(i,1), ber(i,2), '-bo');
hold on;
grid on;

title('BER Performance of MLSE Equalizer in ISI Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
%axis([0 10 1e-6 1]);

legend('SOVA Algorithm', 'CIR: \{0.407,0.815,0.407\}');
%print -djpeg100 MLSE_EQ.jpg;
