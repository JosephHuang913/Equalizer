close;
clear;

ber1 = fopen('ber.log', 'r');
ber2 = fscanf(ber1, '%e');
i=1:1:11;
Eb_No(i) = ber2(2*i-1);
err_rate(i) = ber2(2*i);
semilogy(Eb_No, err_rate, '-b^');
fclose(ber1);

grid on;
title('BER Performance of MAP Equalizer in ISI Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');

legend('BCJR Algorithm', 'CIR: \{0.407,0.815,0.407\}');
%print -djpeg100 MAP_EQ.jpg;
