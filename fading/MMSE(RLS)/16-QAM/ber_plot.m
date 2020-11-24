close; 
clear;

load ber_MMSE_16QAM.log;
i=1:1:31;
semilogy(ber_MMSE_16QAM(i,1), ber_MMSE_16QAM(i,2), '-ro');
hold on;
grid on;

load ber_16QAM.log;
i=1:1:16;
semilogy(ber_16QAM(i,1), ber_16QAM(i,2), '-b*');

title('BER Performance of Adaptive MMSE Equalizer with 16-QAM in Multipath Static Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
axis([0 15 1e-6 1]);

legend('16-QAM, Proakis A Channel', '16-QAM, AWGN', 3);
%print -djpeg100 MMSE_EQ_RLS_16QAM.jpg;
