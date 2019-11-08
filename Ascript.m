clear all
close all
clc

%A.1
%given data 
T=10^(-2); 
over=10;
Ts=T/over; 
A=4; 
%roll-off factor a=0,0.5,1 
[phi1, t1] = srrc_pulse(T,Ts,A,0); 
[phi2, t2] = srrc_pulse(T,Ts,A,0.5); 
[phi3, t3] = srrc_pulse(T,Ts,A,1); 
figure; 
plot(t1,phi1,'r') %hold the current graph 
hold on; 
plot(t2,phi2,'g') 
plot(t3,phi3,'b') 
legend('a = 0', 'a = 0.5', 'a = 1'); 
xlabel('Time (sec)'); 
ylabel('SRRC pulses'); 
title('A1 - Depict SRRC pulses for various rollof factor values (a)')
grid on;
hold off; 
%---------------------------------------------------------------
%A.2

%given data 
Nf=1024; 
Fs=1/Ts; 
fft_phi1 = fftshift(fft(phi1,Nf)*Ts); 
fft_phi2 = fftshift(fft(phi2,Nf)*Ts); 
fft_phi3 = fftshift(fft(phi3,Nf)*Ts); 
F=-Fs/2:Fs/Nf:Fs/2-Fs/Nf; 
spec_phi1 = abs(fft_phi1).^2; 
spec_phi2 = abs(fft_phi2).^2; 
spec_phi3 = abs(fft_phi3).^2; 
%plot energy spectrum 
figure; 
plot(F,spec_phi1,'r')
title('A2 - Common plot of Energy Spectrum')
hold on;

plot(F,spec_phi2,'g') 
plot(F,spec_phi3,'b') 
legend('a = 0', 'a = 0.5', 'a = 1'); 
xlabel('Frequency (Hz)'); 
ylabel('|F(F)|^2'); 
grid on;
%semilogarithmic plot 
figure; 
semilogy(F,spec_phi1,'r') 
hold on; 
semilogy(F,spec_phi2,'g') 
semilogy(F,spec_phi3,'b') 
legend('a = 0', 'a = 0.5', 'a = 1'); 
xlabel('Frequency (Hz)'); 
ylabel('|F(F)|^2');
title('A2 - Common Semilogarithmic plot of Energy Spectrum')
grid on;
hold off;
%------------------------------------------------------------
%A.3
B1=1/(2*T);
B2=1.5/(2*T);
B3=1/T;


%given data 

c1=T/(10.^3); 
c2=T/(10.^5); 
%semilogarithmic plot 
figure; 
semilogy(F,spec_phi1,'r') 
hold on;
semilogy(F,spec_phi2,'g')
semilogy(F,spec_phi3,'b') 
semilogy([F(1) F(end)],[c1 c1],'k') 
semilogy([F(1) F(end)],[c2 c2],'k') 
legend('a = 0', 'a = 0.5', 'a = 1'); 
xlabel('Frequency (Hz)'); 
ylabel('Energy spectrum of SRRC pulses'); 
title('A3 - Depict and define BW of pulses for various a rollof factors');
hold off;






