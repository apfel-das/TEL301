clear all
close all
clc

%A.1

%define values under given assumptions 



T = 10^(-3); 
over = 10;
Ts = T/over; 
A = 4; 

%various roll-off factor values, will be simulated.

[phi1, t1] = srrc_pulse(T,Ts,A,0); 
[phi2, t2] = srrc_pulse(T,Ts,A,0.5); 
[phi3, t3] = srrc_pulse(T,Ts,A,1); 

%open a figure, plot some stuff then hold until done plotting.

figure; 
plot(t1,phi1,'r') 
hold on; 
plot(t2,phi2,'g') 
plot(t3,phi3,'b') 
legend('a = 0', 'a = 0.5', 'a = 1'); 
xlabel('Time (sec)', 'Interpreter', 'latex'); 
ylabel('SRRC pulses', 'Interpreter', 'latex'); 


title('SRRC pulses for various rollof factor values (a $$ \in {0, 0.5 , 1} $$)', 'Interpreter', 'latex')
grid on;
hold off; 



%A.2

%given data
Nf=1024; 
Fs=1/Ts;

%shifted signals

fft_phi1 = fftshift(fft(phi1,Nf)*Ts); 
fft_phi2 = fftshift(fft(phi2,Nf)*Ts); 
fft_phi3 = fftshift(fft(phi3,Nf)*Ts); 

%defining the frequency axis - vector.

F=-Fs/2:Fs/Nf:Fs/2-Fs/Nf; 

%spectrum of shifted signals

spec_phi1 = abs(fft_phi1).^2; 
spec_phi2 = abs(fft_phi2).^2; 
spec_phi3 = abs(fft_phi3).^2;


%plotting the energy spectrum 
figure; 
plot(F,spec_phi1,'r')
title('Common plot of Energy Spectrum', 'Interpreter', 'latex');
hold on;

plot(F,spec_phi2,'g') 
plot(F,spec_phi3,'b') 
legend('a = 0', 'a = 0.5', 'a = 1'); 
xlabel('Frequency (Hz)', 'Interpreter', 'latex'); 
ylabel('${\left|F(F)\right|}^{2}$ ', 'Interpreter', 'latex');
grid on;


%semilogarithmic plot on a new figure

figure; 
semilogy(F,spec_phi1,'r') 
hold on; 
semilogy(F,spec_phi2,'g') 
semilogy(F,spec_phi3,'b') 
legend('a = 0', 'a = 0.5', 'a = 1'); 
xlabel('Frequency (Hz)', 'Interpreter', 'latex'); 
ylabel('${\left|F(F)\right|}^{2}$ ', 'Interpreter', 'latex');
title('Common Semilogarithmic plot of Energy Spectrum', 'Interpreter', 'latex')
grid on;
hold off;


%A.3


%defining the theoritical bandwidth [BW = (1+a)/2T]

B1=1/(2*T);
B2=1.5/(2*T);
B3=1/T;


%cut-off lines 

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
xlabel('Frequency (Hz)', 'Interpreter', 'latex'); 
ylabel('Energy spectrum of SRRC pulses', 'Interpreter', 'latex'); 
title('Depict and define BW of pulses for various a rollof factors', 'Interpreter', 'latex');
hold off;



%B - remove clearance if u have enough horsepower in ur GPU,RAM [I don't].

close all
clear all
clc


%redefining stuff.

T = 10.^(-3);
over = 10;
A = 5;
Ts = T/over;
k=0;
a=[0 0.5 1]; %roll-off factors gathered in a vector



%original signals construction

[phi1, t1] = srrc_pulse(T, Ts, A, a(1));
[phi2, t2] = srrc_pulse(T, Ts, A, a(2));
[phi3, t3] = srrc_pulse(T, Ts, A, a(3));

figure;
plot(t1,phi1,'r')
hold on;
plot(t2,phi2,'g');
hold on;
plot(t3,phi3,'b');
legend('a = 0', 'a=0.5', 'a = 1');
grid on;
title ('Original Signal for various a (roll - of factor) values', 'Interpreter', 'latex');
hold off;


%useful vectors


kVector=0:2*A; 


% initialization

integr1=zeros(1,length(kVector)); 
integr2=zeros(1,length(kVector)); 
integr3=zeros(1,length(kVector));

for j=1:length(a)
    
        for k=0:2*A

            %zero-padding and concatenate   
            phi1_kT=[zeros(1,(1/Ts)*k*T) phi1(1:end-(1/Ts)*k*T)];
            phi2_kT=[zeros(1,(1/Ts)*k*T) phi2(1:end-(1/Ts)*k*T)];
            phi3_kT=[zeros(1,(1/Ts)*k*T) phi3(1:end-(1/Ts)*k*T)];    

            %products
            prod1=phi1.*phi1_kT;  
            prod2=phi2.*phi2_kT;      
            prod3=phi3.*phi3_kT; 

            %intergrals

            integr1(k+1)=sum(prod1)*Ts;    
            integr2(k+1)=sum(prod2)*Ts;        
            integr3(k+1)=sum(prod3)*Ts;     

            %print job for respective rollof factor (a) values and specific
            %k-Values
         

            if(( k ==0) || (k ==1) || ( k == 2) || (k == 3))
                
                
                
              if( j == 1 )  %case a = 0
                
                %draw orig. and delayed signals
                capt=sprintf('Phi(t), Phi(t-kT) [ k = %d  , a = %.1f ] ',k,a(j));  
                figure;  
                plot(t1,phi1,t1,phi1_kT);
                legend('Phi(t)','Phit(t-kT)');
                grid on;
                title(capt, 'Interpreter', 'latex');
                xlabel('Time (t)');
                ylabel('Amplitude of Signals');
                
                %draw product of vectors
                capt=sprintf('Phi(t)Phi(t-kT) [ k = %d  , a = %.1f ] ',k,a(j));  
                figure;  
                plot(t1,prod1);
                grid on;
                title(capt, 'Interpreter', 'latex');
                xlabel('Time (t)');
                ylabel('Amplitude of Product');
             
              elseif( j == 2 ) %case a = 0.5
                  
                %draw orig. and delayed signals  
                capt=sprintf('Phi(t), Phi(t-kT) [ k = %d  , a = %.2f ] ',k,a(j));  
                figure;  
                plot(t2,phi2,t2,phi2_kT);
                legend('Phi(t)','Phit(t-kT)');
                grid on;
                title(capt, 'Interpreter', 'latex');
                xlabel('Time (t)');
                ylabel('Amplitude of Signals');
                
                 %draw product of vectors
                capt=sprintf('Phi(t)Phi(t-kT) [ k = %d  , a = %.1f ] ',k,a(j));  
                figure;  
                plot(t2,prod2);
                grid on;
                title(capt, 'Interpreter', 'latex');
                xlabel('Time (t)');
                ylabel('Amplitude of Product');
                
              else    %case a = 1
                
                 
                  
                %draw orig. and delayed signals  
                capt=sprintf('Phi(t), Phi(t-kT) [ k = %d  , a = %d ] ',k,a(j));
                figure;
                plot(t3,phi3,t3,phi3_kT);
                legend('Phi(t)','Phit(t-kT)');
                grid on;
                title(capt, 'Interpreter', 'latex');
                xlabel('Time (t)');
                ylabel('Amplitude of Signals');
                
                
                %draw product of vectors 
                capt=sprintf('Phi(t)Phi(t-kT) [ k = %d  , a = %.1f ] ',k,a(j));  
                
                figure;  
                plot(t3,prod3);
                grid on;
                title(capt, 'Interpreter', 'latex');
                xlabel('Time (t)');
                ylabel('Amplitude of Product');
                  
              end  


            end
        end
end


%display the integrals, low values will be floored to zero by default

disp('Integral of product (a = 0) , K in [0,10] : '); 
disp(integr1)
disp('Integral of product  (a = 0.5), K in [0,10] : '); 
disp(integr2) 
disp('Integral of product  (a = 1), K in [0,10] : ');
disp(integr3)


%C - comment clearance if available horsepower is present.
clear all 
close all
clc

a = 0.5;
A = 5;
T = 0.1;
over = 10;
Ts = T/over;
N = 100;

%[c1]
%construct N random bits

b = (sign(randn(N, 1)) + 1)/2;


%[c2]
%transform the N-bits Vector created in 2PAM
X = bits_to_2PAM(b);

%simulate X(delta(t))=sum(Xk*delta(t-kT)) [by default]
xDelta=1/Ts*upsample(X,over);

%time Vector
%adds over-1 zeros between symbols respectivelly
tVector=0:Ts:(N+N*(over-1)-1)*Ts;

%drawing xDelta signal
figure;
plot(tVector,xDelta);
grid on;
xlabel('Time(t)');
ylabel('Amplitude');
title('X\delta(t)'); 


%Phi(t) signal pulse generation
[phi,t_phi] = srrc_pulse(T,Ts,A,a);

%time Vector for convolution

convTime =tVector(1)+t_phi(1):Ts:tVector(end)+t_phi(end);


%implement the convolution between the two signals 
x=conv(xDelta,phi)*Ts;

%draw the result

figure;
plot(convTime,x);
xlabel('t');
title('X(t) = X\delta(t) * \Phi(t)');


%constructing second signal,Phi(-t), by inverting time Axis and Values

phi_Inv=phi(end:-1:1);
t_Inv= -t_phi(end:-1:1);

%convoluting inverted signal with Xdelta

z=conv(x,phi_Inv)*Ts;

%new time Vector for the second convolution
convTime_Inv=convTime(1)+t_Inv(1):Ts:convTime(end)+t_Inv(end);

%draw the result
figure;
plot(convTime_Inv,z);
xlabel('t');
title('Z(t) = X(t) * \Phi(-t)');


figure;
plot(convTime_Inv,z);
hold on;
stem([0:N-1]*T,X,'r');
xlabel('t');
legend('Z(t)=X(t)*fi(-t)','X{k}')
