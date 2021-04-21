
%Abd Elelah Arafah - arafaha - 400197623

%% square wave generator which produces graph of time domain input 
clc
clear all
hold off

i = sqrt(-1);%complex
f0=10000; %fundamental freq of input square wave
T0 = 1/f0;  %period of input wave, 0.1 ms peroid
tstep = 0.005*T0;
tt = 0:tstep:3*T0;
no_sample = 3*T0/tstep + 1; %num of samples for - 3*T0
no_sample1 = T0/tstep + 1; %num of samples for - T0


gp1 = square(2*pi*f0*tt,50); %input - square wave in the per
gp_in = [gp1 gp1(2:no_sample1-1) gp1]; %3 cycles of the triangular wave
figure(1)
Hp1 = plot(tt,gp1);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
title('input - time domain')
xlabel(' Time in Seconds ');
ylabel('Amplitude in Voltage');
pause

      
%% Fourier series representation of signal (Amplitude Spectrum and Phase Spectrum)
      
K=1/(pi);
N=100; %num of harmonics
nvec = -N:N;
c_in = zeros(size(nvec)); % amplitude calculation
for n = nvec
    m = n+N+1;
    q = mod(n,2);
    c_in(m) = 1*i*K/n; 
    
    if (n == 0)
      c_in(m) = 0.0;
    end
    if (q == 0)
      c_in(m) = 0.0;
    end
end
f = nvec*f0; %frequency vector
figure(2)
Hp1=stem(f,abs(c_in));
axis([-8*f0 8*f0 0 max(abs(c_in))])
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
title('magnitude spectrum of input')
pause
figure(3)
Hp1=plot(f,angle(c_in)) %phase spectrum
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
axis([-0.1e4 0.1e4 -pi pi])
title('phase spectrum of input')
pause



%% Designing the 2nd order Butterworth filter

R=133; %resistance found in the lab doc (in ohms)
C=0.1e-6; %capacitance found in the lab doc (in F)
fc=12000 %cutoff freq of filter, found in the lab doc (in Hz)

Hf = 1 ./((1+(1.414*i*f/fc)+((1*i*f/fc).^2))) ;%filter transfer function for 2nd order butterworth
c_out = c_in .* Hf; %Fourier coefficients of the filter output

figure(4)
stem(f,abs(c_in),'r','LineWidth',2);
hold on
stem(f,abs(c_out),'b','LineWidth',2);
hold off
axis([-8*f0 8*f0 0 max(abs(c_in))])
Ha = gca;
set(Ha,'Fontsize',16)
title('magnitude spectrum of filter output and input')
Ha = gca;
set(Ha,'Fontsize',16)
legend('input','output')
pause

%% Construct the output signal from the Cout Fourier coefficients

A = zeros(2*N+1,ceil(no_sample));
for n = nvec
    m=n+N+1;
    A(m,:) = c_out(m) .* exp(1i*2*pi*n*f0*tt);
end
gp_out = sum(A);

figure(5)
Hp1 = plot(tt,real(gp_out),'b',tt,gp1,'r');
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
title('filter input and output-time domain')
set(Ha,'Fontsize',16)
legend('output','input')
%% Transfer function plotting 
%plot of H

figure(6)
plot(((f/2*pi)),(((20*log(abs(Hf))))));
xlabel('Normalized Angular Frequncy ');
ylabel('20*log10(H(f))');
title('Normalized 2nd order butterworth filter');
%% Transfer function plotting multiplying by the cutoff frequncy 
%plot of H

figure(7)
o = ((f)*12000);
plot(((o/2*pi)),(((20*log(abs(Hf))))));
xlabel('Frequncy (KHz');
ylabel('20*log10(H(f))');
title('Amplitude Response of the 2nd order butterworth filter with frequncy axis multiplied by cutoff frequncy');