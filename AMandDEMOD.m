%% INTRODUCTION
% Amplitude Modulation
% Abd Elelah Arafah - arafaha - 400197623
%% Sampling Process
clear
hold off
format long e

%Sample Calculation:

NUM = 2^16; %# of fourier samples
Sample = 400000; %Sample rate in Hz
Step = 1/Sample;
Max_t = NUM*Step/2;
Min_t = -Max_t;
tt = Min_t:Step:Max_t-Step;


%frequency sampling:

Max_f = Sample/2; 
Min_f = -Max_f;
Step_f = (Max_f-Min_f)/NUM;
f = Min_f:Step_f:Max_f-Step_f; %Frequency

%% Modulation Process

%The Carrier Signal wave:

fc = 20000;
Ac = 1;
ct = Ac*cos(2*pi*fc*tt);

%The Message Signal wave:

Am = -2;
fm = 2000;
mt = Am*sinc(tt*fm);%sinc function which is equal to sinc(x) = sin(pi*x)/(pi*x)

%Maximum of the message signal:

maxmt = max(abs(mt)); %Retrieves the maximum absolute value of the message signal

%For 200% modulation

ka = 2/maxmt; %2 for 200% modulation

%S(t) signal for Amplitude Modulation:

st = (1+ka*mt).*ct;


% The plot for Carrier Signal in Time Domain

figure(1)
Plot = plot(tt,ct,'b');
set(Plot,'LineWidth',2)
A = gca;
set(A,'Fontsize',14)
Axis = xlabel('Time (in sec) ');
set(Axis,'FontWeight','bold','Fontsize',14)
Axis = ylabel('Carrier c(t)  (in V)');
set(Axis,'FontWeight','bold','Fontsize',14)
title('Carrier Signal in Time domain');
axis([-1e-3 1e-3 min(ct) max(ct)])
pause(1)


% The plot for Message Signal in Time Domain

figure(2)
Plot = plot(tt,mt,'m');
set(Plot,'LineWidth',2)
A = gca;
set(A,'Fontsize',14)
Axis = xlabel('Time (sec) ');
set(Axis,'FontWeight','bold','Fontsize',14)
Axis = ylabel('message  m(t)  (Volt)');
set(Axis,'FontWeight','bold','Fontsize',14)
title('message signal : Time domain');
axis([-2e-3 2e-3 min(mt) max(mt)])
pause(1)

% The plot for Modulated Signal in Time Domain

figure(3)
Plot = plot(tt,st);
set(Plot,'LineWidth',2)
A = gca;
set(A,'Fontsize',16)
Axis = xlabel('Time (sec) ');
set(Axis,'FontWeight','bold','Fontsize',16)
Axis = ylabel('s(t)  (Volt)');
set(Axis,'FontWeight','bold','Fontsize',16)
title('modulated wave : Time domain');
axis([-2e-3 2e-3 min(st) max(st)])
pause(1)

% The plot for Message Signal in Freq Domain

Mf = abs(fftshift(fft(mt)))*Step;
modulated_freq1 = Sample*(-length(st)/2:(length(st)/2)-1)/length(st);

figure(4)
Plot = plot(modulated_freq1,abs(Mf),'m');
set(Plot,'LineWidth',2)
A = gca;
set(A,'Fontsize',14)
Axis = xlabel('Frequency (in Hz) ');
set(Axis,'FontWeight','bold','Fontsize',14)
Axis = ylabel('|M(f)|');
set(Axis,'FontWeight','bold','Fontsize',14)
title('Magnitude Spectrum of the message signal');
axis([-5e3 5e3 0 max(abs(Mf))])


% The plot for Modulated Signal in Freq Domain



Sf = abs(fftshift(fft(st)))*Step;
modulated_freq = Sample*(-length(st)/2:(length(st)/2)-1)/length(st);

figure(5)
Plot = plot(modulated_freq,abs(Sf));
set(Plot,'LineWidth',2)
A = gca;
set(A,'Fontsize',16)
Axis = xlabel('Frequency (Hz) ');
set(Axis,'FontWeight','bold','Fontsize',16)
Axis = ylabel('|S(f)|');
set(Axis,'FontWeight','bold','Fontsize',16)
title('Spectrum of the modulated wave');
axis([-25e3 25e3 0 max(abs(Sf))])


%% Demodulation Process

% Chosen R_LC value for 200% modulation part:

RC = 0.5*((1/fm)+(1/fc));

% Envelope detector code

Output_Env = st;
Val_Compare = 1;


for Time = tt
    if(Val_Compare > 1)
        if(Output_Env(Val_Compare-1) > st(Val_Compare))
            yt0 = Output_Env(Val_Compare-1);
            %time when C starts discharging
            tc = tt(Val_Compare-1);
            Output_Env(Val_Compare) = yt0*exp(-(Time-tc)/RC);
        end
    end
    Val_Compare = Val_Compare+1;
end

Output_Env(1) = Output_Env(2);

% The plot for the Envelope detector output in Time Domain

figure(6)
Plot = plot(tt,Output_Env);
set(Plot,'LineWidth',2)
A = gca;
set(A,'Fontsize',16)
Axis = xlabel('Time (sec) ');
set(Axis,'FontWeight','bold','Fontsize',16)
Axis = ylabel('y(t)  (Volt)');
set(Axis,'FontWeight','bold','Fontsize',16)
title('After the envelope detector');
axis([-2e-3 2e-3 0 max(Output_Env)])
pause(1)


%The plot for the DC Removal output in Time Domain
figure(7)
Output_Env_Final = (Output_Env - 1)/ka;
Plot = plot(tt,Output_Env_Final,'r',tt,mt,'k');
legend('after DC removal','message signal')
set(Plot,'LineWidth',2)
A = gca;
set(A,'Fontsize',16)
Axis = xlabel('Time (sec) ');
set(Axis,'FontWeight','bold','Fontsize',16)
Axis = ylabel('y1(t)  (Volt)');
set(Axis,'FontWeight','bold','Fontsize',16)
title('After the DC removal');
axis([-2e-3 2e-3 min(mt) max(mt)])


CutOff_f = 1.1*fm;
[b,a] = butter(6,CutOff_f/(Sample/2));
message = filter(b,a,Output_Env_Final);
figure(8)
Plot = plot(tt,message,'b',tt,mt,'r');
legend('LPF Output','Message signal')
set(Plot,'LineWidth',3)
Ha = gca;
set(Ha,'Fontsize',14)
Axis=xlabel('Time (sec)');
set(Axis,'FontWeight','bold','Fontsize',14)
Axis=ylabel('m1(t)  (Volt)');
set(Axis,'FontWeight','bold','Fontsize',14)
title('LPF Output');
axis([-2e-3 2e-3 min(mt) max(mt)])