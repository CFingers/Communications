% EE 382
% Lab 5

%%
%Part 1
ts=1.e-4;

t=-0.04:ts:0.04;
Ta=0.01;
m_sig=triangle((t+0.01)/Ta)-triangle((t-0.01)/Ta);
Lfft=length(t); Lfft=2^ceil(log2(Lfft));
M_fre=fftshift(fft(m_sig,Lfft));
freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
B_m=100;    %Bandwidth of the signal is B_m Hz
% Design a simple lowpass filter with bandwidth B_m Hz.
h=fir1(80,[B_m*ts]);
%
kf=80;
m_intg=kf*ts*cumsum(m_sig); 
s_fm=cos(2*pi*300*t+m_intg);
s_pm=cos(2*pi*300*t+pi*m_sig);
Lfft=length(t); Lfft=2^ceil(log2(Lfft)+1);
S_fm=fftshift(fft(s_fm,Lfft));
S_pm=fftshift(fft(s_pm,Lfft));
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);

s_fmdem=diff([s_fm(1) s_fm])/ts/kf;
s_fmrec=s_fmdem.*(s_fmdem>0);
s_dec=filter(h,1,s_fmrec);

% Demodulation
% Using an ideal LPF with bandwidth 200 Hz

Trange1=[-0.04 0.04 -1.2 1.2];

figure(1)
subplot(211); m1=plot(t,m_sig);
axis(Trange1); set(m1, 'Linewidth',2);
xlabel('{\it t} (sec)'); ylabel('{\it m}({\it t})');
title('Message signal');
subplot(212); m2=plot(t,s_dec);
set(m2, 'Linewidth',2);
xlabel('{\it t} (sec)'); ylabel('{\it m}_d({\it t})');
title('Demodulated FM Signal');

figure(2)
subplot(211); td1=plot(t,s_fm);
axis(Trange1); set(td1, 'Linewidth',2);
xlabel('{\it t} (sec)'); ylabel('{\it s}_{\rm FM}({\it t})');
title('FM Signal');
subplot(212); td2=plot(t,s_pm);
axis(Trange1); set(td2, 'Linewidth',2);
xlabel('{\it t} (sec)'); ylabel('{\it m}_{\rm PM})({\it t})');
title('PM Signal');

figure(3)
subplot(211); fp1=plot(t,s_fmdem);
set(fp1, 'Linewidth',2);
xlabel('{\it t} (sec)'); ylabel('{\it d s}_{\rm FM}({\it t})/dt');
title('FM Derivative');
subplot(212); fp2=plot(t,s_fmrec);
set(fp2, 'Linewidth',2);
xlabel('{\it t} (sec)');
title('Rectified FM Derivative');

Frange=[-600 600 0 300];
figure(4)
subplot(211); fd1=plot(freqs,abs(S_fm));
axis(Frange); set(fd1, 'Linewidth',2);
xlabel('{\it f} (Hz)'); ylabel('{\it S}_{\rm FM}({\it f})');
title('FM Amplitude Spectrum');
subplot(212); fd2=plot(freqs,abs(S_pm));
axis(Frange); set(fd2,'Linewidth', 2);
xlabel('{\it f} (Hz)'); ylabel('{\it S}_{\rm PM}({\it f})');
title('PM Amplitude Derivative');
%%
%Part 2
ts=1.e-4;

t=-0.04:ts:0.04;
Ta=0.01;
fc=400;
sig_1=(100.*t+4).^3.*(ustep(t+0.04)-ustep(t+0.02));
sig_2=(100.*t).^3.*(ustep(t+0.02)-ustep(t-0.02));
sig_3=(100.*t-4).^3.*(ustep(t-0.02)-ustep(t-0.04));
m=(-0.125*sig_1)+(0.125*sig_2)+(-0.125*sig_3);
m_sig=m;
Lfft=length(t); Lfft=2^ceil(log2(Lfft));
M_fre=fftshift(fft(m_sig,Lfft));
freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
B_m=100;    %Bandwidth of the signal is B_m Hz
% Design a simple lowpass filter with bandwidth B_m Hz.
h=fir1(80,[B_m*ts]);
%
kf=50*pi;
m_intg=kf*ts*cumsum(m_sig);
s_fm=cos(2*pi*fc*t+m_intg);
s_pm=cos(2*pi*fc*t+pi*m_sig);
Lfft=length(t); Lfft=2^ceil(log2(Lfft)+1);
S_fm=fftshift(fft(s_fm,Lfft));
S_pm=fftshift(fft(s_pm,Lfft));
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);

s_fmdem=diff([s_fm(1) s_fm])/ts/kf;
s_fmrec=s_fmdem.*(s_fmdem>0);
s_dec=filter(h,1,s_fmrec);

% Demodulation
% Using an ideal LPF with bandwidth 200 Hz

Trange1=[-0.04 0.04 -1.2 1.2];

figure(1)
subplot(211); m1=plot(t,m_sig);
axis(Trange1); set(m1, 'Linewidth',2);
xlabel('{\it t} (sec)'); ylabel('{\it m}({\it t})');
title('Message signal');
subplot(212); m2=plot(t,s_dec);
set(m2, 'Linewidth',2);
xlabel('{\it t} (sec)'); ylabel('{\it m}_d({\it t})');
title('Demodulated FM Signal');

figure(2)
subplot(211); td1=plot(t,s_fm);
axis(Trange1); set(td1, 'Linewidth',2);
xlabel('{\it t} (sec)'); ylabel('{\it s}_{\rm FM}({\it t})');
title('FM Signal');
subplot(212); td2=plot(t,s_pm);
axis(Trange1); set(td2, 'Linewidth',2);
xlabel('{\it t} (sec)'); ylabel('{\it m}_{\rm PM})({\it t})');
title('PM Signal');

figure(3)
subplot(211); fp1=plot(t,s_fmdem);
set(fp1, 'Linewidth',2);
xlabel('{\it t} (sec)'); ylabel('{\it d s}_{\rm FM}({\it t})/dt');
title('FM Derivative');
subplot(212); fp2=plot(t,s_fmrec);
set(fp2, 'Linewidth',2);
xlabel('{\it t} (sec)');
title('Rectified FM Derivative');

Frange=[-600 600 0 300];
figure(4)
subplot(211); fd1=plot(freqs,abs(S_fm));
axis(Frange); set(fd1, 'Linewidth',2);
xlabel('{\it f} (Hz)'); ylabel('{\it S}_{\rm FM}({\it f})');
title('FM Amplitude Spectrum');
subplot(212); fd2=plot(freqs,abs(S_pm));
axis(Frange); set(fd2,'Linewidth', 2);
xlabel('{\it f} (Hz)'); ylabel('{\it S}_{\rm PM}({\it f})');
title('PM Amplitude Derivative');
%%
%Part 3
ts=1.e-4;

t=-0.04:ts:0.04;
Ta=0.01;
fc=400;
sig_1=(100.*t+4).^3.*(ustep(t+0.04)-ustep(t+0.02));
sig_2=(100.*t).^3.*(ustep(t+0.02)-ustep(t-0.02));
sig_3=(100.*t-4).^3.*(ustep(t-0.02)-ustep(t-0.04));
m=(-0.125*sig_1)+(0.125*sig_2)+(-0.125*sig_3);
m_sig=m;
Lfft=length(t); Lfft=2^ceil(log2(Lfft));
M_fre=fftshift(fft(m_sig,Lfft));
freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
B_m=100;    %Bandwidth of the signal is B_m Hz
% Design a simple lowpass filter with bandwidth B_m Hz.
h=fir1(80,[B_m*ts]);
%
kf=50*pi;
m_intg=kf*ts*cumsum(m_sig);
s_fm=cos(2*pi*fc*t+m_intg);
s_pm=cos(2*pi*fc*t+0.5*pi*m_sig);
Lfft=length(t); Lfft=2^ceil(log2(Lfft)+1);
S_fm=fftshift(fft(s_fm,Lfft));
S_pm=fftshift(fft(s_pm,Lfft));
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);

s_fmdem=diff([s_fm(1) s_fm])/ts/kf;
s_fmrec=s_fmdem.*(s_fmdem>0);
s_dec=filter(h,1,s_fmrec);

% Demodulation
% Using an ideal LPF with bandwidth 200 Hz

Trange1=[-0.04 0.04 -1.2 1.2];

figure(1)
subplot(211); m1=plot(t,m_sig);
axis(Trange1); set(m1, 'Linewidth',2);
xlabel('{\it t} (sec)'); ylabel('{\it m}({\it t})');
title('Message signal');
subplot(212); m2=plot(t,s_dec);
set(m2, 'Linewidth',2);
xlabel('{\it t} (sec)'); ylabel('{\it m}_d({\it t})');
title('Demodulated PM Signal');

figure(2)
subplot(211); td1=plot(t,s_fm);
axis(Trange1); set(td1, 'Linewidth',2);
xlabel('{\it t} (sec)'); ylabel('{\it s}_{\rm FM}({\it t})');
title('FM Signal');
subplot(212); td2=plot(t,s_pm);
axis(Trange1); set(td2, 'Linewidth',2);
xlabel('{\it t} (sec)'); ylabel('{\it m}_{\rm PM})({\it t})');
title('PM Signal');

figure(3)
subplot(211); fp1=plot(t,s_fmdem);
set(fp1, 'Linewidth',2);
xlabel('{\it t} (sec)'); ylabel('{\it d s}_{\rm FM}({\it t})/dt');
title('PM Derivative');
subplot(212); fp2=plot(t,s_fmrec);
set(fp2, 'Linewidth',2);
xlabel('{\it t} (sec)');
title('Rectified PM Derivative');

Frange=[-600 600 0 300];
figure(4)
subplot(211); fd1=plot(freqs,abs(S_fm));
axis(Frange); set(fd1, 'Linewidth',2);
xlabel('{\it f} (Hz)'); ylabel('{\it S}_{\rm FM}({\it f})');
title('FM Amplitude Spectrum');
subplot(212); fd2=plot(freqs,abs(S_pm));
axis(Frange); set(fd2,'Linewidth', 2);
xlabel('{\it f} (Hz)'); ylabel('{\it S}_{\rm PM}({\it f})');
title('PM Amplitude Derivative');

