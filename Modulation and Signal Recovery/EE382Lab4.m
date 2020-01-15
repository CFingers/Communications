%% EE 382 Lab 4 
%  by. Christopher W. Fingers
%  Professor Morshad
%  Tuesday/Thursday 3:30 pm - 5:45 pm


%% EE 382 Lab 4 Part 1

ts = 1.e-4
t=-0.04:ts:0.04;
Ta = 0.01;
m_sig=triplesinc(t,Ta);
Lfft=length(t); Lfft=2^ceil(log2(Lfft));
M_fre=fftshift(fft(m_sig,Lfft));
freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts);

s_dsb=m_sig.*cos(2*pi*300*t);
Lfft=length(t); Lfft=2^ceil(log2(Lfft)+1);
S_dsb=fftshift(fft(s_dsb,Lfft));
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);

Trange=[-0.03 0.03 -2 2]
figure(1)
subplot(221);td1=plot(t,m_sig);
axis(Trange); set(td1,'Linewidth',2);
xlabel('{\it t} (sec)'); ylabel('{\it m}({\it t})')

subplot(223);td2=plot(t,s_dsb);
axis(Trange); set(td2,'Linewidth',2);
xlabel('{\it t} (sec)'); ylabel('{\it m}({\it t})')

Frange=[-600 600 0 200]
subplot(222);fd1=plot(freqm,abs(M_fre));
axis(Frange); set(fd1,'Linewidth',2);
xlabel('{\it f} (Hz)'); ylabel('{\it m}({\it t})')
subplot(224);fd2=plot(freqs,abs(S_dsb));
axis(Frange); set(fd2,'Linewidth',2);
xlabel('{\it f} (Hz)'); ylabel('{\it m}({\it t})')

%% EE 386 Lab 4 Part 2

ts=1.e-4

t=-0.04:ts:0.04;
Ta=0.01;
m_sig=triangle((t+0.01)/0.01)-triangle((t-0.01)/0.01);
Lm_sig=length(m_sig);
Lfft=length(t); Lfft=2^ceil(log2(Lfft));
M_fre=fftshift(fft(m_sig,Lfft));
freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
B_m=150;    
h=fir1(40,[B_m*ts]);

t=-0.04:ts:0.04;
Ta=0.01;fc=300;
s_dsb=m_sig.*cos(2*pi*fc*t);
Lfft=length(t);  Lfft=2^ceil(log2(Lfft)+1);
S_dsb=fftshift(fft(s_dsb,Lfft));
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);

s_dem=s_dsb.*cos(2*pi*fc*t)*2;
S_dem=fftshift(fft(s_dem,Lfft));

s_rec=filter(h,1,s_dem);
S_rec=fftshift(fft(s_rec,Lfft));

Trange=[-0.025 0.025 -2 2];
figure(1)
subplot(221); td1=plot(t,m_sig);
axis(Trange); set(td1,'Linewidth',1.5);
xlabel('{\it t} (sec)'); ylabel('{\it m}({\it t})')
title('message signal');
subplot(222); td2=plot(t,s_dsb);
axis(Trange); set(td2,'Linewidth',1.5);
xlabel('{\it t} (sec)'); ylabel('{\it t}_{\rm DSB}({\it t})')
title('DSB-SC modulated signal');
subplot(223); td3=plot(t,s_dem);
axis(Trange); set(td3, 'Linewidth',1.5);
xlabel('{\it t} (sec)'); ylabel('{\it e}({\it t})');
title('{\it e}({\it t})');
subplot(224);td4=plot(t,s_rec);
axis(Trange); set(td4,'Linewidth',1.5);
xlabel('{\it t} (sec)'); ylabel('{\it m}_d({\it t})')
title('Recovered signal');

Frange=[-700 700 0 200]
figure(2)
subplot(221); fd1=plot(freqm,abs(M_fre));
axis(Frange); set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('message spectrum');
subplot(222); fd2=plot(freqs,abs(S_dsb));
axis(Frange); set(fd2,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it E}({\it f})');
title('DSB-SC spectrum');

subplot(223);fd3=plot(freqs,abs(S_dem));
axis(Frange);set(fd3,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it E}({\it f})')
title('spectrum of {\it e}({\it t})');
subplot(224); fd4=plot(freqs,abs(S_rec));
axis(Frange); set(fd4,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{]it M}_d({\it f})');
title('recovered spectrum');

Correlation = sum(xcorr(m_sig,s_rec))

%% EE 382 Lab 4 Part 3

ts=1.e-4;
t=-0.04:ts:0.04;
Ta=0.01; fc=300;
m_sig=triangle((t+0.01)/0.01)-triangle((t-0.01)/0.01);
Lm_sig=length(m_sig);
Lfft=length(t);  Lfft=2^ceil(log2(Lfft));
M_fre=fftshift(fft(m_sig,Lfft));
freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
B_m=150;
h=fir1(40,[B_m*ts]);

s_am=(1+m_sig).*cos(2*pi*fc*t);
Lfft=length(t); Lfft=2^ceil(log2(Lfft)+1);
S_am=fftshift(fft(s_am,Lfft));
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);

s_dem=s_am.*(s_am>0);
S_dem=fftshift(fft(s_dem,Lfft));

s_rec=filter(h,1,s_dem);
S_rec=fftshift(fft(s_rec,Lfft));

Trange=[-0.025 0.025 -2 2];
figure(1)
subplot(221); td1=plot(t,m_sig);
axis(Trange); set(td1,'Linewidth',1.5);
xlabel('{\it t} (sec)'); ylabel('{\it m}({\it t})')
title('message signal');
subplot(222); td2=plot(t,s_am);
axis(Trange); set(td2,'Linewidth',1.5);
xlabel('{\it t} (sec)'); ylabel('{\it t}_{\rm DSB}({\it t})')
title('AM modulated signal');

subplot(223); td3=plot(t,s_dem);
axis(Trange); set(td3, 'Linewidth',1.5);
xlabel('{\it t} (sec)'); ylabel('{\it e}({\it t})');
title('rectified signal without local carrier');
subplot(224);td4=plot(t,s_rec);
Trangelow=[-0.025 0.025 -0.5 1];
axis(Trangelow); set(td4,'Linewidth',1.5);
xlabel('{\it t} (sec)'); ylabel('{\it m}_d({\it t})')
title('Detected signal');

Frange=[-700 700 0 200]
figure(2)
subplot(221); fd1=plot(freqm,abs(M_fre));
axis(Frange); set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('message spectrum');
subplot(222); fd2=plot(freqs,abs(S_am));
axis(Frange); set(fd2,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it E}({\it f})');
title('AM spectrum');

subplot(223);fd3=plot(freqs,abs(S_dem));
axis(Frange);set(fd3,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it E}({\it f})')
title('rectified spectrum');
subplot(224); fd4=plot(freqs,abs(S_rec));
axis(Frange); set(fd4,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{]it M}_d({\it f})');
title('recovered spectrum');

% The distortion of the recovered signal is caused by the carrier signals
% frequency.  The carrier signals frequency is too low and the main signal
% acted as a carrier to the carrier.  A carrier signal's frequency needs to
% be much greater than the messanger signal.  
%% EE 382 Lab 4 part 4
ts=1.e-4;
t=-0.04:ts:0.04;
Ta=0.01; fc=300;
m_sig=triangle((t+0.01)/0.01)-triangle((t-0.01)/0.01);
Lm_sig=length(m_sig);
Lfft=length(t);  Lfft=2^ceil(log2(Lfft));
M_fre=fftshift(fft(m_sig,Lfft));
freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
B_m=150;
h=fir1(40,[B_m*ts]);

s_dsb=m_sig.*cos(2*pi*fc*t);
Lfft=length(t); Lfft=2^ceil(log2(Lfft)+1);
S_dsb=fftshift(fft(s_dsb,Lfft));
L_lsb=floor(fc*ts*Lfft);
SSBfilt=ones(1,Lfft);
SSBfilt(Lfft/2-L_lsb+1:Lfft/2+L_lsb)=zeros(1,2*L_lsb);
S_ssb=S_dsb.*SSBfilt;
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
s_ssb=real(ifft(fftshift(S_ssb)));
s_ssb=s_ssb(1:Lm_sig);

s_dem=s_ssb.*cos(2*pi*fc*t)*2;
S_dem=fftshift(fft(s_dem,Lfft));

s_rec=filter(h,1,s_dem);
S_rec=fftshift(fft(s_rec,Lfft));

Trange=[-0.025 0.025 -1 1];
figure(1)
subplot(221); td1=plot(t,m_sig);
axis(Trange); set(td1,'Linewidth',1.5);
xlabel('{\it t} (sec)'); ylabel('{\it m}({\it t})')
title('message signal');
subplot(222); td2=plot(t,s_ssb);
axis(Trange); set(td2,'Linewidth',1.5);
xlabel('{\it t} (sec)'); ylabel('{\it s}_{\rm SSB}({\it t})')
title('SSB-SC modulated signal');

subplot(223); td3=plot(t,s_dem);
axis(Trange); set(td3, 'Linewidth',1.5);
xlabel('{\it t} (sec)'); ylabel('{\it e}({\it t})');
title('after multiplying local carrier');

subplot(224);td4=plot(t,s_rec);
axis(Trange); set(td4,'Linewidth',1.5);
xlabel('{\it t} (sec)'); ylabel('{\it m}_d({\it t})')
title('Recovered signal');

Frange=[-700 700 0 200]
figure(2)
subplot(221); fd1=plot(freqm,abs(M_fre));
axis(Frange); set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('message spectrum');
subplot(222); fd2=plot(freqs,abs(S_ssb));
axis(Frange); set(fd2,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it S}_{rm DSB}({\it f})');
title('AM spectrum');

subplot(223);fd3=plot(freqs,abs(S_dem));
axis(Frange);set(fd3,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it E}({\it f})')
title('detector spectrum');
subplot(224); fd4=plot(freqs,abs(S_rec));
axis(Frange); set(fd4,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{]it M}_d({\it f})');
title('recovered spectrum');


%% EE 382 Lab 4 Part 5 
ts=1.e-4;
t=-0.04:ts:0.04;
Ta=0.01; fc=500;
m_sig=-0.125.*(100.*t+4).^3.*[ustep(t+0.04)-ustep(t+0.02)]+0.125.*(100.*t).^3.*(ustep(t+0.02)-ustep(t-0.02))-0.125.*(100.*t-4).^3.*(ustep(t-0.04)); ;
Lm_sig=length(m_sig);
Lfft=length(t);  Lfft=2^ceil(log2(Lfft));
M_fre=fftshift(fft(m_sig,Lfft));
freqm=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
B_m=150;
h=fir1(40,[B_m*ts]);

s_dsb=m_sig.*cos(2*pi*fc*t);
Lfft=length(t); Lfft=2^ceil(log2(Lfft)+1);
S_dsb=fftshift(fft(s_dsb,Lfft));
L_lsb=floor(fc*ts*Lfft);
SSBfilt=ones(1,Lfft);
SSBfilt(Lfft/2-L_lsb+1:Lfft/2+L_lsb)=zeros(1,2*L_lsb);
S_ssb=S_dsb.*SSBfilt;
freqs=(-Lfft/2:Lfft/2-1)/(Lfft*ts);
s_ssb=real(ifft(fftshift(S_ssb)));
s_ssb=s_ssb(1:Lm_sig);

s_dem=s_ssb.*cos(2*pi*fc*t)*2;
S_dem=fftshift(fft(s_dem,Lfft));

s_rec=filter(h,1,s_dem);
S_rec=fftshift(fft(s_rec,Lfft));

Trange=[-0.025 0.025 -1 1];
figure(1)
subplot(221); td1=plot(t,m_sig);
axis(Trange); set(td1,'Linewidth',1.5);
xlabel('{\it t} (sec)'); ylabel('{\it m}({\it t})')
title('message signal');
subplot(222); td2=plot(t,s_ssb);
axis(Trange); set(td2,'Linewidth',1.5);
xlabel('{\it t} (sec)'); ylabel('{\it s}_{\rm SSB}({\it t})')
title('SSB-SC modulated signal');

subplot(223); td3=plot(t,s_dem);
axis(Trange); set(td3, 'Linewidth',1.5);
xlabel('{\it t} (sec)'); ylabel('{\it e}({\it t})');
title('after multiplying local carrier');

subplot(224);td4=plot(t,s_rec);
axis(Trange); set(td4,'Linewidth',1.5);
xlabel('{\it t} (sec)'); ylabel('{\it m}_d({\it t})')
title('Recovered signal');

Frange=[-700 700 0 200]
figure(2)
subplot(221); fd1=plot(freqm,abs(M_fre));
axis(Frange); set(fd1,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it M}({\it f})');
title('message spectrum');
subplot(222); fd2=plot(freqs,abs(S_ssb));
axis(Frange); set(fd2,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it S}_{rm DSB}({\it f})');
title('AM spectrum');

subplot(223);fd3=plot(freqs,abs(S_dem));
axis(Frange);set(fd3,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{\it E}({\it f})')
title('detector spectrum');
subplot(224); fd4=plot(freqs,abs(S_rec));
axis(Frange); set(fd4,'Linewidth',1.5);
xlabel('{\it f} (Hz)'); ylabel('{]it M}_d({\it f})');
title('recovered spectrum');

