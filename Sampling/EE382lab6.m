%% EE 382 Lab 6
% By. Christopher W. Fingers,
% Tuesday/Thursday 3:30 pm - 5:45pm

% Part 1
clear;clf;
td=0.002;
t=[0:td:1.];
xsig=sin(2*pi*t)-sin(6*pi*t);
Lsig=length(xsig);

ts = 0.02;
Nfactor=ts/td;

[s_out, sg_out,sgh_out,Delta,SQNR] =sampandquant(xsig,16,td,ts);

Lfft=2^ceil(log2(Lsig)+1);
Fmax=1/(2*td);
Faxis=linspace(-Fmax,Fmax,Lfft);
Xsig=fftshift(fft(xsig,Lfft));
S_out=fftshift(fft(s_out,Lfft));

figure(1);
subplot(311); sfig1a=plot(t,xsig,'k');
hold on; sfig1b=plot(t,s_out(1:Lsig),'b'); hold off;
set(sfig1a,'Linewidth',2); set(sfig1b,'Linewidth',2);
xlabel('time');
title('Signal {\it g}({\if t}) and its uniform samples');
subplot(312); sfig1c=plot(Faxis,abs(Xsig));
xlabel('frequency');
axis([-150 150 0 300])
set(sfig1c,'Linewidth',1); title('Spectrum of {\it g}({\if t})');
subplot(313); sfig1d=plot(Faxis,abs(S_out));
xlabel('Frequency');
axis([-150 150 0 300/Nfactor]);

set(sfig1c,'Linewidth',1); title('Spectrum of {\it g}_T({\if t})');

BW=10;
H_lpf=zeros(1,Lfft);H_lpf(Lfft/2-BW:Lfft/2+BW-1)=1;
S_recv=Nfactor*S_out.*H_lpf;
s_recv=real(ifft(fftshift(S_recv)));
s_recv=s_recv(1:Lsig);

figure(2)
subplot(211); sfig2a=plot(Faxis,abs(S_recv));
xlabel('Frequency');
axis([-150 150 0 300]);
title('Spectrum of ideal filtering (reconstruction)');
subplot(212); sfig2b=plot(t,xsig,'k-.',t,s_recv(1:Lsig),'b');
legend('original signal','reconstructed signal');
xlabel('time');
title('original signal versus ideally reconstructed signal');
set(sfig2b,'Linewidth',2);

ZOH=ones(1,Nfactor);
s_ni=kron(downsample(s_out,Nfactor),ZOH);
S_ni=fftshift(fft(s_ni,Lfft));
S_recv2=S_ni.*H_lpf;
s_recv2=real(ifft(fftshift(S_recv2)));
s_recv2=s_recv2(1:Lsig);

figure(3)
subplot(211); sfig3a=plot(t,xsig,'b',t,s_ni(1:Lsig),'b');
xlabel('time');
title('original signal versus flat top reconstruction');
subplot(212); sfig3b=plot(t,xsig,'b',t,s_recv2(1:Lsig),'b--');
legend('original signal','LPF reconstruction');
xlabel('time');
set(sfig3a,'Linewidth',2); set(sfig3b,'Linewidth',2);
title('original flat top reconstruction after LPF');

%% Part 2

clear;clf;
td=0.002;
t=[0:td:1.];
xsig=sin(2*pi*t)-sin(6*pi*t);
Lsig=length(xsig);
Lfft=2^ceil(log2(Lsig)+1);
Xsig=fftshift(fft(xsig,Lfft));
Fmax=1/(2*td);
Faxis=linspace(-Fmax,Fmax,Lfft);
ts=0.02;
Nfact=ts/td;

[s_out,sq_out,sqh_out1,Delta,SQNR]=sampandquant(xsig,16,td,ts);

figure(1);
subplot(211); sfig1=plot(t,xsig,'k',t,sqh_out1(1:Lsig),'b');
set(sfig1,'Linewidth',2);
title('Signal {\it g}({\it t}) and its 16 level PCM signal')
xlabel('time');
[s_out,sq_out,sqh_out2,Delta,SQNR]=sampandquant(xsig,4,td,ts);

subplot(212);sfig2=plot(t,xsig,'k',t,sqh_out2(1:Lsig),'b');
set(sfig2,'Linewidth',2);
title('Signal {\it g}({\it t}) and its 4 level PCM signal')
xlabel('time');

Lfft=2^ceil(log2(Lsig)+1);
Fmax=1/(2*td);
Faxis=linspace(-Fmax,Fmax,Lfft);
SQH1=fftshift(fft(sqh_out1,Lfft));
SQH2=fftshift(fft(sqh_out2,Lfft));

BW=10;
H_lpf=zeros(1,Lfft);H_lpf(Lfft/2-BW:Lfft/2+BW-1);
S1_recv=SQH1.*H_lpf;
s_recv1=real(ifft(fftshift(S1_recv)));
s_recv1=s_recv1(1:Lsig);
S2_recv=SQH2.*H_lpf;
s_recv2=real(ifft(fftshift(S2_recv)));
s_recv2=s_recv2(1:Lsig);

figure(2)
subplot(211);sfig3=plot(t,xsig,'b-',t,s_recv1,'b-.');
legend('original','recovered');
set(sfig3,'Linewidth',2);
title('Signal {\it g}({\it t}) and filtered 16-level PCM signal')
xlabel('time');
subplot(212);sfig4=plot(t,xsig,'b-',t,s_recv2,'b-.');
legend('original','recovered');
set(sfig4,'Linewidth',2);
title('Signal {\it g}({\it t}) and filtered 4-level PCM signal')
xlabel('time');

%% Part 3

clear;clf;
td=0.002;
t=[0:td:1.];
xsig=sin(2*pi*t)-sin(6*pi*t);
Lsig=length(xsig);
ts=0.02;
Nfact=ts/td;
Delta1=0.2;
s_DMout1=deltamod(xsig,Delta1,td,ts);
figure(1);
subplot(311); sfig1=plot(t,xsig,'k',t,s_DMout1(1:Lsig),'b');
set(sfig1,'Linewidth',2);
title('Signal {\it g}({\it t}) and DM signal')
xlabel('time');axis([0 1 -2.2 2.2]);

Delta2=Delta1*2;
s_DMout2=deltamod(xsig,Delta2,td,ts);
subplot(312); sfig2=plot(t,xsig,'k',t,s_DMout2(1:Lsig),'b');
set(sfig2,'Linewidth',2);
title('Signal {\it g}({\it t}) and DM signal')
xlabel('time');axis([0 1 -2.2 2.2]);

Delta3=Delta2*2;
s_DMout3=deltamod(xsig,Delta3,td,ts);
subplot(313); sfig3=plot(t,xsig,'k',t,s_DMout3(1:Lsig),'b');
set(sfig3,'Linewidth',2);
title('Signal {\it g}({\it t}) and DM signal')
xlabel('time');axis([0 1 -2.2 2.2]);