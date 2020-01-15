%% Lab 3 Fourier Series
% By Chrisotpher W. Fingers
% EE 382 Tuesday/Thursday 3:30 pm - 5:45 Pm

%% Part 1 Compact Fourier

% Define variables and vector
count =1;
y =0;
n = 3;
cn = [ 0.504 .244 .125 .084 .063 .0504 .042 .036];
t= 0:0.001:pi;

% Create while loop to add counted values to y vector
while count ~=n
    x = (2/(sqrt(1+16*(count)^2))*cos(2.*count.*t - atan(4*count)));
    y = y +x;   
    count = count+1;
end

% Add A0 value to compact fourier equation
 y = cn(1) + cn(1).*y;
 
hold on
plot(t,y)
xlabel('Time(sec)');
ylabel('Amplitude');
title ('Amplitude vs Time 2, 5 and 7 Harmonics');
 
y = 0;
n = 6;
count =1;
while count ~=n
    x = 2/(sqrt(1+16*(count)^2))*cos(2.*count.*t - atan(4*count));
    y = y +x;
    count = count+1;
end
    y = cn(1) + cn(1).*y;

hold on
plot(t,y)

y =0;
n = 8;
count =1;
while count ~=n
    x = (2/(sqrt(1+16*(count)^2))*cos(2.*count.*t - atan(4*count)));
    y = y +x;
    count = count+1;
    phi = atand(4*count);
end
 y = cn(1) + cn(1).*y;

hold on
plot(t,y)

% As the different harmonics are plotted against eachother the shape of the
% graph becomes more defined to how it should look, an exponential decay.

%% Part B

clc,clear;
% Define variables
t = 0:0.0001:pi;
A0 = 1/2;
w= 0;
n = 1;
count =1;
% Create while loop for the compact fourier equation.
while count ~= n+1
    
    An = (2./(count.*pi)).*sin(count.*pi./2);
    Omega = An.*cos(t.*count);
    w = w+Omega;
    count = count+1;
end

hold on
Endw = 1/2 + w;
plot(t,Endw)
xlabel('time');
ylabel('Amplitude');
title('Amplitude vs Time for Harmonics 1, 3,6,50 and 120');

w = 0;
Endw = 0;
n = 3;
count =1;
while count ~= n+1
    
    An = (2./(count.*pi)).*sin(count.*pi./2);
    Omega = An.*cos(t.*count);
    w = w+Omega;
    count = count+1;
end

hold on
Endw = 1/2 + w;
plot(t,Endw)

w = 0;
Endw = 0;
n = 6;
count =1;
while count ~= n+1
    
    An = (2./(count.*pi)).*sin(count.*pi./2);
    Omega = An.*cos(t.*count);
    w = w+Omega;
    count = count+1;
end

hold on
Endw = 1/2 + w;
plot(t,Endw)

w = 0;
Endw = 0;
n = 20;
count =1;
while count ~= n+1
    
    An = (2./(count.*pi)).*sin(count.*pi./2);
    Omega = An.*cos(t.*count);
    w = w+Omega;
    count = count+1;
end

hold on
Endw = 1/2 + w;
plot(t,Endw)

w = 0;
Endw = 0;
n = 50;
count =1;
while count ~= n+1
    
    An = (2./(count.*pi)).*sin(count.*pi./2);
    Omega = An.*cos(t.*count);
    w = w+Omega;
    count = count+1;
end

hold on
Endw = 1/2 + w;
plot(t,Endw)

w = 0;
Endw = 0;
n = 120;
count =1;
while count ~= n+1
    
    An = (2./(count.*pi)).*sin(count.*pi./2);
    Omega = An.*cos(t.*count);
    w = w+Omega;
    count = count+1;
end

hold on
Endw = 1/2 + w;
plot(t,Endw)

%As the more and more Harmonics are increased the shape of the square wave
%begins to take form.

%%
clc,clear
n = 3;
t = 0:1/12000:0.05;
A=1;
w = 2*pi/(1/60);
count = 1;
xx = 0;
End = 0;

% Use the fourier transform equation for the half wave rectified sine.
for count = 1:n;
   xx = xx + cos(2*count*w*t) /(4*count^2 -1);
end
    End = A/pi + (A/2)*sin(w*t) - (2*A/pi)*xx;
    
    figure
    plot(t,End)
    xlabel('time');
    ylabel('Amplitude');
    title('Amplitude vs Time for a 3rd harmonic half wave rectified sine wave');
    
n = 5;
t = 0:1/12000:0.05;
A=1;
w = 2*pi/(1/60);
count = 1;
xx = 0;
End = 0;

% Use the fourier transform equation for the half wave rectified sine.
for count = 1:n;
   xx = xx + cos(2*count*w*t) /(4*count^2 -1);
end
    End = A/pi + (A/2)*sin(w*t) - (2*A/pi)*xx;
    
    figure
    plot(t,End)
    xlabel('time');
    ylabel('Amplitude');
    title('Amplitude vs Time for a 5th harmonic half wave rectified sine wave');
    
n = 7;
t = 0:1/12000:0.05;
A=1;
w = 2*pi/(1/60);
count = 1;
xx = 0;
End = 0;

% Use the fourier transform equation for the half wave rectified sine.
for count = 1:n;
   xx = xx + cos(2*count*w*t) /(4*count^2 -1);
end
    End = A/pi + (A/2)*sin(w*t) - (2*A/pi)*xx;
    
    figure
    plot(t,End)
    xlabel('time');
    ylabel('Amplitude');
    title('Amplitude vs Time for a 7th harmonic half wave rectified sine wave');
   
 
%% Part 4 Dft/FFt

clc,clear;
Ts=1/64;
T0=4;
N0=T0/Ts;
t=0:Ts:Ts*(N0-1);
t=t';
g=Ts*exp(-2*t);
g(1)=Ts*0.5;
G=fft(g);
[Gp,Gm]=cart2pol(real(G),imag(G));
k=0:N0-1;
k=k';
w=2*pi*k/T0;
subplot(211), stem(w(1:32),Gm(1:32));
subplot(212),stem(w(1:32),Gp(1:32))

%% Part 5 FFt

clc,clear
B=4;
f0=1/4;
Ts=1/(2*B);
T0=1/f0;
N0=T0/Ts;
k=0:N0-1;
k=k';
for m =1:length(k)
    if k(m)>=0 & k(m)<=3, gk(m)=1; end
    if k(m)==4 & k(m)==28 gk(m)=0.5; end
    if k(m)>=5 & k(m)<=27, gk(m)=0; end
    if k(m)>=29 & k(m)<=31, gk(m)=1; end
end
gk=gk';
Gr=fft(gk);
subplot(211),stem(k,gk)
subplot(212),stem(k,Gr)


%% Part 7 Filtering

clc,clear
B=4;
f0=1/4;
Ts=1/(2*B);
T0=1/f0;
N0=T0/Ts;
k=0:N0-1;
k=k';
for m =1:length(k)
    if k(m)>=0 & k(m)<=3, gk(m)=1; end
    if k(m)==4 & k(m)==28 gk(m)=0.5; end
    if k(m)>=5 & k(m)<=27, gk(m)=0; end
    if k(m)>=29 & k(m)<=31, gk(m)=1; end
end
gk=gk';
Gr=fft(gk);
subplot(211),stem(k,gk)
subplot(212),stem(k,Gr)
r=0:32; r=r';
for m = 1:length(r)
    if r(m) >=0 & r(m)<=7, Hr(m)=1; end
    if r(m) >=25 & r(m) <=31, Hr(m) = 1; end
    if r(m) >= 9 & r(m) <=23, Hr(m) =0; end
    if r(m) ==8 & r(m) == 24, Hr(m) = 0.5; end
end
Hr = Hr';
Yr=Gr.*Hr;
yk = ifft(Yr);
clf, stem(k,yk)