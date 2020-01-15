%% Lab 1-2 
%  EE 382 Monday/Wednesday 3:30-5:45pm
%  By Christopher W. Fingers
%  Professor Morshed
%  California State University of Long Beach
%% Question 4i
t = 0:.001:5;
y = exp(-t).*sin(10.*pi.*t).*ustep(t+2);
plot(t,y);
%% Question 4ii
t = 0:.0001:5;
y = exp(2.*t).*cos(20.*t).*rect(t-1);
plot(t,y);
%% Question 4iii
t = 0:.001:5;
y = exp(abs(-t)).*triangle(2.*t);
plot(t,y);
%% Question 5
% Time value is set up to 5 seconds.
t = 0:.001:5;

%First graph shows how the exponential function will normally react with
%not time shift values between a time of 0 and 3.
subplot(5,1,1)
y = exp(-abs(t)./3).*(ustep(t)-ustep(t-3));
plot(t,y);

% The graphs bellow show the effects of various time shifts when applied to
% an exponential decline graph that is set between t = 0 and t = 3.

subplot(5,1,2)
t1 = 4.*t;
y = exp(-abs(t1)./3).*(ustep(t1)-ustep(t1-3));
plot(t1,y);

subplot(5,1,3)
t2 = t-2;
y = exp(-abs(t2)./3).*(ustep(t2)-ustep(t2-3));
plot(t2,y);

subplot(5,1,4)
t3 = t+3;  
y = exp(-abs(t3)./3).*(ustep(t3)-ustep(t3-3));
plot(t3,y);

subplot(5,1,5)
t4 = 3-t;
y = exp(-abs(t4)./3).*(ustep(t4)-ustep(t4-3));
plot(t4,y);
%% Question 6/7

% Set a t value that will show how the x and g graphs will look after a
% certain amount of time has passed.
t= 0:.001:10;

% Set up our x and g functions.
x = ustep(t)-ustep(t-5);
g1= 0.5.*x;
g2 = -1.*x;
g3 = exp(-t./5).*x;
g4 = exp(-t).*x;
g5 = sin(2.*pi.*t).*x;

%Create energy values for x and g values to use in correlation equation
E0 = sum(x.*conj(x))*.001;
E1 = sum(g1.*conj(g1))*.001;
E2 = sum(g2.*conj(g2))*.001;
E3 = sum(g3.*conj(g3))*.001;
E4 = sum(g4.*conj(g4))*.001;
E5 = sum(g5.*conj(g5))*.001;

% Using the Auto Correlation equation, the correlation values for x and g
% values are calculated and displayed in a table at the end of the code.
C0 = sum(x.*conj(x))*.001/(sqrt(E0*E0));
C1 = sum(x.*conj(g1))*.001/sqrt(E1*E1);
C2 = sum(x.*conj(g2))*.001/sqrt(E2*E2);
C3 = sum(x.*conj(g3))*.001/sqrt(E3*E3);
C4 = sum(x.*conj(g4))*.001/sqrt(E4*E4);
C5 = sum(x.*conj(g5))*.001/sqrt(E5*E5);


%Using the correlation equation, the correlation between g1 and x can be
%found.  This gives the signal Cross correlation value.
CC1=sum(x.*conj(g1))*.001/(sqrt(E0*E1));
CC2=sum(x.*conj(g2))*.001/(sqrt(E0*E2));
CC3=sum(x.*conj(g3))*.001/(sqrt(E0*E3));
CC4=sum(x.*conj(g4))*.001/(sqrt(E0*E4));
CC5=sum(x.*conj(g5))*.001/(sqrt(E0*E5));

Correlation = table(C0, C1, C2, C3, C4, C5)
CrossCorrelation = table(CC1, CC2, CC3, CC4, CC5)

%Create plot graphs of t vs x and the g values
subplot(6,1,1)
plot(t,x)
subplot(6,1,2)
plot(t,g1)
subplot(6,1,3)
plot(t,g2)
subplot(6,1,4)
plot(t,g3)
subplot(6,1,5)
plot(t,g4)
subplot(6,1,6)
plot(t,g5)