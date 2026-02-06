%% Exercise 1

N = 4;
wc = 1;
w = linspace(0,2*pi,512);

[b a] = butter(N,wc,'s');
[H w] = freqs(b,a,w);

figure(1)
subplot(211)
plot(w,abs(H))
xlabel('w')
ylabel('|H(jw)|')
subplot(212)
plot(w/(2*pi),20*log10(abs(H)/max(abs(H))))
xlabel('w')
ylabel('|H(jw)| in dB')

%% Exercise 2

Rp = 1;
omegap = 0.4*pi;

Rs = 40;
omegas = 0.5*pi;

wp = tan(omegap/2);
ws = tan(omegas/2);

[N, wa] = buttord(wp,ws,Rp,Rs,'s');
[Ba Aa] = butter(N,wa,'s');

T = 2;
[b, a] = bilinear(Ba, Aa, 1/T);

[HD, w] = freqz(b,a);

figure(2)
plot(w/pi,20*log10(abs(HD)))
axis([0 1 -65 5])
xlabel('w/pi')
ylabel('|H(jw)|')
grid()

%% Exercise 3

Rp = 1;
Rs = 40;
omegap = 0.4*pi;
omegas = 0.5*pi;
wp = tan(omegap/2);
ws = tan(omegas/2);

[N wa] = buttord(wp,ws,Rp,Rs,'s');
[ba aa] = butter(N, wa,'low','s');

T = 2;

[b a] = bilinear(ba,aa,1/T);
[H w] = freqz(b,a);

figure(3)
subplot(211)
plot(w/pi,abs(H))
xlabel('w/\pi')
ylabel('|H(jw)|')
grid()

%% Exercise 4

O = 6;
omega = pi/2;
wc = omega/pi;
Rp = 1;
Rs = 40;

[b1 a1] = butter(O,wc);
[b2 a2] = cheby1(O,Rp,wc);
[b3 a3] = cheby2(O, Rs,wc);
[b4 a4] = ellip(O,Rp,Rs,wc);

[H1 w1] = freqz(b1,a1);
[H2 w2] = freqz(b2,a2);
[H3 w3] = freqz(b3,a3);
[H4 w4] = freqz(b4,a4);

figure(4)
plot(w1/pi,20*log10(abs(H1)))
axis([0 1 -100 5])
hold on
plot(w2/pi,20*log10(abs(H2)))
hold on
plot(w3/pi,20*log10(abs(H3)))
hold on
plot(w4/pi,20*log10(abs(H4)))
hold off

%% Exercise 5

[x Fs] = audioread('speechN.wav');
[Xk w] = freqz(x,1,2^15);
H = freqz(hpfilter,1,2^15);

figure(5)
subplot(311)
plot([0:length(x)-1],x)
subplot(312)
plot(Fs*w(1:500)/(2*pi),abs(Xk(1:500)))
subplot(313)
y = filter(hpfilter,1,x);
plot([0:length(y)-1],y)

%% Exercise 6

[x Fs] = audioread('C:\Users\samuel\Desktop\NuclearLaunchCodes\DSP_trento\kontmateriale\berglN.wav');
[X w] = freqz(x,1,2^15);

figure(6)
subplot(311)
plot(Fs*w/(2*pi),abs(X));
subplot(312)
plot([0:length(x)-1],x)
subplot(313)
y = filter(lpfilter,1,x);
plot([0:length(y)-1],y)

