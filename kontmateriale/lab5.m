%% Exercise 1

h1 = [1 2 3 -2 5 -2 3 2 1];
h2 = [1 2 3 -2 -2 3 2 1];
h3 = [1 2 3 -2 0 2 -3 -2 -1];
h4 = [1 2 3 -2 2 -3 -2 -1];

w = linspace(0,2*pi,512);
ws = linspace(-pi,pi,512);
wn = w/pi;

H1 = freqz(h1,1,'whole');
H2 = freqz(h2,1,'whole');
H3 = freqz(h3,1,'whole');
H4 = freqz(h4,1,'whole');

figure(1)
sgtitle('Amplitude and phase responses')
subplot(211)
title('Amplitude responses of H_i(jw)')
plot(wn,abs(H1))
xlabel('w/pi')
ylabel('|H_i(jw)|')
hold on
plot(wn,abs(H2))
plot(wn,abs(H3))
plot(wn,abs(H4))
hold off
legend('H1','H2','H3','H4')

subplot(212)
title('Phase responses of H_i(jw)')
plot(wn,angle(H1))
xlabel('w/pi')
ylabel('angle(H_i(jw))')
hold on
plot(wn,angle(H2))
plot(wn,angle(H3))
plot(wn,angle(H4))
legend('H1','H2','H3','H4')
hold off

figure(2)
sgtitle('Pole zero plots')
subplot(221)
zplane(h1)
legend('h1')
subplot(222)
zplane(h2)
legend('h2')
subplot(223)
zplane(h3)
legend('h3')
subplot(224)
zplane(h4)
legend('h4')

%% Exercise 2
h = [-4 1 -1 -2 5 6 6 5 -2 -1 1 -4];
H = freqz(h, 1,"whole");

figure(3)
sgtitle('Phase and amplitude response of H(jw)')
subplot(311)
plot(wn,abs(H))
xlabel('w/pi')
ylabel('|H(jw)|')
subplot(312)
plot(wn,angle(H))
xlabel('w/pi')
ylabel('angle(H)')
subplot(313)
zplane(h)

%% Exercise 3
h0 = 2; h1 = 1.5; h2 = -3.2; h3 = -5.2; h4 = 6.4;
h = [h0 h1 h2 h3 h4];

% Type 1
h1 = [h fliplr(h(1:end-1))];

% Type 2
h2 = [h fliplr(h)];

% Type 3
h3 = [h 0 -fliplr(h)];

% Type 4
h4 = [h -fliplr(h)];

figure(4)
subplot(221)
zplane(h1)
legend('h1')
subplot(222)
zplane(h2)
legend('h2')
subplot(223)
zplane(h3)
legend('h3')
subplot(224)
zplane(h4)
legend('h4')

%% Exercise 4
omega_c = pi/2;
wc = omega_c/pi;

N = 50;
n = -N/2:N/2;
h = wc*sinc(wc*n);
H = freqz(h,1,ws);

figure(5)
sgtitle('Truncated filter')
subplot(211)
title('h[n]')
stem(n+N/2,h)
grid()
xlabel('n')
ylabel('h[n]')
subplot(212)
plot(wn,abs(H))
xlabel('w/pi')
ylabel('|H(jw)|')


%% Exercise 5
N = 50;
n = -N/2:N/2;
wc = 1/2;

h1 = ones(1,N+1)';
h2 = hamming(N+1);
h3 = blackman(N+1);

Hones =freqz(h1,1,ws);
Hhamm =freqz(h2,1,ws);
Hblack=freqz(h3,1,ws);

hlp = wc*sinc(wc*n);

figure(6)
subplot(211)
plot(wn,abs(Hones))
hold on
plot(wn,abs(Hhamm))
plot(wn,abs(Hblack))
hold off
legend('Hones','Hhamm','Hblack')
subplot(212)
plot(wn,20*log10(abs(Hones)))
hold on
plot(wn,20*log10(abs(Hhamm)))
plot(wn,20*log10(abs(Hblack)))
hold off
legend('Hones','Hhamm','Hblack')


%% Exercise 6
%window the ideal lp:
h1t = h1.'.*hlp;
h2t = h2.'.*hlp;
h3t = h3.'.*hlp;

H1 = freqz(h1t,1,ws);
H2 = freqz(h2t,1,ws);
H3 = freqz(h3t,1,ws);

figure(7)
subplot(211)
plot(wn,abs(H1))
hold on
plot(wn,abs(H2))
plot(wn,abs(H3))
hold off
legend('H1','H2','H3')
subplot(212)
plot(wn,20*log10(abs(H1)))
hold on
plot(wn,20*log10(abs(H2)))
plot(wn,20*log10(abs(H3)))
hold off
legend('H1','H2','H3')

%% Exercise 7

N = 50;
n = -N/2:N/2;
wc = 1/2;
wc2 = 1/3;

hlp = wc*sinc(wc*n);
hlp2 = wc2*sinc(wc2*n);

Hlp = freqz(hlp,1,ws);
Hlp2 = freqz(hlp2,1,ws);

% hp filter
hhp = [n==0]-hlp;
Hhp = freqz(hhp,1,ws);

% bp filter
hbp = hlp-hlp2;
Hbp = freqz(hbp,1,ws);

% bs filter 1 - Hbp
hbs = [n==0] - hbp;
Hbs = freqz(hbs,1,ws);


figure(8)
sgtitle('Filters obtained by complementary and difference of filters')
subplot(221)
title('Low pass')
plot(ws/pi,abs(Hlp))
legend('Hlp')

subplot(222)
title('High pass = delta - hlp')
plot(ws/pi,abs(Hhp))
legend('Hhp')

subplot(223)
title('Bandpass = hlp1 - hlp2')
plot(ws/pi,abs(Hbp))
legend('Hbp')

subplot(224)
title('Bandstop = delta - hbp')
plot(ws/pi,abs(Hbs))
legend('Hbs')

%% Exercise 8

omegac = pi/2;
wc = omegac/pi;
O = 50;

[b1 a1] = fir1(O,wc);
[b2 a2] = fir1(O,wc,blackman(O+1));

H1 = freqz(b1,a1);
H2 = freqz(b2,a2);

figure(9)
title('Lowpass filters obtained by fir1')
plot(wn,20*log10(abs(H1)))
hold on
plot(wn,20*log10(abs(H2)))
hold off
xlabel('w/pi')
ylabel('|H_i(jw)|')
grid()
legend('|H1(jw)|','|H2(jw)|')

%% Exercise 9

omegac = pi/2;
wc = omegac/pi;
O = 50;
delta = 0.05;

[b a] = fir2(O,[0 wc-delta wc+delta 1], [1 1 0 0]);
H = freqz(b,a);

figure(10)
title('Lowpass filter obtained from fir2')
plot(wn/2,abs(H))
grid()
xlabel('w/pi')
ylabel('|H(jw)|')

%% Exercise 10
O = 35;
f = [0 0.2 0.3 0.7 0.8 1];
m = [0 0 1 1 0 0];
[b a] = firpm(O,f,m);

H = freqz(b,a);
figure(11)
subplot(211)
sgtitle('Passband FIR filter obtained by firpm')
plot(wn/2,abs(H))
xlabel('w/pi')
ylabel('|H(jw)|')
subplot(212)
plot(wn/2,20*log10(abs(H)))
xlabel('w/pi')
ylabel('|H(jw)| in dB')

%% Exercise 11

omega_c = pi/3;
wc = omega_c/pi;
O = 10;
delta = 0.01;

[b1 a1] = fir1(O,wc,"low");
[b2 a2] = fir2(O,[0 wc-delta wc+delta 1],[1 1 0 0]);
[b3 a3] = firpm(O,[0 wc-delta wc+delta 1],[1 1 0 0]);


[H1,w] = freqz(b1,a1);
H2 = freqz(b2,a2);
H3 = freqz(b3,a3);

figure(12)
sgtitle('Amplitude response for fir1 fir2 and firpm fiters')
subplot(211)
plot(w/pi,abs(H1)/max(abs(H1)))
hold on
plot(w/pi,abs(H2)/max(abs(H2)))
hold on
plot(w/pi,abs(H3)/max(abs(H3)))
hold off
xlabel('w/pi')
ylabel('Amplitude')
legend('H1','H2','H3')
subplot(212)
plot(w/pi,20*log10(abs(H1)/max(abs(H1))))
hold on
plot(w/pi,20*log10(abs(H2)/max(abs(H2))))
hold on
plot(w/pi,20*log10(abs(H3)/max(abs(H3))))
hold off
xlabel('w/pi')
ylabel('Amplitude in dB')
legend('H1','H2','H3')


