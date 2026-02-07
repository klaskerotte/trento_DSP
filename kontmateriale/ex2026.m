%19:14

%% Exercise 1
h0 = 0.0024; h1 = 0.0363; h2 = -0.2327; h3 = 0.3827;

h = [h0 h1 h2 h3];

% Type 1
hI = [h fliplr(h(1:end-1))];

% Type 2
hII = [h fliplr(h)];

% Type 3
hIII = [h 0 -fliplr(h)];

% Type 4
hIV = [h -fliplr(h)];

%% Exercise 2
[HI w] = freqz(hI,1,linspace(0,pi,512));
HII = freqz(hII,1);
HIII = freqz(hIII,1);
HIV = freqz(hIV,1);

figure(1)

subplot(221)
plot(w/pi,20*log10(abs(HI)/max(abs(HI))))
title('H_I magnitude in dB')
xline(0.12)
xlabel('w/pi')
ylabel('dB')
grid()

subplot(222)
plot(w/pi,20*log10(abs(HII)/max(abs(HII))))
title('H_{II} magnitude in dB')
xline(0.12)
xlabel('w/pi')
ylabel('dB')
grid()

subplot(223)
plot(w/pi,20*log10(abs(HIII)/max(abs(HIII))))
title('H_{III} magnitude in dB')
xline(0.12)
xlabel('w/pi')
ylabel('dB')
grid()

subplot(224)
plot(w/pi,20*log10(abs(HIV)/max(abs(HIV))))
title('H_{IV} magnitude in dB')
xline(0.12)
xlabel('w/pi')
ylabel('dB')
grid()

%% Exercise 3
N = 400;
m = 0:N-1;
x = cos(0.12*pi*m);
Xk = freqz(x,1,w);

y = filter(hI,1,x);

figure(2)
subplot(211)
stem(m,x)
title('Input x[m] plotted vs. m')
xlabel('m')
ylabel('x[m]')
grid()
subplot(212)
stem(m,y)
title('Input y[m] plotted vs. m')
xlabel('m')
ylabel('y[m]')
grid()

%% Exercise 4
wc = 0.12;
delta = 0.1;
O = 200;

[b a] = fir2(O, [0 wc-delta wc+delta 1], [1 0 0 1]);
Hsb = freqz(b,a,w);

figure(3)
plot(w/pi,20*log10(abs(Hsb)/max(abs(Hsb))))
title('Stopbandfilter H_{sb} attenuating x by at least 60dB')
xline(0.12)
yline(-60)
grid()

%% Exercise 5
y2 = filter(b,a,x);

figure(4)
stem(m,y2)
title('Input y2[m] plotted vs. m')
xlabel('m')
ylabel('y2[m]')
grid()

