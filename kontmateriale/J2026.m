%% Exercise 1

h0 = 0.0024; h1 = 0.0363; h2 = -0.2327; h3 = 0.3827;
h = [h0 h1 h2 h3];

% Type I to IV
hI = [h fliplr(h(1:end-1))];
hII = [h fliplr(h)];
hIII = [h 0 -fliplr(h)];
hIV = [h -fliplr(h)];

figure(1)
subplot(221)
stem([0:length(hI)-1], hI)
subplot(222)
stem([0:length(hII)-1], hII)
subplot(223)
stem([0:length(hIII)-1], hIII)
subplot(224)
stem([0:length(hIV)-1], hIV)

%% Exercise 2

[HI w] = freqz(hI,1);
HII = freqz(hII,1);
HIII = freqz(hIII,1);
HIV = freqz(hIV,1);

figure(2)
subplot(221)
plot(w/pi,20*log10(abs(HI)/max(abs(HI))));
xline(0.12)
title('|H_I(e^{jw})|')
xlabel('w/pi')
ylabel('dB')
grid()
subplot(222)
plot(w/pi,20*log10(abs(HII)/max(abs(HII))));
xline(0.12)
title('|H_{II}(e^{jw})|')
xlabel('w/pi')
ylabel('dB')
grid()
subplot(223)
plot(w/pi,20*log10(abs(HIII)/max(abs(HIII))));
xline(0.12)
title('|H_{III}(e^{jw})|')
xlabel('w/pi')
ylabel('dB')
grid()
subplot(224)
plot(w/pi,20*log10(abs(HIV)/max(abs(HIV))));
xline(0.12)
title('|H_{IV}(e^{jw})|')
xlabel('w/pi')
ylabel('dB')
grid()

%% Exercise 3
N = 400;
m = 0:N-1;

x1 = cos(0.12*pi*m);
% X1(jw) is a spike located at w/pi = 0.12
% This means that only the type I filter is valid
% (checked visually by looking at the dB value for each
% filter type in w/pi = 0.12

y1 = filter(hI,1,x1);

figure(3)
subplot(211)
stem(m,x1)
title('Input x_1[m]')
xlabel('m')
ylabel('x_1[m]')
grid()
subplot(212)
stem(m,y1)
title('Output y_1[m] using filter of type I')
xlabel('m')
ylabel('y_1[m]')
grid()

%% Exercise 4


% Using filterdesigner to create bd
% Settings:
% Responsetype: Bandstop
% Design method: FIR equiripple
% wpass1 = 0.10
% wstop1 = 0.11
% wpass2 = 0.13
% wstop2 = 0.14
% Apass1 = 1
% Astop = 60
% Apass2 = 1

H = freqz(bd,1);

figure(4)
plot(w/pi,20*log10(abs(H)/max(abs(H))));
xline(0.12)
xlabel('w/\pi')
ylabel('dB')
title('Magnitude of filter')
legend('H(e^{jw})|')
grid()

%% Exercise 5
y2 = filter(bd,1,x1);

figure(5)
stem(m,x1)
hold on
stem(m,y2)
hold off
title('filter output y_2[m] using x_1[m] as input')
xlabel('m')
ylabel('y_2[m]')
legend('x_1[m]','y_2[m]')
grid()