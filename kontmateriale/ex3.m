%% Exercise 1

h0 = 0.0024; h1 = 0.0363; h2 = -0.2327; h3 = 0.3827;
h = [h0 h1 h2 h3];

figure(1)
% Type 1
h1 = [h fliplr(h(1:end-1))];
subplot(221)
stem([0:6],h1);

% Type 2
h2 = [h fliplr(h)];
subplot(222)
stem([0:7],h2)

% Type 3
h3 = [h 0 -fliplr(h)];
subplot(223)
stem([0:8],h3)

% Type 4
h4 = [h -fliplr(h)];
subplot(224)
stem([0:7],h4)

%% Exercise 2

[H1,w] = freqz(h1,1);
[H2,w] = freqz(h2,1);
[H3,w] = freqz(h3,1);
[H4,w] = freqz(h4,1);

figure(2)
subplot(221)
title('Type I frequency response')
plot(w/pi,20*log10(abs(H1)/max(abs(H1))))
xlabel('w/pi')
yline(-35)
xline(0.12)
ylabel('dB')
grid()
legend('Normalized frequency response type I')

subplot(222)
title('Type II frequency response')
plot(w/pi,20*log10(abs(H2)/max(abs(H2))))
xlabel('w/pi')
ylabel('dB')
yline(-35)
xline(0.12)
grid()
legend('Normalized frequency response type II')

subplot(223)
title('Type III frequency response')
plot(w/pi,20*log10(abs(H3)/max(abs(H3))))
xlabel('w/pi')
ylabel('dB')
yline(-35)
xline(0.12)
grid()
legend('Normalized frequency response type III')

subplot(224)
title('Type IV frequency response')
plot(w/pi,20*log10(abs(H4)/max(abs(H4))))
xlabel('w/pi')
ylabel('dB')
yline(-35)
xline(0.12)
grid()
legend('Normalized frequency response type IV')


%% Exercise 3
N = 400;
m = 0:N-1;

x1 = cos(0.12*pi*m);
[X1 w] = freqz(x1,1);

% figure(3)
% plot(w/pi,abs(X1))
% hold on
% plot(w/pi,20*log10(abs(H1)/max(abs(H1))))
% plot(w/pi,20*log10(abs(H2)/max(abs(H2))))
% plot(w/pi,20*log10(abs(H3)/max(abs(H3))))
% plot(w/pi,20*log10(abs(H4)/max(abs(H4))))
% yline(-35)
% hold off
% legend('X1','H1','H2','H3','H4')

y1 = filter(h1,1,x1);

figure(3)
sgtitle('Output y1[m] and x1[m]')
subplot(211)
stem(m,y1)
xlabel('m')
ylabel('y1[m]')
subplot(212)
stem(m,x1)
xlabel('m')
ylabel('x1[m]')

%% Exercise 4

O = 500;
wc = 0.12;
delta = 0.05;
[b1 a1] = fir1(O,[wc-delta wc+delta],"stop");
[B1,w] = freqz(b1,1);

figure(4)
title('Stopbandfilter')
plot(w/pi,20*log10(abs(B1)/max(abs(B1))))
xline(0.12)
yline(-60)

%% Exercise 5

y2 = filter(b1,1,x1);

figure(5)
title('Output y2[m] using the stopbandfilter')
stem(m,y2)
hold on
stem(m,x1)
hold off
xlabel('m')
ylabel('y2[m]')
grid()


