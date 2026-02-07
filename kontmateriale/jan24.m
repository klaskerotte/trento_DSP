%% Exercise 1
N = 200;
n = 0:N-1;

x = zeros(1,N);
for r = 1:25;
    x = x + [(n == 8*r-1)];
end

figure(1)
stem(n,x)
title('x[n] plotted vs. n')
xlabel('n')
ylabel('x[n]')
legend('x[n]')
grid()

%% Exercise 2

[Xk w] = freqz(x,1,512);

figure(2)
plot(w/pi,abs(Xk))
title('X[w] plotted vs. w')
xlabel('w')
ylabel('X[w]')
legend('X[w]')
grid()

%% Exercise 3
O = 200;
wc = 0.52;
delta = 0.05;
[b a] = firpm(O,[0 wc wc+delta 1], [0 0 1 1]);
H = freqz(b,a,w);

figure(3)
plot(w/pi,20*log10(abs(H)/max(abs(H))));
xline(0.5)
yline(-65)
title('Filter attenuation')
xlabel('w/pi')
ylabel('dB')
grid()

%% Exercise 4

y = filter(b,a,x);
figure(4)
stem(n,y)
title('y[n] plotted vs. n')
xlabel('n')
ylabel('y[n]')
legend('y[n]')
grid()


%% Exercise 5

ws = 0.55;
delta = 0.10;
wp = ws+delta;
Rp = 3;
Rs = 65;

[O wc] = buttord(wp,ws,Rp,Rs);
[b a] = butter(O,wc,"high");

H = freqz(b,a,w);

figure(5)
plot(w/pi,20*log10(abs(H)/max(abs(H))))
yline(-65)
xline(0.5)
title('H[w] plotted vs. w')
xlabel('w')
ylabel('H[w]')
legend('H[w]')
grid()

