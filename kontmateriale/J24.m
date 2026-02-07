%% Exercise 1
N = 200;
n = 0:N-1;

x = zeros(1,N);

for r = 1:25;
    x = x + [(n==8*r-1)];
end

%% Exercise 2

[Xk w] = freqz(x,1);

figure(1)
plot(w/pi,abs(Xk))
title('Magnitude of X[w]')
legend('|H_x(e^{jw})|')
xlabel('w/\pi')
ylabel('abs(X(e^{jw}))')
grid()


%% Exercise 3

wc = 0.55;
delta = 0.1;
O = 200;

[b a] = fir2(O,[0 wc wc+delta 1], [0 0 1 1]);

H = freqz(b,a);

figure(2)
plot(w/pi,20*log10(abs(H)/max(abs(H))));
title('Magnitude of filter in dB')
xline(0.5)
yline(-65)
xlabel('w/\pi')
ylabel('dB')
legend('|H(e^{jw})|')

grid()

%% Exercise 4
y = filter(b,a,x);

figure(3)
stem(n,y)
title('Filter output with x[n] as input')
legend('y')
xlabel('n')
ylabel('y[n]')
grid()


%% Exercise 5
delta = 0.05;
ws = 0.52;
wp = ws + delta;
Rp = 3;
Rs = 65;


[O wc] = buttord(wp,ws,Rp,Rs);
[b a] = butter(O,wc,"high");

H = freqz(b,a);

figure(4)
plot(w/pi,20*log10(abs(H)/max(abs(H))));
title('Magnitude of butterworthfilter in dB')
xline(0.5)
yline(-65)
xlabel('w/\pi')
ylabel('dB')
legend('|H_{butter}(e^{jw})|')
grid()

%% Exercise 6
y2 = filter(b,a,x);

figure(5)
stem(n,y2)
title('Filter output with x[n] as input')
legend('y2')
xlabel('n')
ylabel('y_2[n]')
grid()
