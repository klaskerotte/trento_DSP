%% Exercise 1
N = 200;
n = 0:N-1;

x = zeros(1,N);

for r = 1:25;
    x = x + 1*[(n == 8*r -1)];
end

figure(1)
sgtitle('x[n] sequence plotted over n')
stem(n,x)
xlabel('n')
ylabel('x[n]')
legend('x[n]')
grid()

%% Exercise 2

[Xk w] = freqz(x,1,linspace(0,pi,1024));

figure(2)
sgtitle('Magnitude of X[k] plotted over normalized frequency')
plot(w/pi, abs(Xk))
xlabel('w/pi')
ylabel('X[k]')
legend('X[k]')
grid()

%% Exercise 3

delta = 0.05;
wc = 0.55;
f = [0 wc wc+delta 1];
a = [0 0 1 1];
m = 300;
[b1 a1] = fir2(m,f,a);
[B1 w] = freqz(b1,a1);

figure(3)
sgtitle('Filter attenuation')
plot(w/pi, 20*log10(abs(B1)/max(abs(B1))))
yline(-65)
xlabel('w/pi')
ylabel('Filter magnitude')
legend('Filter')
grid()

%% Exercise 4

y = filter(b1,a1,x);

figure(4)
sgtitle('Filter output with x[n] as input')
stem(n,y)
xlabel('n')
ylabel('y[n]')
legend('y[n]')
grid()

%% Exercise 5
wc = 0.575;
[M wc] = buttord(wc, wc+0.05,3,65);
[b2 a2] = butter(M,wc,"high");
[B2 w] = freqz(b2,a2);


figure(5)
sgtitle('BW-filter attenuation')
plot(w/pi, 20*log10(abs(B2)/max(abs(B2))))
axis([0 1 -300 5])
yline(-65)
xlabel('w/pi')
ylabel('Butterworthfilter magnitude')
legend('BW-filter')
grid()

%% Exercise 6

y2 = filter(b2,a2,x);

figure(6)
sgtitle('BW-filter output with x[n] as input')
stem(n,y)
xlabel('n')
ylabel('y2[n]')
legend('y2[n]')
grid()

