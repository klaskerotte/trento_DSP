%% Exercise 1
N = 400;
n = 0:N-1;

x1 = 1/50*sin(pi/25*n).*cos(pi/50*n).*sin(pi/2*n);
x2 = 2*cos(4/5*pi*n).*x1;

[X1 w] = freqz(x1,1);
[X2 w] = freqz(x2,1);

figure(1)
sgtitle('Magnitudes of X1 and X2')
plot(w/pi,abs(X1))
legend('X1')
xlabel('w/pi')
ylabel('|X_1(jw)|')
grid()

figure(2)
sgtitle('Magnitudes of X1 and X2')
plot(w/pi,abs(X2))
legend('X2')
xlabel('w/pi')
ylabel('|X_2(jw)|')
grid()


%% Exercise 2
M = 400;
delta = 0.05;
wc1 = 0.4;
wc2 = 0.6;
[b1 a1] = fir2(M,[0 wc1-delta wc1 wc2 wc2+delta 1],[1 1 0 0 1 1]);
[B1 wb] = freqz(b1,a1);

figure(3)
sgtitle('Filter magnitude in dB')
plot(wb/pi, 20*log10(abs(B1)/max(abs(B1))));
xline([0.4 0.6])
xlabel('w/pi')
ylabel('dB')
grid()

%% Exercise 3

y1 = filter(b1,a1,x1);
y2 = filter(b1,a1,x2);

figure(4)
sgtitle('Output y1[n]')
stem(n,y1)
xlabel('n')
ylabel('y1[n]')
legend('y1')
grid()

figure(5)
sgtitle('Output y2[n]')
stem(n,y2)
xlabel('n')
ylabel('y2[n]')
legend('y2')
grid()
