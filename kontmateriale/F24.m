%% Exercise 1
N = 400;
n = 0:N-1;

x1 = 1/50*sin(pi/25*n).*cos(pi/50*n).*sin(pi/2*n);
x2 = 2*cos(4/5*pi*n).*x1;

[X1 w] = freqz(x1,1);
X2 = freqz(x2,1);

figure(1)
plot(w/pi,abs(X1))
title('Magnitude of x1 and x2')
hold on
plot(w/pi,abs(X2))
hold off
xlabel('w/\pi')
ylabel('|X_i(e^{jw})|')
legend('|X_1(e^{jw})|','|X_2(e^{jw})|')
grid()

%% Exercise 2

f = [0 0.38 0.39 0.61 0.62 1];
m = [1 1 0 0 1 1];
O = 400;

[b a] = fir2(O,f,m);
H = freqz(b,a);

figure(2)
plot(w/pi,20*log10(abs(H)/max(abs(H))))
title('Magnitude of x1 and x2')
xlabel('w/\pi')
ylabel('|H_i(e^{jw})|')
legend('|H_1(e^{jw})|')
grid()

%% Exercise 3
y1 = filter(b,a,x1);
y2 = filter(b,a,x2);

figure(3)
stem(n,y1);
title('Filter output with x_1[n] and x_2[n] as input')
hold on
stem(n,y2);
hold off
legend('y_1[n]','y_2[n]')
xlabel('n')
ylabel('y_i[n]')
grid()

