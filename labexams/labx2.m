%% Task 1
% Compute X1(jw) and X2(jw) and plot their magnitudes

N = 400;
n = 0:N-1;

x1 = 1/50*sin(pi/25*n).*cos(pi/50*n).*sin(pi/2*n);
x2 = 2*cos(4/5*pi*n).*x1;

[X1 w] = freqz(x1,1,512);
[X2 w] = freqz(x2,1,512);

subplot(2,1,1)
title('|X1(jw)| and |X2(jw)|')
xlabel('w')
ylabel('Magnitude')
plot(w/pi,abs(X1))
hold on
plot(w/pi,abs(X2))
legend('|X1(jw)|', '|X2(jw)|')
hold off

%% Task 2
% Design a filter to attenuate x1 by 45dB
% that does not attenuate x2


% Mp = 150;
% transitionbands = [0.36 0.4 0.6 0.65];
% fp = [0 transitionbands 1];
% amps = [1 1 0 0 1 1];
% bp = firpm(Mp,fp,amps);

% [BP, w] = freqz(bp,1,1024);
% subplot(2,1,2)
% plot(w/pi,abs(BP))

M = 200;
[b,a] = fir1(M,[0.4 0.6],"stop");
[B,w] = freqz(b,a,1024);
subplot(2,1,2)

plot(w/pi,20*log10(abs(B)))
yline(-45)

% y1 = filter(b,a,x1);
% y2 = filter(b,a,x2);

y1 = conv(b,x1);
y2 = conv(b,x2);

