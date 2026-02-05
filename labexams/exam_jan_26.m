%% Task 1
h = [0.0024 0.0363 -0.2327 0.3827];

h_I = [h fliplr(h(1:(end-1)))];
h_II = [h fliplr(h)];
h_III = [h 0 -fliplr(h)];
h_IV = [h -fliplr(h)];

%% Task 2
[H_I, w1] = freqz(h_I);
[H_II, w2] = freqz(h_II);
[H_III, w3] = freqz(h_III);
[H_IV, w4] = freqz(h_IV);


subplot(2,1,1);
hold on
plot(w1/pi,20*log10(abs(H_I)));
plot(w1/pi,20*log10(abs(H_II)));
plot(w1/pi,20*log10(abs(H_III)));
plot(w1/pi,20*log10(abs(H_IV)));
legend('type I', 'type II', 'type III', 'type IV')
title('Filter magnitudes')
xlabel('\omega / pi')
ylabel('H_i')
hold off


%% Task 3
N = 400;
m = 0:N-1;
x = cos(0.12*pi*m);

[X,wx] = freqz(x);

subplot(2,1,2)
title('Magnitude of X')
xlabel('\omega / pi')
ylabel('|X(e^(j\omega))')
plot(wx/pi,abs(X))
