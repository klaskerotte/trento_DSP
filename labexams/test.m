N=40;
transitionbands = [0.3 0.4, 0.5 0.6, 0.8 0.9];
f=[0 transitionbands 1];
m=[0 0 1 0 1 0 1 4];
w = [1 2 1 1];
b=firpm(N,f,m,w);
[H,w]=freqz(b,1);
subplot(211)
plot(w/pi,abs(H),f,m);
xline([transitionbands])
subplot(212)
plot(w/pi,unwrap(angle(H))/pi);