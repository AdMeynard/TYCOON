t = linspace(0,15,15*50);

am1 = 1+0.2*atan(t-10);
x1 = am1.*cos(2*pi*(2*t+1.4*cos(t)));
x1(1:256) = 0;

am2 = 1 + sqrt(1+0.1*cos(t));
x2 = am2.*cos(2*pi*(5*t+0.5*t.^(1.5)+0.5*cos(t)));

x3 = cos(2*pi*(10*t-0.2*t.^2));

noise = 0.542*randn(size(am1));
x = x1 + x2 + x3+noise;

Fe=50;
s=x;
save('signal_test','t','s','Fe')
plot(t,s)
