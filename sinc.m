function y = sinc(x,B,n)
fc=1/(2*B);
y = sin(2*pi*fc*x)./x/pi/fc/2;
y(x==0) = 1;
