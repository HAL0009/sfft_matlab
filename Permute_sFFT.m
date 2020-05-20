function y = Permute_sFFT(x,sigma,tau)
%sigma为n的模逆元
n=length(x);
i = mod((0:n-1)*sigma+tau,n)+1;
y=x(i);