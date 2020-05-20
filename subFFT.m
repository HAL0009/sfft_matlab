function f=subFFT(x,B)
n=length(x);
if mod(n,B)~=0
    error('n不能被B整除');
end
temp=reshape(x,B,n/B);
y=sum(temp');
f=fft(y)/(length(y));
