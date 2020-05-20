function Xe=EstimationLOOP(L,x,B,Win,d,k,I,sigma,tau)
n=length(x);
G=fft(Win)/(length(Win));
w=exp(2*pi*1i/n);
Xe=zeros(L,n);
for r=1:L
%     sigma=randi([1,n/2])*2-1;
%     tau=randi([0,n-1]);
    y=Permute_sFFT(x,sigma(r),tau(r)).*Win;
    Z=subFFT(y,B);
    for i=1:length(I)   
        h=mod(round(sigma(r)*I(i)*B/n),B);
        o = sigma(r)*I(i)-round(sigma(r)*I(i)*B/n)*n/B;
%         o=mod(o,n/B);
        if o<0 
            o=n+o;
        end       
        Xe(r,I(i)+1)=Z(h+1)*w^(tau(r)*I(i))/G(o+1);
    end
end
