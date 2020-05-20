function [X]=sfft(x,k,d,B,Lc)
n = length(x);
Beta = k/(2^(floor(log2(n)-4)));
if nargin < 5
    Lc = 5;
end
if nargin < 4
    B = k/Beta;
end
if nargin < 3
    d = 1;
end

L = Lc*floor(log2(n));  
Win=window(B,length(x));
sigma=randi([1,n/2],L,1)*2-1;
tau=randi([0,n-1],L,1);
I=LocationLOOP(L,x,B,Win',d,k,sigma,tau);
Isum=[];
for i=1:L
    Isum=[Isum I{i}];
end
Isum=unique(Isum);
s=zeros(1,length(Isum));
for i =1:length(Isum)
    for r=1:L
        if ismember(Isum(i),I{r})
            s(i)=s(i)+1;
        end
    end
end
Isum_new=Isum(s>=L/2);

% Xe=EstimationLOOP(L,x,B,Win',d,k,Isum_new,sigma,tau);
% Xe_new=complex(median(real(Xe)),median(imag(Xe)));
% Xe_new=median(real(Xe))+1i*median(imag(Xe));

X=zeros(1,length(x));
for i = 1 : length(Isum_new)
    X(Isum_new(i)+1)=x*exp(-2*pi*1i*Isum_new(i)*(0:n-1)/n)';
end