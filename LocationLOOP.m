function I=LocationLOOP(L,x,B,Win,d,k,sigma,tau)
n=length(x);
for r=1:L
%     sigma=randi([1,n/2])*2-1;
%     tau=randi([0,n-1]);
    p=Permute_sFFT(x,sigma(r),tau(r));
    y=p.*Win;
    Z=subFFT(y,B);
    [~,indx]=sort(Z);
    J=indx(end-d*k+1:end);
    h=mod(round(sigma(r)*(0:n-1)*B/n),B)+1;
%      h(h==0)=B;
    Lia=ismember(h,J);
    I{r}=find(Lia==1)-1;
%     I{r}=unique(Lob(Lia==1));
end