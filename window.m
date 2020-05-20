  
function [g] = window(B, n)


% w = floor(B*log(n/(2*10^(-9)))/pi)
w=floor(B*log2(n))
if w>n
    w=n;
end
sigma = B*sqrt(log2(n));

% gauss = gausswin(w, (w-1)/sigma/2);
gauss = gausswin(w);
gauss = [gauss; zeros(n-w,1)];
b = [sinc(-w/2:w/2-1,B)'; zeros(n-w,1)];
g = gauss.*b;

end
