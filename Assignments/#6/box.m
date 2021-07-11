function [b]=box(K1,K2)

% K1=3;
% K2=5;
f1=1e9;
f2=0.9e9;
w1=2*pi*f1;
w2=2*pi*f2;



j=zeros(K1+1,(2*K2)+1);


p=K2+1;

%j(1,2)=2

for k2=-K2:1:K2
    for k1=1:1:K1
        m=k2+p;
     j(k1+1,m)=abs(k1*w1+k2*w2);
end 
end

for jj=p+1:(2*K2)+1
    j(1,jj)= abs((jj-p)*w2);
end
for i=0:1:K1
    b(i+1,:)=j((end-i),:);
end
end

