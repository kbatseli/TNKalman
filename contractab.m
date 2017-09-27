function c=contractab(a,b,k)
% c=contractab(a,b,k)
% -------------------
% Contracts the cores of Tensor Network a along mode k(1) with cores of
% Tensor Network b along mode k(2).
%
% c,a,b     =   Tensor Network,
%
% k 		=	vector, k(1) is mode of a and k(2) is mode of b.
%
% Reference
% ---------
%
% A Tensor Network Kalman filter with an application in recursive MIMO Volterra system identification
%
% 2016, Kim Batselier, Zhongming Chen, Ngai Wong

[da,na]=size(a.n);
[db,nb]=size(b.n);
c.n=ones(da,na+nb-4);

for i=1:da
    tempa=reshape(a.core{i},a.n(i,:));
    tempb=reshape(b.core{i},b.n(i,:));
    Ia=1:na;
    Ib=1:nb;
    Ia(k(1))=[];
    Ia=[Ia k(1)];
    Ib(k(2))=[];
    Ib=[k(2) Ib];
    tempa=reshape(permute(tempa,Ia),[prod(a.n(i,Ia(1:end-1))),a.n(i,Ia(end))]);
    tempb=reshape(permute(tempb,Ib),[b.n(i,Ib(1)),prod(b.n(i,Ib(2:end)))]);
    temp=tempa*tempb;
    newsize=[a.n(i,Ia(1:end-1)) b.n(i,Ib(2:end))];
    temp=reshape(temp,newsize);    
    temp=permute(temp,[1 na 2:k(1)-1 na+1:na+nb-3 k(1):na-1 na+nb-2]);
    newsize=newsize([1 na 2:k(1)-1 na+1:na+nb-3 k(1):na-1 na+nb-2]);
    c.core{i}=reshape(temp,[newsize(1)*newsize(2),(newsize(3:end-2)),newsize(end-1)*newsize(end)]);
    c.n(i,:)=[newsize(1)*newsize(2),newsize(3:end-2),newsize(end-1)*newsize(end)];
end


end
