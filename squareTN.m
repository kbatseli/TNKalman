function c=squareTN(a,b)
% c=squareTN(a,b)
% ---------------
% Computes the Tensor Network c of the column-wise outer product of the
% two matrices a and b in the Tensor Network format.
%
% c,a,b     =   Tensor Networks.
%
% Reference
% ---------
%
% A Tensor Network Kalman filter with an application in recursive MIMO Volterra system identification
%
% 2016, Kim Batselier, Zhongming Chen, Ngai Wong

% column-wise outer product of 2 matrices in TN format

[d,n]=size(a.n); 	% n has to be 4!! otherwise not a matrix

c.core=cell(1,d);
tempa=reshape(permute(a.core{1},[2:n 1]),[a.n(1,2) a.n(1,1)*prod(a.n(1,3:end))])';
tempb=reshape(permute(b.core{1},[2:n 1]),[b.n(1,2) b.n(1,1)*prod(b.n(1,3:end))])';

temp=khatri(tempa,tempb);
temp=reshape(temp,[b.n(1,3) b.n(1,end) b.n(1,1) a.n(1,3) a.n(1,end) a.n(1,1) a.n(1,2)]);
c.core{1}=reshape(permute(temp,[3 6 7 1 4 2 5]),[a.n(1,1).*b.n(1,1) a.n(1,2) b.n(1,3)*a.n(1,3) a.n(1,4).*b.n(1,4)]);
for i=2:d
   c.core{i}=tkron(a.core{i},b.core{i});
end
c.n=zeros(d,n+1);
c.n(:,1)=a.n(:,1).*b.n(:,1);
c.n(:,2)=b.n(:,2);
c.n(:,3)=a.n(:,3);
c.n(:,4)=b.n(:,3);
c.n(:,5)=a.n(:,4).*b.n(:,4);

end
