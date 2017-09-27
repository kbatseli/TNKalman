function TN=mkr2tn(U,varargin)
% TN=mkr2tn(U) or TN=mkr2tn(U,d)
% ------------------------------
% When given a single cell argument U, then the Tensor Network corresponding
% with the matrix formed from the row-wise Kronecker product of
% U{1},U{2},...,U{d} is formed. If U is a matrix then the Tensor Network of
% the d-times repeated row-wise Kronecker product of U is formed.
%
% TN 	=   Tensor Network structure,
%
% U 	=   cell/matrix,
%
% d 	=	scalar, argument only required for matrix U to indicated how many times
%           its row-wise Kronecker product needs to be taken.
%
% Reference
% ---------
%
% 11/2016, Kim Batselier

if iscell(U)
    d=length(U);    
else
    d=varargin{1};
    temp=U;
    U=cell(1,d);
    for i=1:d
        U{i}=temp;
    end
end

for i=1:d
    [N,n(i)]=size(U{i});    % first dimension has to be the same over all factor matrices
end
TN.core=cell(1,d);
TN.n=ones(d,4);

% initialize last core
TN.core{d}=U{d};
TN.n(d,:)=[1 N n(d) 1];

for i=d:-1:2
    temp=reshape(TN.core{i},[prod(TN.n(i,1:2)),prod(TN.n(i,3:4))]);
    temp=dotkron(U{i-1},temp);
    temp=reshape(temp,[N*n(i-1),n(i)*TN.n(i,end)]);
    [U1,S1,V1]=svd(temp);
    s=diag(S1);
    tol=eps(s(1))*max([N,n(i),n(i-1)]);
    r=sum(s>tol);
    TN.core{i-1}=U1(:,1:r);
    TN.core{i}=S1(1:r,1:r)*V1(:,1:r)';
    TN.n(i-1,2:end)=[N n(i-1) r];
    TN.n(i,1:end-1)=[r 1 n(i)];    
end
end