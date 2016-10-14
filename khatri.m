function C=khatri(A,B)
% C=khatri(A,B)
% -------------
% Naive implementation of the column-wise Kronecker (Khatri-Rao) product
% of two matrices A and B.
%
% C         =   matrix, each column of this m*n x l matrix is the Kronecker
% 				product of the corresponding columns of A with B,
%
% A         =   m x l matrix,
%
% B			=	n x l matrix.
%
% Reference
% ---------
%
% A Tensor Network Kalman filter with an application in recursive MIMO Volterra system identification
%
% 2016, Kim Batselier, Zhongming Chen, Ngai Wong

% Khatri-Rao product, columnwise Kronecker product kron(A(:,k),B(:,k)) over
% all k columns.

if size(A,2)~=size(B,2)
    error('Matrices A,B should have the same number of columns.')
end

C=zeros(size(A,1)*size(B,1),size(A,2));
for i=1:size(A,2)
   C(:,i)=kron(A(:,i),B(:,i)); 
end

end
