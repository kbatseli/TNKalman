function C=mkhatri(A,d)
% C=mkhatri(A,d)
% --------------
% Calls the khatri.m function d-times on the matrix A.
%
% C         =   matrix, each column of this m^d x l matrix is the Kronecker
% 				product of the corresponding columns of A d-times with itself,
%
% A         =   m x l matrix,
%
% d			=	scalar.
%
% Reference
% ---------
%
% A Tensor Network Kalman filter with an application in recursive MIMO Volterra system identification
%
% 2016, Kim Batselier, Zhongming Chen, Ngai Wong

% Khatri-Rao product, columnwise Kronecker product kron(A(:,k),B(:,k)) over
% all k columns.

C=zeros(size(A,1)^d,size(A,2));
for i=1:size(A,2)
   C(:,i)=mkron(A(:,i),d); 
end

end
