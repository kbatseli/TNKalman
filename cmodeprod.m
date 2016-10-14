function a=cmodeprod(a,u,k,l)
% a=cmodeprod(a,u,k,l)
% --------------------
% Computes the k-mode product of the lth core of a Tensor Network a with
% the matrix u. Typically used for scaling slices of a tensor.
%
% c         =   Tensor Network, result TN,
%
% a         =   Tensor Network, original TN,
%
% u 		=	matrix, matrix to be used in the mode-product,
%
% k 		=	scalar, indicates over which mode the mode-product occurs,
%
% l         =   scalar, indicates on which TN-core the mode-product occurs.
%
% Reference
% ---------
%
% A Tensor Network Kalman filter with an application in recursive MIMO Volterra system identification
%
% 2016, Kim Batselier, Zhongming Chen, Ngai Wong

[d,n]=size(a.n);

temp=reshape(a.core{l},a.n(l,:));
for j=1:length(k)
    I=1:size(a.n,2);
    I(k(j))=[];
    I=[k(j) I];
    temp=reshape(permute(temp,I),[a.n(l,I(1)) prod(a.n(l,I(2:end)))]);
    temp=reshape(u*temp,[size(u,1) a.n(l,I(2:end))]);
    a.n(l,k(j))=size(u,1);
    temp=permute(temp,[2:k(j) 1 k(j)+1:n]);
end
a.core{l}=temp;
end
