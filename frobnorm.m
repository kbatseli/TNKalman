function s=frobnorm(a)
% s=frobnorm(a)
% -------------
% Computes the Frobenius norm of a tensor a in the Tensor Network format.
%
% s         =   scalar, Frobenius norm of tensor a,
%
% a         =   Tensor Network.
%
% Reference
% ---------
%
% A Tensor Network Kalman filter with an application in recursive MIMO Volterra system identification
%
% 2016, Kim Batselier, Zhongming Chen, Ngai Wong

[d,n]=size(a.n);
s=1;
for i=1:d
   Z=s*reshape(a.core{i},[a.n(i,1),prod(a.n(i,2:end))]);   
   s=reshape(a.core{i},[prod(a.n(i,1:end-1)),a.n(i,end)])'*reshape(Z,[prod(a.n(i,1:end-1)),a.n(i,end)]);
end
s=sqrt(s);


end
