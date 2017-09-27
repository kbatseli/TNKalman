function P=initP(sigmas,n,d)
% P=initP(sigmas,n,d)
% -------------------
% Constructs a Tensor Network P corresponding with a n^d x n^d x l tensor
% of diagonal matrix slices. Use this function to initialize the tensor
% of diagonal covariance matrices P(0) of the linear state model.
%
% P         =   Tensor Network, n^d x n^d x l tensor,
%
% sigmas 	=	vector, vector of length l such that sigmas(k) contains
% 				the scaling factor of the kth diagonal slice of P,
%
% n,d     =   scalars.
%
% Reference
% ---------
%
% A Tensor Network Kalman filter with an application in recursive MIMO Volterra system identification
%
% 2016, Kim Batselier, Zhongming Chen, Ngai Wong

% we assume the following ordering of the dimensions
% auxiliary, number outputs, rows, columns, auxiliary

k=length(sigmas);
P.n=[ones(d,1) n*ones(d,2) [k;ones(d-1,1)] ones(d,1)];
for i=1:k
	P.core{1}(:,:,i)=sigmas(i)*eye(n);
for i=2:d
	P.core{i}=eye(n);
end

end
