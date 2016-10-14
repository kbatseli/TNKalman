function m=initm(l,n,d)
% m=initm(l,n,d)
% --------------
% Constructs a Tensor Network *m* corresponding with a n^d x l zero
% matrix. Use this function to initialize the matrix of means M(0) of
% the linear state model.
%
% m         =   Tensor Network, n^d x l zero matrix,
%
% l,n,d     =   scalars.
%
% Reference
% ---------
%
% A Tensor Network Kalman filter with an application in recursive MIMO Volterra system identification
%
% 2016, Kim Batselier, Zhongming Chen, Ngai Wong

% we assume the following ordering of the dimensions
% auxiliary, number outputs, rows, auxiliary

m.n=[ones(d,1) [l;ones(d-1,1)] n*ones(d,1) ones(d,1)];
m.core{1}=zeros(1,n,l);
for i=2:d
	m.core{i}=zeros(1,n,1);
end

end
