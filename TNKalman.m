function [m,P]=TNKalman(m,P,A,Q,r,y,u,tol)
% [m,P]=TNKalman(m,P,A,Q,r,y,u,tol)
% ---------------------------------
% Computes the predict and update step for both the mean vectors m and 
% covariance matrices P in the Tensor Network format specifically for
% the recursive identification of Volterra systems.
%
% m,P 		=	Tensor Networks for matrix of mean vectors M and tensor
% 				of covariance matrices P,
%
% A	 		=	Tensor Network, represents the A matrix in the linear
% 				state space model,
%
% Q 		=	Tensor Network representation of tensor containing
% 				covariance matrices of the Gaussian process noise,
%
% r			=	vector, variances of the Gaussian measurement noise,
%
% y 		=	vector of l measured outputs,
%
% u			=	vector, regressor vector containing input values for 
%				MIMO Volterra system,
%
% tol	    =   scalar, relative approximation error in the TN-rounding.
%
% Reference
% ---------
%
% A Tensor Network Kalman filter with an application in recursive MIMO Volterra system identification
%
% 2016, Kim Batselier, Zhongming Chen, Ngai Wong

% make sure u,r,y are row,col,col vectors, respectively
u=u(:)';
r=r(:);
y=y(:);

%% Predict
%  -------

% M = M x_1 A

% P = P x_1 A x_2 A + Q
if ~isempty(Q)
    P=addTN(P,Q);
    P=roundTN(P,tol);
end

%% Update
%  -------

% v = y- M x_1 u^T
v=y-contract(modeprod(m,u,3));

% s = P x_1 u^T x_2 u^T + r
s=contract(modeprod(P,u,[3 4]))+r;

% k = P x_2 u^T x_3 S^{-1}
k=cmodeprod(modeprod(P,u,4),diag(1./s),2,1);
k=squeeze(k);
k=roundTN(k,tol);

% m = m + k x_2 diag(v)
m=addTN(m,cmodeprod(k,diag(v),2,1));
m=roundTN(m,tol);

% P = P - (k \square k) x_3 S
ks=squareTN(k,k);
ks=roundTN(ks,tol);
P=addTN(P,cmodeprod(ks,diag(-s),2,1));
P=roundTN(P,tol);

end
