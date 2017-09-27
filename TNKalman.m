function [mTN,PTN]=TNKalman(mTN,PTN,A,QTN,r,Y,U,tol)
% [mTN,PTN]=TNKalman(mTN,PTN,A,QTN,r,y,u,tol)
% -------------------------------------------
% Computes the predict and update step for both the mean vectors m and 
% covariance matrices P in the Tensor Network format specifically for
% the recursive identification of Volterra systems.
%
% mTN,PTN   =	Tensor Networks for matrix of mean vectors M and tensor
% 				of covariance matrices P,
%
% A	 		=	Tensor Network, represents the A matrix in the linear
% 				state space model,
%
% QTN 		=	Tensor Network representation of tensor containing
% 				covariance matrices of the Gaussian process noise,
%
% r			=	vector, variances of the Gaussian measurement noise,
%
% Y 		=	matrix, each row corresponds with l measured outputs,
%
% U			=	matrix, regressor vector containing input values for 
%				MIMO Volterra system, each row corresponds with an
%				experiment,
%
% tol	    =   scalar, relative approximation error in the TN-rounding.
%
% Reference
% ---------
%
% A Tensor Network Kalman filter with an application in recursive MIMO Volterra system identification
%
% 2016, Kim Batselier, Zhongming Chen, Ngai Wong

% The assumption is that m,l << n^d

% U is a m x n matrix, each row corresponds with an experiment
% n = p*M+1
[m,n]=size(U);
    
if m==1
    % only 1 row of U, no repeated Kronecker product required
    [mTN,PTN]=m1TNKalman(mTN,PTN,A,QTN,r,Y,U,tol);
else

% Y is m x l matrix, where l is the number of outputs
l=size(Y,2);

% construct TN of m x n^d C matrix, row-wise Kroncker product of u_t
% vectors
d=size(mTN.n,1);
C=mkr2tn(U,d);

%% Predict
%  -------

% M = M x_1 A

% P = P x_1 A x_2 A + Q
if ~isempty(QTN)
    PTN=addTN(PTN,QTN);
    PTN=roundTN(PTN,tol);
end

%% Update
%  -------

% V = Y - C*M, m x l matrix
V=Y-contract(contractab(mTN,C,[2 3]));

% s = C*P*C^T + r, m x m matrix
temp=contractab(PTN,C,[2 3]);
temp=contractab(temp,C,[3 3]);
S=contract(temp)+r;
clear temp

% compute inverse slices of S
Sinv=zeros(size(S));
for i=1:l
    Sinv(:,:,i)=inv(S(:,:,i));
end

% K = P*C^T*S^{-1} + R, TN Kalman gain, n^d x m x l tensor
KTN=contractab(PTN,C,[3 3]);
tempcore=permute(KTN.core{1},[5 1 2 3 4]);
temp=zeros(KTN.n(1,:));
for i=1:l,
    temp(:,:,:,i,:)=permute(reshape(reshape(tempcore(:,:,:,:,i),[KTN.n(1,end)*prod(KTN.n(1,1:2)) KTN.n(1,3)])*Sinv(:,:,i),[KTN.n(1,end),KTN.n(1,1:2),KTN.n(1,3)]),[2 3 4 5 1]);
end
KTN.core{1}=temp;
KTN=roundTN(KTN,tol);
clear temp tempcore

% M = M + K * V, n^d x l matrix
KVTN=KTN;
tempcore=permute(reshape(KTN.core{1},KTN.n(1,:)),[4 5 1 2 3]);
temp=zeros(KTN.n(1,[1:2 4:5]));
for i=1:l
    temp(:,:,i,:)=permute(reshape(reshape(tempcore(i,:,:,:,:),[KTN.n(1,5)*prod(KTN.n(1,1:2)),KTN.n(1,3)])*V(:,i),[1,KTN.n(1,5),KTN.n(1,1:2)]),[3 4 1 2]);
end
KVTN.core{1}=temp;
KVTN.n(1,:)=[KTN.n(1,1:2) 1 KTN.n(1,4:5)];
KVTN=squeeze(KVTN);
mTN=addTN(mTN,KVTN);
mTN=roundTN(mTN,tol);

% P = P - (k \square k) x_3 S
% 2nd core up to last core are tensor kronecker products of K-cores
KK.core=cell(1,d);
KK.n=zeros(d,5);
for i=2:d
    KK.core{i}=reshape(KTN.core{i},KTN.n(i,:));
    KK.core{i}=tkron(KTN.core{i},KTN.core{i});
    KK.n(i,:)=[KTN.n(i,1).^2 n*ones(1,2) KTN.n(i,end-1:end).^2];
end
% first core of KK through slice-wise contraction
firstcore=zeros(1,n,n,l,KK.n(2,1));
for i=1:l
    temp=reshape(KTN.core{1},KTN.n(1,:));
    temp=reshape(permute(temp(:,:,:,i,:),[1 2 4 5 3]),[KTN.n(1,2)*KTN.n(1,5),KTN.n(1,3)]);
    temp=-temp*S(:,:,i)*temp';
    temp=reshape(temp,[1 n KTN.n(1,5) 1 n KTN.n(1,5)]);
    temp=permute(temp,[1 4 2 5 3 6]);
    firstcore(:,:,:,i,:)=reshape(temp,[1 n n 1 KK.n(2,1)]);    
end
KK.core{1}=firstcore;
KK.n(1,:)=[1 n n l KK.n(2,1)]; 
PTN=addTN(PTN,KK);
PTN=roundTN(PTN,tol);

end
end
