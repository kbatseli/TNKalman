function yhat=sim_volterraTN(u,TN)
% yhat=sim_volterraTN(u,TN)
% --------------------------
% Simulates a truncated MIMO Volterra series in the Tensor Network (TN) format for
% given inputs u(:,1),u(:,2),....
%
% yhat      =   matrix, y(:,k) contains the kth simulated output,
%
% u         =   matrix, u(:,k) contains the kth input,
%
% TN        =   Tensor Network.
%
% Reference
% ---------
%
% A Tensor Network Kalman filter with an application in recursive MIMO Volterra system identification
%
% 2016, Kim Batselier, Zhongming Chen, Ngai Wong

% number of inputs
p=size(u,2);       
% number of cores in the TN
d=size(TN.n,1);
% p*M+1
N=TN.n(1,3);
% number of outputs
l=TN.n(1,2);
% compute the Memory from N=p*M+1
M=(N-1)/p;
yhat=zeros(size(u,1),l);
for j=M:size(u,1)
    % The convention for the ordering of the u samples is
    % [1 u_1(t) u_2(t) ... u_p(t) u_1(t-1) u_2(t-1) .... u_p(t-1) ...
    % u_1(t-M+1) ..... u_p(t-M+1) ]
    uj=[1 reshape(u(j:-1:j-M+1,:)',[1,p*M])];
    yhat(j,:)=contract(modeprod(TN,uj,3));
end
end
