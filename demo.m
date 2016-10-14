%% Example SISO Volterra system
% -----------------------------
clear all
d=4; 						% order of Volterra system
M=4;						% Memory of Volterra system
p=1;						% single-input
l=1;						% single-output
sigma_r=1e-2*ones(1,l);		% measurement noise covariance
N=1000;						% number of samples
u=randn(N,p);				% input
y=zeros(N,l);				
h=randn(p*M+1,l);			
Href=mkhatri(h,d);			% Volterra kernel coefficients
for i=M:N
    % The convention for the ordering of the u samples is
    % [1 u_1(t) u_2(t) ... u_p(t) u_1(t-1) u_2(t-1) .... u_p(t-1) ...
    % u_1(t-M+1) ..... u_p(t-M+1) ]    
	uk=[1 reshape(u(i:-1:i-M+1,:)',[1,p*M])];
	y(i,:)=(uk*h).^d + sqrt(sigma_r).*randn(1,l);  % generate output sample
end

tol=1e-1;							% tolerance in TN-rounding
e=zeros(N-M+1,1);
t=zeros(N-M+1,1);
m=initm(l,(p*M+1),d);				% initalize zero mean vector as a TN
P=initP(1e3*ones(1,l),(p*M+1),d);	% initalize covariance matrix as a TN
for i=M:N
        uk=[1 reshape(u(i:-1:i-M+1,:)',[1,p*M])];
        tic;
        [m,P]=TNKalman(m,P,[],[],sigma_r,y(i,:),uk,tol);
        t(i-M+1)=toc;
        temp=Href-contract(m)';
        e(i-M+1)=norm(temp(:))/norm(Href);
end
figure
semilogy(e,'-o');
