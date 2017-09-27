%% Multi-row SISO
clear all
M=4;
d=4;
l=1;
p=1;
n=p*M+1;
sigma_r=1e-8*ones(1,l);
% sigma_r=zeros(1,l);
N=100;
m=5;
u=randn((m+1)*N,p);
y=zeros((m+1)*N,l);
h=randn(p*M+1,l);
Href=mkhatri(h,d);
for i=M:m*N
    % The convention for the ordering of the u samples is
    % [1 u_1(t) u_2(t) ... u_p(t) u_1(t-1) u_2(t-1) .... u_p(t-1) ...
    % u_1(t-M+1) ..... u_p(t-M+1) ]    
	uk=[1 reshape(u(i:-1:i-M+1,:)',[1,p*M])];
	y(i,:)=(uk*h).^d + sqrt(sigma_r).*randn(1,l);
end
clear uk
tol=1e-10;							% tolerance in TN-rounding
e=zeros(N-M+1,m);
t=zeros(N-M+1,m);
for k=1:m
    mTN=initm(l,n,d);                     % initalize zero mean vector as a TN
    PTN=initP(1e1*ones(1,l),n,d);     	% initalize covariance matrix as a TN
    R=sigma_r*eye(k);
%     fi=animatedline;
    for i=M:N
        % cut up dataset into m slices and prepare IO samples
        U=zeros(k,n);
        Y=zeros(k,l);
        for j=1:k
            U(j,:)=[1 reshape(u((j-1)*N+i:-1:(j-1)*N+i-M+1,:)',[1,p*M])];
            Y(j,:)=y((j-1)*N+i,:);
        end
        tic;
        [mTN,PTN]=TNKalman(mTN,PTN,[],[],R,Y,U,tol);
        t(i-M+1,k)=toc;
        temp=Href-contract(mTN);
        e(i-M+1,k)=norm(temp(:))/norm(Href);
%         addpoints(fi,i-M+1,log10(e(i-M+1,k)));
%         drawnow limitrate
    end
%     drawnow
end
semilogy(e);grid on
