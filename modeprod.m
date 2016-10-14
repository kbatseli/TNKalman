function a=modeprod(a,u,k)
% a=modeprod(a,u,k)
% -----------------
% Computes the k-mode product of all cores of the Tensor Network a with 
% the matrix u.
%
% a         =   Tensor Network,
%
% u 		=	matrix, matrix to be used in the mode-product,
%
% k 		=	scalar, indicates over which mode the mode-product occurs.
%
% Reference
% ---------
%
% A Tensor Network Kalman filter with an application in recursive MIMO Volterra system identification
%
% 2016, Kim Batselier, Zhongming Chen, Ngai Wong

[d,n]=size(a.n);

for i=1:d
    temp=reshape(a.core{i},a.n(i,:));
    for j=1:length(k)
        I=1:size(a.n,2);
        I(k(j))=[];
        I=[k(j) I];
        temp=reshape(permute(temp,I),[a.n(i,I(1)) prod(a.n(i,I(2:end)))]);
        temp=reshape(u*temp,[size(u,1) a.n(i,I(2:end))]);
        a.n(i,k(j))=size(u,1);
        temp=permute(temp,[2:k(j) 1 k(j)+1:n]);
    end
    a.core{i}=temp;
end
% % if a mode consists of all ones after the mode product, then we remove
% % that mode
% a.n(:,sum(a.n,1)==size(a.n,1))=[];
% for i=1:d
%     a.core{i}=reshape(a.core{i},a.n(i,:));
% end

end
