function c=squareTN(a,b)
% c=squareTN(a,b)
% ---------------
% Computes the Tensor Network c of the column-wise outer product of the
% two matrices a and b in the Tensor Network format.
%
% c,a,b     =   Tensor Networks.
%
% Reference
% ---------
%
% A Tensor Network Kalman filter with an application in recursive MIMO Volterra system identification
%
% 2016, Kim Batselier, Zhongming Chen, Ngai Wong

% column-wise outer product of 2 matrices in TN format

% l==1
l=a.n(1,3);
d=size(a.n,1);
if l==1
    for i=1:d
        c.core{i}=tkron(a.core{i},b.core{i});
    end
else
    firstcore=zeros(1,a.n(1,2),b.n(1,2),l,a.n(1,4)*b.n(1,4));
    for i=1:l
        tempa=reshape(a.core{1}(:,:,i,:),[1,a.n(1,2),1,a.n(1,4)]);
        tempb=reshape(b.core{1}(:,:,i,:),[1,b.n(1,2),1,b.n(1,4)]);
        tempa=reshape(permute(tempa,[1 2 4 3]),[a.n(1,4)*prod(a.n(1,1:2)) 1]);
        tempb=reshape(permute(tempb,[1 2 4 3]),[b.n(1,4)*prod(b.n(1,1:2)) 1]);
        temp=reshape(tempa*tempb',[a.n(1,1:2) a.n(1,4) b.n(1,1:2) b.n(1,4)]);
        temp=permute(temp,[1 4 2 5 3 6]);
        firstcore(:,:,:,i,:)=reshape(temp,[1,a.n(1,2),b.n(1,2),1,a.n(1,4)*b.n(1,4)]);
    end
    c.core{1}=firstcore;
    for i=2:d
        c.core{i}=tkron(b.core{i},a.core{i});
    end    
end
c.n=zeros(d,5);
c.n(:,1)=a.n(:,1).^2;
c.n(:,2)=a.n(:,2);
c.n(:,3)=a.n(:,2);
c.n(:,4)=a.n(:,3);
c.n(:,5)=a.n(:,4).^2;

% d=size(a.n,1); 	% number of columns of a.n has to be 4!! otherwise not a matrix
% 
% c.core=cell(1,d);
% tempa=reshape(permute(a.core{1},[4 1:3]),[a.n(1,4)*prod(a.n(1,1:2)) a.n(1,3)])';
% tempb=reshape(permute(b.core{1},[4 1:3]),[b.n(1,4)*prod(b.n(1,1:2)) b.n(1,3)])';
% 
% temp=dotkron(tempa,tempb);
% temp=reshape(temp,[a.n(1,3) a.n(1,4) a.n(1,1:2) b.n(1,4) b.n(1,1:2)]);
% c.core{1}=reshape(permute(temp,[3 6 4 7 1 2 5]),[a.n(1,1).*b.n(1,1) a.n(1,2) b.n(1,2) a.n(1,3) a.n(1,4).*b.n(1,4)]);
% for i=2:d
%     c.core{i}=tkron(reshape(b.core{i},b.n(i,:)),reshape(a.core{i},a.n(i,:)));
% end
% c.n=zeros(d,5);
% c.n(:,1)=a.n(:,1).*b.n(:,1);
% c.n(:,2)=a.n(:,2);
% c.n(:,3)=b.n(:,2);
% c.n(:,4)=a.n(:,3);
% c.n(:,5)=a.n(:,4).*b.n(:,4);

end
