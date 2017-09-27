function a=squeeze(a)
% a=squeeze(a)
% ------------
% Removes redundant singleton dimensions from a Tensor Network a.
%
% a     =   Tensor Networks.
%
% Reference
% ---------
%
% A Tensor Network Kalman filter with an application in recursive MIMO Volterra system identification
%
% 2016, Kim Batselier, Zhongming Chen, Ngai Wong

% % if a mode consists of all ones after the mode product, then we remove
% % that mode
[d,n]=size(a.n);
I=sum(a.n(:,2:end-1),1)==size(a.n(:,2:end-1),1);
a.n(:,[false I false])=[];

% if we are left with 2 auxiliaries and 1 free index, then we have a vector,
% however, in the kalman setting we want to explicitly keep the unitary
% column index, so we need to artificially add it
if size(a.n,2)==3           
    % dealing with a vector but need to ones for consistency
    a.n=[a.n(:,1:2) ones(d,1) a.n(:,end)];
end
for i=1:d
    a.core{i}=reshape(a.core{i},a.n(i,:));
end
end
