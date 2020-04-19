%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% At_fWH
%
% Written by: Chengbo Li
% CAAM, Rice University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x= At_fWH1(b, OMEGA, P,siz)
% function X= At_fWH(b, indx, indx2,siz)

N = length(P);
fx = zeros(N,1);
% fx(OMEGA) = b*sqrt(N);

fx(OMEGA) = b/sqrt(N);

x = zeros(N,1);
tmp = ifWHtrans(fx);
x(P) = tmp(1:N);
x=reshape(x,siz);



% y1 = zeros(prod(siz), 1);
% y=y1(:);
% y(indx) = b;
% save y
% 
% X2 =ifWHtrans(y)/sqrt(length(y));
% save X2
% X = 0*y;
% X(indx2) = X2(1:length(y));
% X = reshape(X, siz);