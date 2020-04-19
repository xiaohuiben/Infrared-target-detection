%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A_fWH
%
% Written by: Chengbo Li
% CAAM, Rice University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = A_fWH(x, OMEGA, P)  %%带入参数为x,picks,perm 分别为要处理的图像、M个随机数、N个随机数


N = length(P);

x = x(:);
% fx = fWHtrans(x(P))*sqrt(N);
 fx = fWHtrans(x(P))/sqrt(N);                             %%%%%%  computes the (real) fast discrete Walsh-Hadamard transform with sequency order 
b = fx(OMEGA);