

%Li Li, Hui Li*, Tian Li, Feng Gao. Infrared small target detection in compressive domain[J].

% Electronics Letters,  2014, 50(7):510-512 
%By hui li @BeihangU 2014 3.20 

close all;
clear all;
clc;


T=10% 设定T的值
for i=1:T
        
path(path,genpath(pwd));

if exist('fWHtrans','file') ~= 3
    cd Fast_Walsh_Hadamard_Transform;
    mex -O fWHtrans.cpp
    cd ..;
    fprintf('Finished compiling the C++ code for fast Walsh-Hadamard transform!\n');
end
fprintf('Finished adding paths! Welcome to use  fast Walsh-Hadamard transform as sensing matrix.\n');
% addpath(genpath('./algs/'))
% addpath(genpath('./operators/'))
% addpath(genpath('./PROPACK/'))
% addpath(genpath('./utility/'))

%  I=imread('11fu.jpg');
 I=imread('16.jpg');

if size(I,3)==3
I=rgb2gray(I);
end

I=imresize(I,[128,128]);
% figure(22),imshow(I);

Im=double(I);
[N1,N2]=size(Im);
R =20;  %低秩的选取应用熵来估计，在此暂且取经验值16-20
sratio= 0.0004;%The rank of L
K = ceil(sratio*N1*N2);  

                      %The sparsity level of S
M = floor(0.6*N1*N2);  %The number of measurements to acquire



tp = randperm(N1*N2);
picks = tp(1:M);
for ii = 1:M
    if picks(ii) == 1
        picks(ii) = tp(M+1);
        break;
    end
end
perm = randperm(N1*N2); % column permutations allowable

        %This solves the SpaRCS problem
% 		idx = randperm(N1*r); idx = idx(1:p);
% 		idx2 = randperm(r*r);
% 		A = @(z) Anoiselet(z,idx, idx2);
% 		At = @(z) Atnoiselet(z,idx, idx2, size(M));
                A =@(z)A_fWH(z,picks,perm);
                At= @(z)At_fWH1(z,picks,perm,size(Im));
                
%                 idx=randperm()
%                 AA = @(z) Anoiselet(z,idx, idx2);
% 		Att = @(z) Atnoiselet(z,idx, idx2, size(M));

% b=[];
% for ii=1:((128-r)/r1+1).^2
%     temp1=A(Im(:,ii));
%     b=[b,temp1];
% end    
% [N1,N2]=size(b);
  
 b=A(Im(:));
% b = A(Im(:));
% b = A(Im(:),1); 
tic
[A E err] = sparcs(b, R, K, A, At, 'propack',5e-4, 10, 1); 
toc
e1=E;

A1= A- min(A(:));
mid=max(max(A1));
A=A1/mid.*255;
% A=uint8(A2);
%***********************************************
E1= E - min(E(:));
midE=max(max(E1));
E=E1/midE.*255;
% S11=uint8(S1);


% f1=zeros(128,128);
% c1=zeros(128,128);
% f2=zeros(128,128);
% c2=zeros(128,128);


%  for i=1:r1:128-r+1
%        for j=1:r1:128-r+1
%           e1=reshape(Aa(:,(ceil(i/r1)-1)*((128-r)/r1+1)+ceil(j/r1)),r,r);
% %           f1(i:i+49,j:j+49)=min(f1(i:i+49,j:j+49),e1);
%            f1(i:i+r-1,j:j+r-1)=f1(i:i+r-1,j:j+r-1)+e1;
%           c1(i:i+r-1,j:j+r-1)= c1(i:i+r-1,j:j+r-1)+1;
%        end
%  end
%  A=f1./c1;
% %   A=f1;
% 
% for i=1:r1:128-r+1
%        for j=1:r1:128-r+1
%             e2=reshape(Ea(:,(ceil(i/r1)-1)*((128-r)/r1+1)+ceil(j/r1)),r,r);
% %           e2=median(e2);
% %             f2(i:i+49,j:j+49)=min(f2(i:i+49,j:j+49),e2);
%             f2(i:i+r-1,j:j+r-1)=f2(i:i+r-1,j:j+r-1)+e2;
%           c2(i:i+r-1,j:j+r-1)= c2(i:i+r-1,j:j+r-1)+1;
%        end
%  end
%  E=f2./c2;
%   E=f2;
A_fshow=uint8(A);


%  H = fspecial('average',[3 3]);   %%滤波算子
%  E = imfilter(E,H); 
E_fshow(:,:,i)=E;



end


sum1=0;
for i=1:T
        sum1=sum1+E_fshow(:,:,i);
        mean1=sum1/T;
end


E1=uint8(mean1);

E11=double(E1);
M=max(max(E11));
ro=0.6;
T=E11>=ro*M;
tar=T.*E11;

%  H = fspecial('average',[4 4]);   %%滤波算子
%  estIm8 = imfilter(estIm8,H); 

% threshold1=qiuyuzhi(E_fshow);
% 
% tar= im2bw(E_fshow,threshold1);





figure(1),subplot(2,2,1),imshow(I);title('原图像');
          subplot(2,2,2),imshow(A_fshow);title('低秩矩阵恢复的图像');
          subplot(2,2,3),imshow(E1); title('稀疏矩阵恢复的图像');
               
          subplot(2,2,4),imshow(tar);title('postprocessing');
%           figure(5),imshow(E1),imwrite(E1,'sparse.tif','tif'); 
          
[x,y]=size(I);
[X,Y]=meshgrid(1:y,1:x);
% % grid on;

% [x1,y1]=size(I);
% [X1,Y1]=meshgrid(1:y1,1:x1);

% x=1:150;
% y=1:150;
% [X,Y]=meshgrid(x,y);
% 
% figure(2),subplot(1,3,1),mesh(X,Y,double(I)/255);
%        
% subplot(1,3,2),mesh(X,Y,double(E_fshow)/255);
%           subplot(1,3,3),mesh(X,Y,double(tar)/255);

figure(2),plot3(X,Y,double(I)/255);grid on
    xlim([1 128]);   ylim([1 128]);
         figure(3),plot3(X,Y,double(E1)/255),grid on
         xlim([1 128]);   ylim([1 128]);
  figure(4),plot3(X,Y,double(tar)/255);grid on
    xlim([1 128]);   ylim([1 128]);
