%Please cite the paper "Xiaohong Jia,Tao Lei,Xiaogang Du,Shigang Liu,Hongying Meng,and Asoke K. Nandi,
%Robust Self-Sparse Fuzzy Clustering for Image Segmentation" 
%has been accepted for publication in IEEE Access.

%The code was written by Xiaohong Jia and Tao Lei in 2020.

%We would like to thank Dr. Xu for providing codes related to the paper "Robust
%and Sparse Fuzzy K-Means Clustering"

%references
% Robust and sparse fuzzy k-means clustering
% A new simplex sparse learning model to measure data similarity for clustering

% Welcome to our Research Group website:https://aimv.sust.edu.cn/lwcg.htm

close all
clear all
f_ori=imread('12003.jpg');
[rows,cols,dim]=size(f_ori);
if dim==3
    F=rgb2lab(f_ori);
else
    F=double(f_ori);
end
figure,imshow(f_ori)
Data=reshape(F,rows*cols,dim);
cluster_n=2;
gama=0.2;
maxIter=50;
%% Robust Self-Sparse Fuzzy Clustering for Image Segmentation
[outA,outB,outObj,outNumIter]=RSSFCA(Data',gama,maxIter,cluster_n);
[~,L]=max(outA);
Label=reshape(L,rows,cols);
figure,imshow(Label,[])

%%  connected-component filtering based on area density balance strategy
Th=cluster_n;  %Tune the parameter according to your needs.
RL=CCF_ADB(Label,Th);
R9=Label_image(f_ori,RL);
figure,imshow(R9);