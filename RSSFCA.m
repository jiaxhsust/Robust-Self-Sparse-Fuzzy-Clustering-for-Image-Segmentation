function [outA,outB,outObj,outNumIter] = RSSFCA(X,gama,maxIter,cluster_n)
% A (cxn) is the membership matrix.
% B (mxc) is the cluster center matrix.
c=cluster_n;
n=size(X,2);
thresh = 10^-5;
obj = zeros(maxIter,1);
[dim,data_n]=size(X);
options = [2;	% exponent for the partition matrix U
    100;	% max. number of iteration
    1e-5;	% min. amount of improvement
    0];	% info display during iteration
[center, A] = fcm(X', c, options); % Initialize parameters
B=center';
GMM_WIDTH=1;
%% compute variance
for i = 1:c
    diffs = X'-(ones(data_n, 1) * B(:,i)');
    diffs = diffs.*(sqrt(A(i,:))'*ones(1, dim));
    covars(:,:,i)=(diffs'*diffs)/sum(A(i,:));
    try
        if rank(covars(:,:,i)) < dim
            covars(:,:,i) = covars(:,:,i) + GMM_WIDTH.*eye(dim);
        end
    catch
        covars(:,:,i) = covars(:,:,i) + GMM_WIDTH.*eye(dim);
    end
end

for i = 1:c
    diffs = X'-(ones(data_n, 1) * B(:,i)');
    mahalanobis2 = sum((diffs*inv(covars(:,:,i))) .* diffs,2);  
    D(:,i) =1/2.*(mahalanobis2+log(det(covars(:,:,i)))+dim.*log(2*pi));  
end


MIN_COVAR =1e-3;
init_covars=covars;
for t = 1:maxIter
    t
    %-------- update A when fixed B--------%
    A = updateA_capped_robust(X,B,gama,D); %A denotes degree of membership--c*n,
    
    %-------- update B when fixed A and lamda--------%
    B = updateB_capped_robust(X,A,B); %B denotes centers
    
    %------- calculate the covariance --------%
    for i = 1:c
        diffs = X'-(ones(data_n, 1) * B(:,i)');
        diffs = diffs.*(sqrt(A(i,:))'*ones(1, dim));
        covars(:,:,i)=(diffs'*diffs)/(sum(A(i,:))+eps);
    end
    
    P=sum(covars,3);
    Nan=sum(isnan(P(:)));
    if isreal(covars)==0 || Nan,break; end,
    for j = 1:c
        if min(svd(covars(:,:,j))) < MIN_COVAR
            covars(:,:,j) = init_covars(:,:,j);
        end
    end
    
    Pi=sum(A,2)./data_n;
    for i = 1:c
        diffs = X'-(ones(data_n, 1) * B(:,i)');
        mahalanobis2 = sum((diffs*inv(covars(:,:,i))) .* diffs,2);
        Temp(:,i) =1/2.*(mahalanobis2+log(det(covars(:,:,i)))+dim.*log(2*pi));  %
    end
    D=Temp;
    
    
    for ii=1:n
        AI=D(ii,:);
        if sum(AI<0)>1 ||sum(AI<0)==1
            TempT(ii,:)=AI-min(AI)+eps;
        else
            TempT(ii,:)=AI;
        end
    end
    obj(t)=sum(sum(TempT'.*A + A.^2*gama));
    
    CC{t}=B; % You can also use objective funtion as the convergence condition
    if(t > 1)
        diff = max(max(CC{t-1}-CC{t}));
        if(diff < thresh)
        end
    end
end

% figure,mesh(reshape(E,rows,cols))
outNumIter = t;
outObj = obj(1:t);
outA = A;
outB = B;