function A = updateA_capped_robust(X,B,gama,D)

[~,n] = size(X);  %n denotes the number of samples
[~,c] = size(B);  %c denotes the number of clusters
A = zeros(c,n);   % reserve memory

for i = 1:n
    AI=D(i,:);
    if sum(AI<0)>1 ||sum(AI<0)==1
        Temp=AI-min(AI)+eps;
    else
        Temp=AI;
    end
    D1=sqrt(Temp);
    dnew = - D1/(sqrt(2*gama));
    [anew,~] = EProjSimplex_new(dnew,sqrt(2*gama));
    A(:,i) = (anew./sqrt(2*gama));
end
