function B_new = updateB_capped_robust(X,A,B)

[m,~] = size(B);
Bup = X * A';           % m*c
Bdown = sum(A,2);       % c*1
B_new = Bup./repmat(Bdown',m,1); % m*c
