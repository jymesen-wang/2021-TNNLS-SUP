function [P,F] = GOSUP(X,S,d1,k,lambda)


[d,n] = size(X);
H  = diag(ones(1,n)) - (1/n)*ones(n); 
St =  X*H*X';                         
P  = orth(rand(d,d1));
b  = 10;
K  = 15;


for iter = 1:30
    D = diag(sum(S,2));
    L = D - S;
    
    
    B = lambda*inv(L+lambda*eye(n));
    F = B*X'*P;
    
    
    C     = X*(B'*L*B)*X'+lambda*(X*X'-2*X*B*X'+X*(B'*B)*X');
    a     = trace(P'*C*P)/(trace(P'*St*P)+eps);
    A     = -C + a*St + b*eye(d);
    Eva_A = min(eig(A));
    while Eva_A < 0              
        b     = b*10;
        A     = -C + a*St + b*eye(d);
        Eva_A = min(eig(A));
    end
    A       = (A' + A)/2;          
    P       = A*(P*pinv(P'*A*P)*P')*A;
    Tr_P    = diag(P);             
    [~,Idx] = sort(Tr_P);         
    
    Idx_1   = Idx(d-k+1:d);        
    if k ~= d
        Idx_0      = Idx(1:d-k);   
        P(Idx_0,:) = 0;            
    end
    
    M       = A*P;                
    M1      = M(Idx_1,:);          
    M2      = orth(M1);           
    P       = zeros(d,d1);
    if size(M2,2) > d1
        M2 = M2(:,1:d1);
    elseif size(M2,2) < d1
        M2 = [M2 zeros(k,d1-size(M2,2))];
    end
    for i = 1:k
        P(Idx_1(i),:) = M2(i,:);   
    end
    distFF = Eu2_distance(F',F');
    distFF_sort = sort(distFF,2);  
    g = (K*sum(distFF_sort(:,K+1)) - sum(sum(distFF_sort(:,1:K),1),2))/n + eps;
    for j = 1:n
        ad = -distFF(j,:)/g;
        S(j,:) = EProjSimplex_new(ad);  
    end
end
end