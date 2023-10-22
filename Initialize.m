function [Q, H, S] = Initialize(V, paras,P, MAXITER)

N = size(V,1); 
G = size(V,2);
K = paras(1);
r = paras(2);
la1 = paras(3);  

rand('seed',0)

Q = 1+rand(N, K);
H = 1+rand(K, G);

    for i=1:MAXITER
   
        H = H .* (Q'*(P.*V))./(Q'*(P.*(Q*H))+ la1*H +eps);
        Q = Q .* (P.*V*H')./((P.*(Q*H))*H' + la1*Q +eps);

    end

S = squareform(1-pdist(Q,'cosine'));

B=zeros(size(S));

    for j=1:size(S,1)
        [val,index] = maxk(S(j,:),r);
    
        B(j,index(1:r))=val(1:r);
    end
S=B;
S = max(S,S');
end
