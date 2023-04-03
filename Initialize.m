function [W, H, S] = Initialize(V, paras,P, MAXITER)
%cosine distance
F = size(V,1);
T = size(V,2);
K = paras(1);
r = paras(2);
la1 = paras(3);  

rand('seed',0)

W = 1+rand(F, K);

H = 1+rand(K, T);

%P = zeros(F, T);

%P(V > 0) = 1;

for i=1:MAXITER
    H = H .* (W'*(P.*V))./(W'*(P.*(W*H))+ la1*H +eps);%NMF with miss data
    W = W .* (P.*V*H')./((P.*(W*H))*H' + la1*W +eps);

end
S = squareform(1-pdist(W,'cosine'));
%S = squareform(1./(1+pdist(W,'seuclidean')));
    
B=zeros(size(S));

for j=1:size(S,1)
    [val,index] = maxk(S(j,:),r);
    B(j,index(1:r))=val(1:r);
    %B(j,index(1:r))=ones(1,r);
end

S=B;    
S = max(S,S');
end