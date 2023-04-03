function MP = findMP(data,cn,threshold)
MP = zeros(size(data,1),size(data,2));

S = squareform(1-pdist(data,'cosine'));
%S = squareform(1./(1+pdist(data,'seuclidean')));

D = diag(sum(S));
D1 = diag(sum(S).^-0.5);
[F,~] = eigs(D1*(D-S)*D1,cn,'sm');

idx = kmeans(F,cn);
%idx = kmeans(data,cn);

for i = 1:cn
    id_i = find(idx==i);
    X_i = data(id_i,:);
    P = zeros(size(X_i,1),size(X_i,2));
    P(X_i>0) = 1;
    p_sum = sum(P);
    gene_idx = find(p_sum>(size(P,1)*threshold));
    gene_idx_full = (1:size(X_i,2));
    b = gene_idx_full(ismember(gene_idx_full,gene_idx)==0);
    P(:,b) = 1;
    MP(id_i,:) = P;
end
    
end