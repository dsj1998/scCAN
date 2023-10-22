function [MP] = findMP(data,cn,threshold)

MP = zeros(size(data,1),size(data,2));


S = squareform(1-pdist(data,'cosine'));
D = diag(sum(S));
D1 = diag(sum(S).^-0.5);
[F,~] = eigs(D1*(D-S)*D1,cn,'smallestreal');

idx = kmeans(F,cn);

for i = 1:cn

    id_i = find(idx==i);
    X_i = data(id_i,:);

    P = zeros(size(X_i,1),size(X_i,2));

    P(X_i>0) = 1;
    p_sum = sum(P);

    gene_zero_rate=1-(p_sum/size(X_i,1));
    gene_mean=mean(X_i,1);
    gene_var = var(X_i, [], 1);
    

    confidence_score = (1-gene_zero_rate).*gene_mean./ (((1-gene_zero_rate).*gene_mean)+gene_zero_rate.*gene_var + eps); 

 
    confidence_score_idx=confidence_score<threshold;
    P(:,confidence_score_idx) = 1;
    MP(id_i,:) = P;
end
    
end
