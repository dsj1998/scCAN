function [X] = adaptive_nmf(Init, paras)
%Cosine distance

MAXITER = 1000;
e = 1e-8;
%---------------init--------------------------

V = Init{1};
W = Init{2};
H = Init{3};
S = Init{4};
D = diag(sum(S));
% P = zeros(size(V,1), size(V,2));
% P(V > 0) = 1;
P = Init{5};

r = paras(2);
la1 = paras(3);
la2 = paras(4);
%---------------------------------------------


D1 = diag(sum(S).^-0.5);

L0 = norm(P.*(V-W*H),'fro').^2 + la1.*(norm(W,'fro').^2 + norm(H,'fro').^2)+ la2.*trace(W'*(D1*(D-S)*D1)*W);
  
for i=1:MAXITER
    W = W .* (P.*V*H' + la2*D1*S*D1*W)./((P.*(W*H))*H' + la1*W + la2*D1*D*D1*W + eps);
%     for j=1:size(W,1)
%         W(j,:) = W(j,:)./norm(W(j,:));
%     end
    H = H .* (W'*(P.*V))./(W'*(P.*(W*H))+ la1*H +eps);
         
    L1 = norm(P.*(V-W*H),'fro').^2 + la1.*(norm(W,'fro').^2 + norm(H,'fro').^2)+ la2.*trace(W'*(D1*(D-S)*D1)*W);
    LL = abs(L0-L1)/L1;
            if LL<e
                break
            end
    L0=L1;
    Z(i) = L1;
    
    if i > (MAXITER/2)
        
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
        D = diag(sum(S));
        D1 = diag(sum(S).^-0.5);
    end

end

X = (1-P).*(W*H) + V;

plot(Z); 
end