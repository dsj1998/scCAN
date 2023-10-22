function [X, Q, H, S] = adaptive_nmf(Init, paras)

MAXITER = 1000;
err = 1e-5;


%---------------init--------------------------
V = Init{1};
Q = Init{2};
H = Init{3};
S = Init{4};
D = diag(sum(S));

P = Init{5};

r = paras(2);
la1 = paras(3);
la2 = paras(4);
%---------------------------------------------

D1 = diag(sum(S).^-0.5);
L0 = norm(P.*(V-Q*H),'fro').^2 + la1.*(norm(Q,'fro').^2 + norm(H,'fro').^2)+ la2.*trace(Q'*(D1*(D-S)*D1)*Q);
for i=1:MAXITER
    Q = Q .* (P.*V*H' + la2*D1*S*D1*Q)./((P.*(Q*H))*H' + la1*Q + la2*D1*D*D1*Q + eps);

    H = H .* (Q'*(P.*V))./(Q'*(P.*(Q*H))+ la1*H +eps);
    
    error = mean(mean(abs(V-Q*H)))/mean(mean(V));

     if error< err
            break;
    end     
    if i > (MAXITER/2)
        S = squareform(1-pdist(Q,'cosine'));
        B=zeros(size(S));
        
       for j=1:size(S,1)%
            [val,index] = maxk(S(j,:),r);
           
            B(j,index(1:r))=val(1:r);
          
        end
        S=B;    
        S = max(S,S');
        D = diag(sum(S));
        D1 = diag(sum(S).^-0.5);
    end

end

X = (1-P).*(Q*H) + V;

end
