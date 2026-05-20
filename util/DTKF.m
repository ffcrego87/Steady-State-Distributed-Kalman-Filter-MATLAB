function [xhat,Pnext] = DTKF(xpred,Gamma_ij,ones_i,eta_i,y,yidxs,P,A,Cglob,V,Q)
NAg = length(eta_i);
nx = length(xpred)/NAg;

T = zeros(nx*NAg);
Diag_P = zeros(nx*NAg);

xhat=zeros(NAg*nx,1);
for i=1:NAg
    Cloc = Cglob(yidxs{i},:);
    Vloc = V(yidxs{i},yidxs{i});
    Si = Cloc'*Vloc*Cloc;
    P_aux = ones_i{i}'*pinv(eta_i{i}*P*eta_i{i}');
    AP = A*(P_aux*ones_i{i}+Si)^-1;
    
    idxsx = (i-1)*nx+(1:nx);

    Diag_P(idxsx,idxsx)=AP*Si*AP';

    qloc = Cloc'*Vloc*y(yidxs{i});
    for j=1:NAg
        idxsxn = (j-1)*nx+(1:nx);
        qloc=qloc+P_aux*Gamma_ij{i,j}*xpred(idxsxn);
        
        T(idxsx,idxsxn)=AP*P_aux*Gamma_ij{i,j};
    end
    xhat(idxsx) = AP*qloc;
end

Pnext = T*P*T'+Diag_P+kron(ones(NAg),Q);

end