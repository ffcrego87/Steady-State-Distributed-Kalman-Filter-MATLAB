function [xhat,Omeganextcell] = DKF(xpred,Pi,y,yidxs,Omegacell,Cglob,V)
NAg = length(Omegacell);
nx = length(xpred)/NAg;
xhat = xpred;
Omeganextcell = Omegacell;
for i=1:NAg
    Cloc = Cglob(yidxs{i},:);
    Vloc = V(yidxs{i},yidxs{i});
    Omeganext = Cloc'*Vloc*Cloc;
    idxsx = (i-1)*nx+(1:nx);
    qloc = Cloc'*Vloc*y(yidxs{i});
    for j=1:NAg
        idxsxn = (j-1)*nx+(1:nx);
        qloc=qloc+Pi(i,j)*Omegacell{j}*xpred(idxsxn);
        Omeganext=Omeganext+Pi(i,j)*Omegacell{j};
    end
    xhat(idxsx) = Omeganext^(-1)*qloc;
    Omeganextcell{i}=Omeganext;
end
end