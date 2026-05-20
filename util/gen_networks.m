 
%% Network

Adj = eye(NAg);
for i=1:Circ_order
    Adj = Adj+diag(ones(NAg-i,1),i)+diag(ones(NAg-i,1),-i)+diag(ones(i,1),NAg-i)+diag(ones(i,1),i-NAg);
end
Pi = metropolis(Adj);
eigPi = sort(abs(eig(Pi)),'descend');
sigma2 = eigPi(2);