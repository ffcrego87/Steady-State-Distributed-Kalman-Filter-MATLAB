function Pi = metropolis( Adj )
%metropolis Computes a consensus matrix with metropolis weights with
%sparsity given by the symmetric adjacency matrix Adj.

n = size(Adj); % number of nodes
d = sum(Adj);  % node degree

Pi = zeros(n);

for i = 1:n
    for j=(i+1):n
        if Adj(i,j)
            Pi(i,j)=1/(1+max(d(i),d(j)));
            Pi(j,i)=Pi(i,j);
        end
    end
    Pi(i,i)=1-sum(Pi(i,:));
end

end

