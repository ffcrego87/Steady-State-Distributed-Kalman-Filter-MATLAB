function stats = mc_convergence_P(M, tol, maxit, seed)
% Monte Carlo para estimar iterações até convergir da covariância global P
% Usa DTKF.m e a mesma construção de matrizes do setup_script.m.
%
% Inputs:
%   M     - nº de ensaios Monte Carlo
%   tol   - tolerância para convergência (norma Frobenius)
%   maxit - nº máximo de iterações
%   seed  - seed base (opcional)
%
% Output:
%   stats - struct com resultados e estatísticas

if nargin < 1, M = 50; end
if nargin < 2, tol = 1e-8; end
if nargin < 3, maxit = 5000; end
if nargin < 4, seed = 1; end

rng(seed);

iters = zeros(M,1);
finalDiff = zeros(M,1);

for run = 1:M
    disp(run)
    % --- Carrega parâmetros base (do teu problem.m)
    problem; %#ok<*NASGU>
    % (problem.m define NAg, dim, Amag, Cmag, x0obs, Circ_order, etc.)

    % --- Randomização do "setup" (isto é o que torna o MC relevante)
    % Se quiseres desligar randomização: comenta as 2 linhas seguintes.
    Aloc = eye(dim) + Amag*randn(dim,dim);
    Cloc = eye(dim) + Cmag*randn(dim,dim);

    % --- Dimensões globais
    nx = NAg*dim;
    ny = dim;

    % --- Gera rede circular (equivalente ao gen_networks.m, mas sem depender de metropolis externo)
    Adj = eye(NAg);
    for k = 1:Circ_order
        Adj = Adj + diag(ones(NAg-k,1), k) + diag(ones(NAg-k,1), -k) ...
                  + diag(ones(k,1), NAg-k) + diag(ones(k,1), k-NAg);
    end
    Pi = metropolis_weights(Adj);

    % --- Covariâncias (como no setup_script.m)
    W = eye(nx);                               % informação do processo (no setup: W=eye)
    V = eye((2*NAg-1)*ny);                     % informação da medida (no setup: V=eye)

    % --- Constrói A global
    Aglob = kron(eye(NAg), Aloc);

    % --- Constrói Cglob (como no setup_script.m)
    % Cglob: (2*NAg-1)*ny x nx
    C = cell(NAg,1);

    G = zeros(2,NAg);
    G(1,1)=1; G(2,1)=1; G(1,2)=-1; G(2,NAg)=-1;
    C{1} = kron(G, Cloc);

    for i = 2:(NAg-1)
        G = zeros(2,NAg);
        G(1,i)=1; G(2,i)=1;
        G(1,i+1)=-1;
        G(2,i-1)=-1;
        C{i} = kron(G, Cloc);
    end

    G = zeros(1,NAg);
    G(NAg)=1;
    C{NAg} = kron(G, Cloc);

    Cglob = cell2mat(C);

    % --- yidxs (como no setup_script.m)
    yidxs = cell(NAg,1);
    yidxs{1} = 1:(2*ny);
    for i = 2:(NAg-1)
        yidxs{i} = (i-1)*2*ny + (1:(2*ny));
    end
    yidxs{NAg} = (NAg-1)*2*ny + (1:ny);

    % --- idxs_adj, ones_i, eta_i, Gamma_ij (como no setup_script.m)
    idxs_adj = cell(NAg,1);
    ones_i   = cell(NAg,1);
    eta_i    = cell(NAg,1);
    Gamma_ij = cell(NAg,NAg);

    for i = 1:NAg
        idxs_adj{i} = find(Adj(i,:));
        ones_i{i}   = kron(ones(size(idxs_adj{i}))', eye(nx));

        eta_i{i} = zeros(numel(idxs_adj{i})*nx, NAg*nx);
        for jj = 1:numel(idxs_adj{i})
            e_j = zeros(NAg,1);
            e_j(idxs_adj{i}(jj)) = 1;
            eta_i{i}((jj-1)*nx+(1:nx), :) = kron(e_j', eye(nx));
        end

        for j = 1:NAg
            e_j = zeros(NAg,1);
            e_j(j) = 1;
            Gamma_ij{i,j} = eta_i{i} * kron(e_j, eye(nx));
        end
    end

    % --- Iteração de covariância global (mesma ideia do setup_script)
    P = (1/x0obs)^-2 * eye(NAg*nx);
    P_old = 1e10 * P;

    k = 0;
    y = zeros((2*NAg-1)*ny,1);
    xpred = zeros(nx*NAg,1);

    while (norm(P - P_old, 'fro') > tol) && (k < maxit)
        P_old = P;
        [~, P] = DTKF(xpred, Gamma_ij, ones_i, eta_i, y, yidxs, P, Aglob, Cglob, V, W^-1);
        k = k + 1;
        disp(k)
    end

    iters(run) = k;
    finalDiff(run) = norm(P - P_old, 'fro');
end

% --- Estatísticas
stats.iters = iters;
stats.finalDiff = finalDiff;
stats.mean = mean(iters);
stats.median = median(iters);
stats.p10 = prctile(iters,10);
stats.p90 = prctile(iters,90);
stats.min = min(iters);
stats.max = max(iters);
stats.tol = tol;
stats.maxit = maxit;
stats.M = M;

fprintf('\nMC convergence of global P (tol=%g, maxit=%d, M=%d)\n', tol, maxit, M);
fprintf('iters: mean=%.2f, median=%d, p10=%d, p90=%d, min=%d, max=%d\n', ...
    stats.mean, stats.median, stats.p10, stats.p90, stats.min, stats.max);
fprintf('finalDiff: mean=%g, max=%g\n\n', mean(finalDiff), max(finalDiff));

end

function Pi = metropolis_weights(Adj)
% Metropolis weights para uma matriz de adjacência (inclui self-loops em Adj)
N = size(Adj,1);
deg = sum(Adj>0,2);

Pi = zeros(N,N);
for i = 1:N
    for j = 1:N
        if i ~= j && Adj(i,j) > 0
            Pi(i,j) = 1/(1 + max(deg(i), deg(j)));
        end
    end
    Pi(i,i) = 1 - sum(Pi(i,:));
end
end
