function sim_output = sim_sample(sim_input)

ifinal = sim_input.ifinal;
NAg = sim_input.NAg;
nx = sim_input.nx;
ny = sim_input.ny;
A = sim_input.A;
Pi = sim_input.Pi;
Cglob = sim_input.Cglob;
yidxs = sim_input.yidxs;
wmag = sim_input.wmag;
vmag = sim_input.vmag;
x0mag = sim_input.x0mag;
wobs = sim_input.wobs;
vobs = sim_input.vobs;
x0obs = sim_input.x0obs;
E_cell = sim_input.E_cell;
D_cell = sim_input.D_cell;
Delta_cell = sim_input.Delta_cell;
kfix = sim_input.kfix;
beta = sim_input.beta;
alpha = sim_input.alpha;
check_equivalence = sim_input.check_equivalence;
if check_equivalence
    Omegabar = sim_input.Omegabar;
end
Gamma_ij = sim_input.Gamma_ij;
eta_i = sim_input.eta_i;
ones_i = sim_input.ones_i;
P_s = sim_input.P_s;
    
% initialization
Cglobt = Cglob{1};
xdhat=zeros(NAg*nx,1);
xdghat=zeros(NAg*nx,1);
xdsskfhat=zeros(NAg*nx,1);
xddkfhat=zeros(NAg*nx,1);
xhat=zeros(nx,1);
Omega = (1/x0obs)^2*eye(nx);
W = (1/wobs)^2*eye(nx);
V = (1/vobs)^2*eye(size(Cglob{1},1));
x=x0mag*randn(nx,1);
w=wmag*randn(nx,1);
v=vmag*randn((2*NAg-1)*ny,1);
y=Cglobt*x+v;

P=(1/x0obs)^-2*eye(NAg*nx);

T_cell = cell(NAg,1);
for i=1:NAg
    T_cell{i} = zeros(size(Cglobt,1),nx);
end
O_cell = cell(NAg,1);
for i=1:NAg
    O_cell{i} = zeros(size(Cglob{1},1)*(kfix+1),nx);
end
Omegacell = cell(NAg,1);
for i=1:NAg
    Omegacell{i} = (1/x0obs)^2*eye(nx);
end
OmegacellGramian = cell(NAg,1);
for i=1:NAg
    OmegacellGramian{i} = eye(nx);
end

Chi = alpha*eye(nx);

% log variables
xdhatlog=zeros(NAg*nx,ifinal);
xdghatlog=zeros(NAg*nx,ifinal);
xddkfhatlog=zeros(NAg*nx,ifinal);
xdsskfhatlog=zeros(NAg*nx,ifinal);
xhatlog=zeros(nx,ifinal);
xlog=zeros(nx,ifinal);
wlog=zeros(nx,ifinal);
vlog=zeros((2*NAg-1)*ny,ifinal);
ylog=zeros((2*NAg-1)*ny,ifinal);

xdhatlog(:,1)=xdhat;
xdghatlog(:,1)=xdghat;
xddkfhatlog(:,1)=xddkfhat;
xdsskfhatlog(:,1)=xdsskfhat;
xhatlog(:,1)=xhat;
xlog(:,1)=x;
wlog(:,1)=w;
vlog(:,1)=v;
ylog(:,1)=y;

%cicle
for i=2:ifinal
    At = A{i};
    
    % Distributed observer
    [xdhat,Omegacell]=DKF(xdhat,Pi,y,yidxs,Omegacell,Cglobt,V);
    xdhat=kron(eye(NAg),At)*xdhat;
    for j=1:NAg
        Omegacell{j}=W-W*At*(Omegacell{j}+At'*W*At)^(-1)*At'*W;
    end
    
    % Distributed observer gramian
    Chi = beta*(At^(-1))'*Chi*At^(-1);
    T_prev_cell = T_cell;
    O_prev_cell = O_cell;
    for j = 1:NAg
        % T update
        T_cell{j}=0*T_cell{j};
        for k=1:NAg
            if j==k
                T_cell{j}=T_cell{j}+D_cell{j,k}*Cglobt*At^(-1);
            else
                T_cell{j}=T_cell{j}+D_cell{j,k}*T_prev_cell{k}*At^(-1);
            end
        end
        
        % O update
        O_update = [zeros(size(Cglobt));O_prev_cell{j}]*At^(-1);
        O_cell{j} = O_update(1:(size(Cglob{1},1)*(kfix+1)),:)+E_cell{j}*T_cell{j};
        
        OmegacellGramian{j}=O_prev_cell{j}'*Delta_cell{j}*O_prev_cell{j}+Chi;
        if check_equivalence
            if max(max(abs(OmegacellGramian{j}-Omegabar{j,i})))>10^(-5)
                fprintf('Disagreement at iteration %i, agent %i\n',i,j)
            end
        end
    end
    [xdghat,~]=DKF(xdghat,Pi,y,yidxs,OmegacellGramian,Cglobt,eye(size(Cglobt,1)));
    xdghat=kron(eye(NAg),At)*xdghat;

    % New Distributed observer
    [xddkfhat,P]=DTKF(xddkfhat,Gamma_ij,ones_i,eta_i,y,yidxs,P,At,Cglobt,V,W^-1);

    % New Steady State Distributed observer
    [xdsskfhat,~]=DTKF(xdsskfhat,Gamma_ij,ones_i,eta_i,y,yidxs,P_s,At,Cglobt,V,W^-1);
    disp(norm(P-P_s))

    % Centralized observer
    xhat=(Omega+Cglobt'*V*Cglobt)^(-1)*(Omega*x+Cglobt'*V*y);
    Omega=Omega+Cglobt'*V*Cglobt;
    xhat=At*xhat;
    Omega=W-W*At*(Omega+At'*W*At)^(-1)*At'*W;
    
    %Real system
    x=At*x+w;
    w=wmag*randn(nx,1);
    v=vmag*randn((2*NAg-1)*ny,1);
    Cglobt = Cglob{i};
    y=Cglobt*x+v;
    
    %log variables
    xdhatlog(:,i)=xdhat;
    xdghatlog(:,i)=xdghat;
    xddkfhatlog(:,i)=xddkfhat;
    xdsskfhatlog(:,i)=xdsskfhat;
    xhatlog(:,i)=xhat;
    xlog(:,i)=x;
    wlog(:,i)=w;
    vlog(:,i)=v;
    ylog(:,i)=y;
end

sim_output = struct;

sim_output.xdhatlog=xdhatlog;
sim_output.xdghatlog=xdghatlog;
sim_output.xddkfhatlog=xddkfhatlog;
sim_output.xdsskfhatlog=xdsskfhatlog;
sim_output.xhatlog=xhatlog;
sim_output.xlog=xlog;
sim_output.wlog=wlog;
sim_output.vlog=vlog;
sim_output.ylog=ylog;

end