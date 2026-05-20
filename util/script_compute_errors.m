x_real = logs{1,1}.stateTrajectory;
errors = cell(NAg,1);
nError = zeros(1,size(x_real,2));
nEBound = zeros(1,size(x_real,2));
error = zeros(nx*NAg,size(x_real,2));

for i=1:NAg
    errors{i} = logs{1,1+i}.inputTrajectory((size(C{i},1)+1):end,:)-x_real;
    error((i-1)*nx+(1:nx),:) = errors{i};
end

k1 = 1+alpha^Iter*cf;
k2 = alpha^Iter*k6*(1+k1*Phi_tilde/(beta-bb))*a/(2^nb);
k3 = k1/(1-bt);
k4 = alpha^Iter*k6*(1+k1*Phi_tilde/(1-bb))*b/(2^nb);

for k=1:size(x_real,2)
    nError(k) = norm(error(:,k));
    nEBound(k) = beta^(k-1)*(k1*e_max+k2)+k3*eps/(1-bt)+k4;
end

figure
plot(1:size(x_real,2),nError,1:size(x_real,2),nEBound)
set(gca,'yscale','log')
title('Norm of vector of estimation errors')
legend('Norm of vector of estimation errors','Bound on norm of vector of estimation errors')
xlabel('time')
setNicePlot

disp('number of iterations of consensus')
disp(Iter)
disp(' ')

disp('number of bits transmitted')
disp(nb)
disp(' ')