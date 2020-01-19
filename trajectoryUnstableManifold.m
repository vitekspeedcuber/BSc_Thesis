% по орбите (траектории), точке на орбите (траектории) и 
% времени полёта вдоль многообразия
% возвращает траекторию вдоль неустойчивого многообразия

function xU = trajectoryUnstableManifold(mu, xL, X, T, t0, time)

n = size(X, 1);

i = fix(t0*(n - 1)/T + 1);
alpha = 1 - (t0*(n - 1)/T + 1 - i);

x0 = [alpha*X(rem(i - 1, n) + 1, 1:6)' + (1 - alpha)*X(rem(i, n) + 1, 1:6)';
      reshape(eye(6), [], 1)];
tspan = [t0, T + t0];

opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~, Y] = ode113(@(t, x)func(t, x, mu), tspan, x0, opts);

F = reshape(Y(end, 7:42), 6, 6); % матрица монодромии
[V, D] = eig(F);
[~, index] = max(diag(D));  % индекс максимального собственного значения
maxVector = V(:,index);     % соответствующий этому значению вектор

if (maxVector(1) < 0)
    maxVector = (-1)*maxVector;
end

eps = 1e-8;
xU0 = x0 + [eps*maxVector/norm(maxVector); reshape(zeros(6), [], 1)];
tspanU = [t0, time + t0];

[~, xU] = ode113(@(t, x)func(t, x, mu), tspanU, xU0, opts);

end
