% возвращает приближения для решения уравнения F = 0

function F = initialApproximation(optimParameters, mu, x0, y0, z0, Vx0, Vz0)

Vy0 = optimParameters(1);
T = optimParameters(2);

x0 = [x0; y0; z0; Vx0; Vy0; Vz0];
tspan = [0, T];

opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~, X] = ode113(@(t, x)func3bp(t, x, mu), tspan, x0, opts);

%x = X(end, 1);
y = X(end, 2);
%z = X(end, 3);
Vx = X(end, 4);
%Vy = X(end, 5);
%Vz = X(end, 6);

F = [y; Vx];

end
