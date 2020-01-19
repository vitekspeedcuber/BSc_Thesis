% по параметрам орбиты возвращает начальное значение  фазового вектора
% для интегрирования уравнений движения для плоской орбиты 
% вокруг коллинеарной точки либрации и 
% начальное значение для периода орбиты

function initialX = initialValue(orbit, T)

x = linearApproximationVector(orbit, T);
x0 = x(1);
y0 = x(2);
z0 = x(3);
Vx0 = x(4);
Vy0 = x(5);
Vz0 = x(6);

optimParameters0 = [Vy0, T];
options=optimset;
options.TolFun=1e-10;
optimParameters = fsolve(@(optimParameters)initialApproximation(optimParameters, ... 
    orbit.mu, x0, y0, z0, Vx0, Vz0), optimParameters0, options);

Vy0 = optimParameters(1);
T = 2*optimParameters(2);

initialX = [x0; y0; z0; Vx0; Vy0; Vz0; T];

end
