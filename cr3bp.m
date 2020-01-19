% по параметрам орбиты и оценке половины периода
% возвращает орбиту (траекторию) и период орбиты

function [x, T] = cr3bp(orbit, T)

initialX = initialValue(orbit, T);

mu = orbit.mu;
%xL = orbit.xL;

x0 = initialX(1:6);
T = initialX(7);
tspan = [0, T];

opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~, x] = ode113(@(t, x)func3bp(t, x, mu), tspan, x0,opts);
 
end
