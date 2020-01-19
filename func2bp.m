% функция правых частей уравнений движения в задаче двух тел

function dxdt = func2bp(t, x, mu)

r = x(1:3);
v = x(4:6);

dxdt = [v; -mu*r/(norm(r))^3];

end
