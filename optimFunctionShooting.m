% возвращает функционал: сумма характеристических скоростей

function [J, g] = optimFunctionShooting(s)

dv1 = s(end - 5:end - 3);
dv2 = s(end - 2:end);

J = norm(dv1) + norm(dv2);

N = size(s, 1);
g = zeros(N, 1);

g(end - 5:end - 3) = dv1/norm(dv1);
g(end - 2:end) = dv2/norm(dv2);

end
