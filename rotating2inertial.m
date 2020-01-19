% по значениям фазового вектора (траектории) в ВСК
% возвращает значения фазового вектора (траекторию) в ИСК

function X = rotating2inertial(mu, x, t0, T, InitState)

nRows = size(x, 1);
X = zeros(nRows, 6);

flag = 0;
if (nRows == 1)
    flag = 1;
end

for j = 1:nRows
    if (flag)
        t = t0;
    else
        t = t0 + T*(j - 1)/(nRows - 1);
    end
    
    r = x(j, 1:3)' + [mu; 0; 0]; % начало ИСК в местоположении Солнца
    v = x(j, 4:6)';
    
    A = [cos(t), -sin(t), 0;
         sin(t), cos(t), 0;
         0, 0, 1];
    
    R = A*InitState*r;
    V = A*InitState*(v + [-r(2); r(1); 0]);
    
    X(j, 1:6) = [R; V]';
end

end
