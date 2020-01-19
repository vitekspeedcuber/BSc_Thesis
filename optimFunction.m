% принимает 
% точку на орбите (параметр на отрезке [0; T]) в первой системе 3-ёх тел, 
% время полёта вдоль неустойчивого многообразия, 
% точку на орбите (параметр на отрезке [0; T]) во второй системе 3-ёх тел, 
% время полёта вдоль устойчивого многообразия, 
% время перелёта между точками на многообразиях, 
% плоские орбиты двух систем

% возвращает функционал: сумма характеристических скоростей

function J = optimFunction(parameters, muE, xLE, xE, tE, muV, xLV, xV, tV, InitStateE, InitStateV)

timeCoefficient = 365.25/224.7; % T_Earth / T_Ven
distanceCoefficient = 108200000/149600000; % R_Ven / R_Earth
velocityCoefficient = distanceCoefficient*timeCoefficient;

t01 = parameters(1);
time1 = parameters(2);
um = trajectoryUnstableManifold(muE, xLE, xE, tE, t01, time1);
UM = rotating2inertial(muE, um, t01, time1, InitStateE);

% nRows = size(UM, 1);
% for j = 1:nRows
%     UM(j, 1:3) = (InitStateE*UM(j, 1:3)')';
%     UM(j, 4:6) = (InitStateE*UM(j, 4:6)')';
% end

t02 = parameters(3);
time2 = parameters(4);
sm = trajectoryStableManifold(muV, xLV, xV, tV, t02, time2);
sm = sm(end:-1:1, :);

Tbetween = parameters(5);
SM = rotating2inertial(muV, sm, (t02 - time2)*timeCoefficient,...
    time2*timeCoefficient, InitStateV);

% nRows = size(SM, 1);
% for j = 1:nRows
%     SM(j, 1:3) = distanceCoefficient*(InitStateV*SM(j, 1:3)')';
%     SM(j, 4:6) = velocityCoefficient*(InitStateV*SM(j, 4:6)')';
% end
nRows = size(SM, 1);
for j = 1:nRows
    SM(j, 1:3) = distanceCoefficient*SM(j, 1:3);
    SM(j, 4:6) = velocityCoefficient*SM(j, 4:6);
end

muS = 1 - muE;

r1 = UM(end, 1:3)';
v1 = UM(end, 4:6)';
r2 = SM(1, 1:3)';
v2 = SM(1, 4:6)';

[v1corr, v2corr, ~, ~] = lambert(r1, r2, Tbetween, 0, muS);

% plot(UM(:, 1), UM(:, 2), 'LineWidth', 2, 'color', 'c');
% axis equal;
% grid on;
% hold on;
% quiver(UM(end, 1),UM(end, 2),UM(end, 4),UM(end, 5));
% hold on;
% plot(SM(:, 1), SM(:, 2), 'LineWidth', 2, 'color', 'c');
% axis equal;
% grid on;
% hold on;
% plot(SM(1, 1), SM(1, 2), 'ko', 'LineWidth', 3, 'color', 'y');
% hold on;
% quiver(SM(1, 1),SM(1, 2),SM(1, 4),SM(1, 5));
% hold on;
% plot(SM(end, 1), SM(end, 2), 'ko', 'LineWidth', 3, 'color', 'b');
% hold on;
% quiver(SM(end, 1),SM(end, 2),SM(end, 4),SM(end, 5));
% hold on;
% [~, X] = ode45(@(t, x)func2bp(t, x, muS), [0, Tbetween], [r1; v1corr]);
% plot(X(:, 1), X(:, 2), 'LineWidth', 2, 'color', 'm');
% axis equal;
% grid on;
% hold on;
% quiver(UM(end, 1),UM(end, 2),v1corr(1), v1corr(2));
% hold on;
% quiver(SM(1, 1),SM(1, 2),v2corr(1), v2corr(2));
% hold on;

J = norm(v1corr - v1) + norm(v2corr - v2);

end
