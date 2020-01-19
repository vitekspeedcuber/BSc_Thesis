function trajectory

timeCoefficient = 365.25/224.7; % T_Earth / T_Ven
distanceCoefficient = 108200000/149600000; % R_Ven / R_Earth
velocityCoefficient = distanceCoefficient*timeCoefficient;

DistUnit = 149600000;
muSun = 132712440041.939400; 
VelUnit = (muSun/DistUnit)^(1/2);
TimeUnit = DistUnit/VelUnit/3600/24;

DistUnitV = 108200000;
%muSun = 132712440041.939400; 
VelUnitV = (muSun/DistUnitV)^(1/2);
TimeUnitV = (DistUnitV/VelUnitV/3600/24)*timeCoefficient;


JD0 = juliandate(2044,7,2,0,0,0);

% Используются RCS (вращающаяся система координат - ВСК) и 
% ICS (инерциальная система координат - ИСК)

% Рассматриваются плоские периодические орбиты 
% вблизи коллинеарных точек либрации

%% Система Солнце-Земля

orbitSE.mu = 5.972/(5.972 + 1.989e6);
orbitSE.xL = 1.0100740; % L2

orbitSE.alpha = 0.23;
orbitSE.beta = 0;
orbitSE.gamma = 0;
orbitSE.gamma_dash = 0;
orbitSE.phi1 = 0;
orbitSE.phi2 = 0;

parameters = linearApproximationParameters(orbitSE.mu, orbitSE.xL);
orbitSE.lambda = parameters(1); 
orbitSE.wp = parameters(2); 
orbitSE.wv = parameters(3); 
orbitSE.k1 = parameters(4); 
orbitSE.k2 = parameters(5);
T = parameters(6);

[xE, tE] = cr3bp(orbitSE, T);

figure
Vx = [1.2; 0];
Vz = [0; 1];
Vy = [0.2; 0.8];
quiver(0,0,Vx(1),Vx(2), 'color', 'k', 'LineWidth', 1);
axis equal;
grid on;
hold on;
quiver(0,0,Vz(1),Vz(2), 'color', 'k', 'LineWidth', 1);
hold on;
quiver(0,0,Vy(1),Vy(2), 'color', 'k', 'LineWidth', 1);
hold on;
t = pi/6;
A = [cos(t), -sin(t);
     sin(t), cos(t)];
Vx = A*Vx;
Vy = A*Vy;
quiver(0,0,Vx(1),Vx(2), 'color', 'k', 'LineWidth', 1);
hold on;
quiver(0,0,Vy(1),Vy(2), 'color', 'k', 'LineWidth', 1);
plot(-orbitSE.mu, 0, 'ko', 'LineWidth', 1, 'color', 'y');
axis equal;
grid on;
hold on;
plot(1-orbitSE.mu, 0, 'ko', 'LineWidth', 1, 'color', 'b');
axis equal;
grid on;
hold on;
plot(orbitSE.xL, 0, '.', 'LineWidth', 1, 'color', 'k');
axis equal;
grid on;
hold on;
plot(xE(:,1), xE(:,2), 'color', 'r');
hold on;
xlabel('x, безразмер.');
ylabel('y, безразмер.');

%% Система Солнце-Венера

orbitSV.mu = 4.867/(4.867 + 1.989e6);
orbitSV.xL = 1.0093433; % L2

orbitSV.alpha = 0.002;
orbitSV.beta = 0;
orbitSV.gamma = 0;
orbitSV.gamma_dash = 0;
orbitSV.phi1 = 0;
orbitSV.phi2 = 0;

parameters = linearApproximationParameters(orbitSV.mu, orbitSV.xL);
orbitSV.lambda = parameters(1); 
orbitSV.wp = parameters(2); 
orbitSV.wv = parameters(3); 
orbitSV.k1 = parameters(4); 
orbitSV.k2 = parameters(5);
T = parameters(6);

[xV, tV] = cr3bp(orbitSV, T);
% figure
% plot(orbitSV.xL, 0, '*', 'LineWidth', 2, 'color', 'k');
% axis equal;
% grid on;
% hold on;
% plot(xV(:,1), xV(:,2), 'color', 'b');
% grid on;
% axis equal;
% hold on;
% text(orbitSV.xL, 0.0005, 'L_2');
% sm = trajectoryStableManifold(orbitSV.mu, orbitSV.xL, xV, tV, 0.2, 1.9*tV);
% sm = sm(end:-1:1, :);
% plot(sm(:,1), sm(:,2), 'color', 'c');
% xlabel('x, безразмер.');
% ylabel('y, безразмер.');

%% Построение траектории в простой модели

% поворот к действительному начальному положению Земли
[R1,~] = planetEphemeris(JD0,'Sun','Earth','430');
R1 = R1/DistUnit;
vec1 = [R1(1); R1(2)];
vec2 = [1; 0];
if (vec1(2) > 0)
    phi = acos((vec1/norm(vec1))'*(vec2/norm(vec2)));
else
    phi = -acos((vec1/norm(vec1))'*(vec2/norm(vec2)));
end
InitStateE = [cos(phi), -sin(phi), 0;
              sin(phi), cos(phi), 0;
              0, 0, 1];

% quiver(0,0,vec1(1),vec1(2));
% axis equal;
% grid on;
% hold on;
% поворот к действительному начальному положению Венеры
[R1,~] = planetEphemeris(JD0,'Sun','Venus','430');
R1 = R1/DistUnitV;
vec1 = [R1(1); R1(2)];
% quiver(0,0,vec1(1),vec1(2));
vec2 = [1; 0];
if (vec1(2) > 0)
    phi = acos((vec1/norm(vec1))'*(vec2/norm(vec2)));
else
    phi = -acos((vec1/norm(vec1))'*(vec2/norm(vec2)));
end
InitStateV = [cos(phi), -sin(phi), 0;
              sin(phi), cos(phi), 0;
              0, 0, 1];

t01 = 0.3;
time1 = 0.2;
time2 = 0.2;
t02 = 1.0;
Tbetween = -t01 - time1 + (t02 - time2)*timeCoefficient;
parameters0 = [t01; time1; t02; time2; Tbetween];
options = optimoptions('fmincon','Display','iter','Algorithm', 'sqp',...
    'OptimalityTolerance', 1e-10);

[parameters, J, ~, output] = fmincon(@(parameters)optimFunction(...
    parameters,...
    orbitSE.mu, orbitSE.xL, xE, tE, orbitSV.mu, orbitSV.xL, xV, tV, InitStateE, InitStateV),...
    parameters0, [], [], [-1, -1, 1*timeCoefficient, -1*timeCoefficient, -1], [0],...
    [0.1; 0.1; 0.1; 0.1; 0.1], [tE; tE; tV; tV; 2*pi], [], options);

disp(output.iterations)
t01 = parameters(1);
time1 = parameters(2);
t02 = parameters(3);
time2 = parameters(4);
Tbetween = parameters(5);

%% результаты в простой модели

fprintf('время полёта по плоской орбите в системе Солнце-Земля %f дней\n',...
    t01*TimeUnit);
fprintf('время полёта вдоль неустойчивого многообразия %f дней\n',...
    time1*TimeUnit);
fprintf('время двухимпульсного перелёта %f дней\n',...
    Tbetween*TimeUnit);
fprintf('время полёта вдоль устойчивого многообразия %f дней\n',...
    time2*TimeUnit);
fprintf('время полёта по плоской орбите в системе Солнце-Венера %f дней\n',...
    (tV - t02)*TimeUnit);
fprintf('характеристическая скорость %f км/с\n',...
    J*VelUnit);

%% вывод графиков в простой модели

% Солнце
p1 = plot(0, 0, 'ko', 'LineWidth', 10, 'color', 'y');
axis equal;
grid on;
hold on;

% плоская орбита в системе Солнце-Земля
TE = tE;
XE = rotating2inertial(orbitSE.mu, xE, 0, TE, InitStateE);
% nRows = size(XE, 1);
% for j = 1:nRows
%    XE(j, 1:3) = (InitStateE*XE(j, 1:3)')';
%    XE(j, 4:6) = (InitStateE*XE(j, 4:6)')';
% end

p2 = plot(XE(:, 1), XE(:, 2), '--', 'LineWidth', 2, 'color', 'r');
axis equal;
grid on;
hold on;

% плоская орбита в системе Солнце-Венера
TV = timeCoefficient*tV;
XV = rotating2inertial(orbitSV.mu, xV, 0, TV, InitStateV);
% nRows = size(xV, 1);
% for j = 1:nRows
%    xV(j, 1) = xV(j, 1) + orbitSV.mu;
% end
% tspanDays = linspace(0,tV, nRows)*TimeUnitV;
% XV = rotating2ICCS(xV(:,1:6)', JD0, tspanDays, 'Venus', 'Sun', DistUnitV, VelUnitV)';
% for j = 1:nRows
%    xV(j, 1) = xV(j, 1) - orbitSV.mu;
% end
% nRows = size(XV, 1);
% for j = 1:nRows
%    XV(j, 1:3) = distanceCoefficient*XV(j, 1:3);
%    XV(j, 4:6) = velocityCoefficient*XV(j, 4:6);
% end

% nRows = size(XV, 1);
% for j = 1:nRows
%    XV(j, 1:3) = distanceCoefficient*(InitStateV*XV(j, 1:3)')';
%    XV(j, 4:6) = velocityCoefficient*(InitStateV*XV(j, 4:6)')';
% end
nRows = size(XV, 1);
for j = 1:nRows
   XV(j, 1:3) = distanceCoefficient*XV(j, 1:3);
   XV(j, 4:6) = velocityCoefficient*XV(j, 4:6);
end


p3 = plot(XV(:, 1), XV(:, 2), '--', 'LineWidth', 2, 'color', 'b');
axis equal;
grid on;
hold on;

% неустойчивое многообразие
um = trajectoryUnstableManifold(orbitSE.mu, orbitSE.xL, xE, tE, t01, time1);
UM = rotating2inertial(orbitSE.mu, um, t01, time1, InitStateE);

% nRows = size(UM, 1);
% for j = 1:nRows
%     UM(j, 1:3) = (InitStateE*UM(j, 1:3)')';
%     UM(j, 4:6) = (InitStateE*UM(j, 4:6)')';
% end

p4 = plot(UM(:, 1), UM(:, 2), 'LineWidth', 2, 'color', 'g');
axis equal;
grid on;
hold on;

% устойчивое многообразие
sm = trajectoryStableManifold(orbitSV.mu, orbitSV.xL, xV, tV, t02, time2);
sm = sm(end:-1:1, :);
SM = rotating2inertial(orbitSV.mu, sm, (t02 - time2)*timeCoefficient,...
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

p5 = plot(SM(:, 1), SM(:, 2), 'LineWidth', 2, 'color', 'c');
axis equal;
grid on;
hold on;

% перелёт между двумя многообразиями
muS = 1 - orbitSE.mu;
r1 = UM(end, 1:3)';
v1 = UM(end, 4:6)';
r2 = SM(1, 1:3)';
v2 = SM(1, 4:6)';
[v1corr, v2corr, ~, ~] = lambert(r1, r2, Tbetween, 0, muS);

opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[~, X] = ode113(@(t, x)func2bp(t, x, muS), [0, Tbetween], [r1; v1corr], opts);
p6 = plot(X(:, 1), X(:, 2), 'LineWidth', 2, 'color', 'm');
axis equal;
grid on;
hold on;

% Начальные положения
text(XE(1, 1)-0.1, XE(1, 2) + 0.15, '02.07.2044');
quiver(XE(1, 1), XE(1, 2) + 0.12, 0, -0.12, 'LineWidth', 2.0, 'MaxHeadSize', 1.0);
hold on;

text(XV(1, 1)-0.1, XV(1, 2) + 0.15, '02.07.2044');
quiver(XV(1, 1), XV(1, 2) + 0.12, 0, -0.12, 'LineWidth', 2.0, 'MaxHeadSize', 1.0);
hold on;

text(SM(end, 1)-0.1, SM(end, 2) + 0.15, '15.12.2044');
quiver(SM(end, 1), SM(end, 2) + 0.12, 0, -0.12, 'LineWidth', 2.0, 'MaxHeadSize', 1.0);
hold on;

xlabel('x, безразмер.');
ylabel('y, безразмер.');
legend([p2 p3 p4 p5 p6], 'плоская орбита в системе Солнце-Земля',...
       'плоская орбита в системе Солнце-Венера',...
       'траектория неустойчивого многообразия',...
       'траектория устойчивого многообразия', 'двухимпульсный перелёт',...
       'Location', 'best');

%% начальное приближение для метода параллельной пристрелки

tspanE = linspace(-3*tE, 0, 31);
tspanE = tspanE(1:length(tspanE)-1);
tspanE = [tspanE, linspace(0, t01, 6)];
tspanE = tspanE(1:length(tspanE)-1);

tspanUME = linspace(t01, t01 + time1, 6);
tspanUME = tspanUME(1:length(tspanUME)-1);

tspanT = linspace(t01 + time1, t01 + time1 + Tbetween, 16);
tspanT = tspanT(1:length(tspanT)-1);

tspanSMV = linspace((t02 - time2), t02, 6);
tspanSMV = tspanSMV(1:length(tspanSMV)-1);

tspanV = linspace(t02, tV, 20);

S0 = zeros(length(tspanE) + length(tspanUME) + length(tspanT) + length(tspanSMV) + length(tspanV), 6);

n = size(xE, 1);
for i = 1:length(tspanE)
    t = rem(tspanE(i)+3*tE, tE);
    j = fix(t*(n - 1)/tE + 1);
    alpha = 1 - (t*(n - 1)/tE + 1 - j);
    s0 = alpha*xE(rem(j - 1, n) + 1, 1:6) + (1 - alpha)*xE(rem(j, n) + 1, 1:6);
    S0(i, :) = s0 + [orbitSE.mu, 0, 0, 0, 0, 0];
end

n = size(um, 1);
for i = 1:length(tspanUME)
    t = tspanUME(i) - t01;
    j = fix(t*(n - 1)/time1 + 1);
    alpha = 1 - (t*(n - 1)/time1 + 1 - j);
    s0 = alpha*um(rem(j - 1, n) + 1, 1:6) + (1 - alpha)*um(rem(j, n) + 1, 1:6);
    S0(i+length(tspanE), :) = s0 + [orbitSE.mu, 0, 0, 0, 0, 0];
end

tspanDays = [tspanE, tspanUME]*TimeUnit;
S10 = rotating2ICCS(S0(1:length(tspanE) + length(tspanUME), 1:6)',...
    JD0, tspanDays, 'Earth', 'Sun', DistUnit, VelUnit)';

n = size(X, 1);
for i = 1:length(tspanT)
    t = tspanT(i) - (t01 + time1);
    j = fix(t*(n - 1)/Tbetween + 1);
    alpha = 1 - (t*(n - 1)/Tbetween + 1 - j);
    s0 = alpha*X(rem(j - 1, n) + 1, 1:6) + (1 - alpha)*X(rem(j, n) + 1, 1:6);
    S0(i+length(tspanE)+length(tspanUME), :) = s0;
end

%tspanDays = [tspanT]*TimeUnit;
S20 = S0(length(tspanE) + length(tspanUME) + 1:length(tspanE) + length(tspanUME) + length(tspanT), :);

n = size(sm, 1);
for i = 1:length(tspanSMV)
    t = tspanSMV(i) - (t02 - time2);
    j = fix(t*(n - 1)/(time2) + 1);
    alpha = 1 - (t*(n - 1)/(time2) + 1 - j);
    s0 = alpha*sm(rem(j - 1, n) + 1, 1:6) + (1 - alpha)*sm(rem(j, n) + 1, 1:6);
    S0(i+length(tspanE)+length(tspanUME)+length(tspanT), :) = s0 + [orbitSV.mu, 0, 0, 0, 0, 0];
end

n = size(xV, 1);
for i = 1:length(tspanV)
    t = rem(tspanV(i), tV);
    j = fix(t*(n - 1)/(tV) + 1);
    alpha = 1 - (t*(n - 1)/(tV) + 1 - j);
    s0 = alpha*xV(rem(j - 1, n) + 1, 1:6) + (1 - alpha)*xV(rem(j, n) + 1, 1:6);
    S0(i+length(tspanE)+length(tspanUME)+length(tspanT)+length(tspanSMV), :) = s0 + [orbitSV.mu, 0, 0, 0, 0, 0];
end

tspanDays = [tspanSMV, tspanV]*TimeUnitV;
S30 = rotating2ICCS(S0(length(tspanE) + length(tspanUME) + length(tspanT) + 1:length(tspanE) + length(tspanUME) + length(tspanT)+ length(tspanSMV) + length(tspanV), 1:6)',...
    JD0, tspanDays, 'Venus', 'Sun', DistUnitV, VelUnitV)';
nRows = size(S30, 1);
for j = 1:nRows
   S30(j, 1:3) = distanceCoefficient*S30(j, 1:3);
   S30(j, 4:6) = velocityCoefficient*S30(j, 4:6);
end

tspan = [tspanE, tspanUME, tspanT, tspanSMV, tspanV];
S0 = [S10; S20; S30];

figure
plot(S0(1:35,1), S0(1:35,2), 'color', 'r');
axis equal;
grid on;
hold on;

plot(S0(36:40,1), S0(36:40,2), 'color', 'g');
axis equal;
grid on;
hold on;

plot(S0(41:55,1), S0(41:55,2), 'color', 'm');
axis equal;
grid on;
hold on;

plot(S0(56:60,1), S0(56:60,2), 'color', 'c');
axis equal;
grid on;
hold on;

plot(S0(61:80,1), S0(61:80,2), 'color', 'b');
axis equal;
grid on;

%% построение траектории в эфемеридной модели

splines = PrepareSplines(tspan, JD0, (tE), DistUnit, VelUnit, TimeUnit, DistUnitV, VelUnitV, TimeUnitV);

s0 = [reshape(S0', [], 1); v1 - v1corr; v2corr - v2];
last = S0(end,:)';
% xxx = ine2rot(s0(1:6),JD0,tspan(1)*TimeUnit,'Earth','Sun',DistUnit,VelUnit) - [orbitSE.mu;0;0;0;0;0];
% xxxx = func3bp(tspan(1), xxx, orbitSE.mu);
% yyy = ine2rot(s0(1:6),JD0,tspan(1)*TimeUnit,'Venus','Sun',DistUnit,VelUnit) - [orbitSV.mu;0;0;0;0;0];
% yyyy = func3bp(tspan(1), yyy, orbitSV.mu);
F0 = eye(6);
figure
for i = 1:(length(tspan) - 1)
    x_prev = s0(1 + 6*(i-1):6 + 6*(i-1));
    x0 = [x_prev; reshape(F0, [], 1)];
    tsp = [tspan(i), tspan(i + 1)];
    
    [~, XX] = ode113(@(t, x)funcVariations(t, x, splines), tsp, x0, opts);
    if (i <= 35)
        plot(XX(:,1), XX(:,2), 'color', 'r');
    elseif (i <= 40)
        plot(XX(:,1), XX(:,2), 'color', 'g');
    elseif (i <= 55)
        plot(XX(:,1), XX(:,2), 'color', 'm');
    elseif (i <= 60)
        plot(XX(:,1), XX(:,2), 'color', 'c');
    else
        plot(XX(:,1), XX(:,2), 'color', 'b');
    end
    axis equal;
    grid on;
    hold on;
end

options = optimoptions('fmincon','Display','iter','Algorithm', 'sqp',...
    'CheckGradients', false,...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
    'OptimalityTolerance', 1e-3,'StepTolerance',1e-3,...
    'MaxIterations', 20);

[s, J, ~, output] = fmincon(@(s)optimFunctionShooting(s),s0,...
    [], [], [], [], [], [], [@(s)continuityR(s, last, tspan, splines)], options);

fprintf('характеристическая скорость %f км/с\n',...
    J*VelUnit);

%% построение графика адаптированной траектории

s = s(1:end-6);
F0 = eye(6);

figure
plot(0, 0, 'ko', 'LineWidth', 10, 'color', 'y');
axis equal;
grid on;
hold on;
for i = 1:(length(tspan)-1)
    x_prev = s(1 + 6*(i-1):6 + 6*(i-1));
    
    x0 = [x_prev; reshape(F0, [], 1)]; 
    tsp = [tspan(i), tspan(i + 1)];
    [~, X] = ode113(@(t, x)funcVariations(t, x, splines), tsp, x0, opts);
    if (i <= 35)
        p2=plot(X(:,1), X(:,2), 'color', 'r', 'LineWidth', 2);
    elseif (i <= 40)
        p3=plot(X(:,1), X(:,2), 'color', 'g', 'LineWidth', 2);
    elseif (i <= 55)
        p4=plot(X(:,1), X(:,2), 'color', 'm', 'LineWidth', 2);
    elseif (i <= 60)
        p5=plot(X(:,1), X(:,2), 'color', 'c', 'LineWidth', 2);
    else
        p6=plot(X(:,1), X(:,2), 'color', 'b', 'LineWidth', 2);
    end
    axis equal;
    grid on;
    hold on;
end
xlabel('x, безразмер.');
ylabel('y, безразмер.');
legend([p2 p3 p4 p5 p6], 'плоская орбита в системе Солнце-Земля',...
       'плоская орбита в системе Солнце-Венера',...
       'траектория неустойчивого многообразия',...
       'траектория устойчивого многообразия', 'двухимпульсный перелёт',...
       'Location', 'best');

disp(5);

end
