% функция правых частей уравнений в вариациях

function dxdt = func(t, x, mu)

r = x(1:3);
v = x(4:6);
r1 = sqrt((r(1) + mu)^2 + (r(2))^2 + (r(3))^2);
r2 = sqrt((r(1) - 1 + mu)^2 + (r(2))^2 + (r(3))^2);

Ux = r(1) - (1 - mu)*(r(1) + mu)/r1^3 - mu*(r(1) - 1 + mu)/r2^3;
Uy = r(2) - (1 - mu)*r(2)/r1^3 - mu*r(2)/r2^3;
Uz = -(1 - mu)*r(3)/r1^3 - mu*r(3)/r2^3;

F = reshape(x(7:42), 6, 6);
Omega = [0, -1, 0;
         1, 0, 0;
         0, 0, 0];

a1 = (1 - mu)/r1^3;
a2 = mu/r2^3;
b1 = 3*(1 - mu)/r1^5;
b2 = 3*mu/r2^5;

Uxx = 1 - a1 + b1*(r(1) + mu)^2 - a2 + b2*(r(1) - 1 + mu)^2;
Uyy = 1 - a1 + b1*r(2)^2 - a2 + b2*r(2)^2;
Uzz = -a1 + b1*r(3)^2 - a2 + b2*r(3)^2;
Uxy = b1*(r(1) + mu)*r(2) + b2*(r(1) - 1 + mu)*r(2);
Uyx = Uxy;
Uxz = b1*(r(1) + mu)*r(3) + b2*(r(1) - 1 + mu)*r(3);
Uzx = Uxz;
Uyz = b1*r(2)*r(3) + b2*r(2)*r(3);
Uzy = Uyz;

Urr = [Uxx, Uxy, Uxz;
       Uyx, Uyy, Uyz;
       Uzx, Uzy, Uzz];

% mu_dash = (1 - mu)/(abs(xL + mu))^3 + mu/(abs(xL - 1 + mu))^3;
% 
% Urr = [1 + 2*mu_dash, 0, 0;
%        0, 1 - mu_dash, 0;
%        0, 0, - mu_dash];

A = [zeros(3), eye(3);
     Urr, -2*Omega];

dxdt = [v; (2*v(2) + Ux); (-2*v(1) + Uy); Uz; reshape(A*F, [], 1)];

end
