% по параметрам орбиты вблизи коллинеарной точки либрации и
% оценке половины периода орбиты
% возвращает линейное приближение фазового вектора на половине периода

function X = linearApproximationVector(orbit, T)

xL = orbit.xL;
alpha = orbit.alpha;
beta = orbit.beta;
gamma = orbit.gamma;
gamma_dash = orbit.gamma_dash;
phi1 = orbit.phi1;
phi2 = orbit.phi2;
lambda = orbit.lambda;
wp = orbit.wp;
wv = orbit.wv;
k1 = orbit.k1;
k2 = orbit.k2;

x = xL + alpha*cos(wp*T + phi1) + gamma*exp(lambda*T)...
    + gamma_dash*exp(-lambda*T);
y = -k2*alpha*sin(wp*T + phi1) + k1*(gamma*exp(lambda*T)...
    - gamma_dash*exp(-lambda*T));
z = beta*cos(wv*T + phi2);

Vx = -alpha*wp*sin(wp*T + phi1) + gamma*lambda*exp(lambda*T)...
     - gamma_dash*lambda*exp(-lambda*T);
Vy = -k2*alpha*wp*cos(wp*T + phi1) + k1*(gamma*lambda*exp(lambda*T)...
     + gamma_dash*lambda*exp(-lambda*T));
Vz = -beta*wv*sin(wv*T + phi2);

X = [x; y; z; Vx; Vy; Vz];

end
