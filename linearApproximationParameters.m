% возвращает постоянные интегрирования орбиты в линейном приближении
% вблизи коллинеарной точки либрации и 
% оценку половины периода орбиты

function parameters = linearApproximationParameters(mu, xL)

mu_dash = mu/(abs(xL - 1 + mu))^3 + (1 - mu)/(abs(xL + mu))^3;

lambda = ((mu_dash - 2 + (9*mu_dash^2 - 8*mu_dash)^(1/2))/2)^(1/2);
wp = ((2 - mu_dash + (9*mu_dash^2 - 8*mu_dash)^(1/2))/2)^(1/2);
wv = (mu_dash)^(1/2);

k1 = (lambda^2 - 2*mu_dash - 1)/(2*lambda);
k2 = (wp^2 + 2*mu_dash + 1)/(2*wp);

T = pi/wp; 

parameters = [lambda; wp; wv; k1; k2; T];

end
