% уравнения движения в вариациях относительно Солнца
% в безразмерных единицах

function dydt = funcVariations(t,y,splines)

% Extract variables
rvect = y(1:3);
vvect = y(4:6);
Phi   = reshape(y(7:42),6,6);
% dxdt0 = y(43:48);

% dydt = zeros(6+36+6,1);
dydt = zeros(6+36,1);

dydt(1:3) = vvect;

dfdr = 0;

% Earth gravitational field
muEarth = 5.972/(5.972 + 1.989e6);
[FE,dFEdr] = GravForce(muEarth,splines.SplineEarthEphem,t,rvect);
dydt(4:6) = dydt(4:6) + FE;
dfdr = dfdr + dFEdr;

% Sun gravitational field
muSun = 1 - muEarth;
[FS,dFSdr] = GravForce(muSun,splines.SplineSunEphem,t,rvect);
dydt(4:6) = dydt(4:6) + FS;
dfdr = dfdr + dFSdr;

% Mercury gravitational field
% muMercury = 3.285/(59.72 + 1.989e7);
% [Fm,dFmdr] = GravForce(muMercury,splines.SplineMercuryEphem,t,rvect);
% dydt(4:6) = dydt(4:6) + Fm;
% dfdr = dfdr + dFmdr;

% Venus gravitational field
muVenus = 4.867/(5.972 + 1.989e6);
[FV,dFVdr] = GravForce(muVenus,splines.SplineVenusEphem,t,rvect);
dydt(4:6) = dydt(4:6) + FV;
dfdr = dfdr + dFVdr;

% Mars gravitational field
% muMars = 6.39/(59.72 + 1.989e7);
% [FM,dFMdr] = GravForce(muMars,splines.SplineMarsEphem,t,rvect);
% dydt(4:6) = dydt(4:6) + FM;
% dfdr = dfdr + dFMdr;

% Jupiter gravitational field
% muJupiter = 1.8987/(1.8987 + 1.989e3);
% [FJ,dFJdr] = GravForce(muJupiter,splines.SplineJupiterEphem,t,rvect);
% dydt(4:6) = dydt(4:6) + FJ;
% dfdr = dfdr + dFJdr;

% Saturn gravitational field
% muSaturn = 5.683/(5.683 + 1.989e4);
% [Fs,dFsdr] = GravForce(muSaturn,splines.SplineSaturnEphem,t,rvect);
% dydt(4:6) = dydt(4:6) + Fs;
% dfdr = dfdr + dFsdr;

% Uranus gravitational field
% muUranus = 8.6849/(8.6849 + 1.989e5);
% [FU,dFUdr] = GravForce(muUranus,splines.SplineUranusEphem,t,rvect);
% dydt(4:6) = dydt(4:6) + FU;
% dfdr = dfdr + dFUdr;

% Neptune gravitational field
% muNeptune = 1.024/(1.024 + 1.989e4);
% [FN,dFNdr] = GravForce(muNeptune,splines.SplineNeptuneEphem,t,rvect);
% dydt(4:6) = dydt(4:6) + FN;
% dfdr = dfdr + dFNdr;

% Pluto gravitational field
% muPluto = 1.3/(1.3 + 1.989e8);
% [FP,dFPdr] = GravForce(muPluto,splines.SplinePlutoEphem,t,rvect);
% dydt(4:6) = dydt(4:6) + FP;
% dfdr = dfdr + dFPdr;

% Solar radiation pressure
% [FSRP,dFSRPdr] = SRPForce(splines.SplineSunEphem,splines.SplineEarthEphem,t,rvect);
% dydt(4:6) = dydt(4:6) + FSRP;
% dfdr = dfdr + dFSRPdr;

% Equations for state transition matrix
dfdx = [zeros(3,3),eye(3);dfdr,zeros(3,3)];
dPhidt = dfdx * Phi;
dPhidt = reshape(dPhidt,36,1);
dydt(7:42) = dPhidt;

% Equations for sensitivity vector dxdt0
% dydt(43:48) = dfdx * dxdt0;

end

function [F,dFdrpl] = GravForce(muP,spline,t,rsc)

% Extract phase state of the planet from spline
[breaks,coefs,~,~,d] = unmkpp(spline);
pp = mkpp(breaks,coefs,d);
x = ppval(pp,t);
rpl = x(1:3); % Center --> Planet 
%vpl = x(4:6);

dr      = rsc - rpl; % Planet --> SC
drnorm  = norm(dr);

drc     =     - rpl; %Planet --> Center
drcnorm = norm(drc);

if (drnorm ~= 0)
    sc2pl = -muP*dr/drnorm^3;
else
    sc2pl = [0; 0; 0];
end
if (drcnorm ~= 0)
    c2pl  = -muP*drc/drcnorm^3;
    dFdrpl =  muP*eye(3)/drnorm^3 - 3*muP*(dr*dr.')/drnorm^5 - muP*eye(3)/drcnorm^3 + 3*muP*(drc*drc.')/drcnorm^5;
else
    c2pl = [0; 0; 0];
    dFdrpl =  muP*eye(3)/drnorm^3 - 3*muP*(dr*dr.')/drnorm^5;
end

F = sc2pl - c2pl;
dFdrpl = -dFdrpl;
dFdr   = -muP*eye(3)/drnorm^3 + 3*muP*(dr*dr.')/drnorm^5;

end

% function [F,dFdr] = SRPForce(splineSun,splineEarth,t,rsc)
% 
% % Extract phase state of the Sun from spline
% [breaks,coefs,~,~,d] = unmkpp(splineSun);
% pp = mkpp(breaks,coefs,d);
% x = ppval(pp,t);
% rSun = x(1:3); % Moon --> Sun
% 
% % Extract velocity of the Earth from spline
% [breaks,coefs,~,~,d] = unmkpp(splineEarth);
% pp = mkpp(breaks,coefs,d);
% x = ppval(pp,t);
% rEarth = x(1:3);
% 
% drSun     = rsc - rSun; % Sun --> spacecraft
% drnormSun = norm(drSun);
% 
% % Area to mass ratio (0.005--0.02 for standard s/c without solar sail or big panels)
% PAdm = 0.006;
% 
% % 4.56e-06 N/m^2 in the vicinity of the Earth (i.e. at 1 AU)
% % 2.697489232675493e-03 is acceleration unit (m/s^2, the Earth-Moon system)
% % 3.891673383540797e+02 equals 1 AU
% Coeff = (4.56e-06*PAdm)/2.697489232675493e-03*3.891673383540797e+02^2;
% 
% %[inten,didr] = intensity(rsc,rSun,rEarth);
% 
% F      =  Coeff*drSun/drnormSun^3;
% dFdr   =  (Coeff*eye(3)/drnormSun^3 - 3*Coeff*(drSun*drSun.')/drnormSun^5) + Coeff*(drSun) /drnormSun^3;
% 
% end
