% по значениям фазового вектора (траектории) в ВСК
% возвращает значения фазового вектора (траекторию) в ИСК 
% (с учётом эфемеридного движения)

function Xine = rotating2ICCS(Xrot, JD0, tspanDays, FirstBody, SecondBody, DistUnit, VelUnit)

% RCS basis vectors
JD = JD0 + tspanDays(:);
ncolsTspan = length(tspanDays);
E1 = zeros(3,ncolsTspan);
E2 = zeros(3,ncolsTspan);
E3 = zeros(3,ncolsTspan);
dtheta = zeros(1,ncolsTspan);
[R1,V1] = planetEphemeris(JD,SecondBody,FirstBody,'430');
R1 = R1/DistUnit;
V1 = V1/VelUnit;
for nocolTspan = 1:ncolsTspan
    r = R1(nocolTspan,:).';
    v = V1(nocolTspan,:).';
    E1(:,nocolTspan)   = r/norm(r);
    E3(:,nocolTspan)   = cross(r,v)/norm(cross(r,v));
    E2(:,nocolTspan)   = cross(E3(:,nocolTspan),E1(:,nocolTspan));
    dtheta(nocolTspan) = norm(cross(r,v))/norm(r)^2;
end

% Translate from RCS to ICS
Xine = zeros(6,ncolsTspan);
for nocolXrot = 1:ncolsTspan
    C = [E1(:,nocolXrot),E2(:,nocolXrot),E3(:,nocolXrot)];
    Omega = [dtheta(nocolXrot)*C(1,2),-dtheta(nocolXrot)*C(1,1),0;
             dtheta(nocolXrot)*C(2,2),-dtheta(nocolXrot)*C(2,1),0;
             dtheta(nocolXrot)*C(3,2),-dtheta(nocolXrot)*C(3,1),0];
    Xine(:,nocolXrot) = [C,zeros(3,3);Omega,C]*Xrot(:,nocolXrot);
end

end
