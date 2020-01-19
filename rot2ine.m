function Xine = rot2ine(Xrot,JD0,tspanDays,FirstBody,SecondBody,DistUnit,VelUnit)
%rot2ine   Translate phase states from rotating CS to inertial CS (w.r.t. the 2nd body).
%    
%  Xrot       -- (6xN) phase states in rotating frame with center at the 2nd body
%  JD0        --       reference julian date
%  tspanDays  -- (1xN) time moments relative to JD0 at which states are identified (in days)
%  FirstBody  --       name of the first body in three-body system
%  SecondBody --       name of the second body in the three-body system
%  DistUnit   --       unit of distance ( in km )
%  VelUnit    --       unit of velocity ( in km/s )
%
%  Xine       -- (6xN) phase states in the inertial CS with center at the 2nd body
%
% Author: Maksim Shirobokov
% Date: 30.01.2018

if size(Xrot,1) ~= 6
    error('rot2ine:   size(Xrot,1) ~= 1');
elseif size(Xrot,2) ~= length(tspanDays)
    error('rot2ine:   size(Xrot,2) ~= length(tspanDays)');
end

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
