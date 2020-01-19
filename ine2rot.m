function Xrot = ine2rot(Xine,JD0,tspanDays,FirstBody,SecondBody,DistUnit,VelUnit)
%ine2rot   Translate phase states from inertial CS to the rotating CS (w.r.t. the 2nd body).
%
%  Xine       -- (6xN) phase states in the inertial CS with center at the 2nd body
%  JD0        --       reference julian date
%  tspanDays  -- (1xN) time moments relative to JD0 at which states are identified (in days)
%  FirstBody  --       name of the first body in three-body system
%  SecondBody --       name of the second body in the three-body system
%  DistUnit   --       unit of distance ( in km )
%  VelUnit    --       unit of velocity ( in km/s )
%
%  Xrot       -- (6xN) phase states in rotating frame with center at the 2nd body
%
% Author: Maksim Shirobokov
% Date: 30.01.2018

if size(Xine,1) ~= 6
    error('ine2rot:   size(Xine,1) ~= 6');
elseif size(Xine,2) ~= length(tspanDays)
    error('ine2rot:   size(Xine,2) ~= length(tspanDays)');
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

% Translate from ICS to RCS
Xrot = zeros(6,ncolsTspan);
for nocolXine = 1:size(Xine,2)
    C = [E1(:,nocolXine),E2(:,nocolXine),E3(:,nocolXine)];
    Omega = [dtheta(nocolXine)*C(1,2),-dtheta(nocolXine)*C(1,1),0;
             dtheta(nocolXine)*C(2,2),-dtheta(nocolXine)*C(2,1),0;
             dtheta(nocolXine)*C(3,2),-dtheta(nocolXine)*C(3,1),0];
    Xrot(:,nocolXine) = [C,zeros(3,3);Omega,C]\Xine(:,nocolXine);
end

end