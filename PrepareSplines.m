function splines = PrepareSplines(tspan, JD0, orbPeriod, DistUnit, VelUnit, TimeUnit, DistUnitV, VelUnitV, TimeUnitV)

timeCoefficient = 365.25/224.7; % T_Earth / T_Ven
distanceCoefficient = 108200000/149600000; % R_Ven / R_Earth
velocityCoefficient = distanceCoefficient*timeCoefficient;

dt = 1.0e-01;

t0 = tspan(1) - dt;
tf = tspan(end) + dt;

% t0Days = t0*TimeUnit;
% tfDays = tf*TimeUnit;

Planets = {'Sun';'Venus';'Earth'};%'Mars';'Jupiter';'Saturn';'Uranus'};

splines = struct;

for k = 1:size(Planets,1)
    
    planet = Planets{k};
    
    switch planet
        case 'Sun'
            t0Days = t0*TimeUnit;
            tfDays = tf*TimeUnit;

            N = fix((tf-t0)/orbPeriod*75) + 1;
            JD = JD0 + linspace(t0Days,tfDays,N);
            T  = linspace(t0,tf,N);
    
            % ћассив фазовых состо€ний планеты
            [R,V] = planetEphemeris(JD(:),'Sun',planet,'430');
            R = R/DistUnit;
            V = V/VelUnit;
            X = [R.'; V.'];
%         case 'Mercury'
%             N = fix((tf-t0)/orbPeriod*7) + 2;
        case 'Venus'
            t0Days = t0*TimeUnitV;
            tfDays = tf*TimeUnitV;
            
            N = fix((tf-t0)/orbPeriod*60) + 1;
            JD = JD0 + linspace(t0Days,tfDays,N);
            T  = linspace(t0,tf,N);
    
            % ћассив фазовых состо€ний планеты
            [R,V] = planetEphemeris(JD(:),'Sun',planet,'430');
            R = (R/DistUnitV)*distanceCoefficient;
            V = (V/VelUnitV)*velocityCoefficient;
            X = [R.'; V.'];
        case 'Earth'
            t0Days = t0*TimeUnit;
            tfDays = tf*TimeUnit;
            
            N = fix((tf-t0)/orbPeriod*60) + 1;
            JD = JD0 + linspace(t0Days,tfDays,N);
            T  = linspace(t0,tf,N);
    
            % ћассив фазовых состо€ний планеты
            [R,V] = planetEphemeris(JD(:),'Sun',planet,'430');
            R = R/DistUnit;
            V = V/VelUnit;
            X = [R.'; V.'];
%         case 'Mars'
%             N = fix((tf-t0)/orbPeriod*15) + 1;
%         case 'Jupiter'
%             N = fix((tf-t0)/orbPeriod*15) + 1;
%         case 'Saturn'
%             N = fix((tf-t0)/orbPeriod*20) + 1; % даже 7 можно
%         case 'Uranus'
%             N = fix((tf-t0)/orbPeriod*15) + 1; % даже 7 можно
%         case 'Neptune'
%             N = fix((tf-t0)/orbPeriod*15) + 1; % даже 7 можно
%         case 'Pluto'
%             N = fix((tf-t0)/orbPeriod*15) + 1;
        otherwise
            error('Unknown source of perturbation');
    end
    
%     JD = JD0 + linspace(t0Days,tfDays,N);
%     T  = linspace(t0,tf,N);
%     
%     % ћассив фазовых состо€ний планеты
%     [R,V] = planetEphemeris(JD(:),'Sun',planet,'430');
%     R = R/DistUnit;
%     V = V/VelUnit;
%     X = [R.'; V.'];
    
    Spline = spline(T,X); %#ok<NASGU>
    
    eval(['splines.Spline',planet,'Ephem = Spline;']);
    
end

end
