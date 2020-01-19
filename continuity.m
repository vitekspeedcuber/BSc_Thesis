function [c, ceq, gc, gceq] = continuity(s, first, last, tspan, splines)

N = size(s, 1);
c = [];
gc = [];

ceq = [];
gceq = zeros(N-6, N);

dv1 = s(end - 5:end - 3);
dv2 = s(end - 2:end);

F0 = eye(6);
for i = 1:(length(tspan) - 2)
    x_prev = s(1 + 6*(i-1):6 + 6*(i-1));
    x_next = s(1 + 6*i:6 + 6*i);
    
    if (i == 10)
        ceq = [ceq; x_prev - x_next - [0; 0; 0; dv1]];
        gceq(1 + 6*(i-1):6 + 6*(i-1),1 + 6*(i-1):6 + 6*(i-1)) = eye(6);
        gceq(1 + 6*(i-1):6 + 6*(i-1),1 + 6*i:6 + 6*i) = -eye(6);
        %gceq(1 + 6*(i-1):6 + 6*(i-1),N - 5:N - 3) = [zeros(3);-eye(3)];
    elseif (i == 15)
        ceq = [ceq; x_prev - x_next - [0; 0; 0; dv2]];
        gceq(1 + 6*(i-1):6 + 6*(i-1),1 + 6*(i-1):6 + 6*(i-1)) = eye(6);
        gceq(1 + 6*(i-1):6 + 6*(i-1),1 + 6*i:6 + 6*i) = -eye(6);
        %gceq(1 + 6*(i-1):6 + 6*(i-1),N - 2:N) = [zeros(3);-eye(3)];
    elseif (tspan(i) == tspan(i + 1))
        ceq = [ceq; x_prev - x_next];
        gceq(1 + 6*(i-1):6 + 6*(i-1),1 + 6*(i-1):6 + 6*(i-1)) = eye(6);
        gceq(1 + 6*(i-1):6 + 6*(i-1),1 + 6*i:6 + 6*i) = -eye(6);
    else
        x0 = [x_prev; reshape(F0, [], 1)];
        tsp = [tspan(i), tspan(i + 1)];
        [~, X] = ode113(@(t, x)funcVariations(t, x, splines), tsp, x0);
        ceq = [ceq; X(end,1:6)' - x_next];
        F = reshape(X(end,7:42), 6, 6);
        gceq(1 + 6*(i-1):6 + 6*(i-1),1 + 6*(i-1):6 + 6*(i-1)) = F;
        gceq(1 + 6*(i-1):6 + 6*(i-1),1 + 6*i:6 + 6*i) = -eye(6);
        X = [];
    end
end

i = length(tspan) - 1;
x_prev = s(1 + 6*(i-1):6 + 6*(i-1));

x0 = [x_prev; reshape(F0, [], 1)];
tsp = [tspan(i), tspan(i + 1)];
[~, X] = ode113(@(t, x)funcVariations(t, x, splines), tsp, x0);
ceq = [ceq; X(end,1:6)' - last];
F = reshape(X(end,7:42), 6, 6);
gceq(1 + 6*(i-1):6 + 6*(i-1),1 + 6*(i-1):6 + 6*(i-1)) = F;

%ceq = [ceq; s(1:6) - first];
%gceq(N-5:N,1:6) = eye(6);
gceq = gceq';

end
