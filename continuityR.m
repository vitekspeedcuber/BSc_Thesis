function [c, ceq, gc, gceq] = continuityR(s, last, tspan, splines)

c = [];
gc = [];

N = size(s, 1);

ceq = zeros(N-12,1);
gceq = zeros(N-12, N);

dv1 = s(end - 5:end - 3);
dv2 = s(end - 2:end);

F0 = eye(6);
for i = 1:(length(tspan) - 1)
    x_prev = s(1 + 6*(i-1):6 + 6*(i-1));
    x_next = s(1 + 6*i:6 + 6*i);
    
    x0 = [x_prev; reshape(F0, [], 1)];
    tsp = [tspan(i), tspan(i + 1)];
    opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
    [~, X] = ode113(@(t, x)funcVariations(t, x, splines), tsp, x0,opts);
    F = reshape(X(end,7:42), 6, 6);
    
    if (i == 40)
        ceq(1 + 6*(i-1):6 + 6*(i-1)) = X(end,1:6)' - x_next - [0; 0; 0; dv1];
        gceq(1 + 6*(i-1):6 + 6*(i-1), 1 + 6*(i-1):6 + 6*(i-1)) = F;
        gceq(1 + 6*(i-1):6 + 6*(i-1), 1 + 6*i:6 + 6*i) = -eye(6);
        gceq(4 + 6*(i-1):6 + 6*(i-1),N-5:N-3) = -eye(3);
    elseif (i == 55)
        ceq(1 + 6*(i-1):6 + 6*(i-1)) = X(end,1:6)' - x_next - [0; 0; 0; dv2];
        gceq(1 + 6*(i-1):6 + 6*(i-1), 1 + 6*(i-1):6 + 6*(i-1)) = F;
        gceq(1 + 6*(i-1):6 + 6*(i-1), 1 + 6*i:6 + 6*i) = -eye(6);
        gceq(4 + 6*(i-1):6 + 6*(i-1),N-2:N) = -eye(3);
    else
%         if (i == (length(tspan) - 1))
%             ceq(1 + 6*(i-1):6 + 6*(i-1)) = X(end,1:6)' - last;
%         else
%             ceq(1 + 6*(i-1):6 + 6*(i-1)) = X(end,1:6)' - x_next;
%             gceq(1 + 6*(i-1):6 + 6*(i-1), 1 + 6*i:6 + 6*i) = -eye(6);
%         end
        ceq(1 + 6*(i-1):6 + 6*(i-1)) = X(end,1:6)' - x_next;
        gceq(1 + 6*(i-1):6 + 6*(i-1), 1 + 6*i:6 + 6*i) = -eye(6);
        gceq(1 + 6*(i-1):6 + 6*(i-1), 1 + 6*(i-1):6 + 6*(i-1)) = F;
    end
    X = [];
end

gceq = gceq';

end
