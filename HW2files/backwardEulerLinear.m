function uT = backwardEulerLinear(u0,kappa,A,dt,T)
    N = length(u0);
    K = ceil(T/dt);

    % implement backward Euler here and output just the solution at the
    % final time.
    uT = u0;
    for i=1: K
        u_i = (eye(N) - kappa * dt * A) \ uT;
        uT = u_i;
    end
end

