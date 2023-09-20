function uT = forwardEulerLinear(u0,kappa,A,dt,T)
    N = length(u0);
    K = ceil(T/dt);
    
    % implement forward Euler here and return just the solution at the
    % final time.
    uT = u0;
    for i=1: K
        u_i = uT + kappa * dt * A * uT;
        uT = u_i;
    end
  
end

