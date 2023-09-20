function A = getMatrix(N, L)
    dx = L/N;
    % build A in one line here. (doesn't have to be one line. but it can
    % be. and no loops if you want full pts for efficiency)
    % x = ones(1, N - 1);
    A = (diag(-2 * ones(1, N-1), 0) + diag(ones(1, N-2), 1) + diag(ones(1, N-2), -1)) / (dx^2);
    
end