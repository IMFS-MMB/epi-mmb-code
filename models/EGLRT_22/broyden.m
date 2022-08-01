%% Broyden's Rootfinding Algorithm
function [xstar, fval, iter,distF] = broyden(f,x0,critF,critX,maxiter)
% The function uses the "good" broyden algorithm to solve for the root of a function F.
% X0 is the starting guess. CRITX the precision in variable X, CRITF the precision
% to which the root needs to be found. MAXITER is the maximum number of iterations.
distF= NaN(1,maxiter); distF(1) = 9999; distX = 9999;
iter = 0;
xnow = x0(:); % x needs to be a column vector
Fnow = f(xnow);
Fnow = Fnow(:);  % F needs to be a column vector
Bnow = eye(length(xnow));

while distF(iter+1)>critF && distX>critX && iter<maxiter
    iter = iter+1; % count number of iterations
    Fold = Fnow; % Store last function values
    xold = xnow; % Store last x guess
    xnow = xnow - Bnow*Fnow; % Update x guess
    Fnow = f(xnow); % Evaluate the function
    Fnow = Fnow(:); % If Function is matrix valued then vectorize
    Dx = xnow - xold; % Change in x
    Dy = Fnow - Fold; % Change in f(x)
    % update inverse Jacobian
    Bnow = Bnow + (Dx - Bnow*Dy)*(Dx'*Bnow)/(Dx'*Bnow*Dy);
    distF(iter+1) = max(abs(Fnow)); % Calculate the distance from zero
    distX = max(abs(Dx)); % Calculate the change in variable
end
fval=Fnow; xstar=xnow;
if distF(iter+1)<critF
    %disp('Broyden has found a root of function f')
else
    error('Failed to find a root of functionf')
end