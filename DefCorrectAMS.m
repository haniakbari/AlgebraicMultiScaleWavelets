function [x,erit] = DefCorrectAMS(itol,iter,x,I,R,Ac)
global A L U b;
conv = 0;
r = b - A*x;
nritol = norm(r)*itol;
if (norm(r) < nritol), 
    fprintf('Good Initial Sol, No iteration \n'); erit = 1; return; 
end
erit = zeros(iter,1);
for i=1:iter
p = I * (Ac\(R*r));
x = x + p + U\(L\(r - A*p));
r = b - A*x;
normr = norm(r);
erit(i) = normr;
if (normr < nritol), conv=1; break; end
end
erit = erit(1:i);
Result = sprintf('converged PrecRich:\n \tIters(%d)\trelres(%g)\n',...
    i,itol*normr/nritol);
if ~conv, Result = ['NOT ' Result]; erit(end) = -1; end
display(Result);