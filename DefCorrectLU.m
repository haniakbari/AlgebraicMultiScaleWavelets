function [x,erit] = DefCorrectLU(itol,iter,x)
global A L U b;
conv = 0;
r = b - A*x;
nritol = norm(r)*itol;
erit = zeros(iter,1);
for i=1:iter
x = x + U\(L\r);
r = b - A*x;
normr = norm(r);
erit(i) = normr;
if (normr < nritol), conv=1; break; end
end
Result = sprintf('converged PRM(ILU):\n \tIters(%d)\trelres(%g)\n',...
    i,itol*normr/nritol);
if ~conv, Result = ['NOT ' Result]; end
display(Result);