function r = Mgl(r,Ac,R,I)
global A L U;
p = I*(Ac\(R*r));
r = p + U\(L\(r - A*p));
end
