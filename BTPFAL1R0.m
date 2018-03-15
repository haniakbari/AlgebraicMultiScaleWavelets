function [P,V,A,b]=BTPFAL1R0(K,Q,h)

% Compute transmissibilities by harmonic averaging.
[Nx,Ny] = size(K); 
N  = Nx*Ny;
hx = h(1); 
hy = h(2); 
hz = h(3);
L  = K.^(-1);
tx = 2*hy*hz/hx; 
TX = zeros(Nx+1,Ny);
ty = 2*hx*hz/hy; 
TY = zeros(Nx,Ny+1);
TX(2:Nx,:) = tx./(L(1:Nx-1,:)+L(2:Nx,:));
TY(:,2:Ny) = ty./(L(:,1:Ny-1)+L(:,2:Ny));

% Assemble TPFA discretization matrix.
x1 = reshape(TX(1:Nx,:),N,1); 
x2 = reshape(TX(2:Nx+1,:),N,1);
y1 = reshape(TY(:,1:Ny),N,1); 
y2 = reshape(TY(:,2:Ny+1),N,1);
% 
Pleft = 1; Pright = 0;
dL = zeros(N,1); dL(1:Nx) = tx*K(:,1);      bL = Pleft*dL;
dR = zeros(N,1); dR(N-Nx+1:N) = tx*K(:,Ny); bR = Pright*dR;
%%
DiagVecs = [-y2, -x2, x1+x2+y1+y2+dL+dR, -x1, -y1];
DiagIndx = [-Nx, -1, 0, 1, Nx];
A = spdiags(DiagVecs,DiagIndx,N,N);

q = zeros(N,1);
for i=1:length(Q.ind), q(Q.ind(i)) = Q.val(i); end

b = q + bL + bR;

P = A\b;
P = reshape(P,Nx,Ny);
V.x = zeros(Nx+1,Ny);
V.y = zeros(Nx,Ny+1);
V.x(2:Nx,:) = (P(1:Nx-1,:)-P(2:Nx,:)).*TX(2:Nx,:);
V.y(:,2:Ny) = (P(:,1:Ny-1)-P(:,2:Ny)).*TY(:,2:Ny);
