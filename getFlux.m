function V = getFlux(K,P,h)
% Compute transmissibilities by harmonic averaging.
[Nx,Ny] = size(K(:,:,1)); 
hx = h(1); 
hy = h(2); 
hz = h(3);
L = K.^(-1);
tx = 2*hy*hz/hx; 
TX = zeros(Nx+1,Ny);
ty = 2*hx*hz/hy; 
TY = zeros(Nx,Ny+1);
TX(2:Nx,:) = tx./(L(1:Nx-1,:)+L(2:Nx,:));
TY(:,2:Ny) = ty./(L(:,1:Ny-1)+L(:,2:Ny));

P = reshape(P,Nx,Ny);
V.x = zeros(Nx+1,Ny);
V.y = zeros(Nx,Ny+1);
V.x(2:Nx,:) = (P(1:Nx-1,:)-P(2:Nx,:)).*TX(2:Nx,:);
V.y(:,2:Ny) = (P(:,1:Ny-1)-P(:,2:Ny)).*TY(:,2:Ny);
end