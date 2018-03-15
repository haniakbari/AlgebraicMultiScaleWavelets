function I = localAMS(G,K,h,v)
Q.ind = [1 numel(K)];
Q.val = [1 0];
[~,~,A,~] = BTPFA(K,Q,h);
net = locNet(G,v);
I = reshape(locAMS(A,net),size(K));
end

function NET = locNet(G,v)
v = find(G == -v);
VerI  = find(G<0);
NET.v = find(VerI == v);
nV    = numel(VerI);

EdgeI = find(G==0);
nE    = numel(EdgeI);

nI = numel(G) - nV - nE;
InnerI = zeros(nI,1);
Gp = unique(G(G > 0));
lk = zeros(length(Gp),1);
k = 0;
for i=1:length(Gp)
    lI = find(G==Gp(i));
    lk(i) = length(lI);
    lk0 = lk(i) + k;
    InnerI(k+1:lk0) = lI;
    k = lk0;
end
% NW is the permutation matrix: Natural ordering to wirebasket ordering
NW = [InnerI; EdgeI; VerI];
% WN is the inverse (transpose) of NW: wirebasket to natural ordering
[~,WN] = sort(NW);
NET.NW = NW;
NET.WN = WN;
NET.lk = lk;
NET.nIEV = [nI nE nV];
end

function I = locAMS(A,net)
pA = A(net.NW,:);
pA = pA(:,net.NW);

nI = net.nIEV(1);
n2 = nI+net.nIEV(2);

AII  = pA(1:nI,1:nI);
AIE  = pA(1:nI, nI+1 : n2);

AEI  = pA(nI+1:n2, 1:nI);
tAEE = pA(nI+1:n2, nI+1:n2) + diag(sum(AEI,2));
AEV  = pA(nI+1:n2, n2+1:end);

T2 = tAEE\AEV(:,net.v);
T1 = AII\(AIE*T2);
T3 = zeros(net.nIEV(3),1); T3(net.v) = 1;
I = [T1 ; -T2 ; T3]; %wirebasket ordered
I = I(net.WN); %to natural ordered
end