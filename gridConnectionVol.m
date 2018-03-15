function NET = gridConnectionVol(nx,ny,pos)
N = size(pos,1);
NodesD = zeros(3*N,2);
NetD = zeros(4*N,2);
h = 2.^pos(:,3);
xh = min(nx,pos(:,1)+h);
yh = min(ny,pos(:,2)+h);
k = 1;
for i=1:N
    x = pos(i,1);
    y = pos(i,2);
    NodesD(3*i-2,:) = [x y];
    NodesD(3*i-1,:) = [xh(i) y];
    NodesD(3*i,:)   = [x yh(i)];
    NetD(k,:) = [3*i-2 3*i-1];
    NetD(k+1,:) = [3*i-2 3*i];
    k = k+2;
end
NetD = NetD(1:k-1,:);
[G , ~] = connectBN(nx,ny,NodesD,NetD);
G(1,:) = 0;
G(:,1) = 0;
B = false(nx,ny);
B(sub2ind([nx,ny],NodesD(:,1),NodesD(:,2))) = true;
B(1,1) = 0;
B(1,end) = 0;
B(end,1) = 0;
G = G | B;
B(1,:) = 0; B(:,end) = 0; B(:,1) = 0; B(end,:) = 0;
G = sparse(G+B);
NET = primal(G,find(B));
end

function [B,L] = connectBN(nx,ny,Nod,Net)
B  = zeros(nx,ny);
L  = zeros(nx*ny,1);
iL = 0;
tN = 1;
sbPos = zeros(size(Net,1),4);
for i=1:size(Net,1);    
a = Nod(Net(i,1),:);
b = Nod(Net(i,2),:);
tN = tN+2;
if a(1) > b(1), stx = -1; else stx = 1; end 
if a(2) > b(2), sty = -1; else sty = 1; end 
mx = abs(b(1)-a(1)) - 1;
my = abs(b(2)-a(2)) - 1;
xI = (a(1)+stx:stx:b(1)-stx)';
yI = (a(2)+sty:sty:b(2)-sty)';
midx= floor((a(1)+b(1))/2);
midy= floor((a(2)+b(2))/2);
switch min(mx,my)
    case -1, % in a row
        if my == -1, 
            lxI = length(xI);            
            L(iL+1:iL+lxI) = sub2ind(size(B),xI,a(2)*ones(lxI,1)); iL = iL+lxI;
            B(xI,a(2)) = 1; 
            sbPos(i,:) = [midx midy-1 midx midy+1];
        else
            lyI = length(yI); 
            L(iL+1:iL+lyI) = sub2ind(size(B),a(1)*ones(lyI,1),yI); iL = iL+lyI;
            B(a(1),yI) = 1; 
            sbPos(i,:) = [midx-1 midy midx+1 midy];
        end
    case 0, % two sucssecive row or col
        if my, 
            l2 = floor(length(yI)/2); lyI = length(yI);
            L(iL+1:iL+l2+1) = sub2ind(size(B),a(1)*ones(l2+1,1),yI(1:l2+1)); iL = iL+l2+1;
            B(a(1),yI(1:l2+1)) = 1; 
            
            l3 = lyI-l2;
            L(iL+1:iL+l3) = sub2ind(size(B),b(1)*ones(l3,1),yI(l2+1:end)); iL = iL+l3;
            B(b(1),yI(l2+1:end)) = 1;
            sbPos(i,:) = [midx midy+1 midx+1 midy-1];
        else
            l2 = floor(length(xI)/2); lxI = length(xI);
            L(iL+1:iL+l2+1) = sub2ind(size(B),xI(1:l2+1),a(2)*ones(l2+1,1)); iL = iL+l2+1;
            B(xI(1:l2+1),a(2)) = 1; 
            
            l3 = lxI-l2;
            L(iL+1:iL+l3) = sub2ind(size(B),xI(l2+1:end),b(2)*ones(l3,1)); iL = iL+l3;
            B(xI(l2+1:end),b(2)) = 1;
            sbPos(i,:) = [midx+1 midy midx-1 midy+1];
        end
    otherwise,
        if (mx > my),
            m = floor(mx/my);
            n = mx - m*my;
            r = floor((1:n)*my/n);
            m = m * stx;
            s = a(1);
            k = 1; t=1;
            for j = a(2)+sty:sty:b(2)-sty
                if(n && k==r(t)), jump = stx; t=t+1; else jump = 0; end
                k = k+1;
                ii = (s:stx:s+m+jump)'; lii = length(ii);
                L(iL+1:iL+lii) = sub2ind(size(B),ii,j*ones(lii,1)); iL = iL+lii;
                B(ii,j) = 1;
                s = s+m+jump;
            end
            B(a(1)+stx,a(2)) = 1;
            B(b(1),b(2)-sty) = 1;
            B(b(1)-stx,b(2)) = 1;
            L(iL+1:iL+3) = sub2ind(size(B),[a(1)+stx b(1) b(1)-stx],[a(2) b(2)-sty b(2)]); iL=iL+3;
            sbPos(i,:) = [midx midy+2 midx midy-2];
            m0 = midy+1; while B(midx,m0), m0 = m0+1; end 
            sbPos(i,1:2) = [midx m0];
            m0 = midy-1; while B(midx,m0), m0 = m0-1; end 
            sbPos(i,3:4) = [midx m0];
        else 
            m  = floor(my/mx);
            n = my - m*mx;
            r = floor((1:n)*mx/n);
            m = m * sty;
            s = a(2);
            k = 1; t=1;
            for j = a(1)+stx:stx:b(1)-stx
                if(n && k==r(t)), jump = sty; t=t+1; else jump = 0; end
                k=k+1;
                ii = (s:sty:s+m+jump)'; lii = length(ii);
                L(iL+1:iL+lii) = sub2ind(size(B),j*ones(lii,1),ii); iL = iL+lii;
                B(j,ii) = 1;
                s = s+m+jump;             
            end
            B(a(1),a(2)+sty) = 1;
            B(b(1),b(2)-sty) = 1;
            B(b(1)-stx,b(2)) = 1;
            L(iL+1:iL+3) = sub2ind(size(B),[a(1),b(1),b(1)-stx],[a(2)+sty,b(2)-sty,b(2)]); iL=iL+3;
            m0 = midx+1; while B(m0,midy), m0 = m0+1; end 
            sbPos(i,1:2) = [m0 midy];
            m0 = midx-1; while B(m0,midy), m0 = m0-1; end 
            sbPos(i,3:4) = [m0 midy];
        end
end
end
B = logical(sparse(B));
L = L(1:iL);
end

% primal or control volume of upscaling as dual mesh
function NET = primal(G,iB)
[nx,ny] = size(G);
NET.EdgeI = G==1;
NET.VerI  = G==2;
%% blockAII
G(NET.EdgeI) = -2;
G(NET.VerI)  = -1;
v0 = 0;
i0 = 1;
j0 = 1;
flag0 = 1;
while flag0,
    flag0 = 0;
    for i= i0:nx
        for j = j0:ny, if ~G(i,j), flag0 = 1; break; end; end
        if flag0, i0 = i; j0 = j; break; end
        j0 = 1;
    end
    if flag0,
        u=0;
        % check down if exists!
        if (j>1 && G(i,j-1)>0), u = G(i,j-1); end
        if (i>1 && G(i-1,j)>0), u = G(i-1,j); end
        %check right if exists!
        [~,r,v] = find(G(i,j:end),1);
        if (~isempty(v) && v>0), u = v; end
        %set the value
        if ~u, v0 = v0 + 1; u = v0; end
        %fill right and top of i,j with u         
        if isempty(r), r=ny-j+2; end %right
        G(i,j:j+r-2) = u;
        % fill top
        r = find(G(i:end,j)==-2,1);
        if isempty(r), r=nx-i+2; end
        G(i:i+r-2,j) = u;          
    end
end
for i=1:v0, if length(find(G==i)) < 5, G(G==i) = -2; G(G>i) = G(G>i) - 1; end; end
G(iB+1)=-2; G(iB-1)=-2; G(iB+nx)=-2; G(iB-nx)=-2;
NET.EdgeI = find(G==-2);
NET.VerI  = find(NET.VerI);
G(NET.EdgeI) = 0;  %Edges
G(NET.VerI)  = -1; %Vertices
NET.G = G;
maxG = max(G(:));
%% Ordering
VerI = NET.VerI;
nV   = numel(VerI);
EdgeI = NET.EdgeI;
nE    = numel(EdgeI);

nI = numel(G) - nV - nE;
InnerI = zeros(nI,1);
lk = zeros(maxG,1);
k = 0;
for i=1:maxG
    lI = find(G==i);
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
