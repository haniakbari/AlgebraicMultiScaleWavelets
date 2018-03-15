clc; clear all; close all; %#ok
%%
load K.mat K;
h = [10*.3048 10*.3048 2*.3048];
layer = 76; % 75 76 80 85
Kn = K(:,:,layer); %myPlot(log10(Kn),'b');
test = 0;
uGrid = upsGrid(Kn,h,2,2,test);
save uGrid.mat uGrid;
K = uGrid.K; %myPlot(log10(K),'b');
pureNeum = 1; 
makeFineData(K,h,pureNeum); % save fineData.mat K h A q u v;
%%
global A L U u b nx ny;
load fineData.mat K h A Q u v b;
%%
% Modify A to enforce source terms (Pure Neumann Problem)
%if modA, Ao = A; for i=1:length(Q.ind), a = Q.ind(i); A(a,:) = 0; A(a,a) = 1; end; end
lutol = 1e-6;
itol = 1e-6;
iter = 150;
x0   = zeros(length(b),1);
fprintf('ILU(%g) decomposition\n',lutol);
tic; [L,U] = ilu(A,struct('type','ilutp','droptol',lutol)); toc;
display(nnz(U)); whos U;
%tic; [L,U] = ilu(A,struct('droptol',lutol)); toc;
%%
tic;
[ugmr,erit] = DefCorrectLU(itol,iter,x0); 
%ugmr = gmres(A,q,[],1e-5,30,L,U);
vgmr = getFlux(K,ugmr,h);
toc;
fprintf('\n\tFlux Error: Rech:ILU(%g):\n\t\t x:%.2e y:%.2e\n',...
             lutol,norm(vgmr.x-v.x)/norm(v.x),norm(vgmr.y-v.y)/norm(v.y));
%%
% STRUCTURED AMS
%
nxc = 4; nyc = 12;
[nx,ny]  = size(K);
NetB = setNetBStruct(nx,ny,nxc,nyc);
%%
%Interpolation of ORIGINAL A ??
fprintf('\n*STRUCTURED #Icols %d\n',NetB.nIEV(3));
tic; Is = TI_AMS(A,NetB); Rs  = Is'; Acs = Rs*A*Is; toc;
fprintf('\n STRUCT : AMS+ILU(%g):\n',lutol);
%%
tic;
%[usamslu,eritAMS] = DefCorrectAMS(itol,iter,x0,Is,Rs,Acs);
usamslu = gmres(A,b,[],itol,iter,@(r) Mgl(r,Acs,Rs,Is));
%usamslu = gmres(A,b,[],itol,iter,L,U);
vsamslu = getFlux(K,usamslu,h);
toc;
fprintf('\n\tFlux Error: AMS+ILU:\n\t\t x:%.2e y:%.2e\n',...
    norm(vsamslu.x-v.x)/norm(v.x),norm(vsamslu.y-v.y)/norm(v.y));
%%
% UNSTRUCTURED MESH
load uGrid.mat uGrid;
NET = uSetAllNetMain(uGrid);
%%
fprintf('\n*UNSTRUCTURED: Duals are centered connected\n');
Iuu = TI_AMS(A,NET); 
%%
Iu = Iuu;%.*NET.Supp;
Ru  = Iu'; Acu = Ru*A*Iu;
fprintf('\n UNSTRUCT : AMS+ILU(%g): #Icols %d\n',lutol,size(Iu,2));
%%
tic;
[uuamslu,eritUAMS] = DefCorrectAMS(itol,iter,x0,Iu,Ru,Acu);
%usamslu = gmres(A,b,[],itol,iter,@(r) Mgl(r,Acu,Ru,Iu));
vuamslu = getFlux(K,uuamslu,h);
toc;
fprintf('\n\tFlux Error: AMS+ILU:\n\t\t x:%.2e y:%.2e\n',...
    norm(vuamslu.x-v.x)/norm(v.x),norm(vuamslu.y-v.y)/norm(v.y));
%%
% load uGrid.mat uGrid;
load uGrid.mat uGrid;
[nx,ny]  = size(K);
volNet = gridConnectionVolShifted(nx,ny,uGrid.pos);
%%
tic; Iv = TI_AMS(A,volNet); Rv  = Iv'; Acv = Rv*A*Iv; toc;
fprintf('\n VOL : AMS+ILU(%g): #Icols %d\n',lutol,size(Iv,2));
%%
tic;
[volu,eritUvAMS] = DefCorrectAMS(itol,iter,x0,Iv,Rv,Acv);
volv = getFlux(K,volu,h);
toc;
fprintf('\n\tFlux Error VOL: AMS+ILU:\n\t\t x:%.2e y:%.2e\n',...
    norm(volv.x-v.x)/norm(v.x),norm(volv.y-v.y)/norm(v.y));
%%
%{
close all;
figure(); colormap 'jet';
maxK = abs(log10(min(K(:))));
logK = log10(K)/maxK;
subplot 241; imagesc(logK); title('N*log10(K)');
subplot 245; imagesc(u);  title('fine P');
subplot 242; imagesc((B==1)+logK); title('Struct');
subplot 246; imagesc(abs(u-reshape(usamslu,nx,ny))/norm(u(:)));title('RSU');
subplot 243; imagesc((C==1)+logK); title('UnStructV');
subplot 247; imagesc(abs(u-reshape(uuamslu,nx,ny))/norm(u(:))); title('RUnUV');
subplot 244; imagesc((D==1)+logK); title('UnStructD');
subplot 248; imagesc(abs(u-reshape(vuuamslu,nx,ny))/norm(u(:))); title('RUnUD');
%}
