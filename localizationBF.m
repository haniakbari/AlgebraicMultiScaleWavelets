% Localize propageted basis functions
%%
load IuNet.mat Iu NET; 
 G = NET.G;
 [nx,ny] = size(G);
% puu = reshape(sum(Iu,2),size(G)); figure(); imagesc(puu); axis xy;
%%
nbf = size(Iu,2);
[xD,yD,~] = find(~G);
bfset=[7,10,14,27];
cols = 2;
rows = 2;
k = 0;
figure();
for i=bfset
    bf = reshape(Iu(:,i),[nx,ny]);
    k=k+1; subplot (rows,cols,k);
    imagesc(bf); colormap 'jet'; hold on;
    plot(yD,xD,'w.'); axis xy;
    title(num2str(i));
end