function makeFineData(K,h,pureNeum)
% [nx,ny] = size(K);  Q.ind = [sub2ind(size(K),4,4) (sub2ind(size(K),nx-4,ny-4))];
Q.ind = [1 numel(K)];
if pureNeum
    fprintf('Pure Neumann problem with source term\n');
    Q.val = [1 -1];
    [u,v,A,b] = BTPFA(K,Q,h);
else
    fprintf('Pressure is 1 and 0 on left and right with ZERO source term\n');
    Q.val = [0 0];
    [u,v,A,b] = BTPFAL1R0(K,Q,h);
end
save fineData.mat K h A Q u v b;
end