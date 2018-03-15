B = NET.G <= 0;
I = Iu;
[x,y,~] = find(B);
nbf = size(I,2);
%%
figure();
for i=1:nbf-1
a = I(:,i);
bf = reshape(a,[nx,ny]);
subplot (floor(nbf/3),3,i);  
imagesc(bf); title(num2str(i)); axis off; hold on;
plot(y,x,'w.');
end
%%
 % partition of unity
 figure();
 a = reshape(sum(I,2),[nx,ny]); surf(a); shading flat; axis off; hold on;
 r = sub2ind([nx,ny],x,y);
 title('partition of unity');
 plot3(y,x,a(r),'w.');  drawnow;
 %%
i = 15;
a = I(:,i); bf = reshape(a,[nx,ny]);
figure(); subplot 211,
imagesc(bf); title(num2str(i)); axis off; hold on;
plot(y,x,'w.');
subplot 212, 
i = 17;
a = I(:,i); bf = reshape(a,[nx,ny]);
imagesc(bf); title(num2str(i)); axis off; hold on;
plot(y,x,'w.');
