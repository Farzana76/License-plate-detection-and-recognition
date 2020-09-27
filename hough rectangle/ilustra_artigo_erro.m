load lixo
figure,
axes('position',[.1 .25 .55 .55]), 
imshow(x);
ss=size(retangulos);ss=[ss 1]; % aumenta o valor de ss para evitar erros no caso de um retangulo
hold on;for i=1:ss(3);
    xx=retangulos(:,:,i);
    centro=centros(i,:);
    line(xx(2,:),xx(1,:),'LineWidth',2,'Color','k'); axis ij; hold on;
    plot(centro(2),centro(1),'k+');
end
title('(a)');
h=h/max(h(:));
h(h>.5)=1;

axes('position',[.6 .25 .25 .25]), imshow(1-h);
line([1 1 94 94 1],[1 94 94 1 1],'Color','k');
title('(c)');
axes('position',[.6 .55 .25 .25]), imshow(not(WIN));
line([1 1 71 71 1],[1 71 71 1 1],'Color','k');
title('(b)');
