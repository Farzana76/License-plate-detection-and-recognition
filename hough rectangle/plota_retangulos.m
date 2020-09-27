function void=plota_retangulos(x, retangulos, centros);
%
% plota_retangulos plota os retangulos sobrepostos a imagem x
%  
retangulos_new=[];
centros_new=[];
ss=size(retangulos);ss=[ss 1]; % aumenta o valor de ss para evitar erros no caso de um retangulo
figure,imshow(x);
hold on;for i=1:ss(3);
    xx=retangulos(:,:,i);
    centro=centros(i,:);
     set(gcf,'Color',[.99 .99 .99]);plot(centro(2),centro(1),'+','Color',[.99 .99 .99]);
%    plot(xx(2,:),xx(1,:),'w'); 
    line(xx(2,:),xx(1,:),'Color',[.99 .99 .99],'LineWidth',2);axis ij; 
%    line(xx(2,:),xx(1,:),'Color','k','LineWidth',1);axis ij; 
    
end