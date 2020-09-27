function [retangulos_new, centros_new]=remove_duplicated_paa(img, retangulos, centros, indice, simetry);
%
% [retangulos_new, centros_new]=remove_duplicated_rectangles2(img, retangulos, centros, indice, simetry)
% remove retangulos duplicados, pegando aquele que retorna o menor erro de acordo com
% os valores de simetry
%  
retangulos_new=[];
centros_new=[];
%
% aplica fechamento morfologico para unir centros proximos
%
bw=imdilate((indice>0),ones(5));
%figure,imshow(bw);
%
% busca os clusters de centros
%
[bw,num]=bwlabel(bw);
%
% para cada cluster, acha o erro minimo
%
for i=1:num,
    lista=indice(bw==i);lista=lista(lista>0);
    erro=[];
    for j=1:length(lista);
        erro(j)=norm([1 1 2 2 1].*simetry(lista(j),:));
    end,
    [minimo, onde]=min(erro);onde=onde(1);
    retangulos_new(:,:,i)=retangulos(:,:,lista(onde));
    centros_new(i,:)=centros(lista(onde),:);
end,
