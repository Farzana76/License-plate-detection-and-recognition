function bool=validate_rectangle(x,mag,paralelo,Tperim);
%
% bool=validate_rectangle(x,mag,paralelo,Tperim) retorna o valor 1 se os pontos
% contidos na matriz paralelo correspondem de fato a um paralelogramo.
% Para tal, uma metrica de homogeneidade na imagem original x eh calculada.
%
% Tperim eh um limiar de perimietro, 
%
%
bool=logical(0);
bw=roipoly(x,round(paralelo(2,:)),round(paralelo(1,:)));
perim=bwperim(bw);
edg=imdilate(perim,strel('diamond',2));
erro=-(sum(mag(edg))-sum(perim(perim>0)))/sum(perim(perim>0));
if  erro<Tperim,
    bool=logical(1);
end,
