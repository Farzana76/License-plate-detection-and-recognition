function [retangulos, centros] = detect_rectangles_paa(img, RMin, RMax, amax, D)

% usage img_out = detect_rectangles_paa(img, RMin, RMax, amax, D)
%
% img = edge (binary) image;
% RMin e RMax -> internal and external radii of the ring-like region
% (sliding window)
% amax -> largest possible value for the smallest rectangle size
% D -> binary window containing the centers of the search regions
% (optional)

%
%  artigo submetido ao PAA 2005
%
%

amin=2*RMin;
indice=zeros(size(img));
PONTOS=[];
centros=[];
%RMin=6;
%dT=.56RMax; % discretizacao no angulo
%dS=1/sqrt(2); % discretizacao na distancia
dT=pi/2/(2*RMax-1);
%dS=pi/4;
dS=pi/4;
tic  % inicia o contador de tempo

TAM = 3;
[H_img, W_img] = size(img);
CENTRO = [RMax+1, RMax+1];


index = 0;

% 
% calculo do mapa de bordas
%
%bw = edge(img, 'canny',[.3 .31],.5);figure,imshow(bw);
%bw = edge(img, 'canny',[.1 .2],2);figure,imshow(bw)
%bw = edge(img, 'canny');
%bw = edge(img, 'sobel',.03);
bw=img;
%bw=img;bw(1,:)=1;bw(256,:)=1;bw(:,1)=1;bw(:,256)=1;
%ra=rand(size(bw));
%bw(ra<.01)=1;
%bw(ra>.7)=0;
%save lixo bw;
%figure,imshow(img);hold on
%
%  elimina componentes conexos com poucos elementos e calcula a transformada distancia
%  o espaco de busca serah os maximos locais da TD maiores do que RMin e menores do que RMax
%
%
if nargin==4,
%    mask=[1 0 1 0 1 0 1;0 1 0 1 0 1 0;1 0 1 0 1 0 1;0 1 0 1 0 1 0;1 0 1 0 1 0 1];  % 7 x 7
%    mask=[1 0 1 0 1;0 1 0 1 0;1 0 1 0 1;0 1 0 1 0;1 0 1 0 1];  % 5 x 5
%     mask=[0 0 1 0 0;0 1 0 1 0;1 0 1 0 1;0 1 0 1 0;0 0 1 0 0];  % diamond 5 x 5
    mask = [1 0 1;0 1 0;1 0 1];
%    mask=strel('diamond',2);
    nbw=bw;
    [lab,num]=bwlabel(nbw);
    for i=1:num,
        if sum(nbw(lab==i))<RMin,
            nbw(lab==i)=0;
        end,
    end,
    d=bwdist(nbw);
 %   g=fspecial('gaussian',3);
 g=fspecial('gaussian',7);
    d=filter2(g,d,'same'); % filtro com gaussiana
%    regmax=imdilate(imregionalmax(filter2(ones(3),d,'same'))&(d>=RMin&d<=amax/2),mask);
    regmax=imdilate(imregionalmax(uint8(d))&(d>=RMin&d<=amax/2),mask);
%regmax=imregionalmax(uint8(d))&(d>=RMin&d<=amax/2);
    
clear d lab num
   % bw=nbw;
else,
%    regmax=imdilate(D,strel('diamond',3));
regmax=D;
end
sh=zeros([size(img) 3]);
temp=sh(:,:,3);temp(bw==1)=255;sh(:,:,3)=temp;
temp=sh(:,:,2);temp(regmax==1)=255;sh(:,:,2)=temp;
sh=uint8(sh);
figure,imshow(sh);hold on;
%clear sh

s=size(img);

mat_x=repmat((1:s(1))',[1 s(2)]);
mat_y=repmat(1:s(2),[s(1) 1]);
%if isempty(D),
    %
    % cria um grid com espacamento esp para busca de retangulos
    %
%    esp=1;
%    regmax=(mod(mat_x,esp)==0)&(mod(mat_y,esp)==0);
%end,
    regmax(mat_x-CENTRO(1)<0)=0;regmax(mat_y-CENTRO(2)<0)=0;
    regmax(mat_x+CENTRO(1)-1>s(1))=0;regmax(mat_y+CENTRO(2)-1>s(2))=0;

%
% edgedata contem as informacoes dos centros
%
ix=mat_x(regmax>0);
iy=mat_y(regmax>0);
edgedata=[ix iy]';
%save lixo edgedata

n_pontos = sum(regmax(:));
text = sprintf('Performing windowing... %d janelas!', n_pontos);


%figure,imshow(img);hold on;

%
%  mascara circular para reuzir as bordas 
%
mat_x=repmat((1:(2*RMax+1))',[1 2*RMax+1]);
mat_y=repmat(1:(2*RMax+1),[2*RMax+1 1]);
mask=((mat_x-CENTRO(1)).^2+(mat_y-CENTRO(2)).^2)<=RMax^2&((mat_x-CENTRO(1)).^2+(mat_y-CENTRO(2)).^2)>=RMin^2;

count=1;
%h = waitbar(0,'Please wait...');
for n=1:n_pontos,
    %
    % pega janela retangular RMax x RMax em torno de cada pontod e interesse
    %
  %  i1=max(1,(edgedata(1,n)-RMax));i2=min(s(1),(edgedata(1,n)+RMax));
  %  i3=max(1,(edgedata(2,n)-RMax));i4=min(s(2),(edgedata(2,n)+RMax));
  %  WIN=bw(i1:i2,i3:i4);
    WIN = bw((edgedata(1,n)-RMax):(edgedata(1,n)+RMax), (edgedata(2,n)-RMax):(edgedata(2,n)+RMax));            
    %
    % pega apenas a regiao circular com raio RMax
    %
    WIN(mask==0)=0;
%    save lixo WIN
 %   figure,imshow(WIN)
    %
    %  se o numero de bordas eh maior do que Min_edges, calcula Hough
    %
%    figure,imshow(WIN);pause
 %   Min_edges=20;
    Min_edges=4*(2*RMin+1);
    if sum(WIN(WIN>0))>Min_edges, 
%        save lixo WIN dT dS
        [H,berough, rho, theta, pico, valor_pico]=roda_hough_paa(WIN,dT,dS,RMin);
        if ~isempty(H),
            [PONTOS,simetria] = find_rectangles_paa(img, dT, H, rho, theta, pico, valor_pico, [edgedata(1,n), edgedata(2,n)]);
            if ~isempty(PONTOS),
                retangulos(:,:,count)=PONTOS;
                simetry(count,:)=simetria;
                centros(count,:)=edgedata(:,n)';
                indice(edgedata(1,n),edgedata(2,n))=count;
                count=count+1;
            end,
        end,
   end,
%    waitbar(n/length(edgedata),h)
end;
%close(h);
%img_out = H;
%figure,imshow(bw);hold on;
if ~isempty(retangulos),
    [retangulos,centros]=remove_duplicated_paa(img, retangulos, centros, indice, simetry);
end,
%
% validacao dos retangulos
%
size(retangulos)

Tperim=.3;
s=size(retangulos);
count=0;
if length(s)==2,
    s=[s 1];
end,
figure,imshow(img);
if ~isempty(retangulos);
    for k=1:s(3),
        parcand=[retangulos(1,:,k);retangulos(2,:,k)];
        bool=validate_rectangle(img,bw,parcand,Tperim);
        if bool,
            count=count+1;
            ret_new(:,:,count)=retangulos(:,:,k);
            cent_new(count,:)=centros(k,:);
            % ss=size(retangulos_new);ss=[ss 1]; % aumenta o valor de ss para evitar erros no caso de um retangulo
            hold on;
            xx=ret_new(:,:,count);
            centro=cent_new(count,:);
            line(xx(2,:),xx(1,:),'LineWidth',2,'Color','g'); axis ij; hold on;
            plot(centro(2),centro(1),'r+');
        end,
    end,
    retangulos=ret_new;
    centros=cent_new;
else
    retangulos=[];
    centros=[];
end;

t=toc;  %marca o final do contador de tempo
minutos=floor(t/60);
segundos=(t-minutos*60);
disp(sprintf('tempo decorrido: %2d minutos e %2d segundos',minutos,round(segundos)))
