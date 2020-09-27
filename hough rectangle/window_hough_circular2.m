function [img_out, retangulos, centros] = window_hough_circular2(img, RMin, RMax, D)

% usage img_out = window_hough(img, RMax, nT, nS, Tn, D)
%
% img = imagem em grayscale;
% RMax -> raio da regiao circular (janela)
% D = imagem binaria contendo os centros das janelas

indice=zeros(size(img));
PONTOS=[];
centros=[];
if nargin==3,
    D=[];
end
%RMin=6;
%dT=.56RMax; % discretizacao no angulo
%dS=1/sqrt(2); % discretizacao na distancia
dT=3*pi/4/(2*RMax+1);
dS=3/4;
tic  % inicia o contador de tempo

TAM = 3;
[H_img, W_img] = size(img);
CENTRO = [RMax+1, RMax+1];


index = 0;

% 
% calculo do mapa de bordas
%
%bw = edge(img, 'canny',[.3 .31],.5);figure,imshow(bw);
%bw = edge(img, 'canny',[.1 .2],1);figure,imshow(bw)
%bw = edge(img, 'canny');
%bw = edge(img, 'sobel',.03);
bw=img;
%bw=img;bw(1,:)=1;bw(256,:)=1;bw(:,1)=1;bw(:,256)=1;
%ra=rand(size(bw));
%bw(ra<.01)=1;
%bw(ra>.7)=0;
%save lixo bw;
figure,imshow(bw);hold on;
%figure,imshow(img);hold on
regmax=imdilate(D,ones(3));

s=size(img);

mat_x=repmat((1:s(1))',[1 s(2)]);
mat_y=repmat(1:s(2),[s(1) 1]);



% 
%aumenta a regiao de busca
%
if ~isempty(D),
%    mask=[1 0 1 0 1;0 1 0 1 0;1 0 1 0 1;0 1 0 1 0;1 0 1 0 1];
%    mask=[1 0 1;0 1 0;1 0 1];
   mask=ones(1);
   regmax=imdilate(regmax,mask);
else,
    %
    % cria um grid com espacamento esp para busca de retangulos
    %
    esp=1;
    regmax=(mod(mat_x,esp)==0)&(mod(mat_y,esp)==0);
end,
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
h = waitbar(0,'Please wait...');
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
    if sum(WIN(:))>Min_edges, 
%        save lixo WIN dT dS
        [H,berough, rho, theta, pico, valor_pico]=roda_hough_extended2(WIN,dT,dS,RMin);
        if ~isempty(H),
            [PONTOS,simetria] = find_rectangles2(img, dT, H, rho, theta, pico, valor_pico, [edgedata(1,n), edgedata(2,n)]);
            if ~isempty(PONTOS),
                retangulos(:,:,count)=PONTOS;
                simetry(count,:)=simetria;
                centros(count,:)=edgedata(:,n)';
                indice(edgedata(1,n),edgedata(2,n))=count;
                count=count+1;
            end,
        end,
   end,
    waitbar(n/length(edgedata),h)
end;
close(h);
img_out = H;
if ~isempty(retangulos),
    [retangulos,centros]=remove_duplicated_rectangles2(img, retangulos, centros, indice, simetry);
end,


t=toc;  %marca o final do contador de tempo
minutos=floor(t/60);
segundos=(t-minutos*60);
disp(sprintf('tempo decorrido: %2d minutos e %2d segundos',minutos,round(segundos)))
