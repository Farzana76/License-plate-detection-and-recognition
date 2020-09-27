function [img_out, retangulos,simetry] = window_hough(img, window, nT, nS, D)

% usage img_out = window_hough(img, window, nT, nS, Tn, angulo, D)
%
% img = imagem em grayscale;
% window -> vetor de duas posiç~oes [altura largura]
% nS = discretizaç~ao do rho
% nT = discretizaç~ao do theta
%Tn = limiar para detecç~ao dos pontos na imagem H
Tn=.5;
angulo=90;
radius=(window(1)-1)/2;
if nargin==3,
    D=nT;
    nT=round(pi/atan(1/(2*radius+1)));
    nS=round(window(1)*sqrt(2));
end,
tic  % inicia o contador de tempo

TAM = 3;
[H_img, W_img] = size(img);
%img = double(img);
%img = img/max(img(:));


%tamanho da janela
H_window = window(1);
W_window = window(2);

% janela deve ter dimenç~oes ´impares!
CENTRO = [ceil(H_window/2), ceil(W_window/2)];



index = 0;

% reduz a quantidade pixels para a varredura
bw = edge(img, 'canny',[.1 .2],1.5);figure,imshow(bw)
%bw = edge(img, 'canny',[],2);
%bw = edge(img, 'sobel',.03);
%bw=img;
figure,imshow(bw);hold on;
if (isempty(D)),
    D = dt(bw);
end;


%regmax = imregionalmax(D);
regmax=D;

s=size(regmax);

% 
%aumenta a regiao de busca
%

%mask=[1 0 1 0 1;0 1 0 1 0;1 0 1 0 1;0 1 0 1 0;1 0 1 0 1];
%mask=ones(3);
%regmax=imdilate(regmax,mask);


%   figure, imshow(regmax|bw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mat_x=repmat((1:s(1))',[1 s(2)]);
mat_y=repmat(1:s(2),[s(1) 1]);
%
% cria um grid com espacamento esp para busca de retangulos
%
esp=3;
regmax=(mod(mat_x,esp)==0)&(mod(mat_y,esp)==0);
regmax(mat_x-CENTRO(1)<0)=0;regmax(mat_y-CENTRO(2)<0)=0;
regmax(mat_x+CENTRO(1)-1>s(1))=0;regmax(mat_y+CENTRO(2)-1>s(2))=0;
ix=mat_x(regmax>0);
iy=mat_y(regmax>0);
edgedata=[ix iy]';
%save lixo edgedata

n_pontos = sum(regmax(:));
text = sprintf('Performing windowing... %d janelas!', n_pontos);
%   hw=waitbar(0,text);

%figure,imshow(img);hold on;

count=1;
h = waitbar(0,'Please wait...');
%for n=1:length(edgedata),               
for n=1:n_pontos,
    
   %      if (((edgedata(1,n)-CENTRO(1)+1)>0)&((edgedata(1,n)+CENTRO(1)-1)<s(1))&((edgedata(2,n)-CENTRO(2)+1)>0)&((edgedata(2,n)+CENTRO(2)-1)<s(2)))                      
    %                WIN = img((edgedata(1,n)-CENTRO(1)+1):(edgedata(1,n)+CENTRO(1)-1), (edgedata(2,n)-CENTRO(2)+1):(edgedata(2,n)+CENTRO(2)-1));            

    WIN = bw((edgedata(1,n)-CENTRO(1)+1):(edgedata(1,n)+CENTRO(1)-1), (edgedata(2,n)-CENTRO(2)+1):(edgedata(2,n)+CENTRO(2)-1));            
    %
    %  se o numero de bordas eh maior do que Min_edges, calcula Hough
    %
    Min_edges=60;
    if sum(WIN(:))>Min_edges, 
        [H,berough,m,b, rho, theta, pico]=hough_kittler(WIN,nT,nS,Tn, 1);
        if ~isempty(H),
            [PONTOS,simetria] = find_rectangles(img, WIN, H, berough, rho, theta, nT, nS, b, m, pico, [edgedata(1,n), edgedata(2,n)]);
            if ~isempty(PONTOS),
                retangulos(:,:,count)=PONTOS;
                simetry(count,:)=simetria;
                count=count+1;
            end,
        end,
    end,
    waitbar(n/length(edgedata),h)
end;
close(h);
img_out = H;

toc  %marca o final do contador de tempo