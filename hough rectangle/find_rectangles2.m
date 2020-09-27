function [PONTOS,simetria] = find_rectangles(img_original, dT, H, rho, theta, pico, pic, pos_xy)

%ret = find_simetric_lines(img, H, berough, rho, theta, nT, nS, b, m)
% H  = imagem HOUGH
% berough = imagem HOUGH realçada pela hough_kitler
% rho = rho's retornados pelo kittler
% theta = theta's retornados pelo kittler
% dT eh o espacamento dos angulos theta
% pic eh a altura dos picos
% ret = boolean.  1 se encontrou uma estrutura

%figure,imshow(H/max(H(:)));
angulo = 90;
ret = 0;
ang = 180/pi*theta';
rect=[];
simetria=[];
% 
%   procura retangulos se ha mais do que 3 picos 
%
%
%rho_and_theta = [rho ; 180/pi*theta'],%
if length(theta)>3,
%    for i=1:length(pico);
%        pic(i)=H(pico(i,1),pico(i,2));
%    end,
 %   rho_and_theta = [rho;  180/pi*theta pic],%
    
    %centro = ceil(size(img)/2);
 %   L_theta_intra = 360*dT/pi;  %limiar do angulo para o mesmo par
 %   L_theta_inter = 360*dT/pi/1.5;   %limiar do angulo entre dois pares
 %   L_rho = 2;   %limiar da distancia 
 %   L_pico = 50; %limiar de tamanho dos lados dos retangulos (em %)
    
    L_theta_intra = 3;   %limiar do angulo para o mesmo par
    L_theta_inter = 3;   %limiar do angulo entre dois pares
    L_rho = 3;   %limiar da distancia 
    L_pico = 25; %limiar de tamanho dos lados dos retangulos (em %)

    
    new_pic = []; 
    new_rho = [];
    new_theta = [];
    cont = 0;
    %
    %  procura por pares de picos com mesmo theta, e com rhos simetricos (e aproximadamente mesma altura
    % de pico)
    %
    for index1=1:length(rho);
        for index2=(index1+1):length(rho);
            diferenca_theta = abs(ang(index1)-ang(index2));
            if (diferenca_theta<=L_theta_intra)
                %            ang(index1)
                diferenca_rho = abs(rho(index1)+rho(index2));
                diferenca_pico=abs(H(pico(index1,1),pico(index1,2))-H(pico(index2,1),pico(index2,2)));
                if (diferenca_rho<=L_rho) & (diferenca_pico<=((H(pico(index1,1),pico(index1,2))+ H(pico(index2,1),pico(index2,2)))*L_pico/200)),
                    cont = cont+1;
                    %
                    % pega a m´edia dos angulos e das distancias
                    %
                    new_theta(cont) = (theta(index1)+theta(index2))/2;
                    new_rho(cont) = abs(rho(index1)-rho(index2))/2;
                    new_pic(cont)=(pic(index1)+pic(index2))/2;
                    %
                    % armazena as diferencas de angulos e distancias
                    %
                    dtheta(cont)=diferenca_theta;
                    drho(cont)=diferenca_rho;
                end
            end
        end
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho = new_rho;
    theta = new_theta;
    pic=new_pic;
    %[rho;  180/pi*theta;pic],%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ang = 180/pi*theta';
    indice = 0;
    %
    %   verifica se eh retangulo (ou seja, busca pares de picos com \delta\theta=90) e se
    %   o comprimento de um segmento eh compativel com a distancia entre os segmentos
    %   ortogonais
    %   
    L_comp=.6;
    if (length(ang)>1)
        for n1=1:length(ang),
            n_ang=0; % numero de linhas com inclinacao do (angulo) em referencia a linha (n)
            for n2=(n1+1):length(ang),
                diferenca_ang = abs(ang(n1)-ang(n2))-angulo;
%                diferenca_ang = abs(ang(n2)-ang(n1)-angulo;
                if abs(diferenca_ang)<=L_theta_inter & pic(n1)>L_comp*rho(n2) & pic(n2)>L_comp*rho(n1),
                    indice = indice+1;
                    rect(:,:,indice)=[theta(n1)+diferenca_ang*pi/360 theta(n2)-diferenca_ang*pi/360;rho(n1) rho(n2)]; % armazena rhos e thetas pros retangulos
                    %
                    % medida de simetria de retangulo (um por centro)
                    %
                    simetria=[dtheta(n1) dtheta(n2) drho(n1) drho(n2) abs(diferenca_ang)];
%                    [rho(n1) rho(n2);pic(n1) pic(n2); pic(n1)/rho(n2) pic(n2)/rho(n1)],
                end
            end;
        end;
    end;
    
    
end,
PONTOS = [];

%
% plota o centro da regiao de busca
%


centrox=pos_xy(1);
centroy=pos_xy(2);

%plot(centroy,centrox,'b.'); hold on;

%
% acha as intersecçoes das retas e plota os retangulos
%


if ~isempty(rect),
    ss=size(rect);
    ss=[ss 1];
    for i=1:ss(3),
        
        %  pico,
        ret = 1;
        %        disp(' Retangulo detectado: ');
        %        s=size(img);
        %        centro = ceil(s(1)/2);
        theta=[rect(1,1,i) rect(1,1,i)];theta=[theta rect(1,2,i) rect(1,2,i)];
        rho=[rect(2,1,i) -rect(2,1,i)];rho=[rho rect(2,2,i) -rect(2,2,i)];
        %        [180/pi*theta;rho],%
        cont = 1;
        A=[cos(theta(1)) sin(theta(1));cos(theta(3)) sin(theta(3))];b=[rho(1) rho(3)]';
        P=inv(A)*b;P=P+[centrox centroy]';P1=P;
        A=[cos(theta(2)) sin(theta(2));cos(theta(3)) sin(theta(3))];b=[rho(2) rho(3)]';
        P=inv(A)*b;P=P+[centrox centroy]';P2=P;
        A=[cos(theta(2)) sin(theta(2));cos(theta(4)) sin(theta(4))];b=[rho(2) rho(4)]';
        P=inv(A)*b;P=P+[centrox centroy]';P3=P;
        A=[cos(theta(1)) sin(theta(1));cos(theta(4)) sin(theta(4))];b=[rho(1) rho(4)]';
        P=inv(A)*b;P=P+[centrox centroy]';P4=P;
        x=[P1 P2 P3 P4 P1]; 
        plot(x(2,:),x(1,:),'r'); axis ij; hold on;
%        line(x(2,:),x(1,:),'LineWidth',4,'Color','g'); axis ij; hold on;
        plot(centroy,centrox,'r.'); hold on;
        PONTOS=x;
    end,
end;
