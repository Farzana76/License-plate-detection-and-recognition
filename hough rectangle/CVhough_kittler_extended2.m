function [H,Ts,rho]=CVhough_extended2(edgedata,dT,dS)
%CVhough Hough transform of a binary matrix
%
%function [H,Ts,rho]=CVhough_extended2(edgedata,dT,dS)
%         edgedata a 2-row matrix, with the x and y coordinates of the edges
%         dT is the orientation(thetas) step 
%         dS is the distance step                                
%         H votes histogram (distances - vertical x angles - horizontal)
%         Ts is the orientations vector
%         rho is the distances vector

MAXDIST=1.2;




if nargin<3
   error('wrong number of parameters')
end


row=edgedata(1,:)';
col=edgedata(2,:)';


%defining the range of the orientations of line
Ts=[(-pi/2-dT*5):dT:pi/2+4*dT]';


%cos and sin of all the angles
CsT=cos(Ts);
SnT=sin(Ts);


%solving for distances for all orientations at all nonzero pixels
%size of S is: [length(row) , length(Ts)]
S=row*CsT' + col*SnT';


%mapping:
%         Smin = min(S(:))--> 1
%         Smax = max(S(:))--> nS
%gives (y=mx+b):
%         m=(nS-1)/(Smax-Smin)
%         b=(Smax-nS*Smin)/(Smax-Smin)
%and then round it and get rounded mapped S:rmS


rmS=S/dS;
Smin=min(rmS(:));
Smax=max(rmS(:));

%m =(nS-1)/(Smax-Smin);
%b =(Smax-nS*Smin)/(Smax-Smin);
%rmS=round(m*S + b);




%Note: H is [nT,nS]
%                                 rmS is [nP,nT]  nP:number of edge points

w=1;
H=[];
for k=Smin:Smax,
    isEq=abs(rmS-k)/w;isEq(isEq>1)=0;
    H(k-Smin+1,:)=sum((isEq>0)-2*isEq.^2+isEq.^4);
end
rho=(Smin:Smax)*dS;
