function y= PosEje2(x,c)
%constante
ACx=c(1);       AH=c(2);
HC=c(3);        FB=c(4);
GF=c(5);        AG=c(6);
FE=c(7);        DE=c(8);
DGx=c(9);       DGy=c(10);
alpha=c(11);    th2=c(12);

%incongnitas
ACy=x(1);       th3=x(2);
HB=x(3);        th7=x(4);
thfb=x(5);      thfe=x(6);
th8=x(7);

%ecuaciones
y(1)=AH*cosd(th2)+HC*cosd(th3)-ACx;
y(2)=AH*sind(th2)+HC*sind(th3)-ACy;
y(3)=AG+GF*cosd(th7)+FB*cosd(thfb)-HB*cosd(th3)-AH*cosd(th2);
y(4)=GF*sind(th7)+FB*sind(thfb)-HB*sind(th3)-AH*sind(th2);
y(5)=GF*cosd(th7)+FE*cosd(thfe)-DE*cosd(th8)+DGx;
y(6)=GF*sind(th7)+FE*sind(thfe)-DE*sind(th8)-DGy;
y(7)=thfb-thfe-alpha;

end