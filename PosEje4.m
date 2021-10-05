function y= PosEje4(x,c)
%constante
L2=c(1);    L3=c(2);       
BC=c(3);   L5=c(4);
OO5y=c(5); OO5x=c(6);
CD=c(7);
th2=c(8);%entrada       

%incongnitas
th3=x(1);   th4=x(2);
th5=x(3);   DO5x=x(4);

%ecuaciones
y(1)=L2*cosd(th2)-L3*cosd(th3)+BC*cosd(th4)-L5*cosd(th5)-OO5x;
y(2)=L2*sind(th2)-L3*sind(th3)+BC*sind(th4)-L5*sind(th5)+OO5y;
y(3)=L5*cosd(th5)-CD*cosd(th4)-DO5x;
y(4)=L5*sind(th5)-CD*sind(th4);

end