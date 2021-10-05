function y= PosEje3(x,c)
%constante
AB=c(1);    BC=c(2);       
ADx=c(3);   EC=c(4);
alpha=c(5); CD=c(6);

ADy=c(7);%entrada       

%incongnitas
th2=x(1);   thbc=x(2);
thcd=x(3);  thec=x(4);
DE=x(5);

%ecuaciones
y(1)=AB*cosd(th2)+BC*cosd(thbc)+CD*cosd(thcd)-ADx;
y(2)=AB*sind(th2)+BC*sind(thbc)+CD*sind(thcd)-ADy;
y(3)=EC*cosd(thec)+CD*cosd(thcd);
y(4)=EC*sind(thec)+CD*sind(thcd)-DE;
y(5)=thec-thbc-alpha;

end