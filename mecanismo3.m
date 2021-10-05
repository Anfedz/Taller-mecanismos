%% Condiciones iniciales

%constante
AB=10;   CD=19;
BC=8; ADx=10;
EC=15.8; alpha=19;
k=0;
%entrada
vADy=21:0.1:34.5;
VADy=1;
aADy=0;
%Variables
vth2=0*vADy;  vthbc=vth2;     
vthcd=vth2; vthec=vth2;
vDE=vth2;

vW2=vth2;   vW3=vth2;
vW4=vth2;   vVDE=vth2;

valpha2=vth2;    valpha3=vth2;
valpha4=vth2;    vaDE=vth2;

% Valores iniciales
x_i=[60 90 74 108 33];

%% Calculo
for ADy=vADy
    k=k+1;
    const=[AB BC ADx EC alpha CD ADy];
    xsol=fsolve(@(x) PosEje3(x,const), x_i);

    th2=xsol(1);   thbc=xsol(2);
    thcd=xsol(3);  thec=xsol(4);
    DE=xsol(5);

    vth2(k)=th2;  vthbc(k)=thbc;     
    vthcd(k)=thcd; vthec(k)=thec;
    vDE(k)=DE;
    x_init = xsol; % Componente inercial
    
    %% velocidad
    Av=[AB*cosd(th2+90) BC*cosd(thbc+90) CD*cosd(thcd+90) 0;
        AB*sind(th2+90) BC*sind(thbc+90) CD*sind(thcd+90) 0;
        0 EC*cosd(thec+90) CD*cosd(thcd+90) 0;
        0 EC*sind(thec+90) CD*sind(thcd+90) -1;
        ];
    
    bv=[0;VADy;0;0;];
    
    solv=Av\bv;
    
    W2=solv(1); W3=solv(2); W4=solv(3); VDE=solv(4);
    vW2(k)=W2;   vW3(k)=W3;   vW4(k)=W4;   vVDE(k)=VDE;

    %% aceleracion
    
    Aa=Av;
    
    ba = [AB*W2^2*cosd(th2)+BC*W3^2*cosd(thbc)+CD*W4^2*cosd(thcd);
          AB*W2^2*sind(th2)+BC*W3^2*sind(thbc)+CD*W4^2*sind(thcd)+aADy;
          EC*W3^2*cosd(thec)+CD*W4^2*cosd(thcd);
          EC*W3^2*sind(thec)+CD*W4^2*sind(thcd);
        ];
    sola=Aa\ba;
    alpha2=sola(1); alpha3=sola(2); alpha4=sola(3); aDE=sola(4);
    valpha2(k)=alpha2;    valpha3(k)=alpha3;  valpha4(k)=alpha4;    vaDE(k)=aDE;
    %% Gráfica

    D = [0 0];  A  = [-ADx -ADy];
    B = A + AB*[cosd(th2) sind(th2)];
    C = B + BC*[cosd(thbc) sind(thbc)];
    E = [0 -DE];
    figure(81);
    plot(D(1),D(2),'or'); hold on;
    plot(A(1),A(2),'or');
    plot(B(1),B(2),'or');
    plot(D(1),D(2),'or');
    plot(C(1),C(2),'or');
    plot(E(1),E(2),'or');

    plot([A(1)  B(1)],[A(2) B(2)],'-b');
    plot([C(1)  B(1)],[C(2) B(2)],'-b');
    plot([E(1)  B(1)],[E(2) B(2)],'-b');
    plot([C(1)  D(1)],[C(2) D(2)],'-b');
    plot([C(1)  E(1)],[C(2) E(2)],'-b');
    hold off
    axis equal
    xlim([-40 5])
    ylim([-40 5])
    drawnow
end 
%% gráficas aceleración y velocidad
figure(20),
plot(vADy,vth2,'-d','MarkerIndices',1:10:length(vth2)), hold on, 
plot(vADy,vthbc,'-o','MarkerIndices',1:10:length(vthbc)),plot(vADy,vthcd,'-*','MarkerIndices',1:10:length(vthcd)), 
plot(vADy,vthec,'-+','MarkerIndices',1:10:length(vth7)), 
xlabel('Longitud de entrada [UL]')
ylabel('[grados]')
legend('\theta23','\theta_(BC)','\theta_(CD)','\theta_EC')
title ('Magnitudes angulares del mecanismo VS magnitud de entrada')

figure (21),
plot(vADy,vDE); 
xlabel('Longitud de entrada [UL]')
ylabel('[Unidades de longitud]')
legend('vDE')
title ('Magnitudes lineales del mecanismo VS magnitud de entrada')

figure(22),
plot(vADy,vW2,'-d','MarkerIndices',1:10:length(vW2)), hold on, 
plot(vADy,vW3,'-o','MarkerIndices',1:10:length(vW3)), plot(vADy,vW4,'-*','MarkerIndices',1:10:length(vW4)), 
xlabel('Longitud de entrada [UL]')
ylabel('[rad/s]')
legend('\omega_2','\omega_3','\omega_4')
title ('Variables de velocidad angular del mecanismo VS magnitud de entrada')

figure(23),
plot(vADy,vVDE)
xlabel('Longitud de entrada [UL]')
ylabel('[U.L/s]')
legend('V_DE')
title ('Velocidades relativas lineales del mecanismo VS magnitud de entrada')


figure(24)
plot(vADy,valpha2,'-d','MarkerIndices',1:2:length(valpha2)), hold on, 
plot(vADy,valpha3,'-o','MarkerIndices',1:2:length(valpha3)),plot(vADy,valpha4,'-*','MarkerIndices',1:2:length(valpha4)) , 
xlabel('Longitud de entrada [UL]')
ylabel('[rad/s^2]')
legend('\alpha_2','\alpha_3','\alpha_4')
title ('Variables de aceleración angular del mecanismo VS magnitud de entrada')

figure(25),
plot(vADy,vaDE)
xlabel('Longitud de entrada [UL]')
ylabel('[U.L/s^2]')
legend('a_DE')
title ('Aceleraciones relativas lineales del mecanismo VS angulo de entrada')

