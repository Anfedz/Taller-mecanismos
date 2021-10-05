%% condiciones iniciales
%constante
L2=2.45;    L3=5.42;       
BC=2;   L5=5.4;
OO5y=5.84; OO5x=0.85;
CD=5.311;  k=0;
%entrada    
vth2=152:1:300; W2=1; alpha2=0;  
%Variables
vth3=vth2*0;    vth4=vth2*0;
vth5=vth2*0;    vDO5x=vth2*0;
vW3=vth2*0;     vW4=vth2*0;
vW5=vth2*0;     vVDO5x=vth2*0;
valpha3=vth2*0; valpha4=vth2*0;
valpha5=vth2*0; vaDO5x=vth2*0;
% Valores iniciales
x_i=[130 90 80 1];

%% calculo
for th2=vth2
    k=k+1;
    const=[L2 L3 BC L5 OO5y OO5x CD th2];
    xsol=fsolve(@(x) PosEje4(x,const), x_i);

    th3=xsol(1);   th4=xsol(2);
    th5=xsol(3);   DO5x=xsol(4);

    vth3(k)=th3;  vth4(k)=th4;     
    vth5(k)=th5;  vDO5x(k)=DO5x;
    
    x_i = xsol; % Componente inercial
    
    %% Velocidad
    Av = [-L3*cosd(th3+90) BC*cosd(th4+90) -L5*cosd(th5+90) 0;
          -L3*sind(th3+90) BC*sind(th4+90) -L5*sind(th5+90) 0;
          0 -CD*cosd(th4+90) L5*cosd(th5+90) -1;
          0 -CD*sind(th4+90) L5*sind(th5+90) 0;  
          ];
    bv = [-L2*W2*cosd(th2+90);
          -L2*W2*cosd(th2+90);
          0;
          0;];
    solv=Av\bv;
    
    W3=solv(1); W4=solv(2); W5=solv(3); VDO5x=solv(4);
        
    vW3(k)=W3;     vW4(k)=W4;    vW5(k)=W5;     vVDO5x(k)=VDO5x;
    
    %% Aceleración
    
    Aa=Av;
    ba=[-L2*alpha2*cosd(th2+90)+L2*W2^2*cosd(th2)-L3*W3^2*cosd(th3)+BC*W4^2*cosd(th4)-L5*W5^2*cosd(th5);
        -L2*alpha2*sind(th2+90)+L2*W2^2*sind(th2)-L3*W3^2*sind(th3)+BC*W4^2*sind(th4)-L5*W5^2*sind(th5);
        L5*W5^2*cosd(th5)-CD*W4^2*cosd(th4);
        L5*W5^2*sind(th5)-CD*W4^2*sind(th4);];
    sola=Aa\ba;
    alpha3=sola(1); alpha4=sola(2); alpha5=sola(3); aDO5x=sola(4);
    valpha3(k)=alpha3; valpha4(k)=alpha4; valpha5(k)=alpha5; vaDO5x(k)=aDO5x;
    %% Grafica

    O = [0 0];  O5  = [OO5x -OO5y];
    A = L2*[cosd(th2) sind(th2)];
    B = A - L3*[cosd(th3) sind(th3)];
    C = B + BC*[cosd(th4) sind(th4)];
%     D = C + CD*[cosd(th4) sind(th4)];
    D = O5 + [DO5x 0];
    figure(81);
    plot(O(1),O(2),'or');hold on;
    plot(A(1),A(2),'or'); 
    plot(B(1),B(2),'or');
    plot(C(1),C(2),'or');
    plot(D(1),D(2),'og');
%     plot(D2(1),D2(2),'ob');
    plot(O5(1),O5(2),'or');

    plot([O(1)  A(1)],[O(2) A(2)],'-b');
    plot([A(1)  B(1)],[A(2) B(2)],'-b');
    plot([C(1)  D(1)],[C(2) D(2)],'-b');
    plot([C(1)  O5(1)],[C(2) O5(2)],'-b');
    hold off
    axis equal
    xlim([-5 10])
    ylim([-10 5])
    drawnow
end 

%% gráficas aceleración y velocidad
figure(20),
plot(vth2,vth3,'-d','MarkerIndices',1:10:length(vth3)), hold on, 
plot(vth2,vth4,'-o','MarkerIndices',1:10:length(vth4)),plot(vth2,vth5,'-*','MarkerIndices',1:10:length(vth5)), 
xlabel('Angulo entrada (\theta_2) [grados]')
ylabel('[grados]')
legend('\theta_3','\theta_4','\theta_5')
title ('Magnitudes angulares del mecanismo VS angulo de entrada')

figure (21),
plot(vth2,vDO5x);
xlabel('Angulo entrada (\theta_2) [grados]')
ylabel('[Unidades de longitud]')
legend('DO_5x')
title ('Magnitudes lineales del mecanismo VS angulo de entrada')

figure(22),
plot(vth2,vW3,'-d','MarkerIndices',1:10:length(vW3)), hold on, 
plot(vth2,vW4,'-o','MarkerIndices',1:10:length(vW4)), plot(vth2,vW5,'-*','MarkerIndices',1:10:length(vW5)), 
xlabel('Angulo entrada (\theta_2) [grados]')
ylabel('[rad/s]')
legend('\omega_3','\omega_4','\omega_5')
title ('Variables de velocidad angular del mecanismo VS angulo de entrada')

figure(23),
plot(vth2,vVDO5x),
xlabel('Angulo entrada (\theta_2) [grados]')
ylabel('[U.L/s]')
legend('V_DO5x')
title ('Velocidades relativas lineales del mecanismo VS angulo de entrada')


figure(24)
plot(vth2,valpha3,'-d','MarkerIndices',1:10:length(valpha3)), hold on, 
plot(vth2,valpha4,'-o','MarkerIndices',1:10:length(valpha4)),plot(vth2,valpha5,'-*','MarkerIndices',1:10:length(valpha5)) , 
xlabel('Angulo entrada (\theta_2) [grados]')
ylabel('[rad/s^2]')
legend('\alpha_3','\alpha_4','\alpha_5')
title ('Variables de aceleración angular del mecanismo VS angulo de entrada')

figure(25),
plot(vth2,vaDO5x), 
xlabel('Angulo entrada (\theta_2) [grados]')
ylabel('[U.L/s^2]')
legend('a_DO5x')
title ('Aceleraciones relativas lineales del mecanismo VS angulo de entrada')


