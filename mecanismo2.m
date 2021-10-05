%% Condiciones iniciales

%constante

ACx=1; AH=4; HC=9;
FB=6; GF=7; AG=13;
FE=5; DE=5; DGx=2;
DGy=8; alpha=55; k = 0;

%entrada
vth2=-30:1:65;

%Variables
vACy=vth2*0;    vth3=vth2*0;  vHB=vth2*0;        
vth7=vth2*0;    vthfb=vth2*0; vthfe=vth2*0;
vth8=vth2*0;    vW3=vth2*0; vW6=vth2*0;     
vW7=vth2*0; vW8=vth2*0; vVb=vth2*0;
vVc=vth2*0; valpha3 = vth2*0;   valpha6 = vth2*0; 
valpha7 = vth2*0;   valpha8 = vth2*0; valphab = vth2*0;
valphac = vth2*0; trB = [0*vth2' 0*vth2'];  trC = trB ; trT = trB;

% Valores iniciales
W2 = 10; alpha2 = 0;
x_i=[17 110 8 135 170 100 170];

%% ciclo de calculos
for th2 = vth2
    k =k+1;
    conts = [ACx AH HC FB GF AG FE DE DGx DGy alpha th2];
    
    xsol=fsolve(@(x) PosEje2(x,conts), x_i);
    
    x_init = xsol; % Componente inercial
    
    ACy=xsol(1);       th3=xsol(2);
    HB=xsol(3);        th7=xsol(4);
    thfb=xsol(5);      thfe=xsol(6);
    th8=xsol(7);
    
    vACy(k)=ACy;       vth3(k)=th3;
    vHB(k)=HB;        vth7(k)=th7;
    vthfb(k)=thfb;      vthfe(k)=thfe;
    vth8(k)=th8;
    %% Solucion de velocidad
    Av = [HC*cosd(th3+90) 0 0 0 0 0;
          HC*sind(th3+90) 0 0 0 0 -1;
          -HB*cosd(th3+90) FB*cosd(thfb+90) GF*cosd(th7+90) 0 -cosd(th3) 0;
          -HB*sind(th3+90) FB*sind(thfb+90) GF*sind(th7+90) 0 -sind(th3) 0;
          0 -FE*cosd(thfe+90) -GF*cosd(th7+90) DE*cosd(th8+90) 0 0;
          0 -FE*sind(thfe+90) -GF*sind(th7+90) DE*sind(th8+90) 0 0;];
      
    bv = [-AH*W2*cosd(th2+90);-AH*W2*sind(th2+90);AH*W2*cosd(th2+90);AH*W2*sind(th2+90);0;0;];
    
    solv = Av\bv;
     
    W3 = solv(1); W6 = solv(2); W7 = solv(3); W8 = solv(4); Vb = solv(5); Vc = solv(6);
    
    vW3(k) = W3; vW6(k) = W6; vW7(k) = W7; vW8(k) = W8; vVb(k) = Vb; vVc(k) = Vc;

    %% Solucion de aceleracion 
    
    Aa = Av;
    ba = [-AH*alpha2*cosd(th2+90)+AH*W2^2*cosd(th2)-HC*W3^2*cosd(th3);
          -AH*alpha2*sind(th2+90)+AH*W2^2*sind(th2)-HC*W3^2*sind(th3);
          AH*alpha2*cosd(th2+90)-AH*W2^2*cosd(th2)+2*Vb*W3*cosd(th3+90)-HB*W3^2*cosd(th3)+FB*W6^2*cosd(thfb)+GF*W7^2*cosd(th7);
          AH*alpha2*sind(th2+90)-AH*W2^2*sind(th2)+2*Vb*W3*sind(th3+90)-HB*W3^2*sind(th3)+FB*W6^2*sind(thfb)+GF*W7^2*sind(th7);
          -GF*W7^2*cosd(th7)-FE*W6^2*cosd(thfe)+DE*W8^2*cosd(th8);
          -GF*W7^2*sind(th7)-FE*W6^2*sind(thfe)+DE*W8^2*sind(th8);];
    sola = Av\ba;
    
    alpha3 = sola(1); alpha6 = sola(2); alpha7 = sola(3); alpha8 = sola(4);alphab = sola(5);alphac = sola(6);
    
    valpha3(k) = alpha3; 
    valpha6(k) = alpha6; 
    valpha7(k) = alpha7;
    valpha8(k) = alpha8; 
    valphab(k) = alphab; 
    valphac(k) = alphac;
    
    %% Gráfica mecanismo
    %coordenadas
    A=[0 0];   H=[AH*cosd(th2) AH*sind(th2)];
    G=[AG 0];   B= H + HB*[cosd(th3) sind(th3)];
    C= [ACx ACy];
    D= G + [-DGx DGy];
    F = G + GF*[cosd(th7) sind(th7)]; 
    E = F + FE*[cosd(thfe) sind(thfe)];
    T = (B+E+F)/3;
    trB(k,:) = B;
    trC(k,:) = C;
    trT(k,:) = T;
    %grafica
    figure(80)
    plot(A(1),A(2),'or'), hold on
    plot(H(1),H(2),'or')
    plot(B(1),B(2),'or')
    plot(C(1),C(2),'or')
    plot(E(1),E(2),'or')
    plot(D(1),D(2),'or')
    plot(F(1),F(2),'or')
    plot(G(1),G(2),'or')
    % Elementos
    plot([A(1)  H(1)],[A(2) H(2)],'-b')
    % area()
    plot([B(1) H(1)],[B(2) H(2)],'-b')
    plot([H(1) C(1)],[H(2) C(2)],'-b')
    plot([B(1) E(1)],[B(2) E(2)],'-b')
    plot([B(1) F(1)],[B(2) F(2)],'-b')
    plot([E(1) F(1)],[E(2) F(2)],'-b')
    plot([G(1) F(1)],[G(2) F(2)],'-b')
    plot([D(1) E(1)],[D(2) E(2)],'-b')
    % trayectorias
    plot(trB(1:k,1),trB(1:k,2),'-g')
    plot(trC(1:k,1),trC(1:k,2),'-m')
    plot(trT(1:k,1),trT(1:k,2),'-m')

    
    hold off
    axis equal
    xlim([-5 15])
    ylim([-5 15])
    drawnow
end
%% gráficas posición, aceleración y velocidad
figure(20),
plot(vth2,vth3,'-d','MarkerIndices',1:10:length(vth3)), hold on, 
plot(vth2,vthfe,'-o','MarkerIndices',1:10:length(vthfe)),plot(vth2,vthfb,'-*','MarkerIndices',1:10:length(vthfb)), 
plot(vth2,vth7,'-+','MarkerIndices',1:10:length(vth7)), plot(vth2,vth8,'-s','MarkerIndices',1:10:length(vth8))
xlabel('Angulo entrada (\theta_2) [grados]')
ylabel('[grados]')
legend('\theta_3','\theta_(FE)','\theta_(FB)','\theta_7','\theta_8')
title ('Magnitudes angulares del mecanismo VS angulo de entrada')

figure (21),
plot(vth2,vACy); hold on, plot(vth2,vHB);
xlabel('Angulo entrada (\theta_2) [grados]')
ylabel('[Unidades de longitud]')
legend('AC_y','HB')
title ('Magnitudes lineales del mecanismo VS angulo de entrada')

figure(22),
plot(vth2,vW3,'-d','MarkerIndices',1:10:length(vW3)), hold on, 
plot(vth2,vW6,'-o','MarkerIndices',1:10:length(vW6)), plot(vth2,vW7,'-*','MarkerIndices',1:10:length(vW7)), 
plot(vth2,vW8,'-+','MarkerIndices',1:10:length(vW8))
xlabel('Angulo entrada (\theta_2) [grados]')
ylabel('[rad/s]')
legend('\omega_3','\omega_6','\omega_7','\omega_8')
title ('Variables de velocidad angular del mecanismo VS angulo de entrada')

figure(23),
plot(vth2,vVc), hold on, plot(vth2,vVb)
xlabel('Angulo entrada (\theta_2) [grados]')
ylabel('[U.L/s]')
legend('V_C','V_B')
title ('Velocidades relativas lineales del mecanismo VS angulo de entrada')


figure(24)
plot(vth2,valpha3,'-d','MarkerIndices',1:10:length(valpha3)), hold on, 
plot(vth2,valpha6,'-o','MarkerIndices',1:10:length(valpha6)),plot(vth2,valpha7,'-*','MarkerIndices',1:10:length(valpha7)) , 
plot(vth2,valpha8,'-+','MarkerIndices',1:10:length(valpha8))
xlabel('Angulo entrada (\theta_2) [grados]')
ylabel('[rad/s^2]')
legend('\alpha_3','\alpha_6','\alpha_7','\alpha_8')
title ('Variables de aceleración angular del mecanismo VS angulo de entrada')

figure(25),
plot(vth2,valphac), hold on, plot(vth2,valphab)
xlabel('Angulo entrada (\theta_2) [grados]')
ylabel('[U.L/s^2]')
legend('a_C','a_B')
title ('Aceleraciones relativas lineales del mecanismo VS angulo de entrada')