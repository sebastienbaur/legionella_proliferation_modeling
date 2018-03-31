// ----------------------------------------------------------------------------- 
// ----------------VALEURS NUMERIQUES DES CONSTANTES DU PROBLEME----------------
// -----------------------------------------------------------------------------

T=3000
H=10^(-5)*200;
n=7200
m=200
dz=H/(m+1)
dt=T/(n+1)
D=10^(-10)
k_1 = 10^(-4) ;
k_2 = 2*10^(-4) ;
k_3 = 5/86400 ;
k_4 = 2*10^(-4) ;
k_5 = 2*10^(-4) ;
k_6 = 5/86400 ; 
rho_A = 1;
rho_L = 1;
L_ini=1
A_ini=1
N_ini=1
lambda=1

// ----------------------------------------------------------------------------- 
// ----------------------MATRICES ET VECTEURS DES INCONNUES---------------------
// -----------------------------------------------------------------------------

L=zeros(m+2,n+2)
A=zeros(m+2,n+2)
N=zeros(m+2,n+2)
u=zeros(m+2,n+2)
Ltronq=zeros(m+2,n+2)
Atronq=zeros(m+2,n+2)
Ntronq=zeros(m+2,n+2)
utronq=zeros(m+2,n+2)
e=zeros(1,n+2)
IndiceEpaisseur=zeros(1,n+2)

// ----------------------------------------------------------------------------- 
// -----------------------------CONDITIONS INITIALES----------------------------
// -----------------------------------------------------------------------------

L(1,1) = L_ini;
L(2,1) = L_ini;
A(1,1) = A_ini;
A(2,1) = A_ini;
N(1,1) = N_ini;
N(2,1) = N_ini;
e(1,1) = 2*dz;
IndiceEpaisseur(1,1) = 2;

// ----------------------------------------------------------------------------- 
// -----------------------------BOUCLES DE CALCUL-------------------------------
// -----------------------------------------------------------------------------

for j=0:n //Calcul aux différents instants
    for i=1:m+1 //Calcul aux différentes cotes
        //Pour les légionnelles :
        L(i+1,j+2) = L(i+1,j+1) + dt*(-u(i+1,j+1)*(L(i+2,j+1)-L(i,j+1))/(2*dz) + L(i+1,j+1)*(k_1*N(i+1,j+1)/(k_2+N(i+1,j+1)) + k_3*A(i+1,j+1)/(k_4+A(i+1,j+1)) - ((u(i+2,j+1)-u(i+1,j+1))/dz)));
        L(1,j+2)=L(2,j+2);
//      L(IndiceEpaisseur(1,j+2)+2,j+2) = L(IndiceEpaisseur(1,j+2)+1,j+2);
        //Pour les amibes :
        A(i+1,j+2) = A(i+1,j+1) + dt*(-u(i+1,j+1)*(A(i+2,j+1)-A(i,j+1))/(2*dz) + A(i+1,j+1)*(k_5*N(i+1,j+1)/(k_6+N(i+1,j+1)) - k_3*L(i+1,j+1)/(k_4+A(i+1,j+1))-((u(i+2,j+1)-u(i+1,j+1))/dz)));
        A(1,j+2) = A(2,j+2);
//      A(IndiceEpaisseur(1,j+2)+2,j+2) = A(IndiceEpaisseur(1,j+2)+1,j+2);
        //Pour les nutriments :
        N(i+1,j+2) = N(i+1,j+1) + dt*(D*((N(i+2,j+1)+N(i,j+1)-2*N(i+1,j+1))/(dz*dz)) - N(i+1,j+1)*(k_1*L(i+1,j+1)/(k_2+N(i+1,j+1)) + k_5*A(i+1,j+1)/(k_6+N(i+1,j+1))));
        N(1,j+2)=N(2,j+2);
//      N(IndiceEpaisseur(1,j+2)+2,j+2) = N(IndiceEpaisseur(1,j+2)+1,j+2);
    end
    for i=0:m
        //Pour le champ de vitesses :
        u(i+2,j+2) = u(i+1,j+2) + dz*((k_1*((N(i+1,j+2)+N(i+2,j+2))/2)/(k_2+((N(i+1,j+2)+N(i+2,j+2))/2)) + k_3*((A(i+1,j+2)+A(i+2,j+2))/2)/(k_4+((A(i+1,j+2)+A(i+2,j+2))/2)))*((L(i+1,j+2)+L(i+2,j+2))/2)/rho_L + (k_5*((N(i+1,j+2)+N(i+2,j+2))/2)/(k_6+((N(i+1,j+2)+N(i+2,j+2))/2)) - k_3*((L(i+1,j+2)+L(i+2,j+2))/2)/(k_4+((A(i+1,j+2)+A(i+2,j+2))/2)))*((A(i+1,j+2)+A(i+2,j+2))/2)/rho_A);
    end
    //Détermination de l'épaisseur du biofilm :
    e(1,j+2) = e(1,j+1) + dt*(u(IndiceEpaisseur(1,j+1)+2,j+1) - lambda*e(1,j+1)*e(1,j+1));
    IndiceEpaisseur(1,j+2) = round((e(1,j+2)/dz));
    
    L(IndiceEpaisseur(1,j+2),j+2) = L(IndiceEpaisseur(1,j+2)-1,j+2);
    A(IndiceEpaisseur(1,j+2),j+2) = A(IndiceEpaisseur(1,j+2)-1,j+2);
    N(IndiceEpaisseur(1,j+2),j+2) = N(IndiceEpaisseur(1,j+2)-1,j+2);
end

// ----------------------------------------------------------------------------- 
// ------------------------MISE EN FORME DES RESULTATS--------------------------
// -----------------------------------------------------------------------------

for j=1:n+2
    for i=1:IndiceEpaisseur(1,j)+1
        Ltronq(i,j)=L(i,j);
        Atronq(i,j)=A(i,j);
        Ntronq(i,j)=N(i,j);
        utronq(i,j)=u(i,j);
    end
end

