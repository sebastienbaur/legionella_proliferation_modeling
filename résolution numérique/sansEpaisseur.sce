
// ----------------------------------------------------------------------------- 
// ----------------VALEURS NUMERIQUES DES CONSTANTES DU PROBLEME----------------
// -----------------------------------------------------------------------------

T=18000;
H=5*10^(-6)*200;
n=9*18000;
m=200;
dz=H/(m+1);
dt=T/(n+1);
D=10^(-10);
k_1 = 10^(-4) ;
k_2 = 2*10^(-4) ;
k_3 = 5/86400 ;
k_4 = 2*10^(-4) ;
k_5 = 2*10^(-4) ;
k_6 = 5/86400 ; 
rho_A = 1;
rho_L = 1;
L_ini=1;
A_ini=1;
N_ini=1;


// ----------------------------------------------------------------------------- 
// ----------------------MATRICES ET VECTEURS DES INCONNUES---------------------
// -----------------------------------------------------------------------------

L=zeros(m+2,n+2);
A=zeros(m+2,n+2);
N=zeros(m+2,n+2);
u=zeros(m+2,n+2);
m_L = zeros(1,n+2);
m_A = zeros(1,n+2);
m_N = zeros(1,n+2);
m_L_relatif = zeros(1,n+2);
m_A_relatif = zeros(1,n+2);
m_N_relatif = zeros(1,n+2);


// ----------------------------------------------------------------------------- 
// -----------------------------CONDITIONS INITIALES----------------------------
// -----------------------------------------------------------------------------

L(1,1) = L_ini;
L(2,1) = L_ini;
A(1,1) = A_ini;
A(2,1) = A_ini;
N(1,1) = N_ini;
N(2,1) = N_ini;
// masses de cellules contenues das le biofilm + les quantités relatives associées
m_L(1,1) = dz*10^(-4)*(L(1,1)+L(2,1));
m_A(1,1) = dz*10^(-4)*(A(1,1)+A(2,1));
m_N(1,1) = dz*10^(-4)*(N(1,1)+N(2,1));
m_L_relatif(1,1) = 1;
m_A_relatif(1,1) = 1;
m_N_relatif(1,1) = 1; 
multiplieur = ones(1,m+2);

// ----------------------------------------------------------------------------- 
// -----------------------------BOUCLES DE CALCUL-------------------------------
// -----------------------------------------------------------------------------

for j=0:n //Calcul aux différents instants
    for i=1:m //Calcul aux différentes cotes
        //Pour les légionnelles :
        L(i+1,j+2) = L(i+1,j+1) + dt*(-u(i+1,j+1)*(L(i+2,j+1)-L(i,j+1))/(2*dz) + L(i+1,j+1)*(k_1*N(i+1,j+1)/(k_2+N(i+1,j+1)) + k_3*A(i+1,j+1)/(k_4+A(i+1,j+1)) - ((u(i+2,j+1)-u(i+1,j+1))/dz)));
        L(1,j+2) = L(2,j+2);

        //Pour les amibes :
        A(i+1,j+2) = A(i+1,j+1) + dt*(-u(i+1,j+1)*(A(i+2,j+1)-A(i,j+1))/(2*dz) + A(i+1,j+1)*(k_5*N(i+1,j+1)/(k_6+N(i+1,j+1)) - k_3*L(i+1,j+1)/(k_4+A(i+1,j+1))-((u(i+2,j+1)-u(i+1,j+1))/dz)));
        A(1,j+2) = A(2,j+2);

        //Pour les nutriments :
        N(i+1,j+2) = N(i+1,j+1) + dt*(D*((N(i+2,j+1)+N(i,j+1)-2*N(i+1,j+1))/(dz*dz)) - N(i+1,j+1)*(k_1*L(i+1,j+1)/(k_2+N(i+1,j+1)) + k_5*A(i+1,j+1)/(k_6+N(i+1,j+1))));
        N(1,j+2) = N(2,j+2);

end
        // calcul de la masse totale de légionnelles, amibes et de nutriments
        m_L(1,j+2) = dz*10^(-4)*multiplieur*L(:,j+2);
        m_A(1,j+2) = dz*10^(-4)*multiplieur*A(:,j+2);
        m_N(1,j+2) = dz*10^(-4)*multiplieur*N(:,j+2);
        m_L_relatif(1,j+2) = m_L(1,j+2)/m_L(1,1);
        m_A_relatif(1,j+2) = m_A(1,j+2)/m_A(1,1);
        m_N_relatif(1,j+2) = m_N(1,j+2)/m_N(1,1);
    
    for i=0:m
        //Pour le champ de vitesses :
        u(i+2,j+2) = u(i+1,j+2) + dz*((k_1*((N(i+1,j+2)+N(i+2,j+2))/2)/(k_2+((N(i+1,j+2)+N(i+2,j+2))/2)) + k_3*((A(i+1,j+2)+A(i+2,j+2))/2)/(k_4+((A(i+1,j+2)+A(i+2,j+2))/2)))*((L(i+1,j+2)+L(i+2,j+2))/2)/rho_L + (k_5*((N(i+1,j+2)+N(i+2,j+2))/2)/(k_6+((N(i+1,j+2)+N(i+2,j+2))/2)) - k_3*((L(i+1,j+2)+L(i+2,j+2))/2)/(k_4+((A(i+1,j+2)+A(i+2,j+2))/2)))*((A(i+1,j+2)+A(i+2,j+2))/2)/rho_A);
    end

end



