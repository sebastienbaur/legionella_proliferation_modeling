// ----------------------------------------------------------------------------- 
// VALEURS NUMERIQUES DES CONSTANTES DU PROBLEME
// ----------------------------------------------------------------------------- 

// constantes de Monod : 
k_1 = 5*10^(-5) ; // constante de vitesse légionnelle-nutriment
k_2 = 2*10^(-4) ;
k_3 = 5*10^(-5) ; // constante de vitesse légionnelle-amibe
k_4 = 2*10^(-4) ;
k_5 = 5*10^(-5) ; // constante de vitesse amibe-nutriment
k_6 = 2*10^(-4) ;

// masses volumiques
rho_A = 1;
rho_L = 1;

// coefficient de diffusion
D = 10^(-10);

// coefficient d'arrachement
lambda = 0;

// flux à la limite supérieure de substrat
phi_0 = 1;

// discrétisation du temps et de l'espace
N = 521000;
M = 30;

T = 259200;
L = 10^(-5)*M;



t = linspace(0,T,N+1);
z = linspace(0,L,M+1);

dt = t(2) - t(1);
dz = z(2) - z(1);


// VALEURS INITIALES DES VARIABLES, à t = 0
tCourant = 0;
l0 = zeros(M,1);
l0(1) = 1/2 ;// valeur du vecteur légionnelle à l'instant initial
l0(2)=1/2;
lCourant = l0;

a0 = zeros(M,1); // valeur du vecteur amibe à l'instant initial
a0(1) = 1/2 ;
a0(2)=1/2;
aCourant = a0;
S0 = zeros(M,1); // valeur du vecteur substrat à l'instant initial
S0(1,1) = 1;
S0(2,1) = 1;
SCourant = S0;
//e0 = dz; // épaisseur initiale
//eCourant = e0;
vitesse0 = zeros(M,1);
// vitesse0(1)=10^(-6);
vitesseCourante = vitesse0;

multiplieur = ones(1,M);
m_L = zeros(1,N+1);
m_A = zeros(1,N+1);
m_N = zeros(1,N+1);
m_L_relatif = zeros(1,N+1);
m_A_relatif = zeros(1,N+1);
m_N_relatif = zeros(1,N+1);


// remarque : ci-dessous, la composante de la i_ème ligne, j_ème colonne des matrices l, a, S, et V donne les valeurs de l, a, S et V en z_i, à l'instant t_j
// remarque2 : ci-dessous, le vecteur ligne E contient à la j_ème colonne la valeur de l'épaisseur du biofilm à l'instant t_j
l = [l0];
a = [a0];
S = [S0];
V = [vitesse0];
//E = [e0];
m_L(1,1) = dz*10^(-4)*(l(1,1)+l(2,1));
m_A(1,1) = dz*10^(-4)*(a(1,1)+a(2,1));
m_S(1,1) = dz*10^(-4)*(S(1,1)+S(2,1));
m_L_relatif(1,1) = 1;
m_A_relatif(1,1) = 1;
m_S_relatif(1,1) = 1; 



// ----------------------------------------------------------------------------- 
// QUELQUES MATRICES UTILES DANS LES CALCULS
// ----------------------------------------------------------------------------- 

one = ones(M,1); // vecteur colonne à M lignes qui ne contient que des 1
anotherOne = ones(M-1,1);

identity = diag(one,0); // matrice identité de taille M
surDiag = diag(anotherOne,1); // matrice carrée de taille M dont la surdiagonale ne contient que des 1, le reste est nul
sousDiag = diag(anotherOne,-1); // matrice carrée de taille M dont la sousdiagonale ne contient que des 1, le reste est nul

matriceDeriveePremiere = identity - sousDiag; 
matriceDeriveePremiere(1,1) = 0;
matriceDeriveeSeconde = -2*identity + surDiag + sousDiag; 
matriceDeriveeSeconde(1,1) = 0;
matriceDeriveeSeconde(1,2) = 0;

matriceVitesse = identity;
for i = 1 : M 
    matriceVitesse = matriceVitesse + surDiag^i;
end
matriceVitesse = matriceVitesse';


i=0;

// ----------------------------------------------------------------------------- 
// BOUCLE QUI CALCULE LES DIFFERENTES VALEURS DES LEGIONNELLES, AMIBES, SUBSTRATS, VITESSE, ET EPAISSEUR DU BIOFILM
// ----------------------------------------------------------------------------- 


while (tCourant < T)

    i=i+1;
    // calcul des valeurs des légionnelles, amibes, et substrats à l'instant j pour les différentes abscisses z_i
    lCourant = lCourant + dt * (monodL(SCourant, lCourant, aCourant).*lCourant - (1/dz) * matriceDeriveePremiere * (vitesseCourante .* lCourant));
//  lCourant = lCourant + dt * (monodL(SCourant, lCourant, aCourant) - deriveeVitesse(SCourant, lCourant, aCourant)) .* lCourant - dt * vitesseCourante.*((1/dz)*matriceDeriveePremiere*lCourant);  
    lCourant(1,1) = lCourant(2,1);
    
    aCourant = aCourant + dt * (monodA(SCourant, l(:,$), aCourant).*aCourant - (1/dz) * matriceDeriveePremiere * (vitesseCourante .* aCourant));
//  aCourant = aCourant + dt * (monodA(SCourant,l(:,$),aCourant) - deriveeVitesse(SCourant,l(:,$),aCourant)) .* aCourant - dt*vitesseCourante.*((1/dz)*matriceDeriveePremiere*aCourant);  
    aCourant(1,1) = aCourant(2,1);
    
    SCourant = SCourant + dt * (D * (1/dz^2) * (matriceDeriveeSeconde * SCourant) + consoNutri(SCourant,l(:,$),a(:,$))); //.* SCourant);
    SCourant(1,1) = SCourant(2,1);
    
    // calcul de l'épaisseur grâce à la condition à la limite de/dt = -lambda e^2 + u(t,e(t))
//  indiceEpaisseur = round(eCourant/dz);
//  eCourant = eCourant + dt * (-lambda * eCourant^2 + vitesseCourante(indiceEpaisseur,1) );

    // calcul de la vitesse
    vitesseCourante = dz * matriceVitesse * deriveeVitesse(SCourant,lCourant,aCourant);

    // calcul de la masse totale de légionnelles, amibes et de nutriments

    m_A(1,i+1) = dz*10^(-4)*multiplieur*aCourant;
    m_S(1,i+1) = dz*10^(-4)*multiplieur*SCourant;
    m_L(1,i+1) = dz*10^(-4)*multiplieur*lCourant;
    m_L_relatif(1,i+1) = m_L(1,i+1)/m_L(1,1);
    m_A_relatif(1,i+1) = m_A(1,i+1)/m_A(1,1);
    m_S_relatif(1,i+1) = m_S(1,i+1)/m_S(1,1);
    if m_S(1,i+1) <= 0 then
        S = zeros(size(S,1),1);
    end
    

    // incrémentation du temps
    tCourant = tCourant + dt;

    // troncage des vecteurs au delà de l'épaisseur du biofilm
//  for i = indiceEpaisseur+1:M 
//      lCourant(i,1)=0;
//      aCourant(i,1)=0;
//      SCourant(i,1)=0;
//  end
    
    // conditions limites : flux nul aux interfaces
//  aCourant(1,1)=aCourant(2,1);
//  lCourant(1,1)=lCourant(2,1);
//  SCourant(1,1)=SCourant(2,1);
//  aCourant(indiceEpaisseur,1)=aCourant(indiceEpaisseur-1,1);
//  lCourant(indiceEpaisseur,1)=lCourant(indiceEpaisseur-1,1);
//  SCourant(indiceEpaisseur,1)=SCourant(indiceEpaisseur-1,1);
    
    
    // enregistrement des nouvelles valeurs des grandeurs observées
    a = [a aCourant];
    l = [l lCourant];
    S = [S SCourant];
//  E = [E eCourant];
    V = [V vitesseCourante];


    // permet de tomber pile sur T à la fin
    dt = min(dt, T - tCourant);
    
end



//plot2d(linspace(0,T,N+2)',[m_S(1,:)',m_L(1,:)',m_A(1,:)'],leg = "m_S@m_L@m_A");
//xtitle("avec amibes","t (en s)","m (en kg)") 

plot2d(linspace(0,T,N+1)',[m_S_relatif(1,:)',m_L_relatif(1,:)',m_A_relatif(1,:)'],leg="m_S@m_L@m_A");
xtitle("évolution temporelle des masses relatives avec amibes, k_3 grand","t (en s)","m/m0 (sans unité)");

//plot2d(linspace(0,T,N+2)',[S(1,:)',S(5,:)',S(7,:)',S(10,:)',S(12,:)',S(13,:)'],leg="z=0@z=5@z=7@z=10@z=12@z=13", rect = [0,0,1000,0.4]);
//xtitle("évolution temporelle de S(z,t) à différents z, ","t (en s)","C_s (en g/L)");





