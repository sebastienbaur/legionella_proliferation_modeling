// ----------------------------------------------------------------------------- 
// VALEURS NUMERIQUES DES CONSTANTES DU PROBLEME
// ----------------------------------------------------------------------------- 

// constantes de Monod : 
k_1 = 100*10^(-5) ; // constante de vitesse légionnelle-nutriment
k_2 = 2*10^(-4) ;
k_3 = 100*10^(-5) ; // constante de vitesse légionnelle-amibe
k_4 = 2*10^(-4) ;
k_5 = 100*10^(-5) ; // constante de vitesse amibe-nutriment
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
N = 2000*201;
M = 100;

T = 2000;
L = 10^(-6)*M;



t = linspace(0,T,N+1);
z = linspace(0,L,M+1);

dt = t(2) - t(1);
dz = z(2) - z(1);


// VALEURS INITIALES DES VARIABLES, à t = 0
tCourant = 0;
L0 = zeros(M,1);
L0(1) = 1/2 ;// valeur du vecteur légionnelle à l'instant initial
L0(2)=1/2;
L0(3)=1/2;
L0(4)=1/2;
LCourant = l0;
LAvant = l0;
A0 = zeros(M,1); // valeur du vecteur amibe à l'instant initial
A0(1) = 1/2 ;
A0(2)=1/2;
A0(3)=1/2;
A0(4) = 1/2;
ACourant = A0;
AAvant = A0;
S0 = 100*ones(M,1);
//S0 = zeros(M,1); // valeur du vecteur substrat à l'instant initial
//S0(1,1) = 100;
//S0(2,1) = 100;
SCourant = S0;
//e0 = dz; // épaisseur initiale
//eCourant = e0;
vitesse0 = zeros(M,1);
// vitesse0(1)=10^(-8);
vitesseCourante = vitesse0;

multiplieur = ones(1,M);
m_L = [];
m_A = [];
m_S = [];
m_L_relatif = [];
m_A_relatif = [];
m_S_relatif = [];


// remarque : ci-dessous, la composante de la i_ème ligne, j_ème colonne des matrices l, a, S, et V donne les valeurs de l, a, S et V en z_i, à l'instant t_j
// remarque2 : ci-dessous, le vecteur ligne E contient à la j_ème colonne la valeur de l'épaisseur du biofilm à l'instant t_j
L = [L0];
A = [A0];
S = [S0];
V = [vitesse0];
//E = [e0];
m_L = [m_L dz*10^(-4)*(L(1,1)+L(2,1))];
m_A = [m_A dz*10^(-4)*(A(1,1)+A(2,1))];
m_S = [m_S dz*10^(-4)*(S(1,1)+S(2,1))];
m_L_relatif = [1];
m_A_relatif = [1];
m_S_relatif = [1]; 



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
matriceDeriveeSeconde(M,M) = 0;
matriceDeriveeSeconde(M,M-1) = 0;

matriceVitesse = identity;
for i = 1 : M 
    matriceVitesse = matriceVitesse + surDiag^i;
end
matriceVitesse = matriceVitesse';


i=-1;

// ----------------------------------------------------------------------------- 
// BOUCLE QUI CALCULE LES DIFFERENTES VALEURS DES LEGIONNELLES, AMIBES, SUBSTRATS, VITESSE, ET EPAISSEUR DU BIOFILM
// ----------------------------------------------------------------------------- 


while (tCourant < T)

    i=i+1;
    // calcul des valeurs des légionnelles, amibes, et substrats à l'instant j pour les différentes abscisses z_i
    LCourant = LCourant + dt * (monodL(SCourant, LCourant, ACourant).*LCourant - (1/dz) * matriceDeriveePremiere * (vitesseCourante .* LCourant));
//  lCourant = lCourant + dt * (monodL(SCourant, lCourant, aCourant) - deriveeVitesse(SCourant, lCourant, aCourant)) .* lCourant - dt * vitesseCourante.*((1/dz)*matriceDeriveePremiere*lCourant);  
    LCourant(1,1) = LCourant(2,1);
    
    ACourant = ACourant + dt * (monodA(SCourant, LAvant, ACourant).*ACourant - (1/dz) * matriceDeriveePremiere * (vitesseCourante .* ACourant));
//  aCourant = aCourant + dt * (monodA(SCourant,l(:,$),aCourant) - deriveeVitesse(SCourant,l(:,$),aCourant)) .* aCourant - dt*vitesseCourante.*((1/dz)*matriceDeriveePremiere*aCourant);  
    ACourant(1,1) = ACourant(2,1);
    
    SCourant = SCourant + dt * (D * (1/(dz*dz) ) * (matriceDeriveeSeconde * SCourant) + consoNutri(SCourant,LAvant,AAvant)); //.* SCourant);
    SCourant(1,1) = SCourant(2,1);
    SCourant(M,1) = 100;
    
    // calcul de l'épaisseur grâce à la condition à la limite de/dt = -lambda e^2 + u(t,e(t))
//  indiceEpaisseur = round(eCourant/dz);
//  eCourant = eCourant + dt * (-lambda * eCourant^2 + vitesseCourante(indiceEpaisseur,1) );

    // calcul de la vitesse
    vitesseCourante = dz * matriceVitesse * deriveeVitesse(SCourant,LCourant,ACourant);

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
    if modulo(i,1000) == 0 then
        A = [A aCourant];
        L = [L lCourant];
        S = [S SCourant];
//  E = [E eCourant];
        V = [V vitesseCourante];
        
        m_A = [m_A dz*10^(-4)*multiplieur*ACourant];
        m_S = [m_S dz*10^(-4)*multiplieur*SCourant];
        m_L = [m_L dz*10^(-4)*multiplieur*LCourant];
        m_L_relatif = [m_L_relatif m_L(1,i/1000+2)/m_L(1,1)];
        m_A_relatif = [m_A_relatif m_A(1,i/1000+2)/m_A(1,1)];
        m_S_relatif = [m_S_relatif m_S(1,i/1000+2)/m_S(1,1)];
    end

    AAvant = ACourant;
    LAvant = LCourant;


    // permet de tomber pile sur T à la fin
    dt = min(dt, T - tCourant);
    
end



//plot2d(linspace(0,T,N+2)',[m_S(1,:)',m_L(1,:)',m_A(1,:)'],leg = "m_S@m_L@m_A");
//xtitle("avec amibes","t (en s)","m (en kg)") 

//plot2d(linspace(0,T,N+1)',[m_S_relatif(1,:)',m_L_relatif(1,:)',m_A_relatif(1,:)'],leg="m_S@m_L@m_A");
//xtitle("évolution temporelle des masses relatives avec amibes","t (en s)","m/m0 (sans unité)");

//plot2d(linspace(0,T,N+2)',[S(1,:)',S(5,:)',S(7,:)',S(10,:)',S(12,:)',S(13,:)'],leg="z=0@z=5@z=7@z=10@z=12@z=13", rect = [0,0,1000,0.4]);
//xtitle("évolution temporelle de S(z,t) à différents z, ","t (en s)","C_s (en g/L)");





