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

// discrétisation du temps et de l'espace
N = 10000*201;
M = 100;

T = 10000;
L = 10^(-6)*M;

t = linspace(0,T,N+1);
z = linspace(0,L,M+1);

dt = t(2) - t(1);
dz = z(2) - z(1);


// VALEURS INITIALES DES VARIABLES, à t = 0
tCourant = 0;

L0 = zeros(M,1);
L0(1) = 1/2 ;// valeur du vecteur légionnelle à l'instant initial
L0(2) = 1/2;
L0(3) = 1/2;
L0(4) = 1/2;
LCourant = L0;
LAvant = L0;

A0 = zeros(M,1); // valeur du vecteur amibe à l'instant initial
A0(1) = 1/2 ;
A0(2)=1/2;
A0(3)=1/2;
A0(4) = 1/2;
ACourant = A0;
AAvant = A0;

N0 = zeros(M,1);
N0(1,1) = 10;
N0(2,1) = 10;
NCourant = N0;

e0 = 4*dz; // épaisseur initiale
indiceEpaisseur = round(e0/dz);
eCourant = e0;
eAvant = e0;

vitesse0 = zeros(M,1);
vitesseCourante = vitesse0;

multiplieur = ones(1,M);
m_L = [];
m_A = [];
m_N = [];
m_L_relatif = [];
m_A_relatif = [];
m_N_relatif = [];

// remarque : ci-dessous, la composante de la i_ème ligne, j_ème colonne des matrices l, a, S, et V donne les valeurs de l, a, N et V en z_i, à l'instant t_j
// remarque2 : ci-dessous, le vecteur ligne E contient à la j_ème colonne la valeur de l'épaisseur du biofilm à l'instant t_j
L = [L0];
A = [A0];
Nutriment = [N0];
V = [vitesse0];
E = [e0];
indicesEpaisseur = [indiceEpaisseur]
m_L = [m_L dz*10^(-4)*(L(1,1)+L(2,1))];
m_A = [m_A dz*10^(-4)*(A(1,1)+A(2,1))];
m_N = [m_N dz*10^(-4)*(Nutriment(1,1)+Nutriment(2,1))];
m_L_relatif = [1];
m_A_relatif = [1];
m_N_relatif = [1]; 




// ----------------------------------------------------------------------------- 
// CONDITIONS AUX LIMITES INITIALES
// ----------------------------------------------------------------------------- 

valeurAuBordN = 1;
conditionAuBordN = zeros(M,1);
conditionAuBordN(round(e0/dz),1) = D*dt*valeurAuBordN/(dz*dz);
// remarque : les autres conditions sont des conditions de Neumann, elles sont prises en compte par des modifications sur les matrices




// ----------------------------------------------------------------------------- 
// QUELQUES MATRICES UTILES DANS LES CALCULS
// ----------------------------------------------------------------------------- 

// vecteurs pour construire des matrices
one = ones(M,1); // vecteur colonne à M lignes qui ne contient que des 1
anotherOne = ones(M-1,1);

// matrices pour construire des dérivées
identity = diag(one,0); // matrice identité de taille M
surDiag = diag(anotherOne,1); // matrice carrée de taille M dont la surdiagonale ne contient que des 1, le reste est nul
sousDiag = diag(anotherOne,-1); // matrice carrée de taille M dont la sousdiagonale ne contient que des 1, le reste est nul

// matrices pour calculer des dérivées.
// elles sont un peu modifiées pour prendre en compte des conditions aux limites
matriceDeriveePremiere_ini = surDiag - sousDiag; // pour calculer les dérivées de L,A,N
matriceDeriveePremiere_ini(1,2) = 0; // condition au bord de Neumann (flux nul)
matriceDeriveePremiere = matriceDeriveePremiere_ini;
matriceDeriveePremiere(e0,e0+1) = 0; // condition au bord de Neumann (flux nul)
matriceDeriveeSeconde(e0,e0-1) = 0; // condition au bord de Neumann (flux nul)

matriceDeriveeSeconde_ini = -2*identity + surDiag + sousDiag;
matriceDeriveeSeconde_ini(1,2) = 2; // condition au bord de Neumann (flux nul)
matriceDeriveeSeconde = matriceDeriveeSeconde_ini;
matriceDeriveeSeconde(e0,e0+1) = 0; // condition au bord de Dirichlet (prise en compte avec l'ajout d'un vecteur conditions aux limites)

// matrice utilisée pour calculer la valeur de l'intégrale de la dérivée de la vitesse par la méthode des rectangles à droite
matriceVitesse = identity;
for i = 1 : M 
    matriceVitesse = matriceVitesse + surDiag^i;
end
matriceVitesse = matriceVitesse';




// ----------------------------------------------------------------------------- 
// STOCKAGE DES VALEURS
// ----------------------------------------------------------------------------- 

// les valeurs par instant ne sont pas toutes stockées, on en garde que certaines
i=-1; // on compte chaque itération de la boucle pour savoir quelle valeur est enregistrée
frequenceStockage = 20000; // on en garde 1/fréquenceStockage en fait




// -----------------------------------------------------------------------------------------------------------------
// BOUCLE QUI CALCULE LES DIFFERENTES VALEURS DES LEGIONNELLES, AMIBES, SUBSTRATS, VITESSE, ET EPAISSEUR DU BIOFILM
// -----------------------------------------------------------------------------------------------------------------

while (tCourant < T)

    i=i+1;



    // calcul des valeurs de A,L,N
    LCourant = LCourant + dt * (monodL(NCourant, LCourant, ACourant) - deriveeVitesse(NCourant, LCourant, ACourant)) .* LCourant - dt * vitesseCourante.*((1/(2*dz))*matriceDeriveePremiere*LCourant);  
    ACourant = ACourant + dt * (monodA(NCourant,LAvant,ACourant) - deriveeVitesse(NCourant,LAvant,ACourant)) .* ACourant - dt * vitesseCourante.*((1/(2*dz))*matriceDeriveePremiere*ACourant);
    NCourant = NCourant + dt * (D * (1/(dz*dz) ) * (matriceDeriveeSeconde * NCourant) + consoNutri(NCourant,LAvant,AAvant)) + conditionAuBordN; 

    // calcul de l'épaisseur grâce à la condition à la limite de/dt = -lambda e^2 + u(t,e(t))
    eCourant = eCourant + dt * (-lambda * eCourant^2 + vitesseCourante(indiceEpaisseur,1) );
    indiceEpaisseur = round(eCourant/dz);
    if indiceEpaisseur > round(eAvant/dz) then // lorsque l'épaisseur gagne une couche, on donne comme valeur à L et A la valeur qu'ils avaient dans la couche d'en dessous
        k = round(eAvant/dz);
        while (k <= indiceEpaisseur)
            LCourant(k,1) = LCourant(round(eAvant/dz),1);
            ACourant(k,1) = ACourant(round(eAvant/dz),1);
            k = k+1;
        end
    end

    // calcul de la vitesse
    vitesseCourante = dz * matriceVitesse * deriveeVitesse(NCourant,LCourant,ACourant);



    // troncage des vecteurs au delà de l'épaisseur du biofilm (on peut pour cela se servir de la matrice matriceVitesse dont les lignes sont assez utiles pour faire ça )
    LCourant = LCourant .* (matriceVitesse(indiceEpaisseur,:)');
    ACourant = ACourant .* (matriceVitesse(indiceEpaisseur,:)');
    NCourant = NCourant .* (matriceVitesse(indiceEpaisseur,:)');
    vitesseCourante = vitesseCourante .* (matriceVitesse(indiceEpaisseur,:)');



    // enregistrement des nouvelles valeurs des grandeurs observées
    if modulo(i,frequenceStockage) == 0 then
        A = [A ACourant];
        L = [L LCourant];
        Nutriment = [Nutriment NCourant];
        E = [E eCourant];
        indicesEpaisseur = [indicesEpaisseur indiceEpaisseur]
        V = [V vitesseCourante];
        m_A = [m_A dz*10^(-4)*multiplieur*ACourant];
        m_N = [m_N dz*10^(-4)*multiplieur*NCourant];
        m_L = [m_L dz*10^(-4)*multiplieur*LCourant];
        m_L_relatif = [m_L_relatif m_L(1,i/frequenceStockage+2)/m_L(1,1)];
        m_A_relatif = [m_A_relatif m_A(1,i/frequenceStockage+2)/m_A(1,1)];
        m_N_relatif = [m_N_relatif m_N(1,i/frequenceStockage+2)/m_N(1,1)];
    end

// mise à jour des matrices de dérivée si l'épaisseur a augmenté
    if (round(eAvant) <> round(eCourant)) then
        matriceDeriveePremiere = matriceDeriveePremiere_ini; 
        matriceDeriveePremiere(eCourant,eCourant+1) = 0;
        matriceDeriveeSeconde(eCourant,eCourant-1) = 0;
        matriceDeriveeSeconde = matriceDeriveeSeconde_ini;
        matriceDeriveeSeconde(e0,e0+1) = 0;
        conditionAuBordN = zeros(M,1);
        conditionAuBordN(eCourant,1) = D*dt*valeurAuBordN/(dz*dz);
        end

    // mise à jour de AAvant, LAvant, eAvant
    AAvant = ACourant;
    LAvant = LCourant;
    eAvant = eCourant;


    
    // permet de tomber pile sur T à la fin
    dt = min(dt, T - tCourant);

    // incrémentation du temps
    tCourant = tCourant + dt;
    
end



//plot2d(linspace(0,T,N+2)',[m_N(1,:)',m_L(1,:)',m_A(1,:)'],leg = "m_N@m_L@m_A");
//xtitle("avec amibes","t (en s)","m (en kg)") 

//plot2d(linspace(0,T,N+1)',[m_N_relatif(1,:)',m_L_relatif(1,:)',m_A_relatif(1,:)'],leg="m_N@m_L@m_A");
//xtitle("évolution temporelle des masses relatives avec amibes","t (en s)","m/m0 (sans unité)");

//plot2d(linspace(0,T,N+2)',[S(1,:)',S(5,:)',S(7,:)',S(10,:)',S(12,:)',S(13,:)'],leg="z=0@z=5@z=7@z=10@z=12@z=13", rect = [0,0,1000,0.4]);
//xtitle("évolution temporelle de S(z,t) à différents z, ","t (en s)","C_s (en g/L)");





