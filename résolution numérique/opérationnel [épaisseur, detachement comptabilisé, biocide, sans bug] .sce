// ----------------------------------------------------------------------------- 
// VALEURS NUMERIQUES DES CONSTANTES DU PROBLEME
// ----------------------------------------------------------------------------- 

//Température
Temperature = 37.5;

T_opt_L = 37.5; // température optimale de croissance des légionnelles
T_i_L = 55; // paramètre d'ajustement de la courbe

T_opt_A = 43; // température optimale de croissance des amibes
T_i_A = 55; // paramètre d'ajustement de la courbe

// constantes de Monod : 
k_1 = 150*10^(-5)*exp(-((Temperature - T_opt_L)/(T_i_L - T_opt_L))^2) ; // constante de vitesse légionnelle-nutriment
k_2 = 0.2 ;
k_3 = 200*10^(-5)*exp(-((Temperature - T_opt_L)/(T_i_L - T_opt_L))^2) ; // constante de vitesse légionnelle-amibe
k_4 = 0.2 ;
k_5 = 180*10^(-5)*exp(-((Temperature - T_opt_A)/(T_i_A - T_opt_A))^2) ; // constante de vitesse amibe-nutriment
k_6 = 0.2 ;

// masses volumiques
rho_A = 1;
rho_L = 1;

// coefficients de diffusion
D = 10^(-10);
D_b = 10^(-10); // 3*10^(-10) valeur plus adaptée

// coefficients de mortalité biocide-amibe et biocide-legionnelle
k_m_a = 0.04;
k_m_l = 0.04;

// coefficient d'arrachement
lambda = 15;

// discrétisation du temps et de l'espace
nn = 1000*201;
M = 100;

T = 1000;
H = 10^(-6)*M;

t = linspace(0,T,nn+1);
z = linspace(0,H,M+1);

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
LMortCourant = zeros(M,1);

A0 = zeros(M,1); // valeur du vecteur amibe à l'instant initial
A0(1) = 1/2 ;
A0(2)=1/2;
A0(3)=1/2;
A0(4) = 1/2;
ACourant = A0;
AAvant = A0;
AMortCourant = zeros(M,1);

N0 = zeros(M,1);

B0 = zeros(M,1);

e0 = 4*dz; // épaisseur initiale
indiceEpaisseur = round(e0/dz);
eCourant = e0;
eAvant = e0;

detachement = 0; // par unité de surface de biofilm

vitesse0 = zeros(M,1);
vitesseCourante = vitesse0;

multiplieur = ones(1,M);

// remarque : ci-dessous, la composante de la i_ème ligne, j_ème colonne des matrices l, a, S, et V donne les valeurs de l, a, N et V en z_i, à l'instant t_j
// remarque2 : ci-dessous, le vecteur ligne E contient à la j_ème colonne la valeur de l'épaisseur du biofilm à l'instant t_j
L = [L0];
LMort = [LMortCourant];
AMort = [AMortCourant];
A = [A0];
V = [vitesse0];
E = [e0];
indicesEpaisseur = [indiceEpaisseur]




// ----------------------------------------------------------------------------- 
// CONDITIONS AUX LIMITES DE DIRICHLET POUR LES NUTRIMENTS
// ----------------------------------------------------------------------------- 

valeurAuBordN = 10;
N0(indiceEpaisseur+1,1) = valeurAuBordN; 
NCourant = N0;
Nutriment = [N0];

//valeurAuBordB = 1; // on va plutôt faire un truc progressif sinon le biocide tue toutes les bactéries directement : le biofilm n'a pas le temps de croître
valeurAuBordB = 0;
B0(indiceEpaisseur+1,1) = valeurAuBordB;
BCourant = B0;
Biocide = [B0];
tau = 5*3600;


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
matriceDeriveePremiere(round(e0/dz),round(e0/dz)+1) = 0; // condition au bord de Neumann (flux nul)
matriceDeriveePremiere(round(e0/dz),round(e0/dz)-1) = 0; // condition au bord de Neumann (flux nul)

matriceDeriveeSeconde_ini = -2*identity + surDiag + sousDiag;
matriceDeriveeSeconde_ini(1,2) = 2; // condition au bord de Neumann (flux nul)
matriceDeriveeSeconde = matriceDeriveeSeconde_ini;
//matriceDeriveeSeconde(round(e0/dz),round(e0/dz)+1) = 0; // condition au bord de Dirichlet (prise en compte avec l'ajout d'un vecteur conditions aux limites) // en fait on peut utiliser ça en fixant la valeur de N au niveau de l'épaisseur (côté eau)

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
frequenceStockage = 2000; // on en garde 1/fréquenceStockage en fait




// -----------------------------------------------------------------------------------------------------------------
// BOUCLE QUI CALCULE LES DIFFERENTES VALEURS DES LEGIONNELLES, AMIBES, SUBSTRATS, VITESSE, ET EPAISSEUR DU BIOFILM
// -----------------------------------------------------------------------------------------------------------------

while (tCourant < T)

    i=i+1;

    // calcul des valeurs de A,L,N
    LCourant = LCourant + dt * (-k_m_l*BCourant + monodL(NCourant, LCourant, ACourant) - deriveeVitesse(NCourant, LCourant, ACourant)) .* LCourant - dt * vitesseCourante.*((1/(2*dz))*matriceDeriveePremiere*LCourant);  
    ACourant = ACourant + dt * (-k_m_a*BCourant + monodA(NCourant,LAvant,ACourant) - deriveeVitesse(NCourant,LAvant,ACourant)) .* ACourant - dt * vitesseCourante.*((1/(2*dz))*matriceDeriveePremiere*ACourant);
    LMortCourant = LMortCourant + dt * (k_m_l * BCourant.*LAvant - deriveeVitesse(NCourant,LAvant,AAvant) .* LMortCourant - vitesseCourante.*((1/(2*dz))*matriceDeriveePremiere*LMortCourant));
    AMortCourant = AMortCourant + dt * (k_m_a * BCourant.*AAvant - deriveeVitesse(NCourant,LAvant,AAvant) .* AMortCourant - vitesseCourante.*((1/(2*dz))*matriceDeriveePremiere*AMortCourant));
    NCourant = NCourant + dt * (D * (1/(dz*dz)) * (matriceDeriveeSeconde * NCourant) + consoNutri(NCourant,LAvant,AAvant));// + conditionAuBordN; 
    BCourant = BCourant + dt * (D_b * (1/(dz*dz)) * (matriceDeriveeSeconde * BCourant) - (k_m_l * LAvant + k_m_a * AAvant) .* BCourant);

    // calcul de la quantité de légionnelles détachées
    detachement = detachement + dt*lambda*eCourant*eCourant*LCourant(indiceEpaisseur,1);

    // calcul de l'épaisseur grâce à la condition à la limite de/dt = -lambda e^2 + u(t,e(t))
    eCourant = eCourant + dt * (-lambda*eCourant*eCourant + vitesseCourante(indiceEpaisseur,1) );
    indiceEpaisseur = round(eCourant/dz);

    // calcul de la vitesse
    vitesseCourante = dz * matriceVitesse * deriveeVitesse(NCourant,LCourant,ACourant);



    // troncage des vecteurs au delà de l'épaisseur du biofilm (on peut pour cela se servir de la matrice matriceVitesse dont les lignes sont assez utiles pour faire ça )
    LCourant = LCourant .* (matriceVitesse(indiceEpaisseur,:)');
    ACourant = ACourant .* (matriceVitesse(indiceEpaisseur,:)');
    NCourant = NCourant .* (matriceVitesse(indiceEpaisseur,:)');
    BCourant = BCourant .* (matriceVitesse(indiceEpaisseur,:)');
    LMortCourant = LMortCourant .* (matriceVitesse(indiceEpaisseur,:)');
    AMortCourant = AMortCourant .* (matriceVitesse(indiceEpaisseur,:)');
    vitesseCourante = vitesseCourante .* (matriceVitesse(indiceEpaisseur,:)');
    // mise à jour de la valeur des nutriments à la limite
    NCourant(indiceEpaisseur+1,1) = valeurAuBordN;
    BCourant(indiceEpaisseur+1,1) = 1 - exp(-tCourant/tau);

    // enregistrement des nouvelles valeurs des grandeurs observées
    if modulo(i,frequenceStockage) == 0 then
        A = [A ACourant];
        L = [L LCourant];
        Nutriment = [Nutriment NCourant];
        Biocide = [Biocide BCourant];
        LMort = [LMort LMortCourant];
        AMort = [AMort AMortCourant];
        E = [E eCourant];
        indicesEpaisseur = [indicesEpaisseur indiceEpaisseur]
        V = [V vitesseCourante];
    end

// mise à jour des matrices de dérivée si l'épaisseur a augmenté
    if (round(eAvant/dz) <> indiceEpaisseur) then
        matriceDeriveePremiere = matriceDeriveePremiere_ini; 
        matriceDeriveePremiere(indiceEpaisseur,indiceEpaisseur+1) = 0;
        matriceDeriveePremiere(indiceEpaisseur,indiceEpaisseur-1) = 0;
        //matriceDeriveeSeconde = matriceDeriveeSeconde_ini;
        //matriceDeriveeSeconde(round(e0/dz),round(e0/dz)+1) = 0;
//        NCourant(round(eAvant/dz)+1,1) = 0;
//        BCourant(round(eAvant/dz)+1,1) = 0;
//        if round(eAvant/dz) < indiceEpaisseur then
//            k = round(eAvant/dz) + 1
//            while k <= indiceEpaisseur 
//            NCourant(k,1) = NCourant(round(eAvant/dz),1)
//            BCourant(k,1) = BCourant(round(eAvant/dz),1)
//            k = k + 1
//            end
//            NCourant(indiceEpaisseur+1,1) = valeurAuBordN;
//            BCourant(indiceEpaisseur+1,1) = (1/0.71)*(1 - exp(-tCourant/tau));
//        else
//            k = indiceEpaisseur + 1
//            while k <= round(eAvant/dz) 
//            NCourant(k,1) = 0
//           BCourant(k,1)= 0
//            k = k + 1
//           end
//            NCourant(indiceEpaisseur+1,1) = valeurAuBordN;
//            BCourant(indiceEpaisseur+1,1) = (1/0.71)*(1 - exp(-tCourant/tau));
//        end
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

detachementTotal = detachement * ones(size(detachement),1) // calcul total du détachement

//plot2d(linspace(0,T,N+2)',[m_N(1,:)',m_L(1,:)',m_A(1,:)'],leg = "m_N@m_L@m_A");
//xtitle("avec amibes","t (en s)","m (en kg)") 

//plot2d(linspace(0,T,N+1)',[m_N_relatif(1,:)',m_L_relatif(1,:)',m_A_relatif(1,:)'],leg="m_N@m_L@m_A");
//xtitle("évolution temporelle des masses relatives avec amibes","t (en s)","m/m0 (sans unité)");

//plot2d(linspace(0,T,N+2)',[S(1,:)',S(5,:)',S(7,:)',S(10,:)',S(12,:)',S(13,:)'],leg="z=0@z=5@z=7@z=10@z=12@z=13", rect = [0,0,1000,0.4]);
//xtitle("évolution temporelle de S(z,t) à différents z, ","t (en s)","C_s (en g/L)");





