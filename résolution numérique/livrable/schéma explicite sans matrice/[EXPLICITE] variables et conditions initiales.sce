// ---------------------------------------------------------------------------------------------
//INTRODUCTIONS DE TOUTES LES VARIABLES DU PROBLEME
// ---------------------------------------------------------------------------------------------


// DISCRETISATION DU TEMPS ET DE L'ESPACE
nn = 3600*201;
M = 100;

T = 3600;
H = 10^(-6)*M;

t = linspace(0,T,nn+1);
z = linspace(0,H,M+1);

dt = t(2) - t(1);
dz = z(2) - z(1);


// VALEURS INITIALES DES VARIABLES, à t = 0
tCourant = 0;

L0 = zeros(M,1);
L0(1) = 0.3 ;// valeur du vecteur légionnelle à l'instant initial
L0(2) = 0.2;
L0(3) = 0.1;
L0(4) = 0.05;
LCourant = L0;
LAvant = L0;
LMortCourant = zeros(M,1);

A0 = zeros(M,1); // valeur du vecteur amibe à l'instant initial
A0(1) = 0;
A0(2)=0.1;
A0(3)=0.2;
A0(4) = 0.25;
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

L = [L0];
LMort = [LMortCourant];
AMort = [AMortCourant];
A = [A0];
V = [vitesse0];
E = [e0];
indicesEpaisseur = [indiceEpaisseur]



// CONDITIONS AUX LIMITES DE DIRICHLET POUR LES NUTRIMENTS ET BIOCIDE
valeurAuBordN = 3;
N0(indiceEpaisseur+1,1) = valeurAuBordN; 
NCourant = N0;
Nutriment = [N0];

//valeurAuBordB = 1; // on va plutôt faire un truc progressif sinon le biocide tue toutes les bactéries directement : le biofilm n'a pas le temps de croître
valeurAuBordB = 0;
B0(indiceEpaisseur+1,1) = valeurAuBordB;
BCourant = B0;
Biocide = [B0];
tau = 3600; // on introduit une variation du biocide en 1 - exp(-t/tau) pour ne pas imposer une trop grande quantité lorsque le biofilm est peu développé
