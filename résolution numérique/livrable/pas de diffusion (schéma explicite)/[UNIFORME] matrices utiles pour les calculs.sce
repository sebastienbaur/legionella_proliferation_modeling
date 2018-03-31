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

// matrice utilisée pour calculer la valeur de l'intégrale de la dérivée de la vitesse par la méthode des rectangles à droite
matriceVitesse = identity;
for i = 1 : M 
    matriceVitesse = matriceVitesse + surDiag^i;
end
matriceVitesse = matriceVitesse';
