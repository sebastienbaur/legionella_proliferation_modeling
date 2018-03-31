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
    BCourant(indiceEpaisseur+1,1) = 0.01*(1 - exp(-tCourant/tau));

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
