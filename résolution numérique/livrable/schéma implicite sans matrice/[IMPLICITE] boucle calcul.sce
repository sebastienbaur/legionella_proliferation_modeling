// ----------------------------------------------------------------------------- 
// STOCKAGE DES VALEURS
// ----------------------------------------------------------------------------- 

// les valeurs par instant ne sont pas toutes stockées, on en garde que certaines
i=-1; // on compte chaque itération de la boucle pour savoir quelle valeur est enregistrée
frequenceStockage = 100; // on en garde 1/fréquenceStockage en fait




// -----------------------------------------------------------------------------------------------------------------
// BOUCLE QUI CALCULE LES DIFFERENTES VALEURS DES LEGIONNELLES, AMIBES, SUBSTRATS, VITESSE, ET EPAISSEUR DU BIOFILM
// -----------------------------------------------------------------------------------------------------------------

while (tCourant < T)

    i=i+1;

    // calcul des valeurs de A,L,N
    LCourant = LCourant + dt * (-k_m_l*BCourant + mu0L(NCourant, ACourant) - dudz(NCourant, LCourant, ACourant)) .* LCourant - dt * vitesseCourante.*((1/(2*dz))*matriceDeriveePremiere*LCourant);  
    ACourant = ACourant + dt * (-k_m_a*BCourant + mu0A(NCourant) - dudz(NCourant,LAvant,ACourant)) .* ACourant - dt * vitesseCourante.*((1/(2*dz))*matriceDeriveePremiere*ACourant);
    LMortCourant = LMortCourant + dt * (k_m_l * BCourant.*LAvant - dudz(NCourant,LAvant,AAvant) .* LMortCourant - vitesseCourante.*((1/(2*dz))*matriceDeriveePremiere*LMortCourant));
    AMortCourant = AMortCourant + dt * (k_m_a * BCourant.*AAvant - dudz(NCourant,LAvant,AAvant) .* AMortCourant - vitesseCourante.*((1/(2*dz))*matriceDeriveePremiere*AMortCourant));

    // SCHEMA IMPLICITE EULER
    
//        K = (1+2*dt*D_b/(dz*dz))*identity + diag(dt*(k_m_a*ACourant + k_m_l*LCourant + dudz(NCourant,LCourant,ACourant)),0) + diag((dt/(2*dz))*vitesseCourante(1:$-1,1),1) + diag(-(dt/(2*dz))*vitesseCourante(1:$-1,1),-1) - (dt*D_b/(dz*dz))*(sousDiag+surDiag);
    K = (1+2*dt*D_b/(dz*dz))*identity + diag(dt*(k_m_a*ACourant + k_m_l*LCourant),0) - (dt*D_b/(dz*dz))*(sousDiag+surDiag); // sans convection
    K(1,2) = -2*D_b*dt/(dz*dz);
    K = inv(K);
    K(1,:) = zeros(1,M);
    BCourant = K*BCourant;
    BCourant(1,1) = (BAvant(1,1) + 2*BCourant(2,1)*dt*D_b/(dz*dz))/(1 + dt*(k_m_a*ACourant(1,1) + k_m_l*LCourant(1,1)) + 2*dt*D/(dz*dz));
//    BCourant = BCourant + dt * (D_b * (1/(dz*dz)) * (matriceDeriveeSeconde * BCourant) - (k_m_l * LAvant + k_m_a * AAvant) .* BCourant - BCourant.*dudz(NCourant,LAvant,AAvant) - vitesseCourante .* (matriceDeriveePremiere*BCourant)); // avec convection

    Z = k_1*dt*(ones(M,1)+k_3*ACourant).*LCourant + k_5*dt*ACourant;
//    K = (1+2*dt*D/(dz*dz))*identity - (dt*D/(dz*dz))*(sousDiag+surDiag) + diag(dudz(NCourant,LCourant,ACourant)*dt,0) + diag((dt/(2*dz))*vitesseCourante(1:$-1,1),1) + diag(-(dt/(2*dz))*vitesseCourante(1:$-1,1),-1); // avec convection
    K = (1+2*dt*D/(dz*dz))*identity - (dt*D/(dz*dz))*(sousDiag+surDiag); //sans convection
    K(1,2) = -2*D*dt/(dz*dz);
    K = inv(K);
    K(1,:) = zeros(1,M);
    NCourant = K *(NCourant-Z);
    NCourant(1,1) = (NAvant(1,1) + 2*NCourant(2,1)*D*dt/(dz*dz) - Z(1,1))/(1+2*dt*D/(dz*dz));

//    NCourant = NCourant + dt * (D * (1/(dz*dz)) * (matriceDeriveeSeconde * NCourant) + r_N(NCourant,LAvant,AAvant) - NCourant.*dudz(NCourant,LAvant,AAvant) - vitesseCourante .* (matriceDeriveePremiere*NCourant));// + conditionAuBordN; 


    // calcul de la quantité de légionnelles détachées
    detachement = detachement + dt*lambda*eCourant*eCourant*LCourant(indiceEpaisseur,1);

    // calcul de l'épaisseur grâce à la condition à la limite de/dt = -lambda e^2 + u(t,e(t))
    eCourant = eCourant + dt * (-lambda*eCourant*eCourant + vitesseCourante(indiceEpaisseur,1) );
    indiceEpaisseur = round(eCourant/dz);

    // calcul de la vitesse
    vitesseCourante = dz * matriceVitesse * dudz(NCourant,LCourant,ACourant);



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
    BCourant(indiceEpaisseur+1,1) = 2.2*0.04*(1 - exp(-tCourant/tau));
    if indiceEpaisseur > round(eAvant/dz) then // lorsque l'épaisseur gagne une couche, on donne comme valeur à L et A la valeur qu'ils avaient dans la couche d'en dessous
        k = round(eAvant/dz);
        while (k <= indiceEpaisseur)
           LCourant(k,1) = LCourant(round(eAvant/dz),1);
           ACourant(k,1) = ACourant(round(eAvant/dz),1);
           LMortCourant(k,1) = LMortCourant(round(eAvant/dz),1);
           AMortCourant(k,1) = AMortCourant(round(eAvant/dz),1);
           k = k+1;
        end
    end


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
    BAvant = BCourant;
    NAvant = NCourant;
    
    // permet de tomber pile sur T à la fin
    dt = min(dt, T - tCourant);

    // incrémentation du temps
    tCourant = tCourant + dt;
    
end


testConservationMasse = A + L + AMort + LMort;
