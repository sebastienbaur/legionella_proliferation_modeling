// ----------------------------------------------------------------------------- 
// STOCKAGE DES VALEURS
// ----------------------------------------------------------------------------- 

// les valeurs par instant ne sont pas toutes stockées, on en garde que certaines
i=0; // on compte chaque itération de la boucle pour savoir quelle valeur est enregistrée
frequenceStockage = 1800; // on en garde 1/fréquenceStockage en fait


NCourant = 3*ones(M,1);

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

    LCourant = (LCourant > zeros(M,1)).*LCourant;
    ACourant = (ACourant > zeros(M,1)).*ACourant;
    LMortCourant = (LMortCourant > zeros(M,1)).*LMortCourant;
    AMortCourant = (AMortCourant > zeros(M,1)).*AMortCourant;
    LCourant = epsilon*rho_L*ones(M,1).*(LCourant > epsilon*rho_L*ones(M,1)) + (LCourant <= epsilon*rho_L*ones(M,1)).*LCourant;
    ACourant = epsilon*rho_A*ones(M,1).*(ACourant > epsilon*rho_A*ones(M,1)) + (ACourant <= epsilon*rho_A*ones(M,1)).*ACourant;    
    LMortCourant = epsilon*rho_L*ones(M,1).*(LMortCourant > epsilon*rho_L*ones(M,1)) + (LMortCourant <= epsilon*rho_L*ones(M,1)).*LMortCourant;
    AMortCourant = epsilon*rho_A*ones(M,1).*(AMortCourant > epsilon*rho_A*ones(M,1)) + (AMortCourant <= epsilon*rho_A*ones(M,1)).*AMortCourant;
        
    // calcul de la quantité de légionnelles détachées
    detachement = detachement + dt*lambda*eCourant*eCourant*LCourant(indiceEpaisseur,1);


    // calcul de l'épaisseur grâce à la condition à la limite de/dt = -lambda e^2 + u(t,e(t))
    eCourant = eCourant + dt * (-lambda*eCourant*eCourant + vitesseCourante(indiceEpaisseur,1) );
    indiceEpaisseur = round(eCourant/dz);
    if eCourant > eMax then
        eMax = eCourant;
    end

    // calcul de la vitesse
    vitesseCourante = dz * matriceVitesse * dudz(NCourant,LCourant,ACourant);



    // troncage des vecteurs au delà de l'épaisseur du biofilm (on peut pour cela se servir de la matrice matriceVitesse dont les lignes sont assez utiles pour faire ça )
    LCourant = LCourant .* (matriceVitesse(indiceEpaisseur,:)');
    ACourant = ACourant .* (matriceVitesse(indiceEpaisseur,:)');
    LMortCourant = LMortCourant .* (matriceVitesse(indiceEpaisseur,:)');
    AMortCourant = AMortCourant .* (matriceVitesse(indiceEpaisseur,:)');
    vitesseCourante = vitesseCourante .* (matriceVitesse(indiceEpaisseur,:)');
    BCourant = 0.04*(1 - exp(-tCourant/tau))*ones(M,1);
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
        LMort = [LMort LMortCourant];
        AMort = [AMort AMortCourant];
        E = [E eCourant];
        indicesEpaisseur = [indicesEpaisseur indiceEpaisseur]
        V = [V vitesseCourante];
        evolutionDetachement = [evolutionDetachement detachement];
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


testConservationMasse = A + L + AMort + LMort;

xset('window',0);
plot2d(linspace(0,T,floor(T/frequenceStockage)+1),E);

xset('window',1);
plot2d(linspace(0,M*10^(-6),M),[L(:,1),L(:,round((1/frequenceStockage)*T/5)),L(:,round(2*(1/frequenceStockage)*T/5)),L(:,round(3*(1/frequenceStockage)*T/5)),L(:,round((1/frequenceStockage)*4*T/5)),L(:,round((1/frequenceStockage)*5*T/5))],leg="t=0@t=0.2T@t=0.4T@t=0.6T@t=0.8T@t=T", rect = [0,0,eMax,0.3]);

xset('window',2);
plot2d(linspace(0,M*10^(-6),M),[A(:,1),A(:,round((1/frequenceStockage)*T/5)),A(:,round(2*(1/frequenceStockage)*T/5)),A(:,round(3*(1/frequenceStockage)*T/5)),A(:,round((1/frequenceStockage)*4*T/5)),A(:,round((1/frequenceStockage)*5*T/5))],leg="t=0@t=0.2T@t=0.4T@t=0.6T@t=0.8T@t=T", rect = [0,0,eMax,0.3]);

xset('window',3);
plot2d(linspace(0,M*10^(-6),M),evolutionDetachement');

