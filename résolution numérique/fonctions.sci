// Terme de croissance de Monod qui intervient dans l'équation bilan des légionnelles
function y = monod1 (N,l,a) // N pour nutriment, l légionnelle, et a amibes
    // terme de croissance qui apparaît dans l'équation des légionnelles
    y = ((k_1* ones(size(N,1),1) ./ (k_2*ones(size(N,1),1) + N)) + (k_3*ones(size(N,1),1) ./ (k_4*ones(size(N,1),1) + rho_A * a))) .* l;
endfunction



function y = monod2 (N,l,a)
    // terme de croissance qui apparaît dans l'équation de croissance des amibes
    y = (k_5 * ones(size(N,1),1) ./ (k_6 * ones(size(N,1),1) + N)) .* a - (k_3 * ones(size(N,1),1) ./ (k_4 * ones(size(N,1),1) + rho_A * a)) .* l ;
endfunction



function y = monod3(N,l,a)
    // terme de croissance qui apparaît dans le calcul de l'intégrale (vitesse)
    y = (k_1*ones(size(N,1),1) ./ (k_2*ones(size(N,1),1) + N)) .* l + (k_5*ones(size(N,1),1) ./ (k_6*ones(size(N,1),1) + N)) .*a ;
endfunction



function y = monod4(N,l,a)
    // terme de croissance qui apparaît dans l'équation de diffusion des nutriments
    y = (  (k_1 * ones(size(N,1),1) ./ (k_2*ones(size(N,1),1) + N)) .* (l .* (rho_L*ones(size(N,1),1)))  )  +  (  (k_5*ones(size(N,1),1) ./ (k_6*ones(size(N,1),1) + N)) .* (a .* (rho_A*ones(size(N,1),1)))  );
endfunction



function y = trouveDernierNonNul(xx)
    i = 1;
    while(i <> size(xx,1) & xx($-i+1) == 0  )
        i = i+1;
    end
    y = i;
endfunction
