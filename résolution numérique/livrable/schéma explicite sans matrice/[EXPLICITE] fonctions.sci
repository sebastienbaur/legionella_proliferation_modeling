function y = monodL (N,l,a) // N pour nutriment, l légionnelle, et a amibes
    // terme de croissance qui apparaît dans l'équation des légionnelles
    y = (k_1 * N ./ (k_2*ones(size(N,1),1) + N)) + (k_3*a ./ (k_4*ones(size(N,1),1) + rho_A * a));
endfunction

function y = monodA (N,l,a)
    // terme de croissance qui apparaît dans l'équation de croissance des amibes
    y = k_5 * N ./ (k_6 * ones(size(N,1),1) + N) - k_3 * l  ./ (k_4 * ones(size(N,1),1) + rho_A * a) ;
endfunction


function y = consoNutri (N,l,a)
    y = - k_1*(N.*l) ./ (k_2 * ones(size(N,1),1) + N) - k_5*(N.*a) ./ (k_6 * ones(size(N,1),1) + N) ;
endfunction


function y = deriveeVitesse(N,l,a)
    y = k_1*(1/rho_L)*(N.*l) ./ (k_2 * ones(size(N,1),1) + N) + k_3*((1/rho_L) - 1/rho_A)*(a.*l) ./ (k_4 * ones(size(N,1),1) + a) + k_5*(1/rho_A)*(N.*a) ./ (k_6 * ones(size(N,1),1) + N) ;
endfunction


function y = mu0L(N,a)
    y = (k_1*N./(k_2*ones(size(N,1),1) + N)).*(1+k_3*a);
endfunction


function y = mu0A(N)
    y = k_5*N./(k_6*ones(size(N,1),1)+N);
endfunction


function y = dudz(N,l,a)
    y = (1/(1-epsilon))*((k_1/rho_L)*(ones(size(N,1),1)+k_3*a).*N.*l ./ (k_2*ones(size(N,1),1) + N) + (k_5/rho_A)*N.*a ./ (k_6*ones(size(N,1),1) + N));
endfunction


function y = r_N(N,l,a)
    y = - (k_1*(N.*l).*(1+k_3*a) ./ (k_2 * ones(size(N,1),1) + N)) - (k_5*(N.*a) ./ (k_6 * ones(size(N,1),1) + N))
endfunction
