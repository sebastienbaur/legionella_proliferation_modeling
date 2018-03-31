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
    y = (- k_1*(N.*l) ./ (k_2 * ones(size(N,1),1) + N) - k_5*(N.*a) ./ (k_6 * ones(size(N,1),1) + N))*(1+k_3*a)
endfunction
