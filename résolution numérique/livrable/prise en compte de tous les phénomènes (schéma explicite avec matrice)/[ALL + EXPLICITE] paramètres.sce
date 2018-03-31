//------------------------------------------------------------------------------
// PARAMETRES DU MODELE
//------------------------------------------------------------------------------


//Température
Temperature = 37.5;

T_opt_L = 37.5; // température optimale de croissance des légionnelles
T_i_L = 55; // paramètre d'ajustement de la courbe

T_opt_A = 43; // température optimale de croissance des amibes
T_i_A = 55; // paramètre d'ajustement de la courbe

// constantes de Monod : 
k_1 = 250*10^(-5)*exp(-((Temperature - T_opt_L)/(T_i_L - T_opt_L))^2) ; // constante de vitesse légionnelle-nutriment
k_2 = 0.2 ;
k_3 = 9; // constante de vitesse légionnelle-amibe
//k_4 = 0.2 ;
k_5 = 300*10^(-5)*exp(-((Temperature - T_opt_A)/(T_i_A - T_opt_A))^2) ; // constante de vitesse amibe-nutriment
k_6 = 0.2 ;

// masses volumiques
rho_A = 1;
rho_L = 1;


// coefficients de diffusion
D = 10^(-10);
D_b = 10^(-10); // 3*10^(-10) valeur plus adaptée

// coefficients de mortalité biocide-amibe et biocide-legionnelle
k_m_a = 40*10^(-5);
k_m_l = 40*10^(-5);

// coefficient d'arrachement
lambda = 5;
