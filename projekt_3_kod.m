clear;

nr = 13;
%nr = randi(15,1,1)

% wsp�rz�dne punkt�w i zamiana na radiany
fi_A = 50+(15/60)+(nr*15/60);
fi_A = deg2rad(fi_A);
lambda_A = 20+(45/60);
lambda_A = deg2rad(lambda_A);
 
fi_B = 50+(nr*15/60);
fi_B = deg2rad(fi_B);
lambda_B = 20+(45/60);
lambda_B = deg2rad(lambda_B);
 
fi_C = 50+(15/60)+(nr*15/60);
fi_C = deg2rad(fi_C);
lambda_C = 21+15/60;
lambda_C = deg2rad(lambda_C);
 
fi_D = 50+(nr*15/60);
fi_D = deg2rad(fi_D);
lambda_D = 21+(15/60);
lambda_D = deg2rad(lambda_D);

% obliczenie punktu �redniej szeroko�ci
fi_ss = (fi_B+fi_C)/2;
lambda_ss = (lambda_B+lambda_C)/2;
 
% stale 
a=6378137; 
e2=0.00669437999013;
 
%% algorytm Vincentego

b = a * sqrt(1-e2);                 % b - kr�tsza p�o� elipsoidy
f = 1 - b/a;                        % f - sp�aszczenie elipsoidy
Ub = atan((1-f)*tan(fi_B));         % U - szeroko�� zredukowana(dla danego punktu)
Uc = atan((1-f)*tan(fi_C));

del_lambda = lambda_C - lambda_B;
L2 = del_lambda; % L - r�nica d�ugo�ci na sferze pomocniczej
while 1
L1 = L2;
sin_sigma = sqrt((cos(Uc)*sin(L2))^2 + (cos(Ub)*sin(Uc)-sin(Ub)*cos(Uc)*cos(L2))^2);
cos_sigma = sin(Ub)*sin(Uc)+cos(Ub)*cos(Uc)*cos(L2);
sigma = atan(sin_sigma/cos_sigma); %sigma - odleg�o�� k�towa mi�dzy punktami na sferze
    
% alfa - azymut linii geodezyjnej na r�wniku
sin_alfa = cos(Ub)*cos(Uc)*sin(L2)/sin_sigma;
cos2_alfa = 1-sin_alfa^2;

% sigma_m - odleg�o�� k�towa na sferze od r�wnika do punktu �rodkowego linii geodezyjnej
cos2_sigma_m = cos_sigma - 2*sin(Ub)*sin(Uc)/cos2_alfa;

C = f/16*cos2_alfa*(4+f*(4-3*cos2_alfa));
L2 = del_lambda+(1-C)*f*sin_alfa*(sigma+C*sin_sigma*(cos2_sigma_m+C*cos_sigma*(-1+2*(cos2_sigma_m)^2)));
 dL = abs(rad2deg(L2 - L1));
  if(dL < (0.000001/3600))
     break
  end 
end

u2 = (a^2-b^2)/b^2*cos2_alfa;
A = 1+u2/16384*(4096+u2*(-768+u2*(320-175*u2)));
B = u2/1024*(256+u2*(-128+u2*(74-47*u2)));

del_sigma = B*sin_sigma*(cos2_sigma_m+0.25*B*(cos_sigma*(-1+2*(cos2_sigma_m)^2)-1/6*B*cos2_sigma_m*(-3*4*(sin_sigma)^2)*(-3+4*(cos2_sigma_m)^2)));
    
S = b*A*(sigma-del_sigma);                                                  % S - dlugo�� lini geodezyjnej
Azbc = atan(cos(Uc)*sin(L2)/(cos(Ub)*sin(Uc)-sin(Ub)*cos(Uc)*cos(L2)));     % Azbc - azymut prosty
Azcb = atan(cos(Ub)*sin(L2)/(-sin(Ub)*cos(Uc)+cos(Ub)*sin(Uc)*cos(L2)))+pi; % Azcb - azymut odwrotny

%% algorytm Kivioij
 
S =S/2;
n = round(S/1000);  % n - ilo�� odcink�w
ds = S/n;           % ds - d�ugo�� odcink�w

fi_sr = fi_B;
lambda_sr = lambda_B;
Azbc_sr = Azbc;

mac_fi = 0:(n+1);
mac_lambda = 0:(n+1);
mac_Az = 0:(n+1);

mac_fi(1) = fi_sr;
mac_lambda(1) = lambda_sr;
mac_Az(1) =Azbc_sr;

for i = 1:n
    % 1. obliczmy N i M oraz stala c
    M = a*(1-e2)/sqrt((1-e2*sin(fi_B)^2)^3);
    N = a/(sqrt(1-e2*sin(fi_B)^2));

    % 2. obliczamy dfi, dA
    dfi = cos(Azbc)*ds/M;
    dAzbc = sin(Azbc_sr)*tan(fi_sr)*ds/N;
    dlambda = sin(Azbc)*ds/(N*cos(fi_sr));

    % 3. i 4. obliczamy punkt �rodkowy
    fi_s = fi_sr+0.5*dfi;
    Azbc_s = Azbc_sr+0.5*dAzbc;
    lambda_s = lambda_sr+0.5*dlambda;
    M_s = a*(1-e2)/sqrt((1-e2*sin(fi_s)^2)^3);
    N_s = a/(sqrt(1-e2*sin(fi_s)^2));

    % 5. przyrosty poprawione
    dfi_pop = cos(Azbc_s)*ds/M_s;
    dAzbc_pop = sin(Azbc_s)*tan(fi_s)*ds/N_s;
    dlambda_pop = sin(Azbc_s)*ds/(N_s*cos(fi_s));
    
    % 6. punkty ko�cowe (w ostatniej iteracji obliczenie punktu �rodkowego)
    fi_kon = fi_sr+dfi_pop;
    lambda_kon = lambda_sr+dlambda_pop;
    Azbc_kon = Azbc_sr+dAzbc_pop;
    
    % 7. 
    fi_sr = fi_kon;
    lambda_sr = lambda_kon;
    Azbc_sr = Azbc_kon;
    
    mac_fi(i+1) = fi_sr;
    mac_lambda(i+1) = lambda_sr;
    mac_Az(i+1) = Azbc_sr;
end

mac_fi(n+2) = fi_ss;
mac_lambda(n+2) = lambda_ss;


%% obliczenie pola powierzchni czworokata
% e - pierwiastek z mimo�rodu
P = ((b^2)*(lambda_C-lambda_B)/2)*(((sin(fi_C))/(1-(e2)*(sin(fi_C))^2)+(1/(2*(sqrt(e2))))*log((1+(sqrt(e2))*sin(fi_C))/(1-(sqrt(e2))*sin(fi_C)))) - ((sin(fi_B))/(1-(e2)*(sin(fi_B))^2)+(1/(2*(sqrt(e2))))*log((1+(sqrt(e2))*sin(fi_B))/(1-(sqrt(e2))*sin(fi_B)))));   

%% obliczenie roznicy miedzy punktami �rodkowym, a punktem �redniej szeroko�ci i wyznaczenie azymut�w w tych punktach
% obliczenie: s-dlugo�� lini geodezyjnej, Azbc-azymut prosty, Azcb-azymut odwrotny
% korzystamy z algorytmu Vincentego 
% dane wej�ciowe do wsp�rz�dne punktu �redniej szeroko�ci i wsp�rz�dne punktu �rodkowego

del_lambda_n = lambda_sr - lambda_ss;
Ub_n = atan((1-f)*tan(fi_ss));
Uc_n = atan((1-f)*tan(fi_sr));

% L - r�nica d�ugo�ci na sferze pomocniczej, sigma - odleg�o�� k�towa mi�dzy punktami na sferze
L2_n = del_lambda_n;
while 1
L1_n = L2_n;
sin_sigma_n = sqrt((cos(Uc_n)*sin(L2_n))^2 + (cos(Ub_n)*sin(Uc_n)-sin(Ub_n)*cos(Uc_n)*cos(L2_n))^2);
cos_sigma_n = sin(Ub_n)*sin(Uc_n)+cos(Ub_n)*cos(Uc_n)*cos(L2_n);
sigma_n = atan(sin_sigma_n/cos_sigma_n);
    
% alfa - azymut linii geodezyjnej na r�wniku
sin_alfa_n = cos(Ub_n)*cos(Uc_n)*sin(L2_n)/sin_sigma_n;
cos2_alfa_n = 1-sin_alfa_n^2;

% sigma_m - odleg�o�� k�towa na sferze od r�wnika do punktu �rodkowego linii geodezyjnej
cos2_sigma_m_n = cos_sigma_n - 2*sin(Ub_n)*sin(Uc_n)/cos2_alfa_n;

C_n = f/16*cos2_alfa_n*(4+f*(4-3*cos2_alfa_n));
L2_n = del_lambda_n+(1-C_n)*f*sin_alfa_n*(sigma_n+C_n*sin_sigma_n*(cos2_sigma_m_n+C_n*cos_sigma_n*(-1+2*(cos2_sigma_m_n)^2)));
 dL_n = abs(rad2deg(L2_n - L1_n));
  if(dL_n < (0.000001/3600))
     break
  end 
end
u2_n = (a^2-b^2)/b^2*cos2_alfa_n;
A_n = 1+u2_n/16384*(4096+u2_n*(-768+u2_n*(320-175*u2_n)));
B_n = u2_n/1024*(256+u2_n*(-128+u2_n*(74-47*u2_n)));

del_sigma_n = B_n*sin_sigma_n*(cos2_sigma_m_n+0.25*B_n*(cos_sigma_n*(-1+2*(cos2_sigma_m_n)^2)-1/6*B_n*cos2_sigma_m_n*(-3*4*(sin_sigma_n)^2)*(-3+4*(cos2_sigma_m_n)^2)));
    
% s - dlugo�� lini geodezyjnej, Azbc - azymut prosty, Azcb - azymut odwrotny
S_n = b*A_n*(sigma_n-del_sigma_n);
Azbc_n = atan(cos(Uc_n)*sin(L2_n)/(cos(Ub_n)*sin(Uc_n)-sin(Ub_n)*cos(Uc_n)*cos(L2_n)));
Azcb_n = atan(cos(Ub_n)*sin(L2_n)/(-sin(Ub_n)*cos(Uc_n)+cos(Ub_n)*sin(Uc_n)*cos(L2_n)))+pi;

Azbc = rad2deg(Azbc);
Azcb = rad2deg(Azcb);

%% wizualizacja

%pokazanie wszystkich wyznaczonych punkt�w
%geoscatter(rad2deg(mac_fi),rad2deg(mac_lambda),5, 'ro');

%pokazanie tylko punktu �redniej szeroko�ci i punktu �rodkowego
mac_p_fi = 0:1;
mac_p_fi(1) = fi_ss;
mac_p_fi(2) = fi_sr;
mac_p_lambda = 0:1;
mac_p_lambda(1) = lambda_ss;
mac_p_lambda(2) = lambda_sr;
%geoscatter(rad2deg(mac_p_fi),rad2deg(mac_p_lambda),10, 'ro');

%% Wnioski

%  R��ica mi�dzy punktem �redniej szeroko�ci, a punktem �rodkowym wynosi
%  tylko 56.836[m]. Punkty s� bardzo blisko siebie, wi�c ta warto�� jest
%  ma�a. W przypadku punkt�w du�o bardziej od siebie oddalonych r�nica
%  pomi�dzy punktami znacznie wzro�nie.

%  W przypadku wyznaczania najkr�tszej odleg�o��i pomi�dzy dwoma odleg�ymi
%  punktami trzeba obliczy� t� odleg�o�� na elipsoidzie. Nie mo�na wyznaczy�
%  tej odleg�o�ci za pomoc� linii prostej na mapie. Jest to niezb�dna wiedza
%  przy planowaniu d�ugich trach odbywaj�cych si� "w linii prostej", np.
%  tras statk�w lub samolot�w.
