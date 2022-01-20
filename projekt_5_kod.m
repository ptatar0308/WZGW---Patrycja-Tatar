nr = 13;

fi_A = 50+15/60+nr*15/60;
fi_A = deg2rad(fi_A);
lambda_A = 20+45/60;
lambda_A = deg2rad(lambda_A);
 
fi_B = 50+nr*15/60;
fi_B = deg2rad(fi_B);
lambda_B = 20+45/60;
lambda_B = deg2rad(lambda_B);
 
fi_C = 50+15/60+nr*15/60;
fi_C = deg2rad(fi_C);
lambda_C = 21+15/60;
lambda_C = deg2rad(lambda_C);
 
fi_D = 50+nr*15/60;
fi_D = deg2rad(fi_D);
lambda_D = 21+15/60;
lambda_D = deg2rad(lambda_D);

fi_ss = (fi_B+fi_C)/2;
lambda_ss = (lambda_B+lambda_C)/2;
fi_sr = 0.931574079859961;
lambda_sr = 0.366506367952408;

%fi_sr = 0.931574079859961;
%fi_ss =0.931569488251973;
%lambda_sr =0.366506367952408;
%lambda_ss = 0.366519142918809;

m_fi(1) = fi_A;
m_fi(2) = fi_B;
m_fi(3) = fi_C;
m_fi(4) = fi_D;
m_fi(5) = fi_ss;
m_fi(6) = fi_sr;

m_lambda(1) = lambda_A;
m_lambda(2) = lambda_B;
m_lambda(3) = lambda_C;
m_lambda(4) = lambda_D;
m_lambda(5) = lambda_ss;
m_lambda(6) = lambda_sr;

% h - wysokoœæ wszêdzie taka sama
h = 100;

a = 6378137; 
e2 = 0.00669437999013;
N = a./(sqrt(1-e2.*sin(m_fi).^2));

%% fi, lambda, h GRS80 -> x,y,z w GRS80
xGK = (N+h).*cos(m_fi).*cos(m_lambda);
yGK = (N+h).*cos(m_fi).*sin(m_lambda);
zGK = (N.*(1-e2)+h).*sin(m_fi);

%% x,y,z GRS80 -> x,y,z Krasowskiego
c11 = 0.84076440*10^(-6);
c12 = 4.08960694*10^(-6);
c13 = 0.25613907*10^(-6);
c21 = -4.08960650*10^(-6);
c22 = 0.84076292*10^(-6);
c23 = -1.73888787*10^(-6);
c31 = -0.25614618*10^(-6);
c32 = 1.73888682*10^(-6);
c33 = 0.84077125*10^(-6);
Tx = -33.4297;
Ty = 146.5746;
Tz = 76.2865;

xK = xGK+c11.*xGK+c12.*yGK+c13.*zGK+Tx;
yK = yGK+c21.*xGK+c22.*yGK+c23.*zGK+Ty;
zK = zGK+c31.*xGK+c32.*yGK+c33.*zGK+Tz;

%% algorytm Hirvonena; x,y,z Krasowskiego -> fi,lambda,h Krasowskiego

% 1. r - promien równole¿nika
r = (xK.^2+yK.^2).^(0.5);
% 2. fiK - liczymy pierwsze przybli¿enie fiK
fiK = atan((zK./r).*(1/(1-e2)));
%fiK = deg2rad(fiK);
while 1
% 3. Obliczamy N i h dla aktualnej wartoœci fiK
    NK = a./(sqrt(1-e2.*(sin(fiK)).^2));
    hK = r./(cos(fiK))-NK;
% 4. Liczymy kolejne przybli¿enia fiK+1
    fiK1 = atan((zK./r).*((1-e2.*(NK./(NK+hK))).^(-1)));
% 5. Sprawdzamy, czy spe³niony jest warunek
    epsilon = 0.00005/3600;
    epsilon = deg2rad(epsilon);
    if(fiK1-fiK<epsilon)
        break;
    end
    fiK=fiK1;
end
% 6. ostateczne wartoœci lambda_o,fi_o, h_o
fi_o = fiK1;
lambda_o = atan(yK./xK);
h_o = r./(cos(fi_o))-N_o;
N_o = a./(((1-e2.*((sin(fi_o)).^2))).^(0.5));
% 7.kontrola poprawnoœci obliczeñ
x_o = (N_o+h_o).*(cos(fi_o)).*(cos(lambda_o));
y_o = (N_o+h_o).*(cos(fi_o)).*(sin(lambda_o));
z_o = (N_o.*(1-e2)+h_o).*(sin(fi_o));


