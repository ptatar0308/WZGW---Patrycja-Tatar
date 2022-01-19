nr = 13;

% wspó³rzêdne punktów i zamiana na radiany
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
%fi_ss = 0.931569488251973;
%lambda_sr = 0.366506367952408;
%lambda_ss = 0.366519142918809;

mac_fi = 1:6;
mac_lambda = 1:6;

mac_fi(1) = fi_A;
mac_fi(2) = fi_B;
mac_fi(3) = fi_C;
mac_fi(4) = fi_D;
mac_fi(5) = fi_ss;
mac_fi(6) = fi_sr;

mac_lambda(1) = lambda_A;
mac_lambda(2) = lambda_B;
mac_lambda(3) = lambda_C;
mac_lambda(4) = lambda_D;
mac_lambda(5) = lambda_ss;
mac_lambda(6) = lambda_sr;

mac_fi_s = rad2deg(mac_fi);
mac_lambda_s = rad2deg(mac_lambda);

% stale
a = 6378137; 
e2 = 0.00669437999013;
b = a*sqrt(1-e2);

% drugi mimoœród
e22 = (a^2-b^2)/(b^2);

%% uk³ad Gausa-Krugera dla po³udnika osiowego = 19
% L0 - poludnik osiowy w danym uk³¹dzie; L - lambda punktu
L0GK_19 = 19;
L0GK_19 = deg2rad(L0GK_19);

lGK_19 = mac_lambda-L0GK_19;
t_19=tan(mac_fi);
ni2_19 = e22.*(cos(mac_fi)).^2;

A0_19 = 1-e2/4-3*e2^2/64-5*e2^3/256;
A2_19 = 3/8*(e2+e2^2/4+15*e2^3/128);
A4_19 = 15/256*(e2^2+3*e2^3/4);
A6_19 = 35*e2^3/3072;

sigma_19 = a.*(A0_19.*mac_fi-A2_19.*sin(2.*mac_fi)+A4_19.*sin(4.*mac_fi)-A6_19.*sin(6.*mac_fi));

N = a./(sqrt(1-e2.*sin(mac_fi).^2));
M = a.*(1-e2)./sqrt((1-e2.*sin(mac_fi).^2).^3);
R_19 = sqrt(N.*M);

xGK_19 = sigma_19+((lGK_19.^2)./2).*N.*sin(mac_fi).*cos(mac_fi).*(1+((lGK_19.^2)./12).*(cos(mac_fi).^2).*(5-t_19.^2+9.*ni2_19+4.*(ni2_19.^2))+((lGK_19.^4)./360).*(cos(mac_fi).^4).*(61-58.*(t_19.^2)+(t_19.^4)+270.*ni2_19-330.*ni2_19.*(t_19.^2)));
yGK_19 = lGK_19.*N.*cos(mac_fi).*(1+lGK_19.^2./6.*(cos(mac_fi)).^2.*(1-t_19.^2+ni2_19)+lGK_19.^4/120.*((cos(mac_fi)).^4).*(5-18.*t_19.^2+t_19.^4+14.*ni2_19-58.*ni2_19.*t_19.^2));
    
% pole powierzchni
xGKp_19 = 1:4;
yGKp_19 = 1:4;
xGKp_19(1) = xGK_19(1);
xGKp_19(2) = xGK_19(2);
xGKp_19(3) = xGK_19(4);
xGKp_19(4) = xGK_19(3);

yGKp_19(1) = yGK_19(1);
yGKp_19(2) = yGK_19(2);
yGKp_19(3) = yGK_19(4);
yGKp_19(4) = yGK_19(3);

poleGKp_19 = polyarea(xGKp_19,yGKp_19);
%plot(xGKp,yGKp);

% elementarna skala d³ugoœci
m_19 = 1+(yGK_19.^2)./(2.*R_19.^2)+(yGK_19.^4)./(24.*R_19.^4);
% znieksztalcenie dlugosci
KGK_19 = m_19-1;
% elementarna skala pól
m2_19 = 1+(yGK_19.^2)./(R_19.^2)+(yGK_19.^4)./(3.*R_19.^4);
% znieksztalcenie pól
K2_GK_19 = m2_19-1;


%% uk³ad Gausa-Krugera dla po³udnika osiowego = 21
% L0 - poludnik osiowy w danym uk³¹dzie; L - lambda punktu
L0GK_21 = 21;
L0GK_21 = deg2rad(L0GK_21);

lGK_21 = mac_lambda-L0GK_21;
t_21=tan(mac_fi);
ni2_21 = e22.*(cos(mac_fi)).^2;

A0_21 = 1-e2/4-3*e2^2/64-5*e2^3/256;
A2_21 = 3/8*(e2+e2^2/4+15*e2^3/128);
A4_21 = 15/256*(e2^2+3*e2^3/4);
A6_21 = 35*e2^3/3072;

sigma_21 = a.*(A0_21.*mac_fi-A2_21.*sin(2.*mac_fi)+A4_21.*sin(4.*mac_fi)-A6_21.*sin(6.*mac_fi));

N = a./(sqrt(1-e2.*sin(mac_fi).^2));
M = a.*(1-e2)./sqrt((1-e2.*sin(mac_fi).^2).^3);
R_21 = sqrt(N.*M);

xGK_21 = sigma_21+((lGK_21.^2)./2).*N.*sin(mac_fi).*cos(mac_fi).*(1+((lGK_21.^2)./12).*(cos(mac_fi).^2).*(5-t_21.^2+9.*ni2_21+4.*(ni2_21.^2))+((lGK_21.^4)./360).*(cos(mac_fi).^4).*(61-58.*(t_21.^2)+(t_21.^4)+270.*ni2_21-330.*ni2_21.*(t_21.^2)));
yGK_21 = lGK_21.*N.*cos(mac_fi).*(1+lGK_21.^2./6.*(cos(mac_fi)).^2.*(1-t_21.^2+ni2_21)+lGK_21.^4/120.*((cos(mac_fi)).^4).*(5-18.*t_21.^2+t_21.^4+14.*ni2_21-58.*ni2_21.*t_21.^2));
    
% pole powierzchni
xGKp_21 = 1:4;
yGKp_21 = 1:4;
xGKp_21(1) = xGK_21(1);
xGKp_21(2) = xGK_21(2);
xGKp_21(3) = xGK_21(4);
xGKp_21(4) = xGK_21(3);

yGKp_21(1) = yGK_21(1);
yGKp_21(2) = yGK_21(2);
yGKp_21(3) = yGK_21(4);
yGKp_21(4) = yGK_21(3);

poleGKp_21 = polyarea(xGKp_21,yGKp_21);
%plot(xGKp,yGKp);

% elementarna skala d³ugoœci
m_21 = 1+(yGK_21.^2)./(2.*R_21.^2)+(yGK_21.^4)./(24.*R_21.^4);
% znieksztalcenie dlugosci
KGK_21 = m_21-1;
% elementarna skala pól
m2_21 = 1+(yGK_21.^2)./(R_21.^2)+(yGK_21.^4)./(3.*R_21.^4);
% znieksztalcenie pól
K2_GK_21 = m2_21-1;

%% uklad 1992
L092 = 19;
L092 = deg2rad(L092);
l92 = mac_lambda-L092;
m092 = 0.9993;

x92 = m092.*xGK_19-5300000;
y92 = m092.*yGK_19+500000;

% pole powierzchni
x92p = 1:4;
y92p = 1:4;
x92p(1) = x92(1);
x92p(2) = x92(2);
x92p(3) = x92(4);
x92p(4) = x92(3);

y92p(1) = y92(1);
y92p(2) = y92(2);
y92p(3) = y92(4);
y92p(4) = y92(3);

pole92p = polyarea(x92p,y92p);
%plot(x92p,y92p);

% elementarna skala d³ugoœci
m92 = m_19.*m092;
% zniekszta³cenie d³ugoœci
K92 = m92-1;
% elementarna skala pól
m2_92 = m92.^2;
% zniekszta³cenie pól
K2_92 = m2_92-1;    

%% uk³ad 2000    
L020 = 21;
L020 = deg2rad(L020);
l20 = mac_lambda-L020;
m020 = 0.999923;

x20 = m020.*xGK_21;
y20 = m020.*yGK_21+7*1000000+500000;

% pole powierzchni
x20p = 1:4;
y20p = 1:4;
x20p(1) = x20(1);
x20p(2) = x20(2);
x20p(3) = x20(4);
x20p(4) = x20(3);

y20p(1) = y20(1);
y20p(2) = y20(2);
y20p(3) = y20(4);
y20p(4) = y20(3);

pole20p = polyarea(x20p,y20p);
%plot(x20p,y20p);

% elementarna skala d³ugoœci
m20 = m_21.*m020;
% zniekszta³cenie d³ugoœci
K20 = m20-1;
% elementarna skala pól
m2_20 = m20.^2;
% zniekszta³cenie pól
K2_20 = m2_20-1;



