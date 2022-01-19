clear;

% data
ye = 2021;
m = 11;
d = 24;
h = 0:0.25:23.75;

%obliczanie dnia wed³ug kalendarza juliañskiego 
jd = julday(ye,m,d,0); %[dni]


%gwiazda Regulus (gwiazdozbiór Lwa)

%rektastensja = 10h 09m 31.23s 
re_h = 10;
re_min = 9;
re_s = 31.23;

%deklinacja = 11 stopni 51' 41.8''
de_st = 11;
de_min = 51;
de_s = 41.8;

%1 h - 15 stopni / przeliczenie na dziesiêtne
rekt = re_h*15+(re_min/60)*15+(re_s/3600)*15;
deklin = de_st+(de_min/60)+(de_s/3600);


% 3 miejsca na Ziemi

%wspolrzedne Warszawy(pó³kula pó³nocna)
%(strefa czasowa UTC+01:00)
%lambda_1 = 21.0117800;
%phi_1 = 52.2297700;

%wspolrzedne Wyspy Œwiêtego Tomasza i Ksi¹¿êca (okolice równika)
%(strefa czasowa UTC±0)
%lambda_1 = 6.7273200;
%phi_1 = 0.3365400;

%wspolrzedne Buenos Aires (pó³kula po³udniowa)
%(strefa czasowa UTC-3:00)
lambda_1 = -58.3772300;
phi_1 = -34.6131500;


%Przeliczenie czasu s³onecznego UT na czas gwiazdowy S oraz obliczenie k¹ta godzinnego

%œredni czas gwiazdowy Greenwich
g = GMST(jd); %[stopnie]
% czas uniwersalny UT1
UT1 = h*1.002737909350795; %[godziny] 
%obliczenie czasu gwiazdowego(w stopniach) 
S = UT1*15 + lambda_1 + g; 
%obliczenie k¹ta godzinnego(w stopniach) 
t = S - rekt;    


for i = 1 : 96
    if t(i)<0;t(i) = t(i)+360;end
    if t(i)>360;t(i) = t(i)-360;end 
end


%rozwi¹zanie trójk¹ta paralaktycznego 

%obiczenie odleglosci zenitalnej i azymutu 
Z = acosd(sind(phi_1).*sind(deklin)+cosd(phi_1).*cosd(deklin).*cosd(t));
Az = atand(-cosd(deklin).*sind(t)./(cosd(phi_1).*sind(deklin)-sind(phi_1).*cosd(deklin).*cosd(t)));

%dla Warszaw
%for i = 1:40
%    Az(i) = Az(i)+180;
%end
%for i = 94:96
%    Az(i) = Az(i)+180;
%end

%dla Buenos Aires
for i = 1:11
    Az(i) = Az(i)+180;
end
for i = 69:96
    Az(i) = Az(i)+180;
end

%wysokosc
H = 90 - Z;

%wykres zaleznosci wysokosci od czasu;
%xlabel('h [godziny]');
%ylabel('H [stopnie]');
%title('Wykres zale¿noœci wysokoœci do czasu');
%hold on;
%plot(h,H);
    
%wykres zaleznosci azymutu do czasu
%xlabel('h [godziny]');
%ylabel('Az [stopnie]');
%title('Wykres zale¿noœci azymutu do czasu');
%hold on;
%plot(h,Az);


%transformacja wspolrzednych
r = 1;
x = r.*sind(Z).*cosd(Az);
y = r.*sind(Z).*sind(Az);
z = r.*cosd(Z);

[xh,yh,zh] = sphere(100);
surf(xh,yh,zh,'FaceAlpha',0.6);
hold on
plot3(x,y,z,'ro');
scatter3(x,y,z,44,'r','filled');
