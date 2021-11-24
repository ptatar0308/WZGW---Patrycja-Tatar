clear;
a=6378137;
e2=0.00669437999013;
%wczytanie pliku z danymi punktami polo�enia samolotu w uk�adzie phi,lambda, h
load('C:\moje_pliki\Geoinformatyka - rok II\semestr 3\Wybrane zagadnienia geodezji wy�szej\3 raz\projekt1\matlab\lot21.txt');
%wsp samolotu
phi_samolot=lot21(:,1);
lambda_samolot=lot21(:,2);
h_samolot=lot21(:,3);
%wsp lotniska
phi_lotnisko=55.623830838;
lambda_lotnisko=12.641497434;
h_lotnisko= 5; %m n.p.m

%konwertacja danych z uk�adu xyz na uk�ad neu
%obliczanie pozycji lotniska w uk�adzie xyz
N_lotnisko = a/sqrt(1-(e2*((sind((phi_lotnisko))^2))));
x_lotnisko = (N_lotnisko+h_lotnisko)*cosd(phi_lotnisko)*cosd(lambda_lotnisko);
y_lotnisko = (N_lotnisko+h_lotnisko)*cosd(phi_lotnisko)*sind(lambda_lotnisko);
z_lotnisko = (N_lotnisko*(1-e2)+h_lotnisko)*sind(phi_lotnisko);
%obliczanie pozycji samolotu w danych punktach w uk�adzie xyz
N_samolot = a./((1-e2.*sind(phi_samolot).^2)).^(0.5);
x_samolot = (N_samolot+h_samolot).*cosd(phi_samolot).*cosd(lambda_samolot);
y_samolot = (N_samolot+h_samolot).*cosd(phi_samolot).*sind(lambda_samolot);
z_samolot = (N_samolot.*(1-e2)+h_samolot).*sind(phi_samolot);

%konwertacja danych z uk�adu xyz na uk�ad neu
%obliczanie pozycji samolotu wzgl�dem lotniska
macierz_delt = [ x_lotnisko-x_samolot y_lotnisko-y_samolot z_lotnisko-z_samolot]';
macierz_do_transpozycji = [-sind(phi_lotnisko)*cosd(lambda_lotnisko) -sind(lambda_lotnisko) cosd(phi_lotnisko)*cosd(lambda_lotnisko);
                           -sind(phi_lotnisko)*sind(lambda_lotnisko) cosd(lambda_lotnisko) cosd(phi_lotnisko)*sind(lambda_lotnisko);
                           cosd(phi_lotnisko) 0 sind(phi_lotnisko)]';
macierz_neu = macierz_do_transpozycji * macierz_delt;
%przydzielenie wierszom wsp�rz�dnych neu
n=macierz_neu(1,:);
e=macierz_neu(2,:);
u=macierz_neu(3,:);
%obserwacje w uk�adzie �przestrzennym biegunowym"
azymuty=atand(e./n);
skosna_odleglosc=(n.^2+e.^2+u.^2).^(0.5);
zenitana_odleglosc=acosd(u./skosna_odleglosc);

%znalezienie punktu, gdzie samolot zniknie za lini� horyzontu(brak takiego punktu)

%tworzenie map i wykres�w

% mapa trasu lotu
geoscatter(phi_samolot,lambda_samolot,10, 'xr');
geobasemap landcover
%geoscatter(phi_samolot,lambda_samolot,50,'.r');

% wykres neu
plot3(n,e,u);
grid on
box on
title('Lot w uk�adzie wsp�rz�dnych (n,e,u)');
xlabel('n');
ylabel('e');
zlabel('u');

% wykres A(z)
%plot(zenitana_odleglosc,azymuty);
%grid on
%box on
%xlabel('odleg�o�� zenitalna');
%ylabel('azymut');
%title('zale�no�� A(z)');

% wykres lotu w uk�adzie kartezja�skim
%plot3(x_samolot,y_samolot,z_samolot);
%grid on
%box on
%hold on
%plot3(x_lotnisko,y_lotnisko,z_lotnisko,'gx');
%xlabel('x');
%ylabel('y');
%zlabel('z');
%title('Lot w uk�adzie kartezja�skim (x,y,z)');