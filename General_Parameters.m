%% Proyecto: Turbomáquinas
clear
clc

%% Parametros a analizar:
r = input('Número de velocidades específicas requeridas: ');
Z_max = input('Máximo número de inyectores a analizar: ');


%% Datos de entrada constantes:

g = 9.81;            % [m/s^2]        ACELERACIÓN DE LA GRAVEDAD
k_m = 0.475;         % []             COEFICIENTE DE VELOCIDAD PERIFÉRICA
phi_b = 0.11;        % []             COEFICIENTE DE CARGA VOLUMÉTRICA
Beta_2 = 170;        % [º]            ÁNGULO DE SALIDA DE LA CUCHARA
H = 25.1;            % [m]            ALTURA HIDRÁULICA
Q_total = 185.8e-3;  % [m^3/s]        CAUDAL TOTAL
f = 60;              % [1/s]          FRECUENCIA DE OPERACIÓN DE LA RED EN EL ECUADOR

%% Velocidad específica:

% r = 10;                                   %   NÚMERO DE VELOCIDADES ESPECÍFICAS DENTRO DEL INTERVALO
% d = (0.15 - 0.04)/(r - 1);                %   Criterio Zhang
% d = (2.2285 - 0.2785)/(r - 1);            %   Criterio Csnady
d = (2.2285 - 0.04)/(r - 1);                %   Combinación de Criterios
nq = zeros(1,r);            
nq(1) = 0.04;               
for i = 1:(r - 1)
    nq(i + 1) = nq(i) + d;
end

%% Iteraciones de nq y Z:

% Z_max = 3;                   %      MÁXIMO NÚMERO DE INYECTORES
Co = (2*g*H)^0.5;
Matrix_mother = zeros(8,r,Z_max);
for i = 1:Z_max
    Z = i;
    n_ = nq.*((Z)^0.5*H^0.75/(Q_total^0.5));
    p_ = f./n_;
    p = ceil(p_);
    n = f./p;
    n_rpm = 60.*n;
    nq_ = n.*(Q_total/Z)^0.5/H^0.75;
    Dm = k_m*Co./(pi.*n);
    B = Dm.*nq_./(2.63.*k_m.*phi_b.^0.5);
    Alphao = acos((1+2.*nq_).^(-1));
    Lambda = 0.5./(1-k_m.*(tan(Alphao)./Alphao));
    N = round(pi*(2*Lambda-1)./(k_m.*(nq_.*(1+nq_)).^0.5));
    do_ = ((4*Q_total)/(Z*Co*pi))^0.5;
    do = ones(1,r);
    do = do_.*do;
    N_nozzles = B./do;
    
    
    Matrix_mother(1,:,Z) = n_rpm;
    Matrix_mother(2,:,Z) = nq_;
    Matrix_mother(3,:,Z) = p;
    Matrix_mother(4,:,Z) = Dm;
    Matrix_mother(5,:,Z) = B;
    Matrix_mother(6,:,Z) = N;
    Matrix_mother(7,:,Z) = do;
    Matrix_mother(8,:,Z) = N_nozzles;
    
end

Effht = 2*k_m*(1-k_m )*(1-cos(Beta_2*pi/180));

%% Depuración por el generador (pares de polos):
% Los posibles pares de polos son 4, 6, 8, 10

M_p = zeros(8,1,1);            %    MATRIZ DE ARREGLADA AL NÚMERO DE POLOS
for i = 1:Z_max
    count = 0;
    cond_ = 0;
    for j = 1:r
        cond = Matrix_mother(3, j, i);
        if cond==4 || cond==6 || cond==8
            if cond ~= cond_
                count = count + 1;
                M_p(:,count,i) = Matrix_mother(:,j,i);
                cond_ = M_p(3,count,i);
            end
        end        
    end
end

%% Depuración de forma (Dm/do >7):

M_df = zeros(8, 1, 1);
[m,n,o] = size(M_p);
cond_f = M_p(4,:,:)./M_p(7,:,:);
for i = 1:Z_max
    count = 0;
    for j = 1:n
        if cond_f(1,j,i) > 7
            count = count + 1;
            M_df(:, count,i) = M_p(:,j,i);
        end
    end
end
Velocidad_Angular=Matrix_mother(1,:);
Velocidad_Especifica=Matrix_mother(2,:);
Numero_de_Polos=Matrix_mother(3,:);
Diametro_de_Paso=Matrix_mother(4,:);
Ancho_de_cuchara=Matrix_mother(5,:);
Numero_de_cucharas=Matrix_mother(6,:);
Diametro_de_Chorro=Matrix_mother(7,:);
Numero_de_Inyectores=Matrix_mother(8,:);

plot(Velocidad_Angular,Velocidad_Especifica, 'b-o')
ylabel ('Velocidad específica');
xlabel ('Velocidad angular del rodete [RPM]');
title ('Turbina Pelton; P=45,61 [kW]; Z=3 [Inyectores]');
grid on;
pause;
plot(Velocidad_Angular,Numero_de_Polos, 'g-o')
ylabel ('Numero de Polos');
xlabel ('Velocidad angular del rodete [RPM]');
title ('Turbina Pelton; P=45,61 [kW]; Z=3 [Inyectores]');
grid on;
pause;
plot(Velocidad_Angular,Diametro_de_Paso, 'r-o')
ylabel ('Diametro de Paso [m]');
xlabel ('Velocidad angular del rodete [RPM]');
title ('Turbina Pelton; P=45,61 [kW]; Z=3 [Inyectores]');
grid on;
pause;
plot(Velocidad_Angular,Ancho_de_cuchara, 'k-o')
ylabel ('Ancho de cuchara [m]');
xlabel ('Velocidad angular del rodete [RPM]');
title ('Turbina Pelton; P=45,61 [kW]; Z=3 [Inyectores]');
grid on;
pause;
plot(Velocidad_Angular,Numero_de_cucharas, 'm-o')
ylabel ('Número de cucharas');
title ('Turbina Pelton; P=45,61 [kW]; Z=3 [Inyectores]');
title ('Turbina Pelton; P=45,61 [kW]');
grid on;
pause;
plot(Velocidad_Angular,Diametro_de_Chorro, 'c-o')
ylabel ('Diámetro de Chorro [m]');
xlabel ('Velocidad angular del rodete [RPM]');
title ('Turbina Pelton; P=45,61 [kW]; Z=3 [Inyectores]');
grid on;