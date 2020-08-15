clc
clear all

%% Procedimento para solu��o da equa��o 2D permanente do calor - Gauss Seidel

%Par�metros do problema
Lx = 1;           %comprimento total em x
Ly = 1;           %comprimento total em y
dx = 0.05;        %passo no espa�o em x
dy = 0.05;        %passo no espa�o em y
nx = Lx/dx + 1;   %n�mero de n�s em x
ny = Ly/dy + 1;   %n�mero de n�s em y
tol = 1e-4;       %crit�rio de erro
error = 12;       %inicializando o erro
k = 0;            %inicializando a itera��o

%Discretiza��o do dom�nio
x = linspace(0, Lx, nx);
y = linspace(0, Ly, ny);
[X Y] = meshgrid(x, y);     %criando uma malha 2D x-y 

%Condi��es de contorno
T = zeros(nx,ny);            %inicializando a matriz temperatura
T(ny,:) = sin(pi*x);        %T no topo
T(:,1) = 0;                 %T � esquerda
T(:,nx) = 0;                %T � direita

Told = T;

while error > tol
    for i = 2:nx-1
        T(1,i) = 0.25*(T(1,i+1) + T(1,i-1) + 2*T(2,i)); %T na base
    end %fim do for do x

    for j = 2:ny-1
        for i = 2:nx-1
            T(j,i) = 0.25*(T(j,i+1) + T(j,i-1) + T(j-1,i) + T(j+1,i)); 
        end %fim do for do y
    end %fim do for do x
    
    error_T = max(abs(Told - T));   %calcula o erro de parada do while
    error = max(error_T);
    
    Told = T;   %atualizando o T
    k = k + 1;  %incremento da itera��o
    
    %plotando os resultados - contornos
%     [C, h] = contourf(X,Y,T,40);
%     set(gca,'Xtick',(0:dx:1))
%     set(gca,'Ytick',(0:dy:1))
%     xlabel('coordenada x')
%     ylabel('coordenada y')
%     title('Distribui��o de temperatura');
    
end %fim do while 

