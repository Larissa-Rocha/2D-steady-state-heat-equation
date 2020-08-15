clc
clear all

%% Procedimento para solu��o da equa��o 2D permanente do calor - SOR by lines

%Par�metros do problema
Lx = 1;           %comprimento total em x
Ly = 1;           %comprimento total em y
dx = 0.05;        %passo no espa�o em x
dy = 0.05;        %passo no espa�o em y
nx = Lx/dx + 1;   %n�mero de n�s em x
ny = Ly/dy + 1;   %n�mero de n�s em y
tol = 1e-4;       %crit�rio de erro
error = 10;       %inicializando o erro
k = 0;            %inicializando a itera��o
omega = 1.25;        %relaxation parameter

%Discretiza��o do dom�nio
x = linspace(0, Lx, nx);
y = linspace(0, Ly, ny);
[X Y] = meshgrid(x, y);     %criando uma malha 2D x-y

%Condi��es de contorno
T = zeros(nx,ny);       %inicializando o vetor temperatura                 
T(ny,:) = sin(pi*x);    %T no topo
T(:,1) = 0;             %T � esquerda
T(:,nx) = 0;            %T � direita

Told = T;



while error > tol
    k = k + 1;  %incremento da itera��o
    %Coeficientes da matriz tridiagonal
    b = - omega/4;    %diagonal superior
    d = 1;            %diagonal inferior
    a = - omega/4;    %diagonal principal
    for i = 2 : nx-1
       c(i-1) = (1-omega)*T(1,i) + (omega/2)*T(2,i);
    end
    
    %transformar os escalares a, b e c em vetores linha
    aa = a*ones(1,nx-2);
    bb = aa;
    dd = d*ones(1,nx-2);
    %chamando o TDMA
    UU = TDMAJ(1,nx-2,bb,dd,aa,c);
    %UU = Tridiag(bb,dd,aa,c);
    for i = 2 : nx-1
      T(1,i) = UU(i-1);
    end

    for j = 2 : ny-1
       aa = a*ones(1,nx-2);
       bb = aa;
       dd = d*ones(1,nx-2);
        
       for i = 2 : nx-1
          c(i-1) = (1-omega)*T(j,i) + (omega/4)*(T(j+1,i) + T(j-1,i));
       end
       UU = TDMAJ(1,nx-2,bb,dd,aa,c);
%        UU = Tridiag(bb,dd,aa,c);
       for i = 2 : nx-1
         T(j,i) = UU(i-1);
       end
    end
    
    error_T = max(abs(Told - T));   %calcula o erro de parada do while
    error = max(error_T);
    
    Told = T;   %atualizando o T
    
    
    %plotando os resultados - contornos
%     [C, h] = contourf(X, Y, T,40);
%     set(gca,'Xtick',(0:dx:1))
%     set(gca,'Ytick',(0:dy:1))
%     xlabel('coordenada x')
%     ylabel('coordenada y')
%     title('Distribui��o de temperatura (omega = 1)');
    
end