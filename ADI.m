clc
clear all

%% Procedimento para solu��o da equa��o 2D permanente do calor - ADI com SOR

%Par�metros do problema
Lx = 1;           %comprimento total em x
Ly = 1;           %comprimento total em y
dx = 0.05;        %passo no espa�o em x
dy = 0.05;        %passo no espa�o em y
nx = Lx/dx + 1;   %n�mero de n�s em x
ny = Ly/dy + 1;   %n�mero de n�s em y
tol = 1e-4;       %toler�ncia
error = 10;       %inicializando o erro
itera = 0;        %inicializando a itera��o
omega = 1;        %relaxation parameter

%Discretiza��o do dom�nio
x = linspace(0, Lx, nx);
y = linspace(0, Ly, ny);
[X Y] = meshgrid(x, y);     %criando uma malha 2D x-y

%Condi��es de contorno
T = zeros(ny,nx);       %inicializando a matriz temperatura                 
T(ny,:) = sin(pi*x);    %T no topo
T(:,1) = 0;             %T � esquerda
T(:,nx) = 0;            %T � direita
Told = T;               %pegando os valores iniciais para a primeira itera��o

while error > tol       %loop de itera��es
    
    %VARRENDO POR LINHAS
    
    %Coeficientes da matriz tridiagonal
    b = -1;    %diagonal superior
    d = 4;     %diagonal inferior
    a = -1;    %diagonal principal
    
    for i = 2 : nx-1 %para os n�s da base (cond de contorno)
       c(i-1) = 2*Told(2,i);
    end
    
    %transformar os escalares a, b e c em vetores linha
    aa = a*ones(1,nx-2);
    bb = aa;
    dd = d*ones(1,nx-2);
    
    %chamando o TDMA
    UU = TDMAJ(1,nx-2,bb,dd,aa,c);
    
    for i = 2 : nx-1
      T(1,i) = UU(i-1);
    end
    
    %para os n�s internos
    for j = 2 : ny-1
       aa = a*ones(1,nx-2);
       bb = aa;
       dd = d*ones(1,nx-2);
        
       for i = 2 : nx-1
          c(i-1) = Told(j+1,i) + Told(j-1,i);
       end
       UU = TDMAJ(1,nx-2,bb,dd,aa,c);

       for i = 2 : nx-1
         T(j,i) = UU(i-1);
       end
    end
    
    T12 = T;   %criando a vari�ver T12 que corresponde a itera��o k+1/2
    
    %VARRENDO POR COLUNAS
    for i = 2 : nx-1
        for j = ny-1 : -1 : 1
            T(j,i) = T12(j,i+1) + Told(j,i-1) + Told(j+1,i);
        end
    end
    
    %APLICANDO SOR
    
    Told = T;   %atualizando o T
        
    error_T = max(abs(Told - T));   %calcula o erro de parada do while
    error = max(error_T);
        
    itera = itera + 1;
        
end