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
tol = 1e-4;       %toler�ncia
error = 10;       %inicializando o erro
k = 0;            %inicializando a itera��o
omega = 2;        %relaxation parameter

%Discretiza��o do dom�nio
x = linspace(0, Lx, nx);
y = linspace(0, Ly, ny);
[X Y] = meshgrid(x, y);     %criando uma malha 2D x-y

%Condi��es de contorno
T = ones(nx,ny);       %inicializando a matriz temperatura                 
T(ny,:) = sin(pi*x);    %T no topo
T(:,1) = 0;             %T � esquerda
T(:,nx) = 0;            %T � direita
Told = T;

%Coeficientes da matriz tridiagonal
a = omega/4;        %diagonal superior
b = omega/4;        %diagonal inferior
d = 1;              %diagonal principal
%transformar os escalares a, b e c em vetores linha
aa = a*ones(1,nx-3);
bb = aa;
dd = d*ones(1,nx-2);

while (error > tol)
    
    %para a primeira linha
    %definindo o vetor da direita
    for i = 2 : nx-2
        c(1,i)=(1-omega)*T(1,i) + (omega/2)*T(2,i); %c � um vetor linha (c_1)
    end
    cc = c(1,:);
    %chamando o TDMA para resolver a primeira linha
    B = [0,-bb];
    D = dd;
    A = -aa;
    C = cc;
    U = Tridiag(B,D,A,C);  %solu��o dos n�s internos da primeira linha
    
    %para as linhas seguintes
    for j = 2 : ny-1
        for i = 2 : nx-1
            c(j,i) = (1-omega)*T(j,i) + (omega/4)*T(j+1,i) + T(j-1,i);  %c � um vetor linha (c_j)
            %chamando o TDMA para resolver as linhas seguintes
            B = [0,-bb];
            D = dd;
            A = -aa;
            C = c;
            UU(j-1,:) = Tridiag(B,D,A,C(j-1,:));  %solu��o dos n�s internos das linhas 2 at� ny-1
        end
    end
    %montando a solu��o
    T(1,2:nx-1) = U; %armazenando a primeira linha
    for j = 2:ny-1
        for i = 2:nx-1
           T(j,i) = UU(j-1,i-1);   %armazenado as linhas seguintes
        end
    end

    error_T = max(abs(Told - T));   %calcula o erro de parada do while
    error = max(error_T);
    
    Told = T;   %atualizando o T
    k = k + 1;  %incremento da itera��o
    
end