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
T=ones(nx,ny);
for j=1:nx
    T(ny,j)=T(ny-1,j);  %condi��o da base
    T(1,j)=sin(pi*x(j));    %topo
end
for i=1:ny
    T(i,nx)=0;  %direita
    T(i,1)=0;   %esquerda
end

Tant_meio=ones(npx,npy);
for j=1:npx
    Tant_meio(npy,j)=Tant_meio(npy-1,j);
    Tant_meio(1,j)=sin(pi*x(j));
end
for i=1:npy
    Tant_meio(i,npx)=0;
    Tant_meio(i,1)=0;
end
