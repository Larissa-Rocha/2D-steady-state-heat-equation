clc
clear all

%% Procedimento para solução da equação 2D permanente do calor - ADI com SOR

%Parâmetros do problema
Lx = 1;           %comprimento total em x
Ly = 1;           %comprimento total em y
dx = 0.05;        %passo no espaço em x
dy = 0.05;        %passo no espaço em y
nx = Lx/dx + 1;   %número de nós em x
ny = Ly/dy + 1;   %número de nós em y
tol = 1e-4;       %tolerância
error = 10;       %inicializando o erro
itera = 0;        %inicializando a iteração
omega = 1;        %relaxation parameter

%Discretização do domínio
x = linspace(0, Lx, nx);
y = linspace(0, Ly, ny);
[X Y] = meshgrid(x, y);     %criando uma malha 2D x-y

%Condições de contorno 
T=ones(nx,ny);
for j=1:nx
    T(ny,j)=T(ny-1,j);  %condição da base
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
