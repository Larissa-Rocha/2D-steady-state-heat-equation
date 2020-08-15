clc
clear all

%% Procedimento para solução da equação 2D permanente do calor - SOR by lines

%Parâmetros do problema
Lx = 1;           %comprimento total em x
Ly = 1;           %comprimento total em y
dx = 0.05;        %passo no espaço em x
dy = 0.05;        %passo no espaço em y
m = Lx/dx + 1;   %número de nós (m=nx=ny)
omega = 1;        %relaxation parameter
beta = dx/dy;     %razão de aspecto de malha
a = -omega/(2*(1+beta^2));
iteration = 1000;

%Discretização do domínio
x = linspace(0, Lx, m);
y = linspace(0, Ly, m);
[X Y] = meshgrid(x, y);     %criando uma malha 2D x-y

%Condição inicial
for t = 1:iteration
    for j = 1 : m
        for i = 1 : m
            te(t,j,i)=0;
        end
    end
end

%Condições de contorno
for t = 1 : iteration+1
    for i = 1 : m+1
        te(t,1,i) = sin(pi*x);  %temperatura no topo
        te(t,m+1,i) = 0;    %temperatura na base
    end
end

for t = 1 : iteration+1
    for j = 2 : m+1
        te(t,j,1) = 0;  %temperatura à esquerda
        te(t,j,m+1) = 0;    %temperatura à direita
    end
end
te2(:,:) = te(iteration,:,:);

%Formação da matriz A dos coeficientes
for j = 1 : m
    
    for i = 1 : m
        
        if i==j
            ma(j,i) = 1;
        elseif i==j+1
            ma(j,i) = a;
        elseif i==j-1
            ma(j,i)=a;
        end
        
    end %fim do for do x
    
end %fim do for do y

for t = 1 : iteration
    for j = 2 : m
        for i = 2 : m
            k(i-1) = (1-omega)*te(t,j,i-1) + te(t,j+1,i-1) + te(t+1,j-1,i-1);
        end %fim do for do x 
        
        %Algoritmo de Thomas
        md(1,1) = 1;
        for j1 = 1 : m-1
            for i = 1 : m-1
                if j1~=1
                    if i==j1
                        md(j1,i) = ma(j1,i) - (a^2)/md(j1-1,i);
                    end
                end
                if i==j1+1
                    md(j1,i) = a;
                end
            end
        end
        
        for i = 2 : m-1
            k(i) = k(i) - k(i-1)*a/md(i-1,i-1);
        end
        te(t,j,m)=k(m-1)/md(m-1,m-1);
        
        for i = 1 : m-1
            i = m-i;
            te(t,j,i)=(k(i)-a*te(t,j,i+1))/md(i,i);
        end
            
    end %fim do for do y
end %fim do for da iteração 

te1(:,:)=te(iteration,:,:);
