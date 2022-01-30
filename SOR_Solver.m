function [ U,RES ] = SOR_Solver( FUNC,x,y,Lx,Ly,M,N,MAXERROR,Ws )

dx = Lx/M;
dy = Ly/N;

B = zeros(1,(M+1)*(N+1));
UI = zeros(1,(M+1)*(N+1));
UII = zeros(1,(M+1)*(N+1));

for i = 1:M+1
    for j = 1:N+1
        ip = (j-1)*(M+1)+i;
        B(ip) = FUNC(x(i),y(j));
    end
end

i = 1;
for j = 1:N+1
    ip = (j-1)*(M+1)+i;
    B(ip) = 0.0;
    UI(ip) = 0.0;
    UII(ip) = 0.0;
end
i = M+1;
for j = 1:N+1
    ip = (j-1)*(M+1)+i;
    B(ip) = 0.0;
    UI(ip) = 0.0;
    UII(ip) = 0.0;
end
j = 1;
for i = 2:M
    ip = (j-1)*(M+1)+i;
    B(ip) = 0.0;
    UI(ip) = 0.0;
    UII(ip) = 0.0;
end
j = N+1;
for i = 2:M
    ip = (j-1)*(M+1)+i;
    B(ip) = 0.0;
    UI(ip) = 0.0;
    UII(ip) = 0.0;
end

[ Ap,An,As,Ae,Aw,MAT ] = MatrixCalc( dx,dy,M,N );

ITER = 1;
Residual = 1.0;
RES(ITER,1) = ITER;
RES(ITER,2) = Residual;
fprintf('Iteration : %d \n',ITER);

while (Residual > MAXERROR)
    for i = 2:M
        for j = 2:N
            ip = (j-1)*(M+1)+i;
            Sum = Aw(ip)*UII(ip-1)+As(ip)*UII(ip-M-1)+Ae(ip)*UI(ip+1)+An(ip)*UI(ip+M+1);
            UII(ip) = (B(ip)-Sum)/Ap(ip);
        end
    end
    
    i = 1;
    for j = 1:N+1
        ip = (j-1)*(M+1)+i;
        UI(ip) = 0.0;
        UII(ip) = 0.0;
    end
    i = M+1;
    for j = 1:N+1
        ip = (j-1)*(M+1)+i;
        UI(ip) = 0.0;
        UII(ip) = 0.0;
    end
    j = 1;
    for i = 2:M
        ip = (j-1)*(M+1)+i;
        UI(ip) = 0.0;
        UII(ip) = 0.0;
    end
    j = N+1;
    for i = 2:M
        ip = (j-1)*(M+1)+i;
        UI(ip) = 0.0;
        UII(ip) = 0.0;
    end
    
    Residual = sqrt(sum((UII-UI).^2))/((M+1)*(N+1));
    
    if (ITER == 1)
        NORM = Residual;
    end
    Residual = Residual/NORM;
    RES(ITER+1,1) = ITER+1;
    RES(ITER+1,2) = Residual;
    
    ITER = ITER+1;
    fprintf('Iteration : %d \n',ITER);
    
    UI = Ws*UII+(1.0-Ws)*UI;
end

U = zeros(M+1,N+1);
for i = 1:M+1
    for j = 1:N+1
        ip = (j-1)*(M+1)+i;
        U(i,j) = UII(ip);
    end
end

end

