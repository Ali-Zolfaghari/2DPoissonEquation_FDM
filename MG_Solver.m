function [ U,RES ] = MG_Solver( FUNC,x,y,Lx,Ly,M,N,MaxI,MaxE,MAXERROR,Ws )

M2 = M/2;
N2 = N/2;
M4 = M/4;
N4 = N/4;
M8 = M/8;
N8 = N/8;

dx = Lx/M;
dy = Ly/N;
dx2 = Lx/M2;
dy2 = Ly/N2;
dx4 = Lx/M4;
dy4 = Ly/N4;
dx8 = Lx/M8;
dy8 = Ly/N8;

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
[ Ap2,An2,As2,Ae2,Aw2,MAT2 ] = MatrixCalc( dx2,dy2,M2,N2 );
[ Ap4,An4,As4,Ae4,Aw4,MAT4 ] = MatrixCalc( dx4,dy4,M4,N4 );
[ Ap8,An8,As8,Ae8,Aw8,MAT8 ] = MatrixCalc( dx8,dy8,M8,N8 );

ITER = 1;
Residual = 1.0;
RES(ITER,1) = ITER;
RES(ITER,2) = Residual;
fprintf('Iteration : %d \n',ITER);

while (Residual > MAXERROR)
    
    U2h = zeros(1,(M2+1)*(N2+1));
    U4h = zeros(1,(M4+1)*(N4+1));
    U8h = zeros(1,(M8+1)*(N8+1));
    U2g = zeros(1,(M2+1)*(N2+1));
    U4g = zeros(1,(M4+1)*(N4+1));
    U8g = zeros(1,(M8+1)*(N8+1));
    
    [ Uh,Res,Iter ] = GSSolver( UII,UI,B,Ap,An,As,Ae,Aw,M,N,Ws,MaxI,MaxE,1 );
    
    [ Rh ] = ResidualCalc( Uh,B,Ap,An,As,Ae,Aw,M,N );
    
    [ R2h ] = RestrictCalc( Rh,M,N );
    
    [ U2h,Res,Iter ] = GSSolver( U2h,U2g,R2h,Ap2,An2,As2,Ae2,Aw2,M2,N2,Ws,MaxI,MaxE,1 );
    
    [ R2h ] = ResidualCalc( U2h,R2h,Ap2,An2,As2,Ae2,Aw2,M2,N2 );
    
    [ R4h ] = RestrictCalc( R2h,M2,N2 );
    
    [ U4h,Res,Iter ] = GSSolver( U4h,U4g,R4h,Ap4,An4,As4,Ae4,Aw4,M4,N4,Ws,MaxI,MaxE,1 );
    
    [ R4h ] = ResidualCalc( U4h,R4h,Ap4,An4,As4,Ae4,Aw4,M4,N4 );
    
    [ R8h ] = RestrictCalc( R4h,M4,N4 );
    
    [ U8h,Res,Iter ] = GSSolver( U8h,U8g,R8h,Ap8,An8,As8,Ae8,Aw8,M8,N8,Ws,MaxI,MaxE,2 );
    
    [ I4h ] = ProlongateCalc( U8h,M8,N8 );
    
    U4h = U4h+I4h;
    
    [ I2h ] = ProlongateCalc( U4h,M4,N4 );
    
    U2h = U2h+I2h;
    
    [ Ih ] = ProlongateCalc( U2h,M2,N2 );
    
    UII = Uh+Ih;
    
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
    
    UI = UII;
    
end

U = zeros(M+1,N+1);
for i = 1:M+1
    for j = 1:N+1
        ip = (j-1)*(M+1)+i;
        U(i,j) = UII(ip);
    end
end

end

