function [ U,RES ] = CGP_Solver( FUNC,x,D,Lx,Ly,M,N,MAXERROR )

dx = Lx/M;
dy = Ly/N;

B = zeros((M+1)*(N+1),1);
X = zeros((M+1)*(N+1),1);

for i = 1:M+1
    for j = 1:N+1
        ip = (j-1)*(M+1)+i;
        B(ip,1) = FUNC(x(i),D(j));
    end
end

i = 1;
for j = 1:N+1
    ip = (j-1)*(M+1)+i;
    B(ip,1) = 0.0;
end
i = M+1;
for j = 1:N+1
    ip = (j-1)*(M+1)+i;
    B(ip,1) = 0.0;
end
j = 1;
for i = 2:M
    ip = (j-1)*(M+1)+i;
    B(ip,1) = 0.0;
end
j = N+1;
for i = 2:M
    ip = (j-1)*(M+1)+i;
    B(ip,1) = 0.0;
end

[ AD,An,As,Ae,Aw,MAT ] = MatrixCalc( dx,dy,M,N );

ITER = 1;
Residual = 1.0;
RES(ITER,1) = ITER;
RES(ITER,2) = Residual;
fprintf('Iteration : %d \n',ITER);

A = MAT;

% MA = ILU(A);
setup.type = 'nofill';
setup.milu = 'off';
setup.droptol = 0.1;
[L,U] = ilu(sparse(A),setup);
MA = L*U;

R = B-A*X;
Z = MA\R;
D = Z;
RTZ_old = R'*Z;
NORM = sqrt(abs(RTZ_old));

while (Residual > MAXERROR && ITER < ((M+1)*(N+1)))
    
    AD = A*D;
    Alpha = RTZ_old/(D'*AD);
    X = X+Alpha*D;
    R = R-Alpha*AD;
    Z = MA\R;
    RTZ_new = R'*Z;
    Beta = RTZ_new/RTZ_old;
    D = Z+Beta*D;
    RTZ_old = RTZ_new;

    i = 1;
    for j = 1:N+1
        ip = (j-1)*(M+1)+i;
        X(ip,1) = 0.0;
    end
    i = M+1;
    for j = 1:N+1
        ip = (j-1)*(M+1)+i;
        X(ip,1) = 0.0;
    end
    j = 1;
    for i = 2:M
        ip = (j-1)*(M+1)+i;
        X(ip,1) = 0.0;
    end
    j = N+1;
    for i = 2:M
        ip = (j-1)*(M+1)+i;
        X(ip,1) = 0.0;
    end
    
    Residual = sqrt(abs(RTZ_new))/NORM;
    
    RES(ITER+1,1) = ITER+1;
    RES(ITER+1,2) = Residual;
    
    ITER = ITER+1;
    fprintf('Iteration : %d \n',ITER);
    
end

U = zeros(M+1,N+1);
for i = 1:M+1
    for j = 1:N+1
        ip = (j-1)*(M+1)+i;
        U(i,j) = X(ip,1);
    end
end

end

