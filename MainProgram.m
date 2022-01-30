clear,clc,close all
format compact
format long



% input
M = 16;
N = 16;

Lx = 1.0;
Ly = 1.0;

MAXERROR = 0.000001;

MaxI_MG = 5;
MaxE_MG = 0.000001;
Wsor_MG = 0.85;
Wsor_SOR = 0.85;

func = @(x,y)(sin(pi*x)*sin(pi*y));

XPlot = 0.5*Lx;

% initial
for i = 1:M+1
    for j = 1:N+1
        ip = (j-1)*(M+1)+i;
        x(i) = (i-1)*(Lx/M);
        y(j) = (j-1)*(Ly/N);
    end
end

% solve
fprintf('=========================== SOR =============================\n');

[ U_SOR,Res_SOR ] = SOR_Solver( func,x,y,Lx,Ly,M,N,MAXERROR,Wsor_SOR );

fprintf('=========================== MG =============================\n');

[ U_MG,Res_MG ] = MG_Solver( func,x,y,Lx,Ly,M,N,MaxI_MG,MaxE_MG,MAXERROR,Wsor_MG );

fprintf('=========================== CG =============================\n');

[ U_CG,Res_CG] = CG_Solver( func,x,y,Lx,Ly,M,N,MAXERROR );

fprintf('=========================== CGP =============================\n');

[ U_CGP,Res_CGP ] = CGP_Solver( func,x,y,Lx,Ly,M,N,MAXERROR );

fprintf('============================================================\n');

% plot
PLOT;







