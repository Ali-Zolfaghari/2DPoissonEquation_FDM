%***************************************************************************************************
%*   Solve 2D poisson equation by presented code.
%*   I take no responsibilities for any errors in the code or damage thereby.
%*   Please notify me at zolfaghari1992iut@gmail.com if the code is used in any type of application.
%***************************************************************************************************
%*   Developer   : Ali Zolfaghari Sichani (14-08-2018)
%***************************************************************************************************
%*   References  : 
%*   Computational Fluid Mechanics and Heat Transfer.
%*   by John C. Tannehill (Author), Dale Anderson (Author), Richard H. Pletcher (Author).
%***************************************************************************************************
%*   Poisson Equation in two-dimensional square domain. (solving by centeral difference scheme)   :   
%*   Uxx + Uyy = f(x,y)
%*   1 : SUCCESSIVE OVER RELAXATION
%*   2 : MULTIGRID + GAUSS SEIDEL
%*   3 : CONJUGATE GRADIENT
%*   4 : PRECONDITIONED CONJUGATE GRADIENT
%*   Inputs      :
%*   M           (number of division of domain in x-direction      )
%*   N           (number of division of domain in y-direction      )
%*   Lx          (length of domain in x-direction                  )
%*   Ly          (length of domain in y-direction                  )
%*   Wsor_MG     (relaxation factor in MULTIGRID + GAUSS SEIDEL    )
%*   Wsor_SOR    (relaxation factor in SUCCESSIVE OVER RELAXATION  )
%*   MaxI_MG     (max. allowable itrerations of gauss              )
%*   MaxE_MG     (max. allowable error of gauss                    )
%*   MAXERROR    (max. allowable error                             )
%*   Outputs     :
%*   plot numerical solutions
%***************************************************************************************************


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







