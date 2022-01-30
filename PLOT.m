
[ Mesh1,Mesh2 ] = meshgrid(x,y);

f = find( x >= XPlot , 1 );

s = (U_SOR(f,:)-U_SOR(f+1,:))./(x(f)-x(f+1));
UI_SOR = U_SOR(f,:)+s.*(XPlot-x(f));

s = (U_MG(f,:)-U_MG(f+1,:))./(x(f)-x(f+1));
UI_MG = U_MG(f,:)+s.*(XPlot-x(f));

s = (U_CG(f,:)-U_CG(f+1,:))./(x(f)-x(f+1));
UI_CG = U_CG(f,:)+s.*(XPlot-x(f));

s = (U_CGP(f,:)-U_CGP(f+1,:))./(x(f)-x(f+1));
UI_CGP = U_CGP(f,:)+s.*(XPlot-x(f));

figure(1);
contourf(Mesh1,Mesh2,U_SOR',15);title('SOR','FontWeight','bold','FontSize',14,'FontName','Times New Roman','FontAngle','italic');
colorbar('location','EastOutside');colormap('jet');

figure(2);
contourf(Mesh1,Mesh2,U_MG',15);title('MG','FontWeight','bold','FontSize',14,'FontName','Times New Roman','FontAngle','italic');
colorbar('location','EastOutside');colormap('jet');

figure(3);
contourf(Mesh1,Mesh2,U_CG',15);title('CG','FontWeight','bold','FontSize',14,'FontName','Times New Roman','FontAngle','italic');
colorbar('location','EastOutside');colormap('jet');

figure(4);
contourf(Mesh1,Mesh2,U_CGP',15);title('CGP','FontWeight','bold','FontSize',14,'FontName','Times New Roman','FontAngle','italic');
colorbar('location','EastOutside');colormap('jet');

figure(5);hold on;grid on;
semilogy(Res_SOR(:,1),Res_SOR(:,2),'b-','linewidth',1.5);
semilogy(Res_MG(:,1),Res_MG(:,2),'r-*','linewidth',1.5);
semilogy(Res_CG(:,1),Res_CG(:,2),'k-s','linewidth',1.5);
semilogy(Res_CGP(:,1),Res_CGP(:,2),'g-.','linewidth',1.5);
MAX = max([max(Res_SOR(:,1)),max(Res_MG(:,1)),max(Res_CG(:,1)),max(Res_CGP(:,1))]);
XL = 1:1:(MAX+10);
YL = MAXERROR*ones(1,length(XL));
semilogy(XL,YL,'c--','linewidth',2);
xlabel('Iteration','FontWeight','bold','FontSize',14,'FontName','Times New Roman','FontAngle','italic');
ylabel('Log \epsilon','FontWeight','bold','FontSize',14,'FontName','Times New Roman','FontAngle','italic');
legend('SOR','MG','CG','CGP');

figure(6);
surf(Mesh1,Mesh2,U_SOR');title('SOR','FontWeight','bold','FontSize',14,'FontName','Times New Roman','FontAngle','italic');
colorbar('location','EastOutside');colormap('jet');

figure(7);
surf(Mesh1,Mesh2,U_MG');title('MG','FontWeight','bold','FontSize',14,'FontName','Times New Roman','FontAngle','italic');
colorbar('location','EastOutside');colormap('jet');

figure(8);
surf(Mesh1,Mesh2,U_CG');title('CG','FontWeight','bold','FontSize',14,'FontName','Times New Roman','FontAngle','italic');
colorbar('location','EastOutside');colormap('jet');

figure(9);
surf(Mesh1,Mesh2,U_CGP');title('CGP','FontWeight','bold','FontSize',14,'FontName','Times New Roman','FontAngle','italic');
colorbar('location','EastOutside');colormap('jet');

figure(10);hold on;grid on;
plot(x,UI_SOR,'b-','linewidth',2.0);
plot(x,UI_MG,'r-*','linewidth',2.0);
plot(x,UI_CG,'k-s','linewidth',2.0);
plot(x,UI_CGP,'g-.','linewidth',2.0);
title(['U @ X = ',num2str(XPlot)],'FontWeight','bold','FontSize',14,'FontName','Times New Roman','FontAngle','italic');
xlabel('X','FontWeight','bold','FontSize',14,'FontName','Times New Roman','FontAngle','italic');
ylabel('U','FontWeight','bold','FontSize',14,'FontName','Times New Roman','FontAngle','italic');
legend('SOR','MG','CG','CGP');






