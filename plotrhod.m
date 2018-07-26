function plotrhod

figure;hold on;box on;

load('rhodemtdata-T70.mat');
plot(nplt,rhod,'-k','LineWidth',3); 
hold on;
load('rhodemtdata-T130.mat');
plot(nplt,rhod,'-r','LineWidth',3); 
load('rhodemtdata-T190.mat');
plot(nplt,rhod,'-b','LineWidth',3); 

set(gca,'FontSize',20);
% title('(a)','FontSize', 20, 'Interpreter', 'latex');
% text(-55,23,'(a)','FontSize',30);
text(-55,19,'Inhomogeneous','FontSize',22);
text(-55,17,'momentum drag','FontSize',22);
text(-44,15,'$(\eta <0)$','Interpreter', 'latex','FontSize',22);


text(22,11,'$n_{P}=0$', 'Interpreter', 'latex','FontSize',25);
legend({'$T=70K$' '$T=130K$' '$T=190K$'}, 'Interpreter', 'latex','FontSize',20, 'Location', 'NorthEast')
% xlim([-30 30])
% ylim([0 30])
xlabel('$n_{A} (10^{10} \mathrm{cm}^{-2})$', 'FontSize', 30, 'Interpreter', 'latex');
ylabel('$ \tilde{\rho}^{\mathrm{EMT}}_{D} (\Omega)$', 'FontSize', 30, 'Interpreter', 'latex');

h=gca;
h.XTick=(-60:20:60);
h.XMinorTick='on';
h.YMinorTick='on';
legend({'$T=70K$' '$T=130K$' '$T=190K$' }, 'Interpreter', 'latex','FontSize',20, 'Location','NorthEast','Orientation','Vertical')

print('-dpdf','Fig5.pdf')

end