close all;
clear all;

folder ='DAT1/';
A = importdata(fullfile(folder, 'ECILS_COMP'));
B = importdata(fullfile(folder, 'KNOT_INF'));
C = importdata(fullfile(folder, 'BISECTION'));
D = importdata(fullfile(folder, 'data.txt'));


format long
figure(1)
hold on

if strcmp(folder,'DAT1/')
axis([min(C(:,1)) 1000 min(C(:,2)) max(D(:,2))])
rect = [0.1, 0.65, .25, .25]
elseif strcmp(folder,'DAT2/')
axis([min(B(:,1)) 1000 min(C(:,2)) max(D(:,2))])
rect = [0.1, 0.65, .25, .25]
else
axis([min(B(:,1)) 1000 min(C(:,2)) max(D(:,2))])
rect = [0.1, 0.65, .25, .25]
end


h1 = plot (A(:,1), A(:,2), 'o',  'markers',3,'LineWidth',0.5);
h2 = plot (B(:,1), B(:,2), 'x', 'markers', 2, 'LineWidth',0.5);
h3 = plot (C(:,1), C(:,2), '--',  'markers',2,'LineWidth',2);
h4 = plot (D(:,1), D(:,2), '-',  'markers', 1,'LineWidth',1);
h = legend({'Weighted','Max Point','Mid Point','Data'},'FontSize',16,'Interpreter','latex','Position',rect)
set(gca,'position',[0 0 1 1])
axis off

set(gca,'xtick',[],'ytick',[])
name = fullfile(folder,'knot_placement_plot')
print(name,'-dpdf','-besfit')

