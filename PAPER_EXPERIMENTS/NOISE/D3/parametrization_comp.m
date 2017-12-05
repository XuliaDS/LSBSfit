close all;
clear all;

A = importdata('WEIGHTED');
B = importdata('MAX');
C = importdata('MID');
D = importdata('data.txt');


format long
figure(1)
hold on

rect = [0.1, 0.75, .25, .25]

axis([min(C(:,1)) 1200 min(C(:,2)) max(D(:,2))])

h1 = plot (A(:,1), A(:,2), 'o',  'markers',3,'LineWidth',0.5);
h2 = plot (B(:,1), B(:,2), 'x', 'markers', 2, 'LineWidth',0.5);
h3 = plot (C(:,1), C(:,2), '--',  'markers',2,'LineWidth',2);
h4 = plot (D(:,1), D(:,2), '-',  'markers', 1,'LineWidth',1);
h = legend({'Weighted','Max Point','Mid Point','Data'},'FontSize',16,'Interpreter','latex','Position',rect)
set(gca,'position',[0 0 1 1])
axis off

set(gca,'xtick',[],'ytick',[])
name = ('knot_placement_plot')
print(name,'-dpdf','-besfit')

