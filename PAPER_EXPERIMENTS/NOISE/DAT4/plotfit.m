close all
clear all
A = importdata('ECILS');
B = importdata('data.txt');

figure(1)

hold on
h2 = plot(B(:,1), B(:,2), 'o', 'markers',3)
h1 = plot(A(:,1), A(:,2), '-',  'markers',3,'LineWidth',1);
ax = gca % Get handle to current axes.
ax.XColor = 'w'; % Red
ax.YColor = 'w'; % Blue
hold off

h = legend({'Data','Fit'},'FontSize',16,'Interpreter','latex','Location','northeast')
%set(gca,'position',[0.1 0.1 0.5 0.5])
axis off

set(gca,'xtick',[],'ytick',[])
print('wiggle','-dpdf','-besfit')