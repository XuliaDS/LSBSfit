close all;
clear all;


A = importdata('WEIGHTED');
B = importdata('data.txt');
%C = importdata(fullfile(folder, 'DLS_long'));
%D = importdata(fullfile(folder, 'data.txt'));
hold on
h2 = plot (B(:,1), B(:,2), '-+', 'markers', 2, 'LineWidth',0.5);
h1 = plot (A(:,1), A(:,2), '-o',  'markers',1,'LineWidth',0.5);

ax = gca % Get handle to currenkt axes.
ax.XColor = 'w'; % Red
ax.YColor = 'w'; %
hold off

format long
figure(1)
hold on

h = legend({'Data','Fit'},'FontSize',16,'Interpreter','latex','Location','northeast')


name = 'plot'
print(name,'-dpdf','-besfit')

