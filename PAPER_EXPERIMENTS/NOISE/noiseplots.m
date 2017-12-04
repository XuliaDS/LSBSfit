close all;
clear all;

folder ='DAT3/';

A = importdata(fullfile(folder, 'ECILS'));
B = importdata(fullfile(folder, 'DLS_short'));
C = importdata(fullfile(folder, 'DLS_long'));
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



h2 = plot (B(:,1), B(:,2), 'x', 'markers', 5, 'LineWidth',0.5);
h3 = plot (C(:,1), C(:,2), '--',  'markers',5,'LineWidth',2);
h1 = plot (A(:,1), A(:,2), 'o',  'markers',4,'LineWidth',0.5);
h4 = plot (D(:,1), D(:,2), '-',  'markers', 1,'LineWidth',1);
ax = gca % Get handle to currenkt axes.
ax.XColor = 'w'; % Red
ax.YColor = 'w'; %
hold off
if strcmp(folder,'DAT1/') 
h = legend({ '\begin{tabular}{p{0.8cm}l}DLS&n = 51\end{tabular}',...
'\begin{tabular}{p{0.8cm}r}DLS&n = 200\end{tabular}',...
'\begin{tabular}{p{0.8cm}r}ECILS&n = 51\end{tabular}',...
'\begin{tabular}{p{0.8cm}r}Data&\end{tabular}'}...
,'FontSize',16,'Interpreter','latex','Position',rect)
elseif strcmp(folder,'DAT2/')
h = legend({ '\begin{tabular}{p{0.8cm}l}DLS&n = 47\end{tabular}',...
'\begin{tabular}{p{0.8cm}r}DLS&n = 200\end{tabular}',...
'\begin{tabular}{p{0.8cm}r}ECILS&n = 47\end{tabular}',...
'\begin{tabular}{p{0.8cm}r}Data&\end{tabular}'}...
,'FontSize',16,'Interpreter','latex','Position',rect)
else 
    h = legend({ '\begin{tabular}{p{0.8cm}l}DLS&n = 45\end{tabular}',...
'\begin{tabular}{p{0.8cm}r}DLS&n = 200\end{tabular}',...
'\begin{tabular}{p{0.8cm}r}ECILS&n = 45\end{tabular}',...
'\begin{tabular}{p{0.8cm}r}Data&\end{tabular}'}...
,'FontSize',16,'Interpreter','latex','Position',rect)

end

set(gca,'position',[0 0 1 1])
axis off

set(gca,'xtick',[],'ytick',[])
name = fullfile(folder,'plot')
print(name,'-dpdf','-besfit')

