close all
clear all
plotWiggle = 1
plotNaca   = 0
plotSnake  = 0



if plotWiggle ==1
A = importdata('wiggle_spline');
B = importdata('wiggle_optimal');
C = importdata('ice_matlab.dat');

nPt = 70
OFFSET = 70
AA = ones(nPt,2);
BB = ones(nPt,2);
CC = ones(nPt,2);
for i = 1: nPt
 AA(i,1)= A(i+OFFSET,1);
 AA(i,2)= A(i+OFFSET,2);
 BB(i,1)= B(i+OFFSET,1);
 BB(i,2)= B(i+OFFSET,2);
 CC(i,1)= C(i+OFFSET,1);
 CC(i,2)= C(i+OFFSET,2);
 end
format long
figure(1)

hold on
h1 = plot (AA(:,1), AA(:,2), '.',  'markers', 5);
h2 = plot (BB(:,1), BB(:,2), 'o',  'markers', 5);
h3 = plot (CC(:,1), CC(:,2), '--k', 'markers', 3);
ax = gca % Get handle to current axes.
ax.XColor = 'w'; % Red
ax.YColor = 'w'; % Blue
hold off
h = legend({'Uniform', 'Centripetal',  'Data'},'FontSize',16,'Interpreter','latex','Location','northeastoutside')
set(gca,'position',[0.1 0.1 0.6 0.6])
axis off

set(gca,'xtick',[],'ytick',[])
print('wiggle','-dpdf','-besfit')

end
%%

if plotNaca ==1
A = importdata('nacaM.dat');
B = importdata('nacaSpline');
nPt = 200;

AA = ones(nPt,2);
BB = ones(nPt,2);

divA = floor( size(A,1)/ (nPt-1))
divB = floor( size(B,1)/ (nPt-1))
for i = 1: nPt
    i
 AA(i,1)= A(1 + (i-1)*divA,1);
 AA(i,2)= A(1 + (i-1)*divA,2);
 BB(i,1)= B(1+(i-1)*divB,1);
 BB(i,2)= B(1+(i-1)*divB,2);
 
end

format long
figure(2)

hold on
axis on
h2 = plot (AA(:,1),AA(:,2),'+-','markers',4);
h1 = plot (BB(:,1), BB(:,2), 'o-','markers',4);
hold off
h = legend({'Data','Bspline'},'FontSize',32,'Interpreter','latex','Location','northeastoutside')
set(gca,'xtick',[],'ytick',[])
ax = gca % Get handle to current axes.
print('naca27','-dpdf','-besfit','-fillpage')


end
%%

%%

if plotSnake == 1
A = importdata('snakeM.dat');
B = importdata('snakeSpline');
nPt = 250;


divA = floor( size(A,1)/ (nPt-1))
divB = floor( size(B,1)/ (nPt-1))
AA = ones(nPt ,3);
BB = ones(nPt  ,3);
for i = 1: nPt
 AA(i,1)= A(1 + (i-1)*divA,1);
 AA(i,2)= A(1 + (i-1)*divA,2);
 AA(i,3)= A(1 + (i-1)*divA,3);
 BB(i,1)= B(1 + (i-1)*divB,1);
 BB(i,2)= B(1 + (i-1)*divB,2);
 BB(i,3)= B(1 + (i-1)*divB,3);
 
end
   
format long
figure(3)
hold on
rotate3d on
h1 = plot3 (AA(:,1), AA(:,2),AA(:,3), '+-','markers',4);
h2 = plot3 (BB(:,1), BB(:,2),BB(:,3), 'o-','markers',4);

hold off
h = legend({'Data','Bspline'},'FontSize',32,'Interpreter','latex','Location','northeastOutside')
set(gca,'xtick',[],'ytick',[],'ztick',[])
%axis off
gca = figure(3)
direction = [0 0 1];
rotate([h1,h2],direction,45)
print('snake','-dpdf')
end
%%
