close all
clear all

FONT_CP = 10;
FONT_KNOTS = 12;
FONTSIZE_LEGEND = 12;
A = importdata('KNOTS');
B = importdata('Bspline');
C = importdata('CP');
sizeA = size(A,1);
for i = 1:sizeA 
    A(i,2) = -0.65;
end
format long
figure(1)

hold on
h1 = plot (B(:,1), B(:,2), '.',  'markers', 5);
h2 = plot (C(:,1), C(:,2), '-o', 'markers', 3);
h3 = plot (A(:,1), A(:,2), '-+', 'markers', 3);
text (A(1,1)-0.05,A(1,2)-0.15,'$u_0$','FontSize',FONT_KNOTS,'Interpreter','latex');
text (A(2,1)-0.05,A(1,2)-0.3,'$u_1$','FontSize',FONT_KNOTS,'Interpreter','latex');
text (A(3,1)-0.05,A(1,2)-0.45,'$u_2$','FontSize',FONT_KNOTS,'Interpreter','latex');
text (A(4,1)-0.05,A(1,2)-0.6,'$u_3$','FontSize',FONT_KNOTS,'Interpreter','latex');
text (A(5,1)-0.05,A(5,2)-0.15,'$u_4$','FontSize',FONT_KNOTS,'Interpreter','latex');
text (A(9,1)-0.05,A(9,2)-0.15,'$u_8$','FontSize',FONT_KNOTS,'Interpreter','latex');
text (A(10,1)-0.05,A(13,2)-0.15,'$u_{9}$','FontSize',FONT_KNOTS,'Interpreter','latex');
text (A(11,1)-0.05,A(13,2)-0.3,'$u_{10}$','FontSize',FONT_KNOTS,'Interpreter','latex');
text (A(12,1)-0.05,A(13,2)-0.45,'$u_{11}$','FontSize',FONT_KNOTS,'Interpreter','latex');
text (A(13,1)-0.05,A(13,2)-0.6,'$u_{12}$','FontSize',FONT_KNOTS,'Interpreter','latex');
text (C(1,1) - 0.1,C(1,2) + 0.2,'$\mathbf P_{0}$','FontSize',FONT_CP,'Interpreter','latex');

text (C(3,1) - 0.1,C(3,2) + 0.2,'$\mathbf P_{2}$','FontSize',FONT_CP,'Interpreter','latex');
text (C(9,1) + 0.05,C(9,2) + 0.1,'$\mathbf P_{8}$','FontSize',FONT_CP,'Interpreter','latex');


ax = gca % Get handle to current axes.
ax.XColor = 'w'; % Red
ax.YColor = 'w'; % Blue
hold off

rect = [0.45, 0.6, .0, .05]
h = legend({'B-spline', 'Control Polygon','Knot Sequence'},'FontSize',FONTSIZE_LEGEND,'Interpreter','latex','Orientation','horizontal','Position',rect);%Location','northeastoutside')
set(gca,'position',[0.1 0.1 0.7 0.5])
axis off

set(gca,'xtick',[],'ytick',[])
print('sineBspline','-dpdf','-besfit')
