clear all
blue = [0 0.4470 0.7410];
red = [0.6350 0.0780 0.1840];
purple = [.63,.13,.94];
att = red;
rep = blue;
green = [0.4660 0.6740 0.1880];
syms x y z real
syms alpha beta gamma real
syms d1 d2 d3 real

d = [d1, d2, d3];
vf = [x*(alpha-alpha*x-(alpha+beta+d1)*y+d2*z);
    y*(beta+d1*x-beta*y-(beta+gamma+d3)*z);
    z*(gamma-(gamma+alpha+d2)*x+d3*y-gamma*z)] %system (1.4)
F = [x;y;z;x+y+z-1]
K = F;
for i = 1:size(F,1)
    K(i) = collect(simplify(gradient(F(i),[x,y,z]).'*vf/F(i)),[x,y,z]);
end
K

sp = solve(vf); sp = [sp.x,sp.y,sp.z]
%%
num = [-2,1,1,0,3,0] %parameters

vfi = subs(vf,[alpha,beta,gamma,d],num)

spi = solve(vfi);  spi = [spi.x,spi.y,spi.z] %equilibra

J = jacobian(vfi,[x,y,z]) %Jacobian matrix
temp = [];
for i = 1:size(spi,1)
    if min(spi(i,:)) >= 0
        temp = [temp;i];
    end
end
spi = spi(temp,:)
lam = [];
ji = jacobian(vfi,[x,y,z]);
attP = []; repP = [];
for i = 1:size(spi,1)
    disp('------------------------------------------')
    spi(i,:),jii = subs(ji,[x,y,z],spi(i,:));
    lam = double(expand(eig(jii)))
    if lam < 0
        attP = [attP,i];
    elseif lam > 0
        repP = [repP, i];
    end
end

%% 
% 

%determine the invariant surface
spi = double(spi)
H1 = x*y*z;
L = [0.5-0.5*z; 1-2*z]
p = [1/2,1/6,1/3];
c1 = subs(H1,[x,y,z],p)

A2 = -6*z^2 + 6*z - 1;
resz = solve(A2),double(resz)
zba = resz(resz<2/4)
pba = [subs(L,z,zba)',zba] 
c2 = subs(H1, [x,y,z], pba)

s = spi(end,:) 
%% 

%plot intersections
warning('off');
axisrange = [0,3.5,0,3.5,0,1.3];
figure; hold on;

space = linspace(0,1,100);
plot3(zeros(1,100),space*axisrange(4),zeros(1,100),'k--','linewidth',1.2);
plot3(zeros(1,100),zeros(1,100),space*axisrange(6),'k--','linewidth',1.2);
plot3(space*axisrange(2),zeros(1,100),zeros(1,100),'k--','linewidth',1.2);
text('Interpreter','latex','String','$O$','Position',[0,-0.2,-0.2], 'FontSize',20);
text('Interpreter','latex','String','$x_1$','Position',[axisrange(2)+.2,-.2,0.05],'FontSize',15);
text('Interpreter','latex','String','$x_2$','Position',[-.2,axisrange(4)-.01,0.14],'FontSize',15);
text('Interpreter','latex','String','$x_3$','Position',[0,-0.3,axisrange(6)],'FontSize',15);

ck = (c1+c2)/2
fimplicit3(H1-ck,[0,2,0,2,0,2],'FaceColor',blue,'FaceAlpha',.3,'LineStyle',"none")


res = double(solve(subs(H1^2-ck^2,[x;y],L),z));
HL2 = double(subs([L',z],z,res))
pc0 = zeros(2,3); j = 0;
for i = 1:size(HL2,1)
    if HL2(i,:) > 0
        j = j + 1;
        pc0(j,:) = HL2(i,:);
    end
end
pc0
plot3(pc0(:,1),pc0(:,2),pc0(:,3),'o','markerSize',5,'MarkerFaceColor',green,'MarkerEdgeColor',green)

zi = linspace(-1,min(pc0(:,3)),100);
xy = subs(L,z,zi);
plot3(xy(1,:),xy(2,:),zi,'k','LineWidth',.8)

zi = linspace(max(pc0(:,3)),1,100);
xy = subs(L,z,zi);
plot3(xy(1,:),xy(2,:),zi,'k','LineWidth',.8)

zi = linspace(min(pc0(:,3)),max(pc0(:,3)),100);
xy = subs(L,z,zi);
plot3(xy(1,:),xy(2,:),zi,'k--','LineWidth',.8)

view(60,-10)
axis equal; axis off;
axis(axisrange);

hold off
print('-dpng','-r300','fig2.3(i)_1.png');
%% 
% 

%plot phase portrait
warning('off');
axisrange = [0,3.5,0,3.5,0,1.3];
figure; hold on;

space = linspace(0,1,100);
plot3(zeros(1,100),space*axisrange(4),zeros(1,100),'k--','linewidth',1.2);
plot3(zeros(1,100),zeros(1,100),space*axisrange(6),'k--','linewidth',1.2);
plot3(space*axisrange(2),zeros(1,100),zeros(1,100),'k--','linewidth',1.2);
text('Interpreter','latex','String','$O$','Position',[0,-0.2,-0.2], 'FontSize',20);
text('Interpreter','latex','String','$x_1$','Position',[axisrange(2)+.2,-.2,0.05],'FontSize',15);
text('Interpreter','latex','String','$x_2$','Position',[-.2,axisrange(4)-.01,0.14],'FontSize',15);
text('Interpreter','latex','String','$x_3$','Position',[0,-0.3,axisrange(6)],'FontSize',15);

fill3([1,0,0],[0,1,0],[0,0,1],'--','FaceColor',[0.93,0.91,0.91],...
    'FaceAlpha',.6,'EdgeLighting',"none")%,'EdgeColor','none')

the = linspace(0,pi/2,100); [th1,th2] = meshgrid(the,the);

sp_ = [double(spi)];%eye(3)*2];

xi = linspace(0,1,100);
plot3(xi,1-xi,zeros(1,100),'k','LineWidth',1.2)
plot3(zeros(1,100),xi,1-xi,'k','LineWidth',1.2)
plot3(xi,zeros(1,100),1-xi,'k','LineWidth',1.2)

ck = (c1+c2)/2;
fimplicit3(H1-ck,[0,2,0,2,0,2],'FaceColor',blue,'FaceAlpha',.3,'LineStyle',"none")


res = double(solve(subs(H1^2-ck^2,[x;y],L),z));
HL2 = double(subs([L',z],z,res))
pc0 = zeros(2,3);
j = 0;
for i = 1:size(HL2,1)
    if HL2(i,:) > 0
        j = j + 1;
        pc0(j,:) = HL2(i,:);
    end
end
pc0
plot3(pc0(:,1),pc0(:,2),pc0(:,3),'o','markerSize',5,'MarkerFaceColor',green,'MarkerEdgeColor',green)

p0 = [0.4,0.2,0.2]; 
res = double(solve(subs(H1.^2-ck^2,[x,z],p0([1,3])),y));
p0(2) = res(2);
[~,Y] = ode45(@df,[0,20],p0); i1=1;i2=30;
plot3(Y(i1:i2,1),Y(i1:i2,2),Y(i1:i2,3),'color',red,'linewidth',1.2)

p0 = [1.2,2,.07];
res = double(solve(subs(H1.^2-ck^2,[x,z],p0([1,3])),y));
p0(2) = res(2);
[~,Y] = ode45(@df,[0,0.3],p0);
plot3(Y(:,1),Y(:,2),Y(:,3),'color',red,'linewidth',1.2)
[~,Y] = ode45(@df_pos,[0,8.4],p0);
plot3(Y(:,1),Y(:,2),Y(:,3),'color',red,'linewidth',1.2)

p0 = [0.45,2,0.11];
res = double(solve(subs(H1.^2-ck^2,[x,z],p0([1,3])),y));
p0(2) = res(2);
[~,Y] = ode45(@df,[0,61],p0);
i = 33; 
plot3(Y(1:i,1),Y(1:i,2),Y(1:i,3),'color',red,'linewidth',1.2)
[~,Y] = ode45(@df_pos,[0,100],p0);
i = 56; 
plot3(Y(1:i,1),Y(1:i,2),Y(1:i,3),'color',red,'linewidth',1.2)

p0 = [0.587,2.119,0.05];
res = double(solve(subs(H1.^2-ck^2,[x,z],p0([1,3])),y));
p0(2) = res(2);
[~,Y] = ode45(@df,[0,10],p0);i = 14;i1=1;
plot3(Y(i1:i,1),Y(i1:i,2),Y(i1:i,3),'color',red,'linewidth',1.2)
[~,Y] = ode45(@df_pos,[0,10],p0);i = 8;i1=1;
plot3(Y(i1:i,1),Y(i1:i,2),Y(i1:i,3),'color',red,'linewidth',1.2)

view(80,20)

axis equal; axis off;
axis(axisrange);
hold off
print('-dpng','-r300','fig2.3(i)_2.png');
%%
function vf = df(~,xyz)
x = xyz(1); y = xyz(2); z = xyz(3);
num = [-2,1,1,0,3,0];
alpha = num(1); beta = num(2); gamma = num(3);
d1 = num(4); d2 = num(5); d3 = num(6);
vf = [x*(alpha-alpha*x-(alpha+beta+d1)*y+d2*z);
    y*(beta+d1*x-beta*y-(beta+gamma+d3)*z);
    z*(gamma-(gamma+alpha+d2)*x+d3*y-gamma*z)];
end


function vf = df_pos(~,xyz)
x = xyz(1); y = xyz(2); z = xyz(3);
num = -[-2,1,1,0,3,0];
alpha = num(1); beta = num(2); gamma = num(3);
d1 = num(4); d2 = num(5); d3 = num(6);
vf = [x*(alpha-alpha*x-(alpha+beta+d1)*y+d2*z);
    y*(beta+d1*x-beta*y-(beta+gamma+d3)*z);
    z*(gamma-(gamma+alpha+d2)*x+d3*y-gamma*z)];
end

function vf = df_sp2(~,xyz)
x = xyz(1); y = xyz(2); z = xyz(3);
num = [-2,1,1,0,3,0];
alpha = num(1); beta = num(2); gamma = num(3);
d1 = num(4); d2 = num(5); d3 = num(6);
vf = [x*(alpha-alpha/4*x^2-(alpha+beta+d1)/4*y^2+d2/4*z^2);
    y*(beta+d1/4*x^2-beta/4*y^2-(beta+gamma+d3)/4*z^2);
    z*(gamma-(gamma+alpha+d2)/4*x^2+d3/4*y^2-gamma/4*z^2)];
end