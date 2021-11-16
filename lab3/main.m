%% Physikalische Geodäsie 1
%% Exercise 2 - Gravitation
% Autoren: Jens Kappler, Hannes Nübel
clc 
clear
close all

k = 1;         
load PREM
G = 6.672e-11;

m_E = 5.9736e24;
R_E = 6371e3;

m_M =  7.349e22;
R_M = 1738e3;
r_EM = 384400e3;
X_M = [r_EM*cosd(k*10); r_EM*sind(k*10); 0];

%% Superposition of Earth's and Moon's gravitational fields
% Task 1-4:

x_V = linspace(-4e7,4e8,500);
y_V = linspace(-4e7,2e8,500);
x_a = linspace(-4e7,4e8,30);
y_a = linspace(-4e7,2e8,30);

[X_V,Y_V] = meshgrid(x_V,y_V);
[X_a,Y_a] = meshgrid(x_a,y_a);

V_E = potential(X_V,Y_V,0,0,m_E);
V_M = potential(X_V,Y_V,X_M(1),X_M(2),m_M);

figure; hold on;
% rectangle('Position',[-R_E -R_E 2*R_E 2*R_E],'Curvature',[1 1],'FaceColor',[.5 .5 0])
% rectangle('Position',[X_M(1)-R_M X_M(2)-R_M 2*R_M 2*R_M],'Curvature',[1 1],'FaceColor',[.5 .5 0])

labels = floor(logspace(0,9));
[C,h] = contour(X_V,Y_V,V_E+V_M,labels); axis equal;
clabel(C,h,labels(2:2:end))
% 
% [ax_E,ay_E] = a_sphere(X_a,Y_a,R_E,0,0,rho_E,G);
% [ax_M,ay_M] = a_sphere(X_a,Y_a,R_M,X_M(1),X_M(2),rho_M,G);
% 
% ax_sup = ax_E + ax_M;
% ay_sup = ay_E + ay_M;
% a_norm = sqrt(ax_sup.^2+ay_sup.^2);
% quiver(X_a, Y_a, ax_sup./a_norm, ay_sup./a_norm, 0.5); axis equal;print('-dpng','grav');

% % Gravitational potential and attraction of spherical shells
% Task 7:
% R_C = 3500e3;
% R_M = 6400e3;
% rho_C = 11200;
% rho_M = 4300;
% 
% r = linspace(0,4*R_M,1000);
% 
% V_C = V_sphere(r,zeros(size(r)),0,0,R_C,rho_C,G);
% V_M = V_shell(R_C,R_M,rho_M,r,G);
% V = V_C + V_M;
% 
% a_C= a_sphere_scal(R_C,rho_C,r,G);
% a_M = a_shell(R_C,R_M,rho_M,r,G);
% a = a_C + a_M;
% 
% figure
% subplot(2,1,1);
% plot(r/1000,V);grid on;hold on;
% xlabel('r in [km]');
% ylabel('V in [m?/s²]');
% title('Potential of the Earth [m?/s²]');
% xlim([0 4*R_M*1e-3]);
% 
% subplot(2,1,2);
% plot(r/1000,a);grid on;hold on;
% xlabel('r in [km]');
% ylabel('a in [m/s²]');
% title('Attraction of the Earth [m/s²]');
% xlim([0 4*R_M*1e-3]);print('-dpng','shells');
% 
% % PREM density model of the Earth
% figure; plot(PREM(:,1),PREM(:,2)); title('PREM Modell \rho(r'')'); grid on
% xlabel('Radial coordinate r''[km]');
% ylabel('Density \rho[kg/m^3]');print('-dpng','prem');
% 
% Task 8
% steps = size(PREM,1);
% r = linspace(0,2*max(PREM(:,1))*1e3,steps)';
% 
% a_P = zeros(steps,1);
% V_P = zeros(steps,1);
% 
% for i = 1:steps-1
%     V_P = V_P + V_shell(PREM(i,1)*1e3,PREM(i+1,1)*1e3,PREM(i,2),r,G);
%     a_P = a_P + a_shell(PREM(i,1)*1e3,PREM(i+1,1)*1e3,PREM(i,2),r,G);
% end
% 
% figure
% subplot(2,1,1);
% plot(r/1000,V_P);grid on;hold on;
% xlabel('r in [km]');
% ylabel('V in [m?/s²]');
% title('Potential of the Earth [m?/s²]');
% xlim([0 max(r)*1e-3]);
% 
% subplot(2,1,2);
% plot(r/1000,a_P);grid on;hold on;
% xlabel('r in [km]');
% ylabel('a in [m/s²]');
% title('Attraction of the Earth [m/s²]');
% xlim([0 max(r)*1e-3]);
% 
% Task 9
% V_surf = V_P(end/2-0.5);
% a_surf = a_P(end/2-0.5);
% 
% subplot(2,1,1); s1=scatter(r(end/2-0.5)/1000,V_surf,'filled'); 
% subplot(2,1,2); s2=scatter(r(end/2-0.5)/1000,a_surf,'filled'); 
% 
% Task 10
% [V_max,V_maxi] = max(V_P);
% r_V_maxi = r(V_maxi);
% [a_max,a_maxi] = max(abs(a_P));
% r_a_maxi = r(a_maxi);
% 
% subplot(2,1,1); s3=scatter(r_V_maxi/1000,V_max,'filled');
% legend([s1,s3],'surface','maximum');
% subplot(2,1,2); s4=scatter(r_a_maxi/1000,-a_max,'filled');
% legend([s2,s4],'surface','maximum','Location','southeast');
% print('-dpng','shells_prem')

%% Functions
function [V] = V_sphere(X,Y,R,X_m,Y_m,rho,G)
    r = sqrt((X - X_m).^2 + (Y - Y_m).^2);
    
    V = 4/3*pi * G * rho * R^3 .* 1./(r+eps) .*(r>R);
    V = V + 2*pi * G * rho .* (R^2 - r.^2./3) .*(r<=R);
end

function [a_x,a_y] = a_sphere(X,Y,R,X_m,Y_m,rho,G)
    rx = X - X_m; 
    ry = Y - Y_m;
    r = sqrt((rx).^2 + (ry).^2);
    
    a_x = -G*4/3 * pi * rho * R^3./r.^3 .* rx .*(r>R);
    a_y = -G*4/3 * pi * rho * R^3./r.^3 .* ry .*(r>R);
    
    a_x = a_x - 4/3 * pi * G * rho .* rx .*(r<=R);
    a_y = a_y - 4/3 * pi * G * rho .* ry .*(r<=R);
end

function [a] = a_sphere_scal(R,rho,r,G)
    a_out = -(4/3*pi*G*rho) * R.^3./(r.^2+eps) .* (r>=R);
    a_in = -(4/3*pi*G*rho) .* r .* (r<R);
    
    a = a_out + a_in;
end

function [V] = V_shell(R_in,R_out,rho,r,G)
    V_in = (2*pi*G*rho) * (R_out^2 - R_in^2) .*(r<=R_in);
    V_sh = ((2*pi*G*rho) * (R_out^2 - r.^2./3) - 4/3*pi*G*rho*R_in^3 ./(r+eps)) .*((r>R_in)&(r<R_out));
    V_out = (4/3*pi*G*rho) * (R_out^3 - R_in^3) ./(r+eps) .*(r>=R_out);
    
    V = V_in + V_sh + V_out;     
end

function [a] = a_shell(R_in,R_out,rho,r,G)
    a_out = -(4/3*pi*G*rho) * (R_out^3-R_in^3) ./(r.^2+eps) .*(r>=R_out);
    a_sh = -(4/3*pi*G*rho) * (r.^3-R_in^3) ./(r.^2+eps) .*((r>R_in) & (r<R_out));
    % a_in = 0;
    
    a = a_out + a_sh;
end

function V = potential(x,y,x1,y1,m)
G = 6.672e-11;
V = G*(m./sqrt((x-x1).^2+(y-y1).^2));
end