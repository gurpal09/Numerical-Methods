
% Stability Region for Forward Euler
rho_FE = @(zeta) zeta - 1;
sigma_FE = @(zeta) 1;
z_FE = @(rho_AB4,sigma_AB4) rho_AB4./sigma_AB4;

theta = 0:.01:2*pi;
zeta = exp(i*theta);

rho_out = rho_FE(zeta);
sigma_out = sigma_FE(zeta);

z_out = z_FE(rho_out,sigma_out);

figure
whitebg('white')
plot(z_out)
fill(real(z_out),imag(z_out),'r')
hold on

title('Stability Region for Forward Euler')
xlabel('Re(z)')
ylabel('Im(z)')
grid
axis([-3 1 -2 2])


% Stability Region for Backward Euler
rho_BE = @(zeta) zeta - 1;
sigma_BE = @(zeta) -zeta;
z_BE = @(rho_AB4,sigma_AB4) rho_AB4./sigma_AB4;

rho_out = rho_BE(zeta);
sigma_out = sigma_BE(zeta);

z_out = z_BE(rho_out,sigma_out);


figure
rectangle('Position',[-3 -2 4 4],'FaceColor','red')
hold on
plot(z_out)
fill(real(z_out),imag(z_out),'white')
hold on

title('Stability Region for Backward Euler')
xlabel('Re(z)')
ylabel('Im(z)')
grid on
set(gca,'layer','top');
axis([-3 1 -2 2])



