clc; 
clear; 
close all;

% Input parameter
[infileextn] = uigetfile('*.xlsx','Select File InputInversion.xls...');
itin=xlsread(infileextn);

% Input model sintetik 
[infileextn] = uigetfile('*.dat','Select grav_obs...');
grav=load(infileextn);
f=grav(:,2);
f_space=grav(:,1);
f_space=5:5:85;


for i = 1:length(f)
    if i == length(f)
        break
    end
    aa(i) = (f(i+1)+f(i))/2;
end


for i =1:length(f)+length(aa)
    if mod(i,2)==1
        bb(i)=[f((i/2)+0.5)];
    else
        bb(i)=[aa(i/2)];
    end
end
f=[bb];





% Parameter Model Blok 2D
nlay=itin(1,1);             % Banyak grid
nxg=itin(1,2);            % Banyak cell lateral
nzg=itin(1,3);            % Banyak cell vetikal
dx =itin(1,4);              % Dimensi cell lateral (m)
dh =itin(1,5);              % Dimensi cell vetikal (m)
model=[nxg nzg dx dh];
l(1,1:nlay)=itin(2,1:nlay); % Batas bawah pencarian
u(1,1:nlay)=itin(3,1:nlay); % Batas atas pencarian
l=repelem(l,4);
u=repelem(u,4);
x0=l+(u-l)./2;              % Tebakan model awal
% x0=repelem(x0,2);
kmax=1000;                  % Jumlah Iterasi
q=1;                        % Faktor quenching
TolFun=1e-8;                % Toleransi relatif

% Inversi global simulated annealing
[xx,ff,dff]=sim_anl(f,x0,l,u,kmax,q,TolFun,model);
e_min=min(dff);
ff_min=find(dff==e_min);
ff_sol=ff(:,ff_min);
xsol=reshape(xx',nxg*2*nzg*2,kmax);
xx_sol=xsol(:,ff_min);
VV2=reshape(xx_sol,nxg*2,nzg*2);
VV=VV2';
zSA1 = [0:dh:(nzg-1)*dh]';
zSA = zSA1+10/2;

% for i = 1:10
%     [a,b,c]=sim_anl(f,x0,l,u,kmax,q,TolFun,model);
%     e_min=min(c)
%     ff_min=find(c==e_min);
%     ff_sol=b(:,ff_min);
%     xsol=reshape(a',nlay,kmax);
%     xx_sol=xsol(:,ff_min);
%     VV2=reshape(xx_sol,nxg,nzg);
%     VV=VV2';
%     zSA1 = [0:dh:(nzg-1)*dh]';
%     zSA = zSA1+10/2;
% end


% Hasil inversi global simulated annealing
figure(1)
subplot(2,1,1); 
plot(f_space,f,'b-o','LineWidth',2.0,'MarkerSize',3);hold on
subplot(2,1,1); 
plot(f_space,ff_sol([1:length(ff_sol)-1]),'r-o','LineWidth',2.0,'MarkerSize',3);
title('Penampan Hasil Inversi SA 2D Gravitasi','fontweight','bold','fontsize',8)
ylabel('Anomali Medan Gravitasi [mGal]','fontsize',7);
xlabel('Spasi [m]','fontsize',7);
% ylim([0.05,0.6])
% yticks([0.05,0.1:.1:0.6])
% yticks([0:0.8])
set(gca,'fontsize',7);
hleg= legend('G-obs','Inversi SA','Location','NorthEast');
set(hleg,'FontAngle','italic','FontSize',6);
subplot(2,1,2); imagesc(f_space,zSA,VV);
set(gca,'XAxisLocation','top','fontsize',7,'XMinorTick','on');
ylabel('Kedalaman [m]','fontsize',7);
colorbar('horiz');
colormap('default');
grid on
set(gca, 'GridLineStyle', '-');
set(gca,'XGrid','on','YGrid','on','ZGrid','on')
grid(gca,'minor');
writematrix(VV,"model.xlsx");
