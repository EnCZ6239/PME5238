%% LINEAR Stability ANALYSIS of a constricted channel
% (January 2023, StabFem 3.8 ) 

%% First we set a few global variables needed by StabFem

%{%
tic
close all;
addpath('../../SOURCES_MATLAB');
SF_Start('verbosity',4);
SF_DataBase('create','./WORK/');
%}

%% ##### COMPUTING THE MESH 

%{%
ffmesh = SF_Mesh('Mesh_Channel.edp','Options',{'Xinlet',-6,'Xoutlet',15,'d',0.5,'alpha',15,'n',20});
%SF_Plot(ffmesh,'boundary','on','bdlabels',2,'bdcolors','k');
%{
SF_Plot(ffmesh);
title('Domínio inteiro e malha inicial');

figure();
SF_Plot(ffmesh,'xlim',[-2 6],'ylim',[-0.5 0.5]);
title('Malha inicial próxima da constrição');
%set(gcf,'position',[80,100,1900,250]) %%Plot do mesh gerado

figure();
SF_Plot(ffmesh,'xlim',[-2 6],'ylim',[-0.5 0.5]); %%Plot com zoom do mesh gerado
%}

%% ##### Computing base flow at Re = 1


bf=SF_BaseFlow(ffmesh,'solver','Newton_2D.edp','Re',1);

%{%
figure();
SF_Plot(bf,'ux','xlim',[-2 6],'ylim',[-0.5 0.5],'colormap','redblue');
hold on;

SF_Plot(bf,{'ux','uy'},'xlim',[-2 6],'ylim',[-0.5 0.5]);
set(gcf,'position',[80,200,1900,250]) %%Plot do cambo base para Re = 1, ux e streamlines
hold off;
%}


%% #### Computing base flow and adapting mesh for increasing Re until desired 

%{%
VRe = [5,50,100];

for i = 1:length(VRe)
    fprintf('---------------------------------------------------- RESOLVENDO PARA RE = %d \n',VRe(i));
    Mask = SF_Mask(bf,'rectangle',[-1, 5 -0.5 0.5 0.0125 ]) ;
    bf = SF_Adapt(bf, Mask);
    bf = SF_BaseFlow(bf,'Re',VRe(i));
end
%{
figure();
SF_Plot(bf,'ux','xlim',[-2 6],'ylim',[-0.5 0.5],'colormap','redblue','colorrange',[-1 2],'boundary','on','bdlabels',2,'bdcolors','k');
set(gcf,'position',[80,200,1900,250]) %%Plot do cambo base, pressão e streamlines
hold on;
SF_Plot(bf,{'ux','uy'},'xlim',[-2 6],'ylim',[-0.5 0.5]);
%set(gca,'DataAspectRatio',[2 1 1])
hold off;

title("Escoamento base - campo u_x e perfis de velocidade para Re = "+VRe(i));
%}
BaseTime = toc;



%% #### Computing first eigenmode, adapting mesh, computing nev first eigenmodes

tic

[ev1,em1] = SF_Stability(bf,'solver','Stab_2D.edp','shift',5,'nev',1,'sort','LR','Options',{'type','D'});

%{%
Mask = SF_Mask(bf,'rectangle',[-1, 5 -0.5 0.5 0.0125 ]) ;
bf = SF_Adapt(bf,em1,Mask);
bf = SF_BaseFlow(bf,'Re',VRe(i));
[ev,em] = SF_Stability(bf,'solver','Stab_2D.edp','shift',5,'nev',1,'sort','LR','Options',{'type','D'});
%}
%{
figure();
SF_Plot(bf.mesh,'xlim',[-2 6],'ylim',[-0.5 0.5]); %%Plot com zoom do mesh gerado
title('Malha adaptada para escoamento base e primeiro modo');
%}

[Gi,nGi] = max(abs(imag(ev)));
%{

figure();
plot(ev,'ko')

ylim([-0.05, Gi+0.5])
grid on;
title("Autovalores para Re = "+VRe(i));
%}

%% #### Ploting some eigenmodes
%{
%First mode

%em = SF_LoadFields(em);

figure();
SF_Plot(em(1),'ux.re','xlim',[-1 9],'ylim',[-0.5 0.5],'colormap','redblue','colorrange','cropcentered',...
    'boundary','on','bdlabels',2,'bdcolors','k');
box on;
set(gca,'FontSize', 14); 
title("Primeiro modo, componente ux, Re = "+ VRe(i));

figure();
SF_Plot(em(1),'uy.re','xlim',[-1 9],'ylim',[-0.5 0.5],'colormap','redblue','colorrange','cropcentered',...
    'boundary','on','bdlabels',2,'bdcolors','k');
box on;
set(gca,'FontSize', 14); 
title("Primeiro modo, componente uy, Re = "+ VRe(i));

figure();
SF_Plot(em(1),'p.re','xlim',[-1 9],'ylim',[-0.5 0.5],'colormap','redblue','colorrange','cropcentered',...
    'boundary','on','bdlabels',2,'bdcolors','k');
box on;
set(gca,'FontSize', 14); 
title("Primeiro modo, pressão , Re = "+ VRe(i));

%{%
%Second mode
figure();
SF_Plot(em(2),'ux.re','xlim',[-1 9],'ylim',[-0.5 0.5],'colormap','redblue','colorrange','cropcentered',...
    'boundary','on','bdlabels',2,'bdcolors','k');
box on;
set(gca,'FontSize', 14); 
title("Segundo modo, componente ux, Re = "+ VRe(i));

figure();
SF_Plot(em(2),'uy.re','xlim',[-1 9],'ylim',[-0.5 0.5],'colormap','redblue','colorrange','cropcentered',...
    'boundary','on','bdlabels',2,'bdcolors','k');
box on;
set(gca,'FontSize', 14); 
title("Segundo modo, componente uy, Re = "+ VRe(i));

figure();
SF_Plot(em(2),'p.re','xlim',[-1 9],'ylim',[-0.5 0.5],'colormap','redblue','colorrange','cropcentered',...
    'boundary','on','bdlabels',2,'bdcolors','k');
box on;
set(gca,'FontSize', 14); 
title("Segundo modo, pressão , Re = "+ VRe(i));

%Oscilating mode
figure();
SF_Plot(em(nGi),'ux.re','xlim',[-1 9],'ylim',[-0.5 0.5],'colormap','redblue','colorrange','cropcentered',...
    'boundary','on','bdlabels',2,'bdcolors','k');
box on;
set(gca,'FontSize', 14); 
title("Modo oscilatório, componente ux, Re = "+ VRe(i));

figure();
SF_Plot(em(nGi),'uy.re','xlim',[-1 9],'ylim',[-0.5 0.5],'colormap','redblue','colorrange','cropcentered',...
    'boundary','on','bdlabels',2,'bdcolors','k');
box on;
set(gca,'FontSize', 14); 
title("Modo oscilatório, componente uy, Re = "+ VRe(i));

figure();
SF_Plot(em(nGi),'p.re','xlim',[-1 9],'ylim',[-0.5 0.5],'colormap','redblue','colorrange','cropcentered',...
    'boundary','on','bdlabels',2,'bdcolors','k');
box on;
set(gca,'FontSize', 14); 
title("Modo oscilatório, pressão , Re = "+ VRe(i));
%}
EigTime = toc;


%% #### For a specific mode, computing adjoint eigenmode and structural sensitivity

[evA,emA] = SF_Stability(bf,'solver','Stab_2D.edp','shift',5,'nev',1,'sort','LR','Options',{'type','A'});

emS = SF_Sensitivity(bf,em(1),emA,'solver','Sensitivity2D.edp');

%emA = SF_LoadFields(emA);
%emS = SF_LoadFields(emS);

%{
figure();
SF_Plot(emA,'ux.re','xlim',[-1 9],'ylim',[-0.5 0.5],'colormap','redblue','colorrange','cropcentered',...
    'boundary','on','bdlabels',2,'bdcolors','k');
box on;
set(gca,'FontSize', 14); 
title("Adjunto do primeiro modo, componente ux, Re = "+ VRe(i));

figure();
SF_Plot(emA,'uy.re','xlim',[-1 9],'ylim',[-0.5 0.5],'colormap','redblue','colorrange','cropcentered',...
    'boundary','on','bdlabels',2,'bdcolors','k');
box on;
set(gca,'FontSize', 14); 
title("Adjunto do primeiro modo, componente uy, Re = "+ VRe(i));

figure();
SF_Plot(emA,'p.re','xlim',[-1 9],'ylim',[-0.5 0.5],'colormap','redblue','colorrange','cropcentered',...
    'boundary','on','bdlabels',2,'bdcolors','k');
box on;
set(gca,'FontSize', 14); 
title("Adjunto do primeiro modo, pressão , Re = "+ VRe(i));
%}

%{
figure();
SF_Plot(emS,'S.re','xlim',[-0.5 6],'ylim',[-0.5 0.5],'colormap','ice','boundary','on','bdlabels',2,'bdcolors','k');
hold on; 
SF_Plot(bf,{'ux','uy'},'xlim',[-0.5 6],'ylim',[-0.5 0.5]);

set(gca,'FontSize', 14); 
title("Sensibilidade estrutural do primeiro modo, Re = "+ VRe(i));

%}

%{
figure();
SF_Plot(emA(1),'p','xlim',[-2 6],'ylim',[-0.5 0.5],'colormap','redblue','colorrange','cropcentered','boundary','on','bdlabels',2,'bdcolors','k');
set(gcf,'position',[80,200,1900,250]) %%Plot do cambo base, pressão e streamlines
hold on;
SF_Plot(emA(1),{'ux','uy'},'xlim',[-2 6],'ylim',[-0.5 0.5]);
%set(gca,'DataAspectRatio',[2 1 1])
hold off;
%}




