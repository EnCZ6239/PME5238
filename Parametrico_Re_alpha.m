%% LINEAR Stability ANALYSIS of a constricted channel, parametric analysis varying Re and alpha

% (January 2023, StabFem 3.8 ) 
%

%%
%
% First we set a few global variables needed by StabFem
%
%{%
tic
close all;
addpath('../../SOURCES_MATLAB');
SF_Start('verbosity',4);
SF_DataBase('create','./WORK/');


%% ##### CHAPTER 1 : COMPUTING THE MESH 
% 

%{%
Valpha = [15,25,35,45,55];

VRe = [100,150,200,250,300];

Maskadapt = 0.01;

%}
tic
%{
Valpha = [10,60];

VRe = [200,300];
%}
VReinit = [5,50,100];

for a=1:length(Valpha)
    
ffmesh = SF_Mesh('Mesh_Channel.edp','Options',{'Xinlet',-6,'Xoutlet',14,'d',0.5,'alpha',Valpha(a),'n',18});
clear bf;
bf=SF_BaseFlow(ffmesh,'solver','Newton_2D.edp','Re',1,'ncores',1);



for i = 1:length(VReinit)
    fprintf('---------------------------------------------------- RESOLVENDO PARA RE = %d \n',VReinit(i));
    Mask = SF_Mask(bf,'rectangle',[-1, 5 -0.5 0.5 Maskadapt]) ;
    bf = SF_Adapt(bf, Mask);
    bf = SF_BaseFlow(bf,'Re',VReinit(i));
end


figure();
SF_Plot(bf,'ux','xlim',[-1 9],'ylim',[-0.55 0.55],'colormap','redblue','colorrange',[-0.5 2],'boundary','on','bdlabels',2,'bdcolors','k');
set(gcf,'position',[80,200,1900,250]) %%Plot do cambo base para Re = 250, pressão e streamlines
hold on;
SF_Plot(bf,{'ux','uy'},'xlim',[-2 6],'ylim',[-0.5 0.5]);
hold off;
title("Escoamento base - campo u_x e perfis de velocidade para Re = "+VReinit(i));


%[~,em] = SF_Stability(bf,'solver','Stab_2D.edp','shift',5,'nev',1,'sort','LR','Options',{'type','D'});
%Mask = SF_Mask(bf,'rectangle',[-1, 5 -0.5 0.5 Maskadapt]) ;
%bf = SF_Adapt(bf,em,Mask);

%bf = SF_BaseFlow(bf,'Re',VReinit(i));
[ev,~] = SF_Stability(bf,'solver','Stab_2D.edp','shift',5,'nev',1,'sort','LR','Options',{'type','D'});

MatrizLambda(a,1) = ev;

for j = 2:length(VRe)
    fprintf('---------------------------------------------------- RESOLVENDO PARA RE = %d \n',VRe(j));
    %Mask = SF_Mask(bf,'rectangle',[-1, 5 -0.5 0.5 Maskadapt]) ;
    %bf = SF_Adapt(bf, Mask);
    bf = SF_BaseFlow(bf,'Re',VRe(j));
    
    [ev,em] = SF_Stability(bf,'solver','Stab_2D.edp','shift',5,'nev',1,'sort','LR','Options',{'type','D'});
    MatrizLambda(a,j) = ev;
    
    figure();
    SF_Plot(em,'ux','xlim',[-1 9],'ylim',[-0.5 0.5],'colormap','redblue','colorrange','cropcentered',...
    'boundary','on','bdlabels',2,'bdcolors','k');
    box on;
    set(gca,'FontSize', 14); 
    title("Primeiro modo, componente ux, Re = "+ VRe(j));
    
end

end
TimeTeste = toc;





