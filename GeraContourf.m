Valpha = [15,25,35,45,55];

VRe = [100,150,200,250,300];

[MRe,Malpha] = meshgrid(VRe,Valpha);

K = real(MatrizLambda(:,[1:5]));

%K = flipud(K);

figure(100)
contourf(MRe,Malpha,K,[0 1]);
c.LineWidth = 3;
colormap gray;
title("Curva de estabilidade neutra")
xlabel('\Re')
ylabel('\alpha')


