load('/Volumes/Shield/CMAQ-DDM/data/LCoeff.mat');
load('/Volumes/Shield/CMAQ-DDM/data/QCoeff.mat');

N = [0.01, 1, 3.9, 0.1, 1, 2, 0.1, 3, 1.1, 0.9, 0.9, 1.1, 2.3, 1.3, 0.8];
V = [0.01, 1, 4.9, 1, 0.1, 0.1, 4, 0.1, 1.1, 0.9, 1, 1, 3.6, 1.5, 0.6];

for siteNum = 1: 1
	% find D
	sheets = {'1percent', 'base', '1985', 'n10v100', 'n100v10', 'n200v10', 'n10v400', 'n300v10', 'n110v110', 'n90v90', 'n90v100', 'n110v100', '2001'};
    LC = LCoeff(siteNum, :);
    QC = QCoeff(siteNum, :);
    figure(1)
    n = [0: 0.01: 3.4];
	v = [0: 0.01: 4.9];
	[NOx,VOC] = meshgrid(n,v);
    LO3 = exp(LC(1) + LC(2) * NOx + LC(3) * VOC + LC(4) * NOx.^2 + LC(5) * NOx .* VOC + LC(6) * VOC.^2);
    QO3 = QC(1) + QC(2) * NOx + QC(3) * VOC + QC(4) * NOx.^2 + QC(5) * NOx .* VOC + QC(6) * VOC.^2;
    LO3_CMAQ = exp(LC(1) + LC(2) * N + LC(3) * V + LC(4) * N.^2 + LC(5) * N .* V + LC(6) * V.^2);
    QO3_CMAQ = QC(1) + QC(2) * N + QC(3) * V + QC(4) * N.^2 + QC(5) * N .* V + QC(6) * V.^2;
    Diff = LO3 - QO3;
    contourf(VOC, NOx, Diff), hold on
    scatter(V, N, 'filled', 'r')
    colorbar
    hold off

    figure(2)
    n = [0: 0.1: 3.4];
	v = [0: 0.1: 4.9];
    [NOx,VOC] = meshgrid(n,v);
    LO3 = exp(LC(1) + LC(2) * NOx + LC(3) * VOC + LC(4) * NOx.^2 + LC(5) * NOx .* VOC + LC(6) * VOC.^2);
    LO3_OBV = exp(LC(1) + LC(2) * N + LC(3) * V + LC(4) * N.^2 + LC(5) * N .* V + LC(6) * V.^2);
    QO3 = QC(1) + QC(2) * NOx + QC(3) * VOC + QC(4) * NOx.^2 + QC(5) * NOx .* VOC + QC(6) * VOC.^2;
    QO3_OBV = QC(1) + QC(2) * N + QC(3) * V + QC(4) * N.^2 + QC(5) * N .* V + QC(6) * V.^2;
    surf(VOC, NOx, LO3), hold on,
    surf(VOC, NOx, QO3)
    scatter3(V, N, LO3_OBV, 'filled', 'r')
    scatter3(V, N, QO3_OBV, 'filled', 'r')
    % for k = 1:13
    %     line(N(k), V(k), [LO3_CMAQ(k):QO3_CMAQ(k)]);
    % end
    % view(3)
    hold off
end