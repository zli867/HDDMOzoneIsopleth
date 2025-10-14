T = readtable('/Users/zongrunli/Desktop/CMAQ-DDM/data/Observations.xlsx');
load('/Users/zongrunli/Desktop/CMAQ-DDM/data/LCoeff.mat');
load('/Users/zongrunli/Desktop/CMAQ-DDM/data/QCoeff.mat');
siteName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Rubidoux', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernardino'};
siteTableName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Riverside', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernadino'};
baseNOx = 557.384;
baseVOC = 521.582;
for siteNum = 1: 1
    QO3 = [];
    LO3 = [];
    ODV = [];
    for i = 1: size(T, 1)
        temp = T(i, :);
        if temp.Year ~= 1111 & temp.Year ~= 2222 & temp.Year ~= 3333 & strcmp(temp.SiteName, siteTableName{siteNum})
            NOx = temp.NOx/baseNOx;
            VOC = temp.VOC/baseVOC;
            QO3Temp = QCoeff(siteNum, 1) + QCoeff(siteNum, 2) * NOx + QCoeff(siteNum, 3) * VOC + QCoeff(siteNum, 4) * NOx.^2 + QCoeff(siteNum, 5) * NOx .* VOC + QCoeff(siteNum, 6) * VOC.^2;
            LO3Temp = exp(LCoeff(siteNum, 1) + LCoeff(siteNum, 2) * NOx + LCoeff(siteNum, 3) * VOC + LCoeff(siteNum, 4) * NOx.^2 + LCoeff(siteNum, 5) * NOx .* VOC + LCoeff(siteNum, 6) * VOC.^2);
            ODVTemp = temp.ODV;
            QO3 = [QO3, QO3Temp];
            LO3 = [LO3, LO3Temp];
            ODV = [ODV, ODVTemp];
        end
    end
    % plot
    figure(siteNum)
    maxValue = ceil(max(max(max(ODV), max(LO3)), max(QO3))) + 1;
    scatter(ODV, QO3, 'fill', 'r'), hold on, 
    scatter(ODV, LO3, 'fill', 'b'), xlabel('Observed Ozone (ppb)'), ylabel('Isopleth Estimated Ozone (ppb)')
    plot([0: maxValue], [0: maxValue], '--'), hold on,
    xlim([0, maxValue]), ylim([0, maxValue]);
    title(siteName{siteNum});
    fit = fitlm(ODV, QO3);
    RSquare = fit.Rsquared.Ordinary;
    text1 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    fit2 = fitlm(ODV, LO3);
    RSquare = fit2.Rsquared.Ordinary;
    text2 = strcat('R^2 (Log-Quadratic): ', num2str(RSquare, '%.3f'));
    txt = {text1, text2};
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
    text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
    legend('Quadratic', 'Log-Quadratic', 'Location','southeast')
    axis square;
    hold off;
    % save figures
    saveas(gcf, strcat('/Users/zongrunli/Desktop/CMAQ-DDM/figures/', siteName{siteNum}, ' obv.jpg'))
end