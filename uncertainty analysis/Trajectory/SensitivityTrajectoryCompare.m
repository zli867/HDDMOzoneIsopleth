% regression model
T = readtable('/Users/zongrunli/Desktop/CMAQ-DDM/data/Observations.xlsx');
load('/Users/zongrunli/Desktop/CMAQ-DDM/data/LCoeff.mat');
load('/Users/zongrunli/Desktop/CMAQ-DDM/data/QCoeff.mat');
siteName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Rubidoux', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernardino'};
siteTableName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Riverside', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernadino'};
baseNOx = 395.748;
baseVOC = 401.412;
for siteNum = 1: 1
    y = [];
    X = [];
    NOx = [];
    VOC = [];
    year = [];
    for i = 1: size(T, 1)
        temp = T(i, :);
        if strcmp(temp.SiteName, siteTableName{siteNum})
            NOxtemp = temp.NOx/baseNOx;
            VOCtemp = temp.VOC/baseVOC;
            ODVTemp = temp.ODV;
            if temp.Year ~= 1111 & temp.Year ~= 2222 & temp.Year ~= 3333
                NOx = [NOx, NOxtemp];
                VOC = [VOC, VOCtemp];
                year = [year, temp.Year];
            end
            % NOx = [NOx, NOxtemp];
            % VOC = [VOC, VOCtemp];
            y = [y; ODVTemp];
            x = [1, NOxtemp, VOCtemp, NOxtemp^2, VOCtemp* NOxtemp, VOCtemp^2];
            X = [X; x];
        end
    end
    % pseudo inverse
    b = pinv(X) * y;

    %plot
	% regression sensitive
    % dO3/dNOx = 0
    regress_NOx = b(2) + 2* b(4) .* NOx + b(5) * VOC;
    % dO3/dVOC
    regress_VOC = b(3) + b(5) .* NOx + 2*b(6) .* VOC;

    % Quadratic sensitivity
    % dO3/dNOx
    Q_NOx = QCoeff(siteNum, 2) + 2* QCoeff(siteNum, 4) .* NOx + QCoeff(siteNum, 5) .* VOC;
    % dO3/dVOC
    Q_VOC = QCoeff(siteNum, 3) + QCoeff(siteNum, 5) .* NOx + 2* QCoeff(siteNum, 6) .* VOC;

    % Log-Quadratic sensitivity
    % dO3/dNOx
    LO3 = exp(LCoeff(siteNum, 1) + LCoeff(siteNum, 2) * NOx + LCoeff(siteNum, 3) * VOC + LCoeff(siteNum, 4) * NOx.^2 + LCoeff(siteNum, 5) * NOx .* VOC + LCoeff(siteNum, 6) * VOC.^2);
    L_NOx = (LCoeff(siteNum, 2) + 2* LCoeff(siteNum, 4) .* NOx + LCoeff(siteNum, 5) .* VOC) .* LO3;
    % dO3/dVOC
    L_VOC = (LCoeff(siteNum, 3) + 2* LCoeff(siteNum, 6) .* VOC + LCoeff(siteNum, 5) .* NOx) .* LO3;

    
    % plot dO3/dNOx
    figure(siteNum)
    subplot(1, 2, 1)
    maxValue = ceil(max(max(max(regress_NOx), max(Q_NOx)), max(L_NOx)));
    minValue = floor(min(min(min(regress_NOx), min(Q_NOx)), min(L_NOx)));
    scatter(regress_NOx, Q_NOx,'filled'), hold on
    scatter(regress_NOx, L_NOx,'filled'),
    % h = colorbar;
    % set(get(h, 'label'), 'string', 'year');
    plot([minValue: maxValue], [minValue: maxValue], '--'),
    xlim([minValue, maxValue]), ylim([minValue, maxValue]),
    fit1 = fitlm(regress_NOx, Q_NOx);
    fit2 = fitlm(regress_NOx, L_NOx);
    Rsquare1 = fit1.Rsquared.Ordinary;
    Rsquare2 = fit2.Rsquared.Ordinary;
    txt1 = strcat('R^2 (Quadratic): ', num2str(Rsquare1, '%.3f'));
    txt2 = strcat('R^2 (Log-Quadratic): ', num2str(Rsquare2, '%.3f'));
    txt = {txt1, txt2};
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
    text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left'),
    xlabel('Empirical Sensitivity'), ylabel('Isopleth Sensitivity'),
    title('\partial O_3/\partial NO_x')
    legend('Quadratic', 'Log-Quadratic', 'Location','southeast'), 
    axis square;
    hold off;

    % plot dO3/dVOC
    subplot(1, 2, 2)
    maxValue = ceil(max(max(max(regress_VOC), max(Q_VOC)), max(L_VOC)));
    minValue = floor(min(min(min(regress_VOC), min(Q_VOC)), min(L_VOC)));
    scatter(regress_VOC, Q_VOC, 'filled'), hold on
    scatter(regress_VOC, L_VOC,  'filled'),
    % h = colorbar;
    % set(get(h, 'label'), 'string', 'year');
    plot([minValue: maxValue], [minValue: maxValue], '--'),
    xlim([minValue, maxValue]), ylim([minValue, maxValue]),
    fit1 = fitlm(regress_VOC, Q_VOC);
    fit2 = fitlm(regress_VOC, L_VOC);
    Rsquare1 = fit1.Rsquared.Ordinary;
    Rsquare2 = fit2.Rsquared.Ordinary;
    txt1 = strcat('R^2 (Quadratic): ', num2str(Rsquare1, '%.3f'));
    txt2 = strcat('R^2 (Log-Quadratic): ', num2str(Rsquare2, '%.3f'));
    txt = {txt1, txt2};
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
    text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left'),
    xlabel('Empirical Sensitivity'), ylabel('Isopleth Sensitivity'),
    title('\partial O_3/\partial VOC')
    legend('Quadratic', 'Log-Quadratic', 'Location','southeast')
    axis square;
    hold off;

    sgtitle(siteName{siteNum})

    % save figures
    saveas(gcf, strcat('/Users/zongrunli/Desktop/CMAQ-DDM/figures/', siteName{siteNum}, ' senTraj.jpg'))
end