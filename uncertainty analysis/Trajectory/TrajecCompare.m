% regression model
% load('/Users/zongrunli/Desktop/CMAQ-DDM/data/LCoeffWeightedAll.mat');
% load('/Users/zongrunli/Desktop/CMAQ-DDM/data/QCoeffWeightedAll.mat');
% load('/Users/zongrunli/Desktop/CMAQ-DDM/data/LCoeffPart.mat');
% load('/Users/zongrunli/Desktop/CMAQ-DDM/data/QCoeffPart.mat');
load('/Users/zongrunli/Desktop/CMAQ-DDM/data/LCoeff.mat');
load('/Users/zongrunli/Desktop/CMAQ-DDM/data/QCoeff.mat');
T = readtable('/Users/zongrunli/Desktop/CMAQ-DDM/data/Observations.xlsx');
siteName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Rubidoux', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernardino'};
siteTableName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Riverside', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernadino'};
baseNOx = 395.748;
baseVOC = 401.412;
% save data
%siteList = [];
%ODV_L = [];
%S_N = [];
%S_V = [];

for siteNum = 1: 24
    y = [];
    X = [];
    NOx = [];
    VOC = [];
    year = [];
    QO3 = [];
    LO3 = [];
    ODV = [];
    for i = 1: size(T, 1)
        temp = T(i, :);
        if strcmp(temp.SiteName, siteTableName{siteNum})
            NOxtemp = temp.NOx/baseNOx;
            VOCtemp = temp.VOC/baseVOC;
            ODVTemp = temp.ODV;
            if temp.Year ~= 2222 & temp.Year ~= 3333 & temp.Year ~= 4444
                NOx = [NOx, NOxtemp];
                VOC = [VOC, VOCtemp];
                year = [year, temp.Year];
                ODVTemp = temp.ODV;
                ODV = [ODV, ODVTemp];
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
    % dO3/dNOx
    regress_NOx = b(2) + 2* b(4) .* NOx + b(5) * VOC;
    % dO3/dVOC
    regress_VOC = b(3) + b(5) .* NOx + 2*b(6) .* VOC;

    % Quadratic sensitivity
    % dO3/dNOx
    QO3 = QCoeff(siteNum, 1) + QCoeff(siteNum, 2) * NOx + QCoeff(siteNum, 3) * VOC + QCoeff(siteNum, 4) * NOx.^2 + QCoeff(siteNum, 5) * NOx .* VOC + QCoeff(siteNum, 6) * VOC.^2;
    Q_NOx = QCoeff(siteNum, 2) + 2* QCoeff(siteNum, 4) .* NOx + QCoeff(siteNum, 5) .* VOC;
    % dO3/dVOC
    Q_VOC = QCoeff(siteNum, 3) + QCoeff(siteNum, 5) .* NOx + 2* QCoeff(siteNum, 6) .* VOC;

    % Log-Quadratic sensitivity
    % dO3/dNOx
    LO3 = exp(LCoeff(siteNum, 1) + LCoeff(siteNum, 2) * NOx + LCoeff(siteNum, 3) * VOC + LCoeff(siteNum, 4) * NOx.^2 + LCoeff(siteNum, 5) * NOx .* VOC + LCoeff(siteNum, 6) * VOC.^2);
    L_NOx = (LCoeff(siteNum, 2) + 2* LCoeff(siteNum, 4) .* NOx + LCoeff(siteNum, 5) .* VOC) .* LO3;
    % dO3/dVOC
    L_VOC = (LCoeff(siteNum, 3) + 2* LCoeff(siteNum, 6) .* VOC + LCoeff(siteNum, 5) .* NOx) .* LO3;

    % Quadratic Model
    % plot O3
    figure(siteNum)
    subplot(2, 3, 1)
    maxValue = ceil(max(max(QO3), max(ODV)));
    minValue = floor(min(min(QO3), min(ODV)));
    scatter(ODV, QO3, 30, year, 'filled'), hold on
%     h = colorbar;
%     set(get(h, 'label'), 'string', 'year');
    plot([minValue: maxValue], [minValue: maxValue], '--'),
    xlim([minValue, maxValue]), ylim([minValue, maxValue]),
    fit = fitlm(ODV, QO3);
    Rsquare = fit.Rsquared.Ordinary;
    txt = strcat('R^2: ', num2str(Rsquare, '%.3f'));
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
    text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
    xlabel('Observation'), ylabel('Quadratic Model'),
    title('O_3'),
    axis square;
    hold off;
    % plot dO3/dNOx
    subplot(2, 3, 2)
    maxValue = ceil(max(max(regress_NOx), max(Q_NOx)));
    minValue = floor(min(min(regress_NOx), min(Q_NOx)));
    scatter(regress_NOx, Q_NOx, 30, year, 'filled'), hold on
%     h = colorbar;
%     set(get(h, 'label'), 'string', 'year');
    plot([minValue: maxValue], [minValue: maxValue], '--'),
    xlim([minValue, maxValue]), ylim([minValue, maxValue]),
    fit = fitlm(regress_NOx, Q_NOx);
    Rsquare = fit.Rsquared.Ordinary;
    txt = strcat('R^2: ', num2str(Rsquare, '%.3f'));
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
    text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
    xlabel('Empirical Sensitivity'), ylabel('Quadratic Model Sensitivity'),
    title('\partial O_3/\partial NO_x'),
    axis square;
    hold off;

    % plot dO3/dVOC
    subplot(2, 3, 3)
    maxValue = ceil(max(max(regress_VOC), max(Q_VOC)));
    minValue = floor(min(min(regress_VOC), min(Q_VOC)));
    scatter(regress_VOC, Q_VOC, 30, year, 'filled'), hold on
%     h = colorbar;
%     set(get(h, 'label'), 'string', 'year');
    plot([minValue: maxValue], [minValue: maxValue], '--'),
    xlim([minValue, maxValue]), ylim([minValue, maxValue]),
    fit = fitlm(regress_VOC, Q_VOC);
    Rsquare = fit.Rsquared.Ordinary;
    txt = strcat('R^2: ', num2str(Rsquare, '%.3f'));
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
    text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
    xlabel('Empirical Sensitivity'), ylabel('Quadratic Model Sensitivity'),
    title('\partial O_3/\partial VOC'),
    axis square;
    hold off;

    % plot Log-Quadratic VOC
    % ODV
    subplot(2, 3, 4)
    maxValue = ceil(max(max(LO3), max(ODV)));
    minValue = floor(min(min(LO3), min(ODV)));
    scatter(ODV, LO3, 30, year, 'filled'), hold on
%     h = colorbar;
%     set(get(h, 'label'), 'string', 'year');
    plot([minValue: maxValue], [minValue: maxValue], '--'),
    xlim([minValue, maxValue]), ylim([minValue, maxValue]),
    fit = fitlm(ODV, LO3);
    Rsquare = fit.Rsquared.Ordinary;
    txt = strcat('R^2: ', num2str(Rsquare, '%.3f'));
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
    text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
    xlabel('Observation'), ylabel('Log Quadratic Model'),
    title('O_3'),
    axis square;
    hold off;
    % plot dO3/dNOx
    subplot(2, 3, 5)
    maxValue = ceil(max(max(regress_NOx), max(L_NOx)));
    minValue = floor(min(min(regress_NOx), min(L_NOx)));
    scatter(regress_NOx, L_NOx, 30, year, 'filled'), hold on
%     h = colorbar;
%     set(get(h, 'label'), 'string', 'year');
    plot([minValue: maxValue], [minValue: maxValue], '--'),
    xlim([minValue, maxValue]), ylim([minValue, maxValue]),
    fit = fitlm(regress_NOx, L_NOx);
    Rsquare = fit.Rsquared.Ordinary;
    txt = strcat('R^2: ', num2str(Rsquare, '%.3f'));
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
    text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
    xlabel('Empirical Sensitivity'), ylabel('Log-Quadratic Model Sensitivity'),
    title('\partial O_3/\partial NO_x'),
    axis square;
    hold off;

    % plot dO3/dVOC
    subplot(2, 3, 6)
    maxValue = ceil(max(max(regress_VOC), max(L_VOC)));
    minValue = floor(min(min(regress_VOC), min(L_VOC)));
    scatter(regress_VOC, L_VOC, 30, year, 'filled'), hold on
%     h = colorbar;
%     set(get(h, 'label'), 'string', 'year');
    plot([minValue: maxValue], [minValue: maxValue], '--'),
    xlim([minValue, maxValue]), ylim([minValue, maxValue]),
    fit = fitlm(regress_VOC, L_VOC);
    Rsquare = fit.Rsquared.Ordinary;
    txt = strcat('R^2: ', num2str(Rsquare, '%.3f'));
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
    text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
    xlabel('Empirical Sensitivity'), ylabel('Log-Quadratic Model Sensitivity'),
    title('\partial O_3/\partial VOC'),
    axis square;
    hold off;
    
    hp = get(subplot(2,3,6),'Position');
    h = colorbar('Position', [hp(1)+hp(3)+0.015  hp(2)  0.015  hp(2)+hp(3)*3.2]);
    h.Title.String = "Year";
%     set(get(h, 'label'), 'string', 'year');
    sgtitle(siteName{siteNum})

    % save figures
    % saveas(gcf, strcat('/Users/zongrunli/Desktop/DDM-CMAQ Updates/figures top 4/', siteName{siteNum}, ' sensitive.jpg'))
    % save data
    %for m = 1: size(QO3, 2)
        %siteList = [siteList; convertCharsToStrings(siteName{siteNum})];
    %end
    %ODV_L = [ODV_L; LO3'];
    %S_N = [S_N; L_NOx'];
    %S_V = [S_V; L_VOC'];
end