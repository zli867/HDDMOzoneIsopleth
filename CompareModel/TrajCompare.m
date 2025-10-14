% regression model
load('/Volumes/Shield/CMAQ-DDM/data/LCoeff.mat');
load('/Volumes/Shield/CMAQ-DDM/data/QCoeff.mat');
T = readtable('/Volumes/Shield/CMAQ-DDM/data/Observations.xlsx');
siteName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Rubidoux', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernardino'};
siteTableName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Riverside', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernadino'};
baseNOx = 395.748;
baseVOC = 401.412;
% save data
%siteList = [];
%ODV_L = [];
%S_N = [];
%S_V = [];

for siteNum = 1: 1
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
            NOxtemp = temp.NOx;
            VOCtemp = temp.VOC;
            ODVTemp = temp.ODV;
            if temp.Year ~= 2222 & temp.Year ~= 3333 & temp.Year ~= 4444
                NOx = [NOx, NOxtemp];
                VOC = [VOC, VOCtemp];
                year = [year, temp.Year];
                ODVTemp = temp.ODV;
                ODV = [ODV, ODVTemp];
            end
            y = [y; ODVTemp];
            x = [1, NOxtemp, VOCtemp, NOxtemp^2, VOCtemp* NOxtemp, VOCtemp^2];
            X = [X; x];
        end
    end
    % pseudo inverse
    b = pinv(X) * log10(y);

	% regression sensitive
    regress_NOx = log(10) * ODV .* (b(2) + 2* b(4) .* NOx + b(5) * VOC);
    % dO3/dVOC
    regress_VOC = log(10) * ODV .* (b(3) + b(5) .* NOx + 2*b(6) .* VOC);

    % Quadratic sensitivity
    % dO3/dNOx
    NOx_norm = NOx./baseNOx;
    VOC_norm = VOC./baseVOC;
    QO3 = QCoeff(siteNum, 1) + QCoeff(siteNum, 2) * NOx_norm + QCoeff(siteNum, 3) * VOC_norm + QCoeff(siteNum, 4) * NOx_norm.^2 + QCoeff(siteNum, 5) * NOx_norm .* VOC_norm + QCoeff(siteNum, 6) * VOC_norm.^2;
    Q_NOx = (QCoeff(siteNum, 2) + 2* QCoeff(siteNum, 4) .* NOx_norm + QCoeff(siteNum, 5) .* VOC_norm)/baseNOx;
    % dO3/dVOC
    Q_VOC = (QCoeff(siteNum, 3) + QCoeff(siteNum, 5) .* NOx_norm + 2* QCoeff(siteNum, 6) .* VOC_norm)/baseVOC;

    % Log-Quadratic sensitivity
    % dO3/dNOx
    LO3 = exp(LCoeff(siteNum, 1) + LCoeff(siteNum, 2) * NOx_norm + LCoeff(siteNum, 3) * VOC_norm + LCoeff(siteNum, 4) * NOx_norm.^2 + LCoeff(siteNum, 5) * NOx_norm .* VOC_norm + LCoeff(siteNum, 6) * VOC_norm.^2);
    L_NOx = ((LCoeff(siteNum, 2) + 2* LCoeff(siteNum, 4) .* NOx_norm + LCoeff(siteNum, 5) .* VOC_norm)/baseNOx) .* LO3;
    % dO3/dVOC
    L_VOC = ((LCoeff(siteNum, 3) + LCoeff(siteNum, 5) .* NOx_norm + 2* LCoeff(siteNum, 6) .* VOC_norm)/baseVOC) .* LO3;
    
    labelFontsize = 15;
    % Ozone
    % Quadratic Model
    % plot O3
    f = figure(siteNum);
    f.Position = [0 0 1200 400];
    subplot(1, 3, 1)
    scatter(ODV, QO3, 50, year, 'filled', '^'), hold on
    fit = fitlm(ODV, QO3);
    Rsquare_Q = fit.Rsquared.Ordinary;
    xlabel('Observed O_3', 'fontsize', labelFontsize), ylabel('Model Estimated O_3', 'fontsize', labelFontsize),
    scatter(ODV, LO3, 50, year, 'filled', 's'), hold on
    fit = fitlm(ODV, LO3);
    Rsquare_L = fit.Rsquared.Ordinary;
    maxValue = max(max(max(ODV(:)), max(QO3(:))), max(LO3(:)));
    minValue = min(min(min(ODV(:)), min(QO3(:))), min(LO3(:)));
    plot([minValue:0.01:maxValue], [minValue:0.01:maxValue], '--')
    xlim([minValue, maxValue]), ylim([minValue, maxValue])
    txt = ['R^2 (Q): ', num2str(round(Rsquare_Q, 2), '%.2f'); 'R^2 (L): ', num2str(round(Rsquare_L, 2), '%.2f')];
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.02;
    text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
    xtick = xticklabels;
    ytick = [];
    for tick = 1: length(xtick)
        ytick = [ytick, str2num(xtick{tick})];
    end
    set(gca, 'Ytick', ytick)
%     set(gca,'Xtick',minValue:30:maxValue)
%     set(gca,'Ytick',minValue:30:maxValue)
%     legend(['Quadratic'; 'R^2: ', num2str(round(Rsquare_Q, 2), '%.2f')], ...
%         ['Log Quadratic'; '  R^2: ', num2str(round(Rsquare_L, 2), '%.2f'), '  '])
    hold off;
    daspect([1, 1, 1]);
    
    % plot dO3/dNOx
    subplot(1, 3, 2)
    scatter(regress_NOx, Q_NOx, 50, year, 'filled', '^'), hold on
    fit = fitlm(regress_NOx, Q_NOx);
    Rsquare_Q = fit.Rsquared.Ordinary;
    xlabel('Empirical dO_3/dNO_x', 'fontsize', labelFontsize), ylabel('Model Estimated dO_3/dNO_x', 'fontsize', labelFontsize),
    scatter(regress_NOx, L_NOx, 50, year, 'filled', 's'), hold on
    fit = fitlm(regress_NOx, L_NOx);
    Rsquare_L = fit.Rsquared.Ordinary;
    maxValue = max(max(max(regress_NOx(:)), max(Q_NOx(:))), max(L_NOx(:)));
    minValue = min(min(min(regress_NOx(:)), min(Q_NOx(:))), min(L_NOx(:)));
    plot([minValue:0.01:maxValue], [minValue:0.01:maxValue], '--')
    xlim([minValue, maxValue]), ylim([minValue, maxValue])
    txt = ['R^2 (Q): ', num2str(round(Rsquare_Q, 2), '%.2f'); 'R^2 (L): ', num2str(round(Rsquare_L, 2), '%.2f')];
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.02;
    text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
    xtick = xticklabels;
    ytick = [];
    for tick = 1: length(xtick)
        ytick = [ytick, str2num(xtick{tick})];
    end
    set(gca, 'Ytick', ytick)
    hold off;
    daspect([1, 1, 1]);
    title(siteName{siteNum}, 'fontsize', 18)
    
    % plot dO3/dVOC
    subplot(1, 3, 3)
    scatter(regress_VOC, Q_VOC, 50, year, 'filled', '^'), hold on
    fit = fitlm(regress_VOC, Q_VOC);
    Rsquare_Q = fit.Rsquared.Ordinary;
    xlabel('Empirical dO_3/dVOC', 'fontsize', labelFontsize), ylabel('Model Estimated dO_3/dVOC', 'fontsize', labelFontsize),
    scatter(regress_VOC, L_VOC, 50, year, 'filled', 's'), hold on
    fit = fitlm(regress_VOC, L_VOC);
    Rsquare_L = fit.Rsquared.Ordinary;
    maxValue = max(max(max(regress_VOC(:)), max(Q_VOC(:))), max(L_VOC(:)));
    minValue = min(min(min(regress_VOC(:)), min(Q_VOC(:))), min(L_VOC(:)));
    plot([minValue:0.01:maxValue], [minValue:0.01:maxValue], '--')
    xlim([minValue, maxValue]), ylim([minValue, maxValue])
    txt = ['R^2 (Q): ', num2str(round(Rsquare_Q, 2), '%.2f'); 'R^2 (L): ', num2str(round(Rsquare_L, 2), '%.2f')];
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.02;
    text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
    xtick = xticklabels;
    ytick = [];
    for tick = 1: length(xtick)
        ytick = [ytick, str2num(xtick{tick})];
    end
    set(gca, 'Ytick', ytick)
    hold off;
    daspect([1, 1, 1]);
    
    hp = get(subplot(1,3,3),'Position');
    h = colorbar('Position', [hp(1)+hp(3)+0.015  hp(2)+0.08  0.015  hp(2)+hp(3)*2.6]);
    h.Title.String = "Year";
    h.FontSize = 15;
    
    lgd = legend({'Quadratic model Traj','Log Quadratic Model Traj'},'Orientation','horizontal', 'fontsize', 15);
    set(lgd, 'Position', [0.5, 0.02, 0.07, 0.05])
    
    % % save figures
	% saveas(gcf, strcat('/Volumes/Shield/CMAQ-DDM/Figure/', siteName{siteNum}, '_compare_traj.png'))
end