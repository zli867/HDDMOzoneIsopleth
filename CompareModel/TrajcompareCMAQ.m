load('/Users/zongrunli/Desktop/CMAQ-DDM/data/LCoeff.mat');
load('/Users/zongrunli/Desktop/CMAQ-DDM/data/QCoeff.mat');

N = [0.01, 1, 3.9, 0.1, 1, 2, 0.1, 3, 1.1, 0.9, 0.9, 1.1, 2.3, 1.3, 0.8];
V = [0.01, 1, 4.9, 1, 0.1, 0.1, 4, 0.1, 1.1, 0.9, 1, 1, 3.6, 1.5, 0.6];
x1 = [2: 13: 301];
x2 = [14: 13: 313];
siteName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Rubidoux', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernardino'};
baseNOx = 395.748;
baseVOC = 401.412;
for k = 1: 24
	site{k} = ['D', int2str(x1(k)), ':', 'I', int2str(x2(k))];
end

for siteNum = 1: 24
	% find D
	sheets = {'1percent', 'base', '1985', 'n10v100', 'n100v10', 'n200v10', 'n10v400', 'n300v10', 'n110v110', 'n90v90', 'n90v100', 'n110v100', '2001', '2011', '2028'};
    D = [];
    
    % max ozone data
	for i = 1: 15
		sheet = sheets{i};
		temp = xlsread('/Users/zongrunli/Desktop/CMAQ-DDM/data/sens_summary_4k_DDM_1985up.xlsx', sheet, site{siteNum});
		temp(:, 2) = temp(:, 2)/N(i);
		temp(:, 3) = temp(:, 3)/V(i);
		temp(:, 4) = temp(:, 4)/(N(i) * V(i));
		temp(:, 5) = temp(:, 5)/N(i)^2;
		temp(:, 6) = temp(:, 6)/V(i)^2;
		temp = sort(temp, 1, 'descend');
		tempData  = temp(1, :);
		D = [D, transpose(tempData)];
    end
    
    dataNum = [1:1:15];
    NOx = N;
    VOC = V;
    f = figure(siteNum);
    f.Position = [0 0 1200 400];
    % Ozone
    subplot(1, 3, 1)
    QO3 = QCoeff(i, 1) + QCoeff(i, 2) * NOx + QCoeff(i, 3) * VOC + QCoeff(i, 4) * NOx.^2 + QCoeff(i, 5) * NOx .* VOC + QCoeff(i, 6) * VOC.^2;
    LO3 = exp(LCoeff(i, 1) + LCoeff(i, 2) * NOx + LCoeff(i, 3) * VOC + LCoeff(i, 4) * NOx.^2 + LCoeff(i, 5) * NOx .* VOC + LCoeff(i, 6) * VOC.^2);
    CMAQO3 = D(1, :);
    maxValue = max(max(max(QO3), max(LO3)), max(CMAQO3)) + 1;
    minValue = min(min(min(QO3), min(LO3)), min(CMAQO3)) - 1;
    
%     scatter(ODV, QO3, 50, year, 'filled', '^')
    scatter(CMAQO3, QO3, 50, dataNum, 'fill', '^'), hold on, 
    scatter(CMAQO3, LO3, 50, dataNum, 'fill', 's'), xlabel('CMAQ Estimated O_3 (ppb)', 'fontsize', 15), ylabel('Isopleth Estimated O_3 (ppb)', 'fontsize', 15)
    plot([minValue: 0.01: maxValue], [minValue: 0.01: maxValue], '--'), hold on,
    xlim([minValue, maxValue]), ylim([minValue, maxValue]);
    fit = fitlm(CMAQO3, QO3);
    RSquare = fit.Rsquared.Ordinary;
    txt1 = strcat('R^2 (Quadratic): ', num2str(round(RSquare, 2), '%.2f'));
    fit2 = fitlm(CMAQO3, LO3);
    RSquare = fit2.Rsquared.Ordinary;
    txt2 = strcat('R^2 (Log-Quadratic): ', num2str(round(RSquare, 2), '%.2f'));
    txt = {txt1, txt2};
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
    text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left', 'fontsize', 13)
    daspect([1, 1, 1]);
    % Use same tick for x and y axis
    xtick = xticklabels;
    ytick = [];
    for tick = 1: length(xtick)
        ytick = [ytick, str2num(xtick{tick})];
    end
    set(gca, 'Ytick', ytick)
    hold off;

    % dO3/dNOx
    subplot(1, 3, 2)
    QO3_NOx = (QCoeff(i, 2) + 2*QCoeff(i, 4) * NOx + QCoeff(i, 5) * VOC)/baseNOx;
    LO3_NOx = ((LCoeff(i, 2) + 2*LCoeff(i, 4) * NOx + LCoeff(i, 5) * VOC) .* exp(LCoeff(i, 1) + LCoeff(i, 2) * NOx + LCoeff(i, 3) * VOC + LCoeff(i, 4) * NOx.^2 + LCoeff(i, 5) * NOx .* VOC + LCoeff(i, 6) * VOC.^2))/baseNOx;
    DDM_NOx = D(2, :)/baseNOx;
    maxValue = max(max(max(QO3_NOx), max(LO3_NOx)), max(DDM_NOx));
    minValue = min(min(min(QO3_NOx), min(LO3_NOx)), min(DDM_NOx));

    scatter(DDM_NOx, QO3_NOx, 50, dataNum, 'fill', '^'), hold on, 
    scatter(DDM_NOx, LO3_NOx, 50, dataNum, 'fill', 's'), xlabel('DDM Estimated dO_3/dNO_x', 'fontsize', 15), ylabel('Isopleth Estimated dO_3/dNO_x', 'fontsize', 15)
    plot([minValue: 0.01: maxValue], [minValue: 0.01: maxValue], '--'), hold on,
    xlim([minValue, maxValue]), ylim([minValue, maxValue]);
    fit = fitlm(DDM_NOx, QO3_NOx);
    RSquare = fit.Rsquared.Ordinary;
    txt1 = strcat('R^2 (Quadratic): ', num2str(round(RSquare, 2), '%.2f'));
    fit2 = fitlm(DDM_NOx, LO3_NOx);
    RSquare = fit2.Rsquared.Ordinary;
    txt2 = strcat('R^2 (Log-Quadratic): ', num2str(round(RSquare, 2), '%.2f'));
    txt = {txt1, txt2};
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
    text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left', 'fontsize', 13)
    daspect([1, 1, 1]);
    title(siteName{siteNum}, 'fontsize', 18)
    % Use same tick for x and y axis
    xtick = xticklabels;
    ytick = [];
    for tick = 1: length(xtick)
        ytick = [ytick, str2num(xtick{tick})];
    end
    set(gca, 'Ytick', ytick)
    hold off;
    
    % dO3/dVOC
    subplot(1, 3, 3)
    QO3_VOC = (QCoeff(i, 3) + QCoeff(i, 5) * NOx + 2* QCoeff(i, 6) * VOC)/baseVOC;
    LO3_VOC = ((LCoeff(i, 3) + LCoeff(i, 5) * NOx + 2* LCoeff(i, 6) * VOC) .* exp(LCoeff(i, 1) + LCoeff(i, 2) * NOx + LCoeff(i, 3) * VOC + LCoeff(i, 4) * NOx.^2 + LCoeff(i, 5) * NOx .* VOC + LCoeff(i, 6) * VOC.^2))/baseVOC;
    DDM_VOC = D(3, :)/baseVOC;
    maxValue = max(max(max(QO3_VOC), max(LO3_VOC)), max(DDM_VOC));
    minValue = min(min(min(QO3_VOC), min(LO3_VOC)), min(DDM_VOC));

    scatter(DDM_VOC, QO3_VOC, 50, dataNum, 'fill', '^'), hold on, 
    scatter(DDM_VOC, LO3_VOC, 50, dataNum, 'fill', 's'), xlabel('DDM Estimated dO_3/dVOC', 'fontsize', 15), ylabel('Isopleth Estimated dO_3/dVOC', 'fontsize', 15)
    plot([minValue: 0.01: maxValue], [minValue:0.01: maxValue], '--'), hold on,
    xlim([minValue, maxValue]), ylim([minValue, maxValue]);
    fit = fitlm(DDM_VOC, QO3_VOC);
    RSquare = fit.Rsquared.Ordinary;
    txt1 = strcat('R^2 (Quadratic): ', num2str(round(RSquare, 2), '%.2f'));
    fit2 = fitlm(DDM_VOC, LO3_VOC);
    RSquare = fit2.Rsquared.Ordinary;
    txt2 = strcat('R^2 (Log-Quadratic): ', num2str(round(RSquare, 2), '%.2f'));
    txt = {txt1, txt2};
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
    text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left', 'fontsize', 13)
    legend('Quadratic', 'Log-Quadratic', 'Location','southeast')
    daspect([1, 1, 1]);
    % Use same tick for x and y axis
    xtick = xticklabels;
    ytick = [];
    for tick = 1: length(xtick)
        ytick = [ytick, str2num(xtick{tick})];
    end
    set(gca, 'Ytick', ytick)
    hold off;
    
    hp = get(subplot(1,3,3),'Position');
    h = colorbar('Position', [hp(1)+hp(3)+0.015  hp(2)+0.08  0.015  hp(2)+hp(3)*2.6], ...
                 'XTickLabel',sheets, ...
                 'XTick', 1:1:15);
    h.Title.String = "Data";
    h.FontSize = 15;
        
    lgd = legend({'Quadratic model Traj','Log Quadratic Model Traj'},'Orientation','horizontal', 'fontsize', 15);
    set(lgd, 'Position', [0.5, 0.02, 0.07, 0.05])
    
    % save figures
	% saveas(gcf, strcat('/Users/zongrunli/Desktop/CMAQ-DDM/Figure/', siteName{siteNum}, '_compare_CMAQ_traj.png'))
end