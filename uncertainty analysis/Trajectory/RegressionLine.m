load('/Users/zongrunli/Desktop/CMAQ-DDM/data/LCoeff.mat');
load('/Users/zongrunli/Desktop/CMAQ-DDM/data/QCoeff.mat');
figureNum = 1;

N = [0.01, 1, 3.4, 0.1, 1, 2, 0.1, 3, 1.1, 0.9, 0.9, 1.1, 2.04];
V = [0.01, 1, 4.9, 1, 0.1, 0.1, 4, 0.1, 1.1, 0.9, 1, 1, 1.95];
x1 = [2: 13: 301];
x2 = [14: 13: 313];
site = {};
siteTableName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Riverside', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernadino'};
siteName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Rubidoux', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernardino'};
baseNOx = 557.384;
baseVOC = 521.582;
for k = 1: 24
	site{k} = ['D', int2str(x1(k)), ':', 'I', int2str(x2(k))];
end

for siteNum = 1: 1
	% find D
	sheets = {'1percent', 'base', '1985', 'n10v100', 'n100v10', 'n200v10', 'n10v400', 'n300v10', 'n110v110', 'n90v90', 'n90v100', 'n110v100', '2001'};
    D = [];
    
    % Average the top 4 ozone data
	for i = 1: 13
		sheet = sheets{i};
		temp = xlsread('/Users/zongrunli/Desktop/CMAQ-DDM/data/sens_summary_4k_DDM_1985up.xlsx', sheet, site{siteNum});
		temp(:, 2) = temp(:, 2)/N(i);
		temp(:, 3) = temp(:, 3)/V(i);
		temp(:, 4) = temp(:, 4)/(N(i) * V(i));
		temp(:, 5) = temp(:, 5)/N(i)^2;
		temp(:, 6) = temp(:, 6)/V(i)^2;
		temp = sort(temp, 1, 'descend');
		tempData  = temp(1:4, :);
		tempData = mean(tempData);
		D = [D, transpose(tempData)];
    end
    
    NOx = N;
    VOC = V;
    % data points
    % Ozone
    QO3 = QCoeff(i, 1) + QCoeff(i, 2) * NOx + QCoeff(i, 3) * VOC + QCoeff(i, 4) * NOx.^2 + QCoeff(i, 5) * NOx .* VOC + QCoeff(i, 6) * VOC.^2;
    LO3 = exp(LCoeff(i, 1) + LCoeff(i, 2) * NOx + LCoeff(i, 3) * VOC + LCoeff(i, 4) * NOx.^2 + LCoeff(i, 5) * NOx .* VOC + LCoeff(i, 6) * VOC.^2);
    CMAQO3 = D(1, :);
    % dO3/dNOx
    QO3_NOx = QCoeff(i, 2) + 2*QCoeff(i, 4) * NOx + QCoeff(i, 5) * VOC;
    LO3_NOx = (LCoeff(i, 2) + 2*LCoeff(i, 4) * NOx + LCoeff(i, 5) * VOC) .* LO3;
    DDM_NOx = D(2, :);
    % dO3/dVOC
    QO3_VOC = QCoeff(i, 3) + QCoeff(i, 5) * NOx + 2* QCoeff(i, 6) * VOC;
    LO3_VOC = (LCoeff(i, 3) + LCoeff(i, 5) * NOx + 2* LCoeff(i, 6) * VOC) .* LO3;
    DDM_VOC = D(3, :);
    % dO3/dVOCdNOx
    QO3_NOx_VOC = QCoeff(i, 5) * ones(1, 13);
    LO3_NOx_VOC = LO3_VOC .* (LCoeff(i, 2) + 2*LCoeff(i, 4) * NOx + LCoeff(i, 5) * VOC) + LCoeff(i, 5) * LO3;
    DDM_NOx_VOC = D(4, :);
    %dO3/dNOx^2
    QO3_NOx_2 = 2* QCoeff(i, 4) * ones(1, 13);
    LO3_NOx_2 = LO3_NOx .* (LCoeff(i, 2) + 2*LCoeff(i, 4) * NOx + LCoeff(i, 5) * VOC) + 2* LCoeff(i, 4) * LO3;
    DDM_NOx_2 = D(5, :);
    %dO3/dVOC^2
    QO3_VOC_2 = 2* QCoeff(i, 6) * ones(1, 13);
    LO3_VOC_2 = (LCoeff(i, 3) + LCoeff(i, 5) * NOx + 2* LCoeff(i, 6) * VOC) .* LO3_VOC + 2* LCoeff(i, 6) * LO3;
    DDM_VOC_2 = D(6, :);

    Matrix = [QO3; LO3; CMAQO3; QO3_NOx; LO3_NOx; DDM_NOx; QO3_VOC; LO3_VOC; DDM_VOC; QO3_NOx_VOC; LO3_NOx_VOC; DDM_NOx_VOC; QO3_NOx_2; LO3_NOx_2; DDM_NOx_2; QO3_VOC_2; LO3_VOC_2; DDM_VOC_2];
    maxValue = max(max(Matrix));
    minValue = min(min(Matrix));
    
    % Quadratic
    figure(figureNum)
    scatter(CMAQO3, QO3, 'filled'), hold on
    scatter(DDM_NOx, QO3_NOx, 'filled'),
    scatter(DDM_VOC, QO3_VOC, 'filled'),
    scatter(DDM_NOx_VOC, QO3_NOx_VOC, 'filled'),
    scatter(DDM_NOx_2, QO3_NOx_2, 'filled'),
    scatter(DDM_VOC_2, QO3_VOC_2, 'filled'),
    plot([minValue: maxValue], [minValue: maxValue], '--'),
    xlabel('CMAQ Estimated Ozone (ppb)'), ylabel('Quadratic Model Estimated Ozone (ppb)')
    fit = fitlm(CMAQO3, QO3);
    RSquare = fit.Rsquared.Ordinary;
    txt1 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    fit = fitlm(DDM_NOx, QO3_NOx);
    RSquare = fit.Rsquared.Ordinary;
    txt2 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    fit = fitlm(DDM_VOC, QO3_VOC);
    RSquare = fit.Rsquared.Ordinary;
    txt3 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    fit = fitlm(DDM_NOx_VOC, QO3_NOx_VOC);
    RSquare = fit.Rsquared.Ordinary;
    txt4 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    fit = fitlm(DDM_NOx_2, QO3_NOx_2);
    RSquare = fit.Rsquared.Ordinary;
    txt5 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    fit = fitlm(DDM_VOC_2, QO3_VOC_2);
    RSquare = fit.Rsquared.Ordinary;
    txt6 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    txt = {txt1, txt2, txt3, txt4, txt5, txt6};
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
    % text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
    title(horzcat('Quadratic Model', ' ', siteName{siteNum}))
    legend('O_3', '\partial O_3/\partial NO_x', '\partial O_3/\partial VOC',...
     '\partial^2 O_3/\partial NO_xVOC', '\partial^2 O_3/\partial {NO_x}^2', '\partial^2 O_3/\partial VOC^2', 'Location','northwest')
    axis square;
    hold off
    saveas(gcf, strcat('/Users/zongrunli/Desktop/CMAQ-DDM/figures/', siteName{siteNum}, ' quadraScatter.jpg'))

    figureNum = figureNum + 1;
    % Log-Quadratic
    figure(figureNum)
    scatter(CMAQO3, LO3, 'filled'), hold on
    scatter(DDM_NOx, LO3_NOx, 'filled'),
    scatter(DDM_VOC, LO3_VOC, 'filled'),
    scatter(DDM_NOx_VOC, LO3_NOx_VOC, 'filled'),
    scatter(DDM_NOx_2, LO3_NOx_2, 'filled'),
    scatter(DDM_VOC_2, LO3_VOC_2, 'filled'),
    plot([minValue: maxValue], [minValue: maxValue], '--'),
    xlabel('CMAQ Estimated Ozone (ppb)'), ylabel('Quadratic Model Estimated Ozone (ppb)')
    fit = fitlm(CMAQO3, LO3);
    RSquare = fit.Rsquared.Ordinary;
    txt1 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    fit = fitlm(DDM_NOx, LO3_NOx);
    RSquare = fit.Rsquared.Ordinary;
    txt2 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    fit = fitlm(DDM_VOC, LO3_VOC);
    RSquare = fit.Rsquared.Ordinary;
    txt3 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    fit = fitlm(DDM_NOx_VOC, LO3_NOx_VOC);
    RSquare = fit.Rsquared.Ordinary;
    txt4 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    fit = fitlm(DDM_NOx_2, LO3_NOx_2);
    RSquare = fit.Rsquared.Ordinary;
    txt5 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    fit = fitlm(DDM_VOC_2, LO3_VOC_2);
    RSquare = fit.Rsquared.Ordinary;
    txt6 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    txt = {txt1, txt2, txt3, txt4, txt5, txt6};
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
    % text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
    title(horzcat('Log-Quadratic Model', ' ', siteName{siteNum}))
    legend('O_3', '\partial O_3/\partial NO_x', '\partial O_3/\partial VOC',...
     '\partial^2 O_3/\partial NO_xVOC', '\partial^2 O_3/\partial {NO_x}^2', '\partial^2 O_3/\partial VOC^2', 'Location','northwest')
    axis square;
    hold off
    figureNum = figureNum + 1;
    %

    % save figures
    saveas(gcf, strcat('/Users/zongrunli/Desktop/CMAQ-DDM/figures/', siteName{siteNum}, ' logScatter.jpg'))
end