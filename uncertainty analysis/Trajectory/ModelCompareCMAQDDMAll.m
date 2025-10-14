% load('/Users/zongrunli/Desktop/CMAQ-DDM/data/LCoeffWeightedAll.mat');
% load('/Users/zongrunli/Desktop/CMAQ-DDM/data/QCoeffWeightedAll.mat');
% load('/Users/zongrunli/Desktop/CMAQ-DDM/data/LCoeff.mat');
% load('/Users/zongrunli/Desktop/CMAQ-DDM/data/QCoeff.mat');
load('/Users/zongrunli/Desktop/CMAQ-DDM/data/LCoeffPart.mat');
load('/Users/zongrunli/Desktop/CMAQ-DDM/data/QCoeffPart.mat');

N = [0.01, 1, 3.9, 0.1, 1, 2, 0.1, 3, 1.1, 0.9, 0.9, 1.1, 2.3, 1.3, 0.8];
V = [0.01, 1, 4.9, 1, 0.1, 0.1, 4, 0.1, 1.1, 0.9, 1, 1, 3.6, 1.5, 0.6];
x1 = [2: 13: 301];
x2 = [14: 13: 313];
site = {};
siteTableName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Riverside', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernadino'};
siteName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Rubidoux', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernardino'};
baseNOx = 395.748;
baseVOC = 401.412;
for k = 1: 24
	site{k} = ['D', int2str(x1(k)), ':', 'I', int2str(x2(k))];
end

for siteNum = 1: 1
	% find D
	sheets = {'1percent', 'base', '1985', 'n10v100', 'n100v10', 'n200v10', 'n10v400', 'n300v10', 'n110v110', 'n90v90', 'n90v100', 'n110v100', '2001', '2011', '2028'};
    D = [];
    
    % Average the top ozone data
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
    
    NOx = N;
    VOC = V;
    figure(siteNum)
    % Ozone
    subplot(2, 3, 1)
    QO3 = QCoeff(i, 1) + QCoeff(i, 2) * NOx + QCoeff(i, 3) * VOC + QCoeff(i, 4) * NOx.^2 + QCoeff(i, 5) * NOx .* VOC + QCoeff(i, 6) * VOC.^2;
    LO3 = exp(LCoeff(i, 1) + LCoeff(i, 2) * NOx + LCoeff(i, 3) * VOC + LCoeff(i, 4) * NOx.^2 + LCoeff(i, 5) * NOx .* VOC + LCoeff(i, 6) * VOC.^2);
    CMAQO3 = D(1, :);
    maxValue = ceil(max(max(max(QO3), max(LO3)), max(CMAQO3))) + 1;

    scatter(CMAQO3, QO3, 'fill', 'r'), hold on, 
    scatter(CMAQO3, LO3, 'fill', 'b'), xlabel('CMAQ Estimated Ozone (ppb)'), ylabel('Isopleth Estimated Ozone (ppb)')
    plot([0: maxValue], [0: maxValue], '--'), hold on,
    xlim([0, maxValue]), ylim([0, maxValue]);
    fit = fitlm(CMAQO3, QO3);
    RSquare = fit.Rsquared.Ordinary;
    txt1 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    fit2 = fitlm(CMAQO3, LO3);
    RSquare = fit2.Rsquared.Ordinary;
    txt2 = strcat('R^2 (Log-Quadratic): ', num2str(RSquare, '%.3f'));
    txt = {txt1, txt2};
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
%     text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
    % legend('Quadratic', 'Log-Quadratic', 'Location','southeast')
    axis square;
    hold off;

    % dO3/dNOx
    subplot(2, 3, 2)
    QO3_NOx = QCoeff(i, 2) + 2*QCoeff(i, 4) * NOx + QCoeff(i, 5) * VOC;
    LO3_NOx = (LCoeff(i, 2) + 2*LCoeff(i, 4) * NOx + LCoeff(i, 5) * VOC) .* exp(LCoeff(i, 1) + LCoeff(i, 2) * NOx + LCoeff(i, 3) * VOC + LCoeff(i, 4) * NOx.^2 + LCoeff(i, 5) * NOx .* VOC + LCoeff(i, 6) * VOC.^2);
    DDM_NOx = D(2, :);
    maxValue = ceil(max(max(max(QO3_NOx), max(LO3_NOx)), max(DDM_NOx))) + 1;
    minValue = floor(min(min(min(QO3_NOx), min(LO3_NOx)), min(DDM_NOx)));

    scatter(DDM_NOx, QO3_NOx, 'fill', 'r'), hold on, 
    scatter(DDM_NOx, LO3_NOx, 'fill', 'b'), xlabel('CMAQ Estimated \partial O_3/\partial NO_x'), ylabel('Isopleth Estimated \partial O_3/\partial NO_x')
    plot([minValue: maxValue], [minValue: maxValue], '--'), hold on,
    xlim([minValue, maxValue]), ylim([minValue, maxValue]);
    fit = fitlm(DDM_NOx, QO3_NOx);
    RSquare = fit.Rsquared.Ordinary;
    txt1 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    fit2 = fitlm(DDM_NOx, LO3_NOx);
    RSquare = fit2.Rsquared.Ordinary;
    txt2 = strcat('R^2 (Log-Quadratic): ', num2str(RSquare, '%.3f'));
    txt = {txt1, txt2};
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
%     text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
%     legend('Quadratic', 'Log-Quadratic', 'Location','southeast')
    axis square;
    hold off;
    
    % dO3/dVOC
    subplot(2, 3, 3)
    QO3_VOC = QCoeff(i, 3) + QCoeff(i, 5) * NOx + 2* QCoeff(i, 6) * VOC;
    LO3_VOC = (LCoeff(i, 3) + LCoeff(i, 5) * NOx + 2* LCoeff(i, 6) * VOC) .* exp(LCoeff(i, 1) + LCoeff(i, 2) * NOx + LCoeff(i, 3) * VOC + LCoeff(i, 4) * NOx.^2 + LCoeff(i, 5) * NOx .* VOC + LCoeff(i, 6) * VOC.^2);
    DDM_VOC = D(3, :);
    maxValue = ceil(max(max(max(QO3_VOC), max(LO3_VOC)), max(DDM_VOC))) + 1;
    minValue = floor(min(min(min(QO3_VOC), min(LO3_VOC)), min(DDM_VOC)));

    scatter(DDM_VOC, QO3_VOC, 'fill', 'r'), hold on, 
    scatter(DDM_VOC, LO3_VOC, 'fill', 'b'), xlabel('CMAQ Estimated \partial O_3/\partial VOC'), ylabel('Isopleth Estimated \partial O_3/\partial VOC')
    plot([minValue: maxValue], [minValue: maxValue], '--'), hold on,
    xlim([minValue, maxValue]), ylim([minValue, maxValue]);
    fit = fitlm(DDM_VOC, QO3_VOC);
    RSquare = fit.Rsquared.Ordinary;
    txt1 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    fit2 = fitlm(DDM_VOC, LO3_VOC);
    RSquare = fit2.Rsquared.Ordinary;
    txt2 = strcat('R^2 (Log-Quadratic): ', num2str(RSquare, '%.3f'));
    txt = {txt1, txt2};
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
%     text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
%     legend('Quadratic', 'Log-Quadratic', 'Location','southeast')
    axis square;
    hold off;
    
    % d2O3/dVOCdNOx
    subplot(2, 3, 4)
    QO3_VOC_NOx = QCoeff(i, 5) * ones(1, 15);
    LO3_VOC_NOx = LCoeff(i, 5) * LO3 + (LCoeff(i, 3) + LCoeff(i, 5) * NOx + 2* LCoeff(i, 6) * VOC) .* LO3_NOx;
    DDM_VOC_NOx = D(4, :);
    maxValue = ceil(max(max(max(QO3_VOC_NOx), max(LO3_VOC_NOx)), max(DDM_VOC_NOx))) + 1;
    minValue = floor(min(min(min(QO3_VOC_NOx), min(LO3_VOC_NOx)), min(DDM_VOC_NOx)));

    scatter(DDM_VOC_NOx, QO3_VOC_NOx, 'fill', 'r'), hold on, 
    scatter(DDM_VOC_NOx, LO3_VOC_NOx, 'fill', 'b'), xlabel('CMAQ Estimated \partial^2 O_3/\partial VOC NO_x'), ylabel('Isopleth Estimated \partial^2 O_3/\partial VOC \partial NO_x')
    plot([minValue: maxValue], [minValue: maxValue], '--'), hold on,
    xlim([minValue, maxValue]), ylim([minValue, maxValue]);
    fit = fitlm(DDM_VOC_NOx, QO3_VOC_NOx);
    RSquare = fit.Rsquared.Ordinary;
    txt1 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    fit2 = fitlm(DDM_VOC_NOx, LO3_VOC_NOx);
    RSquare = fit2.Rsquared.Ordinary;
    txt2 = strcat('R^2 (Log-Quadratic): ', num2str(RSquare, '%.3f'));
    txt = {txt1, txt2};
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
%     text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
%     legend('Quadratic', 'Log-Quadratic', 'Location','southeast')
    axis square;
    hold off;
    
    % d2O3/dNOx2
    subplot(2, 3, 5)
    QO3_NOx_2 = 2 * QCoeff(i, 4) * ones(1, 15);
    LO3_NOx_2 = 2 * LCoeff(i, 4) * LO3 + (LCoeff(i, 2) + 2*LCoeff(i, 4) * NOx + LCoeff(i, 5) * VOC) .* LO3_NOx;
    DDM_NOx_2 = D(5, :);
    maxValue = ceil(max(max(max(QO3_NOx_2), max(LO3_NOx_2)), max(DDM_NOx_2))) + 1;
    minValue = floor(min(min(min(QO3_NOx_2), min(LO3_NOx_2)), min(DDM_NOx_2)));

    scatter(DDM_NOx_2, QO3_NOx_2, 'fill', 'r'), hold on, 
    scatter(DDM_NOx_2, LO3_NOx_2, 'fill', 'b'), xlabel('CMAQ Estimated \partial^2 O_3/{\partial NO_x}^2'), ylabel('Isopleth Estimated \partial^2 O_3/{\partial NO_x}^2')
    plot([minValue: maxValue], [minValue: maxValue], '--'), hold on,
    xlim([minValue, maxValue]), ylim([minValue, maxValue]);
    fit = fitlm(DDM_NOx_2, QO3_NOx_2);
    RSquare = fit.Rsquared.Ordinary;
    txt1 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    fit2 = fitlm(DDM_NOx_2, LO3_NOx_2);
    RSquare = fit2.Rsquared.Ordinary;
    txt2 = strcat('R^2 (Log-Quadratic): ', num2str(RSquare, '%.3f'));
    txt = {txt1, txt2};
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
%     text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
%     legend('Quadratic', 'Log-Quadratic', 'Location','southeast')
    axis square;
    hold off;
    
    
    % d2O3/dVOC2
    subplot(2, 3, 6)
    QO3_VOC_2 = 2* QCoeff(i, 6) * ones(1, 15);
    LO3_VOC_2 = 2 * LCoeff(i, 6) * LO3 + (LCoeff(i, 3) + LCoeff(i, 5) * NOx + 2* LCoeff(i, 6) * VOC) .* LO3_VOC;
    DDM_VOC_2 = D(6, :);
    maxValue = ceil(max(max(max(QO3_VOC_2), max(LO3_VOC_2)), max(DDM_VOC_2))) + 1;
    minValue = floor(min(min(min(QO3_VOC_2), min(LO3_VOC_2)), min(DDM_VOC_2)));

    scatter(DDM_VOC_2, QO3_VOC_2, 'fill', 'r'), hold on, 
    scatter(DDM_VOC_2, LO3_VOC_2, 'fill', 'b'), xlabel('CMAQ Estimated \partial^2 O_3/{\partial VOC}^2'), ylabel('Isopleth Estimated \partial^2 O_3/{\partial VOC}^2')
    plot([minValue: maxValue], [minValue: maxValue], '--'), hold on,
    xlim([minValue, maxValue]), ylim([minValue, maxValue]);
    fit = fitlm(DDM_VOC_2, QO3_VOC_2);
    RSquare = fit.Rsquared.Ordinary;
    txt1 = strcat('R^2 (Quadratic): ', num2str(RSquare, '%.3f'));
    fit2 = fitlm(DDM_VOC_2, LO3_VOC_2);
    RSquare = fit2.Rsquared.Ordinary;
    txt2 = strcat('R^2 (Log-Quadratic): ', num2str(RSquare, '%.3f'));
    txt = {txt1, txt2};
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
%     text(NW(1), NW(2), txt, 'VerticalAlignment','top', 'HorizontalAlignment','left')
%     legend('Quadratic', 'Log-Quadratic', 'Location','southeast')
    axis square;
    hold off;
    

    sgtitle(siteName{siteNum})

    % save figures
    % saveas(gcf, strcat('/Users/zongrunli/Desktop/DDM-CMAQ Updates/figures top 4/', siteName{siteNum}, ' CMAQTraj.jpg'))
end