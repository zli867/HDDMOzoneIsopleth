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

    % Matrix = [QO3; LO3; CMAQO3; QO3_NOx; LO3_NOx; DDM_NOx; QO3_VOC; LO3_VOC; DDM_VOC; QO3_NOx_VOC; LO3_NOx_VOC; DDM_NOx_VOC; QO3_NOx_2; LO3_NOx_2; DDM_NOx_2; QO3_VOC_2; LO3_VOC_2; DDM_VOC_2];
    % maxValue = max(max(Matrix));
    % minValue = min(min(Matrix));
    
    %Calculate
    x = 1:156;
    resQO3 = (CMAQO3 - QO3) .^ 2;
    resLO3 = (CMAQO3 - LO3) .^ 2;
    resQNOx = (DDM_NOx - QO3_NOx) .^2;
    resLNOx = (DDM_NOx - LO3_NOx) .^2;
    resQVOC = (DDM_VOC - QO3_VOC) .^2;
    resLVOC = (DDM_VOC - LO3_VOC) .^2;
    resQNOxVOC = (DDM_NOx_VOC - QO3_NOx_VOC) .^2;
    resLNOxVOC = (DDM_NOx_VOC - LO3_NOx_VOC) .^2;
    resQNOx2 = (DDM_NOx_2 - QO3_NOx_2).^2;
    resLNOx2 = (DDM_NOx_2 - LO3_NOx_2).^2;
    resQVOC2 = (DDM_VOC_2 - QO3_VOC_2).^2;
    resLVOC2 = (DDM_VOC_2 - LO3_VOC_2).^2;


    % resQO3 = sort(resQO3, 1, 'descend');
    % resQO3 = resQO3(1:10);
    % resLO3 = sort(resLO3, 1, 'descend');
    % resLO3 = resLO3(1:10);
    % resQNOx = sort(resQNOx, 1, 'descend');
    % resQNOx = resQNOx(1:10);
    % resLNOx = sort(resLNOx, 1, 'descend');
    % resLNOx = resLNOx(1:10);
    % resQVOC = sort(resQVOC, 1, 'descend');
    % resQVOC = resQVOC(1:10);
    % resLVOC = sort(resLVOC, 1, 'descend');
    % resLVOC = resLVOC(1:10);
    % resQNOxVOC = sort(resQNOxVOC, 1, 'descend');
    % resQNOxVOC = resQNOxVOC(1:10);
    % resLNOxVOC = sort(resLNOxVOC, 1, 'descend');
    % resLNOxVOC = resLNOxVOC(1:10);
    % resQNOx2 = sort(resQNOx2, 1, 'descend');
    % resQNOx2 = resQNOx2(1:10);
    % resLNOx2 = sort(resLNOx2, 1, 'descend');
    % resLNOx2 = resLNOx2(1:10);
    % resQVOC2 = sort(resQVOC2, 1, 'descend');
    % resQVOC2 = resQVOC2(1:10);
    % resLVOC2 = sort(resLVOC2, 1, 'descend');
    % resLVOC2 = resLVOC2(1:10);


    bar(resQO3), hold on
    bar(14: 26, resLO3),
    bar(27: 39, resQNOx)
    bar(40: 52, resQNOx)
    bar(53: 65, resQNOx)
    bar(66: 78, resQNOx)
    bar(79: 91, resQNOx)
    bar(92: 104, resQNOx)
    bar(105: 117, resQNOx)
    bar(118: 130, resQNOx)
    bar(131: 143, resQNOx)
    bar(144: 156, resQNOx), hold off;
  
    % bar(resQO3), hold on
    % bar(11: 20, resLO3),
    % bar(21: 30, resQNOx)
    % bar(31: 40, resQNOx)
    % bar(41: 50, resQNOx)
    % bar(51: 60, resQNOx)
    % bar(61: 70, resQNOx)
    % bar(71: 80, resQNOx)
    % bar(81: 90, resQNOx)
    % bar(91: 100, resQNOx)
    % bar(101: 110, resQNOx)
    % bar(111: 120, resQNOx), hold off;
    % save figures
    % saveas(gcf, strcat('/Users/zongrunli/Desktop/CMAQ-DDM/figures/', siteName{siteNum}, ' CMAQTraj.jpg'))
end