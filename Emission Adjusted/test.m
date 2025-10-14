T = readtable('/Users/zongrunli/Desktop/CMAQ-DDM/data/Observations.xlsx');
load('/Users/zongrunli/Desktop/CMAQ-DDM/data/QCoeff.mat');
siteName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Rubidoux', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernardino'};
siteTableName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Riverside', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernadino'};
baseNOx = 557.384;
baseVOC = 521.582;
for siteNum = 1: 1
    ODV = [];
	E_N = [];
	E_V = [];
	t = [];
    for i = 1: size(T, 1)
        temp = T(i, :);
        if temp.Year ~= 1111 & temp.Year ~= 2222 & temp.Year ~= 3333 & strcmp(temp.SiteName, siteTableName{siteNum})
            NOx = temp.NOx/baseNOx;
            VOC = temp.VOC/baseVOC;
            ODVTemp = temp.ODV;
            t = [t; temp.Year - 1980];
            ODV = [ODV;ODVTemp];
			E_N = [E_N; NOx];
			E_V = [E_V; VOC];
        end
    end
    % adjusted_N = (E_N .* (beta(1) + beta(2) .* t))
    % adjusted_V = (E_V .* (beta(3) + beta(4) .* t))
    fun = @(beta)ODV - (QCoeff(siteNum, 1) + QCoeff(siteNum, 2) .* (E_N .* (beta(1) + beta(2) .* t)) + QCoeff(siteNum, 3) .* (E_V .* (beta(3) + beta(4) .* t)) + QCoeff(siteNum, 4) .* (E_N .* (beta(1) + beta(2) .* t)) .^2 + QCoeff(siteNum, 5) .* (E_N .* (beta(1) + beta(2) .* t)) .* (E_V .* (beta(3) + beta(4) .* t)) + QCoeff(siteNum, 6) .* (E_V .* (beta(3) + beta(4) .* t)) .^2);
    x0 = [1, 0, 1, 0];
    for i = 1:100
    	result = lsqnonlin(fun, x0)
    	x0 = result;
    end
    beta = result;
    ODV_adjusted = (QCoeff(siteNum, 1) + QCoeff(siteNum, 2) .* (E_N .* (beta(1) + beta(2) .* t)) + QCoeff(siteNum, 3) .* (E_V .* (beta(3) + beta(4) .* t)) + QCoeff(siteNum, 4) .* (E_N .* (beta(1) + beta(2) .* t)) .^2 + QCoeff(siteNum, 5) .* (E_N .* (beta(1) + beta(2) .* t)) .* (E_V .* (beta(3) + beta(4) .* t)) + QCoeff(siteNum, 6) .* (E_V .* (beta(3) + beta(4) .* t)) .^2);
    figure(1)
    maxValue = max(max(ODV), max(ODV_adjusted)) + 1;
    scatter(ODV, ODV_adjusted), hold on
    plot([0:1:maxValue], [0:1:maxValue], '--')
    xlabel("ODV"), 
    ylabel("Adjusted ODV")
    hold off
    
    adjusted_N = (E_N .* (beta(1) + beta(2) .* t));
    adjusted_V = (E_V .* (beta(3) + beta(4) .* t));

    figure(2)
    maxValue = max(max(E_N), max(adjusted_N)) + 1;
    scatter(E_N, adjusted_N), hold on
    plot([0:1:maxValue], [0:1:maxValue], '--'),
    xlabel("NO_x Emission")
    ylabel("Adjusted NO_x Emission")
    hold off

    figure(3)
    maxValue = max(max(E_V), max(adjusted_V)) + 1;
    scatter(E_V, adjusted_V), hold on
    plot([0:1:maxValue], [0:1:maxValue], '--'),
    xlabel("VOC Emission")
    ylabel("Adjusted VOC Emission")
    hold off
    
end
