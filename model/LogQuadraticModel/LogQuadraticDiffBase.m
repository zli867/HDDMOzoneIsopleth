LCoeff = [];
T = readtable('/Users/zongrunli/Desktop/CMAQ-DDM/data/Observations.xlsx');
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

% Set Base
base = 100;

for siteNum = 1: 24
	% find D
	sheets = {'1percent', 'base', '1985', 'n10v100', 'n100v10', 'n200v10', 'n10v400', 'n300v10', 'n110v110', 'n90v90', 'n90v100', 'n110v100', '2001'};
	D = [];
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
    % D

    % change D to log scale
    for i = 1:13
        O = D(1, i);
        D(1, i) = log(D(1, i))/log(base);
        D(2, i) = ((1/O) * D(2, i))/log(base);
        D(3, i) = ((1/O) * D(3, i))/log(base);
        D(4, i) = ((1/O) * D(4, i) - (1/O^2) * D(2, i) * D(3, i))/log(base);
        D(5, i) = ((1/O) * D(5, i) - (1/O^2) * D(2, i) * D(2, i))/log(base);
        D(6, i) = ((1/O) * D(6, i) - (1/O^2) * D(3, i) * D(3, i))/log(base);
    end
    % D
    % b
    b = [];
    for i = 1: size(D, 2)
        b = [b; D(:, i)];
    end

    % A
    A = [];
    for i = 1: size(D, 2)
        tempA = [1, N(i), V(i), N(i)^2, N(i)*V(i), V(i)^2;
                0, 1, 0, 2*N(i), V(i), 0;
                0, 0, 1, 0, N(i), 2*V(i);
                0, 0, 0, 0, 1, 0;
                0, 0, 0, 2, 0, 0;
                0, 0, 0, 0, 0, 2];
        A = [A; tempA];
    end

    fit = fitlm(A,b,'Intercept',false);
    a0 = fit.Coefficients.Estimate(1);
    a1 = fit.Coefficients.Estimate(2);
    a2 = fit.Coefficients.Estimate(3);
    a3 = fit.Coefficients.Estimate(4);
    a4 = fit.Coefficients.Estimate(5);
    a5 = fit.Coefficients.Estimate(6);
    
	scatterN = [];
    scatterV = [];
	% select scatter data
    for i = 1: size(T, 1)
        temp = T(i, :);
        if temp.Year ~= 1111 & temp.Year ~= 2222 & temp.Year ~= 3333 & strcmp(temp.SiteName, siteTableName{k})
            scatterN = [scatterN, temp.NOx/baseNOx];
            scatterV = [scatterV, temp.VOC/baseVOC];
        end
    end

    % plot the isopleth
	n = [0: 0.01: 3.4];
	v = [0: 0.01: 4.9];
	[NOx,VOC] = meshgrid(n,v);
	O3 = base .^(a0 + a1 * NOx + a2 * VOC + a3 * NOx.^2 + a4 * NOx .* VOC + a5 * VOC.^2);
	LCoeff = [LCoeff; a0, a1, a2, a3, a4, a5];
	y_1 = (-a4/(2*a3)) * v - (a1/(2*a3));
	y_2 = ((a4-2*a5)/(a4-2*a3)) * v + (a2 - a1)/(2 * a3 - a4);
	figure(siteNum)
	% isopleth
	[c, h] = contourf(VOC, NOx, O3, 10, 'ShowText','on');
  	h.LevelList=round(h.LevelList,0);  %rounds levels to 0 decimal place
  	clabel(c,h);
	colormap(jet), 
	colorbar, hold on;
	% sensitive lines
	plot(v, y_1, 'w', 'LineWidth', 1.5)
	plot(v, y_2, 'r', 'LineWidth', 1.5)
	scatter(scatterV, scatterN, '*', 'b')
	xlim([0, 4.9])
	ylim([0, 3.4])
	title(siteName{siteNum});
	xlabel('VOC Emission'), ylabel('NO_x Emission')
	hold off;
	% save figures
	% saveas(gcf, strcat('/Users/zongrunli/Desktop/PycharmProjects/Ozone Interpolation/logModel Figures/', siteName{siteNum}, '.jpg'))
end