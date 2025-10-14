LCoeff = [];
T = readtable('/Users/zongrunli/Desktop/CMAQ-DDM/data/Observations.xlsx');
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

for siteNum = 1: 24
	% find D
	sheets = {'1percent', 'base', '1985', 'n10v100', 'n100v10', 'n200v10', 'n10v400', 'n300v10', 'n110v110', 'n90v90', 'n90v100', 'n110v100', '2001', '2011', '2028'};
	D = [];
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
    % D

    % change D to log scale
    for i = 1:15
        O = D(1, i);
        D(1, i) = log(D(1, i));
        D(2, i) = (1/O) * D(2, i);
        D(3, i) = (1/O) * D(3, i);
        D(4, i) = (1/O) * D(4, i) - (1/O^2) * D(2, i) * D(3, i);
        D(5, i) = (1/O) * D(5, i) - (1/O^2) * D(2, i) * D(2, i);
        D(6, i) = (1/O) * D(6, i) - (1/O^2) * D(3, i) * D(3, i);
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
        if temp.Year ~= 2222 & temp.Year ~= 3333 & temp.Year ~= 4444 & strcmp(temp.SiteName, siteTableName{k})
            scatterN = [scatterN, temp.NOx];
            scatterV = [scatterV, temp.VOC];
        end
    end

    % plot the isopleth
	n = [0: 0.01: 5];
	v = [0: 0.01: 6.25];
	[NOx,VOC] = meshgrid(n,v);
	O3 = exp(a0 + a1 * NOx + a2 * VOC + a3 * NOx.^2 + a4 * NOx .* VOC + a5 * VOC.^2);
	LCoeff = [LCoeff; a0, a1, a2, a3, a4, a5];
    y_1 = (-a4/(2*a3)) * v - (a1/(2*a3));
	y_2 = ((a4-2*a5)/(a4-2*a3)) * v + (a2 - a1)/(2 * a3 - a4);
	figure(siteNum)
	% isopleth
    contourf(VOC * baseVOC, NOx * baseNOx, O3, 100, 'LineStyle', 'none'), hold on
    [c, h] = contour(VOC * baseVOC, NOx * baseNOx, O3, 5,'k', 'LineWidth', 1.5, 'ShowText','on');
	h.LevelList=round(h.LevelList,0);  %rounds levels to 0 decimal place
	clabel(c,h, 'FontSize', 14);
	colormap(jet), 
	cb = colorbar;
    hold on;
    cb.Title.String = "ODV ppb";
    set(cb,'YTick',[0: 50: max(O3(:))])
	title(siteName{siteNum}, 'FontSize', 20);
	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1900])
    set(gca,'FontSize',15);
    % sensitive lines
	plot(v * baseVOC, y_1 * baseNOx, 'w--', 'LineWidth', 1.5)
	plot(v * baseVOC, y_2 * baseNOx, 'k--', 'LineWidth', 1.5)
    scatter(scatterV, scatterN, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    daspect([1, 1, 1]);
	% save figures
	% saveas(gcf, strcat('/Users/zongrunli/Desktop/DDM-CMAQ Updates/figures top 4/', siteName{siteNum}, ' Lmodel.jpg'))
end