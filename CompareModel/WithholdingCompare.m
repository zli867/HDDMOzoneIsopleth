load('/Volumes/Shield/CMAQ-DDM//data/QCoeff.mat');
load('/Volumes/Shield/CMAQ-DDM//data/LCoeff.mat');
N = [0.01, 1, 3.9, 0.1, 1, 2, 0.1, 3, 1.1, 0.9, 0.9, 1.1, 2.3, 1.3, 0.8];
V = [0.01, 1, 4.9, 1, 0.1, 0.1, 4, 0.1, 1.1, 0.9, 1, 1, 3.6, 1.5, 0.6];
x1 = [2: 13: 301];
x2 = [14: 13: 313];
baseNOx = 395.748;
baseVOC = 401.412;
site = {};
siteName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Rubidoux', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernardino'};
for k = 1: 24
	site{k} = ['D', int2str(x1(k)), ':', 'I', int2str(x2(k))];
end

for siteNum = 1: 24
    % Quadratic model
    % find D The original Data
	sheets = {'1percent', 'base', '1985', 'n10v100', 'n100v10', 'n200v10', 'n10v400', 'n300v10', 'n110v110', 'n90v90', 'n90v100', 'n110v100', '2001', '2011', '2028'};
    D = [];
    CoeffWithholding = [];
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

    for removeIndex = 1:15
        data = [1: 15];
        data(removeIndex) = [];
        selectData = data;
    
        % started generate the mathematical model
        % b
        b = [];
        for i = selectData
            b = [b; D(:, i)];
        end

        % A
        A = [];
        for i = selectData
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
        CoeffWithholding = [CoeffWithholding; a0, a1, a2, a3, a4, a5];
    end

    % calculate mean of bias and bias' standard deviation
	n = [0: 0.01: 4.5];
	v = [0: 0.01: 6.25];
    [NOx,VOC] = meshgrid(n,v);
    OriginalO3 = QCoeff(siteNum, 1) + QCoeff(siteNum, 2) * NOx + QCoeff(siteNum, 3) * VOC + QCoeff(siteNum, 4) * NOx.^2 + QCoeff(siteNum, 5) * NOx .* VOC + QCoeff(siteNum, 6) * VOC.^2;
    avg = 0;
    std = 0;
    for t = 1:15
        O3 = CoeffWithholding(t, 1) + CoeffWithholding(t, 2) * NOx + CoeffWithholding(t, 3) * VOC + CoeffWithholding(t, 4) * NOx.^2 + CoeffWithholding(t, 5) * NOx .* VOC + CoeffWithholding(t, 6) * VOC.^2; 
        avg = avg + (O3 - OriginalO3);
    end
    avg = avg/15;

    for t = 1:15
        O3 = CoeffWithholding(t, 1) + CoeffWithholding(t, 2) * NOx + CoeffWithholding(t, 3) * VOC + CoeffWithholding(t, 4) * NOx.^2 + CoeffWithholding(t, 5) * NOx .* VOC + CoeffWithholding(t, 6) * VOC.^2;
        diff = O3 - OriginalO3;
        std = std + (diff - avg).^2;
    end
    std = sqrt(std./15);
    
    Q_avg = avg;
    Q_std = std;
    
    % log quadratic
    D = [];
    CoeffWithholding = [];
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

    % transform D
    for i = 1:15
        O = D(1, i);
        D(1, i) = log(D(1, i));
        D(2, i) = (1/O) * D(2, i);
        D(3, i) = (1/O) * D(3, i);
        D(4, i) = (1/O) * D(4, i) - (1/O^2) * D(2, i) * D(3, i);
        D(5, i) = (1/O) * D(5, i) - (1/O^2) * D(2, i) * D(2, i);
        D(6, i) = (1/O) * D(6, i) - (1/O^2) * D(3, i) * D(3, i);
    end

    for removeIndex = 1:15
        data = [1: 15];
        data(removeIndex) = [];
        selectData = data;
    
        % started generate the mathematical model
        % b
        b = [];
        for i = selectData
            b = [b; D(:, i)];
        end

        % A
        A = [];
        for i = selectData
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
        CoeffWithholding = [CoeffWithholding; a0, a1, a2, a3, a4, a5];
    end

    % calculate mean of bias and bias' standard deviation
	n = [0: 0.01: 4.5];
	v = [0: 0.01: 6.25];
    [NOx,VOC] = meshgrid(n,v);
    OriginalO3 = exp(LCoeff(siteNum, 1) + LCoeff(siteNum, 2) * NOx + LCoeff(siteNum, 3) * VOC + LCoeff(siteNum, 4) * NOx.^2 + LCoeff(siteNum, 5) * NOx .* VOC + LCoeff(siteNum, 6) * VOC.^2);
    avg = 0;
    std = 0;
    for t = 1:15
        O3 = exp(CoeffWithholding(t, 1) + CoeffWithholding(t, 2) * NOx + CoeffWithholding(t, 3) * VOC + CoeffWithholding(t, 4) * NOx.^2 + CoeffWithholding(t, 5) * NOx .* VOC + CoeffWithholding(t, 6) * VOC.^2); 
        avg = avg + (O3 - OriginalO3);
    end
    avg = avg/15;

    for t = 1:15
        O3 = exp(CoeffWithholding(t, 1) + CoeffWithholding(t, 2) * NOx + CoeffWithholding(t, 3) * VOC + CoeffWithholding(t, 4) * NOx.^2 + CoeffWithholding(t, 5) * NOx .* VOC + CoeffWithholding(t, 6) * VOC.^2);
        diff = O3 - OriginalO3;
        std = std + (diff - avg).^2;
    end
    std = sqrt(std./15);
    
    L_avg = avg;
    L_std = std;
    
    mean_range = [min(min(Q_avg(:)), min(L_avg(:))), max(max(Q_avg(:)), max(L_avg(:)))];
    std_range = [min(min(Q_std(:)), min(L_std(:))), max(max(Q_std(:)), max(L_std(:)))];
    % plot figures
    % mean differences
    f = figure(siteNum);
    f.Position = [0 0 900 600];
    subplot(2, 2, 1)
    contourf(VOC * baseVOC, NOx * baseNOx, Q_avg, 100, 'LineStyle', 'none'), hold on
    scatter(V * baseVOC, N * baseNOx, '*', 'w')
    [c, h] = contour(VOC * baseVOC, NOx * baseNOx, Q_avg, 6, 'k', 'LineWidth', 1.5, 'ShowText','on');
    h.LevelList=round(h.LevelList,0);  %rounds levels to 0 decimal place
	clabel(c,h, 'FontSize', 14);
	colormap(jet)
    cb = colorbar;
    cb.Title.String = "Mean Bias";
    cb.FontSize = 15;
    caxis([mean_range(1), mean_range(2)])
	title(['Quadratic Model (', siteName{siteNum}, ')'], 'FontSize', 20);
	ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2100]), set(gca,'YTick',[0 : 500 : 1600])
    set(gca,'FontSize',15);
	hold off;
    daspect([1, 1, 1]);

    % std
    subplot(2, 2, 3)
    contourf(VOC * baseVOC, NOx * baseNOx, Q_std, 100, 'LineStyle', 'none'), hold on
    scatter(V * baseVOC, N * baseNOx, '*', 'w')
    [c, h] = contour(VOC* baseVOC, NOx* baseNOx, Q_std, 6, 'k', 'LineWidth', 1.5, 'ShowText','on');
    h.LevelList=round(h.LevelList,0);  %rounds levels to 0 decimal place
	clabel(c,h, 'FontSize', 14);
	colormap(jet)
    cb = colorbar;
	cb.Title.String = "Std Bias";
    cb.FontSize = 15;
	caxis([std_range(1), std_range(2)])
	title('Quadratic Model', 'FontSize', 20);
	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2100]), set(gca,'YTick',[0 : 500 : 1600])
    set(gca,'FontSize',15);
	hold off;
    daspect([1, 1, 1]);
    
    % mean differences
    figure(siteNum)
    subplot(2, 2, 2)
    contourf(VOC * baseVOC, NOx * baseNOx, L_avg, 100, 'LineStyle', 'none'), hold on
    scatter(V * baseVOC, N * baseNOx, '*', 'w')
    [c, h] = contour(VOC * baseVOC, NOx * baseNOx, L_avg, 6, 'k', 'LineWidth', 1.5, 'ShowText','on');
    h.LevelList=round(h.LevelList,0);  %rounds levels to 0 decimal place
	clabel(c,h, 'FontSize', 14);
	colormap(jet)
    cb = colorbar;
    cb.Title.String = "Mean Bias";
    cb.FontSize = 15;
    caxis([mean_range(1), mean_range(2)])
	title('Log Quadratic Model', 'FontSize', 20);
    set(gca,'XTick',[0 : 500 : 2100]), set(gca,'YTick',[0 : 500 : 1600])
    set(gca,'FontSize',15);
	hold off;
    daspect([1, 1, 1]);

    % std
    subplot(2, 2, 4)
    contourf(VOC * baseVOC, NOx * baseNOx, L_std, 100, 'LineStyle', 'none'), hold on
    scatter(V * baseVOC, N * baseNOx, '*', 'w')
    [c, h] = contour(VOC* baseVOC, NOx* baseNOx, L_std, 6, 'k', 'LineWidth', 1.5, 'ShowText','on');
    h.LevelList=round(h.LevelList,0);  %rounds levels to 0 decimal place
	clabel(c,h, 'FontSize', 14);
	colormap(jet)
    cb = colorbar;
	cb.Title.String = "Std Bias";
    cb.FontSize = 15;
	caxis([std_range(1), std_range(2)])
	title('Log Quadratic Model', 'FontSize', 20);
	xlabel('VOC Emission (tons/day)', 'FontSize', 16)
    set(gca,'XTick',[0 : 500 : 2100]), set(gca,'YTick',[0 : 500 : 1600])
    set(gca,'FontSize',15);
	hold off;
    daspect([1, 1, 1]);
    
    % save figures
	saveas(gcf, strcat('/Volumes/Shield/CMAQ-DDM/Figure/', siteName{siteNum}, '_data_withholding.png'))
end