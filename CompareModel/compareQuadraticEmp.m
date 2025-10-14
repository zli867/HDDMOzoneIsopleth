% regression model
load('/Volumes/Shield/CMAQ-DDM/data/QCoeff.mat');
T = readtable('/Volumes/Shield/CMAQ-DDM/data/Observations.xlsx');
siteName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Rubidoux', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernardino'};
siteTableName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Riverside', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernadino'};
baseNOx = 395.748;
baseVOC = 401.412;

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
    
    qua_coef = QCoeff(siteNum, :);
	n = [0: 10: 1600];
	v = [0: 10: 2100];
	[N,V] = meshgrid(n,v);
    O3_emp = 10.^(b(1) + b(2) * N + b(3) * V + b(4) * N.^2 + b(5) * N.*V + b(6) * V.^2);
    O3_N_emp = O3_emp .* log(10) .* (b(2) + 2* b(4) * N + b(5) * V);
    O3_V_emp = O3_emp .* log(10) .* (b(3) + b(5) * N + 2 * b(6) * V);
    y_1_emp = -(b(5)/(2*b(4))) * v - (b(2)/(2*b(4)));
    y_2_emp = (2* b(6) - b(5))/(2*b(4) - b(5)) * v + (b(3) - b(2))/(2*b(4) - b(5));
    
    NOx_norm = N./baseNOx;
    VOC_norm = V./baseVOC;
	O3 = qua_coef(1) + qua_coef(2) * NOx_norm + qua_coef(3) * VOC_norm + qua_coef(4) * NOx_norm.^2 + qua_coef(5) * NOx_norm .* VOC_norm + qua_coef(6) * VOC_norm.^2;    
    O3_N = (qua_coef(2) + 2 * qua_coef(4) * NOx_norm + qua_coef(5) * VOC_norm)/baseNOx;
    O3_V = (qua_coef(3) + qua_coef(5) * NOx_norm + 2 * qua_coef(6) * VOC_norm)/baseVOC;
	y_1 = -(baseNOx/baseVOC) * (qua_coef(5)/(2*qua_coef(4))) * v - qua_coef(2) * baseNOx/(2*qua_coef(4));
	y_2 = (((2 * qua_coef(6)/baseVOC^2) - (qua_coef(5)/(baseVOC * baseNOx)))* v + (qua_coef(3)/baseVOC - qua_coef(2)/baseNOx)) /(2* qua_coef(4)/baseNOx^2 - qua_coef(5)/(baseNOx*baseVOC));
    
    diff_O3 = O3_emp - O3;
    diff_O3_N = O3_N_emp - O3_N;
    diff_O3_V = O3_V_emp - O3_V;
    
%     O3_range = max(max(O3(:)), max(O3_emp(:)));
%     O3_N_range = [min(min(O3_N(:)), min(O3_N_emp(:))), max(max(O3_N(:)), max(O3_N_emp(:)))];
%     O3_V_range = [min(min(O3_V(:)), min(O3_V_emp(:))), max(max(O3_V(:)), max(O3_V_emp(:)))];
    O3_range = 300;
    O3_N_range = [-0.3, 0.3];
    O3_V_range = [-0.05, 0.2];
    
    % Plot Figures
	f = figure(siteNum);
    f.Position = [0 0 1200 800];
	% isopleth
    % Plot for Empirical model
    subplot(3, 3, 1)
    contourf(V, N, O3_emp, 100, 'LineStyle', 'none'), hold on
%     [c, h] = contour(V, N, O3_emp, 5,'k', 'LineWidth', 1.5, 'ShowText','on');
    [c, h] = contour(V, N, O3_emp, [0: 30: O3_range], 'k', 'LineWidth', 1.5, 'ShowText','on');
	h.LevelList=round(h.LevelList,0);  %rounds levels to 0 decimal place
	clabel(c,h, 'FontSize', 14);
	colormap(jet), 
	cb = colorbar;
    hold on;
    cb.Title.String = "ODV";
    caxis([0 O3_range])
%     set(cb,'YTick',[0: 30: O3_range])
	title(['Empirical Model (', siteName{siteNum}, ')'], 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    % sensitive lines
	plot(v, y_1_emp, 'w--', 'LineWidth', 1.5)
	plot(v, y_2_emp, 'k--', 'LineWidth', 1.5)
    scatter(VOC, NOx, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    daspect([1, 1, 1]);
    % dO3/dNOx
    subplot(3, 3, 4)
    contourf(V, N, O3_N_emp, 100, 'LineStyle', 'none'), hold on
    [c, h] = contour(V, N, O3_N_emp, [O3_N_range(1):0.025:O3_N_range(2)],'k', 'LineWidth', 1.5, 'ShowText','on');
	h.LevelList=round(h.LevelList,3);  %rounds levels to 0 decimal place
	clabel(c,h, 'FontSize', 14);
	colormap(jet), 
	cb = colorbar;
    caxis([O3_N_range(1) O3_N_range(2)])
    hold on;
	title('dO_3/dE_{NO_x} Empirical', 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    scatter(VOC, NOx, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    daspect([1, 1, 1]);
    % dO3/dVOC
    subplot(3, 3, 7)
    contourf(V, N, O3_V_emp, 100, 'LineStyle', 'none'), hold on
    [c, h] = contour(V, N, O3_V_emp, [O3_V_range(1):0.025:O3_V_range(2)],'k', 'LineWidth', 1.5, 'ShowText','on');
	h.LevelList=round(h.LevelList,3);  %rounds levels to 0 decimal place
	clabel(c,h, 'FontSize', 14);
	colormap(jet), 
	cb = colorbar;
    caxis([O3_V_range(1) O3_V_range(2)])
    hold on;
	title('dO_3/dE_{VOC} Empirical', 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    scatter(VOC, NOx, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    daspect([1, 1, 1]);
    
    % Plot for quadratic model
    subplot(3, 3, 2)
    contourf(V, N, O3, 100, 'LineStyle', 'none'), hold on
    [c, h] = contour(V, N, O3, [0: 30: O3_range],'k', 'LineWidth', 1.5, 'ShowText','on');
	h.LevelList=round(h.LevelList,0);  %rounds levels to 0 decimal place
	clabel(c,h, 'FontSize', 14);
	colormap(jet), 
	cb = colorbar;
    hold on;
    cb.Title.String = "Peak O_3";
    caxis([0 O3_range])
%     set(cb,'YTick',[0: 30: O3_range])
	title('Quadratic Model', 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    % sensitive lines
	plot(v, y_1, 'w--', 'LineWidth', 1.5)
	plot(v, y_2, 'k--', 'LineWidth', 1.5)
    scatter(VOC, NOx, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    daspect([1, 1, 1]);
    % dO3/dNOx
    subplot(3, 3, 5)
    contourf(V, N, O3_N, 100, 'LineStyle', 'none'), hold on
    [c, h] = contour(V, N, O3_N, [O3_N_range(1):0.025:O3_N_range(2)],'k', 'LineWidth', 1.5, 'ShowText','on');
	h.LevelList=round(h.LevelList,3);  %rounds levels to 0 decimal place
	clabel(c,h, 'FontSize', 14);
	colormap(jet), 
	cb = colorbar;
    caxis([O3_N_range(1) O3_N_range(2)])
    hold on;
	title('dO_3/dE_{NO_x} Quadratic', 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    scatter(VOC, NOx, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    daspect([1, 1, 1]);
    % dO3/dVOC
    subplot(3, 3, 8)
    contourf(V, N, O3_V, 100, 'LineStyle', 'none'), hold on
    [c, h] = contour(V, N, O3_V, [O3_V_range(1):0.025:O3_V_range(2)],'k', 'LineWidth', 1.5, 'ShowText','on');
	h.LevelList=round(h.LevelList,3);  %rounds levels to 0 decimal place
	clabel(c,h, 'FontSize', 14);
	colormap(jet), 
	cb = colorbar;
    caxis([O3_V_range(1) O3_V_range(2)])
    hold on;
	title('dO_3/dE_{VOC} Quadratic', 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    xlabel('VOC Emission (tons/day)', 'FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    scatter(VOC, NOx, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    daspect([1, 1, 1]);
    
    % Plot for differences
    subplot(3, 3, 3)
    contourf(V, N, diff_O3, 100, 'LineStyle', 'none'), hold on
	colormap(jet), 
	cb = colorbar;
    hold on;
    cb.Title.String = "diff";
	title('Diff (Emp - Quadratic)', 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    scatter(VOC, NOx, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    caxis([-40, 40])
    daspect([1, 1, 1]);
    % dO3/dNOx
    subplot(3, 3, 6)
    contourf(V, N, diff_O3_N, 100, 'LineStyle', 'none'), hold on
	colormap(jet), 
	cb = colorbar;
    hold on;
	title('Diff (Emp - Quadratic)', 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    scatter(VOC, NOx, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    caxis([-0.2, 0.2])
    daspect([1, 1, 1]);
    % dO3/dVOC
    subplot(3, 3, 9)
    contourf(V, N, diff_O3_V, 100, 'LineStyle', 'none'), hold on
	colormap(jet), 
	cb = colorbar;
    hold on;
	title('Diff (Emp - Quadratic)', 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    scatter(VOC, NOx, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    caxis([-0.1, 0.1])
    daspect([1, 1, 1]);
    
    % save figures
	saveas(gcf, strcat('/Users/zongrunli/Desktop/CMAQ-DDM/Figure/', siteName{siteNum}, '_compare_quadratic_emp.png'))
end