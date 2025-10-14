% % Load data
load('/Volumes/Shield/CMAQ-DDM/data/LCoeff.mat')
load('/Volumes/Shield/CMAQ-DDM/data/QCoeff.mat')
T = readtable('/Volumes/Shield/CMAQ-DDM/data/Observations.xlsx');
siteName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Rubidoux', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernardino'};
siteTableName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Riverside', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernadino'};
baseNOx = 395.748;
baseVOC = 401.412;
for siteNum = 1: 24
	scatterN = [];
    scatterV = [];
	% select scatter data
    for i = 1: size(T, 1)
        temp = T(i, :);
        if temp.Year ~= 2222 & temp.Year ~= 3333 & temp.Year ~= 4444 & strcmp(temp.SiteName, siteTableName{siteNum})
            scatterN = [scatterN, temp.NOx];
            scatterV = [scatterV, temp.VOC];
        end
    end
    
    qua_coef = QCoeff(siteNum, :);
    log_coef = LCoeff(siteNum, :);
    
    
    ratio = baseVOC/baseNOx;
	n = [0: 10: 1600];
	v = [0: 10: 2100];
	[NOx,VOC] = meshgrid(n,v);
    NOx_norm = NOx./baseNOx;
    VOC_norm = VOC./baseVOC;
	O3 = qua_coef(1) + qua_coef(2) * NOx_norm + qua_coef(3) * VOC_norm + qua_coef(4) * NOx_norm.^2 + qua_coef(5) * NOx_norm .* VOC_norm + qua_coef(6) * VOC_norm.^2;    
    O3_N = (qua_coef(2) + 2 * qua_coef(4) * NOx_norm + qua_coef(5) * VOC_norm)/baseNOx;
    O3_V = (qua_coef(3) + qua_coef(5) * NOx_norm + 2 * qua_coef(6) * VOC_norm)/baseVOC;
	y_1 = -(baseNOx/baseVOC) * (qua_coef(5)/(2*qua_coef(4))) * v - qua_coef(2) * baseNOx/(2*qua_coef(4));
	y_2 = (((2 * qua_coef(6)/baseVOC^2) - (qua_coef(5)/(baseVOC * baseNOx)))* v + (qua_coef(3)/baseVOC - qua_coef(2)/baseNOx)) /(2* qua_coef(4)/baseNOx^2 - qua_coef(5)/(baseNOx*baseVOC));
    
    O3_L = exp(log_coef(1) + log_coef(2) * NOx_norm + log_coef(3) * VOC_norm + log_coef(4) * NOx_norm.^2 + log_coef(5) * NOx_norm .* VOC_norm + log_coef(6) * VOC_norm.^2);
    O3_N_L = O3_L .* (log_coef(2) + 2 * log_coef(4) * NOx_norm + log_coef(5) * VOC_norm)/baseNOx;
    O3_V_L = O3_L .* (log_coef(3) + log_coef(5) * NOx_norm + 2 * log_coef(6) * VOC_norm)/baseVOC;
	y_1_L = -(baseNOx/baseVOC) * (log_coef(5)/(2*log_coef(4))) * v - log_coef(2) * baseNOx/(2*log_coef(4));
	y_2_L = (((2 * log_coef(6)/baseVOC^2) - (log_coef(5)/(baseVOC * baseNOx)))* v + (log_coef(3)/baseVOC - log_coef(2)/baseNOx)) /(2* log_coef(4)/baseNOx^2 - log_coef(5)/(baseNOx*baseVOC));
    
    diff_O3 = O3-O3_L;
    diff_O3_N = O3_N - O3_N_L;
    diff_O3_V = O3_V - O3_V_L;
    
%     O3_range = max(max(O3(:)), max(O3_L(:)));
%     O3_N_range = [min(min(O3_N(:)), min(O3_N_L(:))), max(max(O3_N(:)), max(O3_N_L(:)))];
%     O3_V_range = [min(min(O3_V(:)), min(O3_V_L(:))), max(max(O3_V(:)), max(O3_V_L(:)))];
    O3_range = 300;
    O3_N_range = [-0.4, 0.4];
    O3_V_range = [-0.05, 0.2];
    
    % Plot Figures
    % Plot for quadratic model
	f = figure(siteNum);
    f.Position = [0 0 1200 800];
	% isopleth
    subplot(3, 3, 1)
    contourf(VOC, NOx, O3, 100, 'LineStyle', 'none'), hold on
    [c, h] = contour(VOC, NOx, O3, [0: 30: O3_range],'k', 'LineWidth', 1.5, 'ShowText','on');
	h.LevelList=round(h.LevelList,0);  %rounds levels to 0 decimal place
	clabel(c,h, 'FontSize', 14);
	colormap(jet), 
	cb = colorbar;
    hold on;
    cb.Title.String = "Peak O_3";
    caxis([0 O3_range])
%     set(cb,'YTick',[0: 30: O3_range])
	title(['Quadratic Model (', siteName{siteNum}, ')'], 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    % sensitive lines
	plot(v, y_1, 'w--', 'LineWidth', 1.5)
	plot(v, y_2, 'k--', 'LineWidth', 1.5)
    scatter(scatterV, scatterN, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    daspect([1, 1, 1]);
    % dO3/dNOx
    subplot(3, 3, 4)
    contourf(VOC, NOx, O3_N, 100, 'LineStyle', 'none'), hold on
    [c, h] = contour(VOC, NOx, O3_N, [O3_N_range(1):0.025:O3_N_range(2)],'k', 'LineWidth', 1.5, 'ShowText','on');
	h.LevelList=round(h.LevelList,3);  %rounds levels to 0 decimal place
	clabel(c,h, 'FontSize', 14);
	colormap(jet), 
	cb = colorbar;
    caxis([O3_N_range(1) O3_N_range(2)])
    hold on;
	title('dO_3/dE_{NO_x} Quadratic', 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    scatter(scatterV, scatterN, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    daspect([1, 1, 1]);
    % dO3/dVOC
    subplot(3, 3, 7)
    contourf(VOC, NOx, O3_V, 100, 'LineStyle', 'none'), hold on
    [c, h] = contour(VOC, NOx, O3_V, [O3_V_range(1):0.025:O3_V_range(2)],'k', 'LineWidth', 1.5, 'ShowText','on');
	h.LevelList=round(h.LevelList,3);  %rounds levels to 0 decimal place
	clabel(c,h, 'FontSize', 14);
	colormap(jet), 
	cb = colorbar;
    caxis([O3_V_range(1) O3_V_range(2)])
    hold on;
	title('dO_3/dE_{VOC} Quadratic', 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    scatter(scatterV, scatterN, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    daspect([1, 1, 1]);
    
    % Plot for log quadratic model
    subplot(3, 3, 2)
    contourf(VOC, NOx, O3_L, 100, 'LineStyle', 'none'), hold on
    [c, h] = contour(VOC, NOx, O3_L,[0: 30: O3_range], 'k', 'LineWidth', 1.5, 'ShowText','on');
	h.LevelList=round(h.LevelList,0);  %rounds levels to 0 decimal place
	clabel(c,h, 'FontSize', 14);
	colormap(jet), 
	cb = colorbar;
    hold on;
    cb.Title.String = "Peak O_3";
    caxis([0 O3_range])
%     set(cb,'YTick',[0: 30: O3_range])
	title('Log Quadratic Model', 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    % sensitive lines
	plot(v, y_1_L, 'w--', 'LineWidth', 1.5)
	plot(v, y_2_L, 'k--', 'LineWidth', 1.5)
    scatter(scatterV, scatterN, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    daspect([1, 1, 1]);
    % dO3/dNOx
    subplot(3, 3, 5)
    contourf(VOC, NOx, O3_N_L, 100, 'LineStyle', 'none'), hold on
    [c, h] = contour(VOC, NOx, O3_N_L, [O3_N_range(1):0.025:O3_N_range(2)],'k', 'LineWidth', 1.5, 'ShowText','on');
	h.LevelList=round(h.LevelList,3);  %rounds levels to 0 decimal place
	clabel(c,h, 'FontSize', 14);
	colormap(jet), 
	cb = colorbar;
    caxis([O3_N_range(1) O3_N_range(2)])
    hold on;
	title('dO_3/dE_{NO_x} Log', 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    scatter(scatterV, scatterN, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    daspect([1, 1, 1]);
    % dO3/dVOC
    subplot(3, 3, 8)
    contourf(VOC, NOx, O3_V_L, 100, 'LineStyle', 'none'), hold on
    [c, h] = contour(VOC, NOx, O3_V_L, [O3_V_range(1):0.025:O3_V_range(2)],'k', 'LineWidth', 1.5, 'ShowText','on');
	h.LevelList=round(h.LevelList,3);  %rounds levels to 0 decimal place
	clabel(c,h, 'FontSize', 14);
	colormap(jet), 
	cb = colorbar;
    caxis([O3_V_range(1) O3_V_range(2)])
    hold on;
	title('dO_3/dE_{VOC} Log', 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    xlabel('VOC Emission (tons/day)', 'FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    scatter(scatterV, scatterN, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    daspect([1, 1, 1]);
    
    % Plot for differences
    subplot(3, 3, 3)
    contourf(VOC, NOx, diff_O3, 100, 'LineStyle', 'none'), hold on
	colormap(jet), 
	cb = colorbar;
    hold on;
    cb.Title.String = "diff";
	title('Diff (Quadratic - Log Quadratic)', 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    scatter(scatterV, scatterN, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    caxis([-80, 10])
    daspect([1, 1, 1]);
    % dO3/dNOx
    subplot(3, 3, 6)
    contourf(VOC, NOx, diff_O3_N, 100, 'LineStyle', 'none'), hold on
	colormap(jet), 
	cb = colorbar;
    hold on;
	title('Diff (Quadratic - Log Quadratic)', 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    scatter(scatterV, scatterN, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    caxis([-0.2, 0.2])
    daspect([1, 1, 1]);
    % dO3/dVOC
    subplot(3, 3, 9)
    contourf(VOC, NOx, diff_O3_V, 100, 'LineStyle', 'none'), hold on
	colormap(jet), 
	cb = colorbar;
    hold on;
	title('Diff (Quadratic - Log Quadratic)', 'FontSize', 20);
% 	xlabel('VOC Emission (tons/day)', 'FontSize', 16), ylabel('NO_x Emission (tons/day)','FontSize', 16),
    set(gca,'XTick',[0 : 500 : 2500]), set(gca,'YTick',[0 : 500 : 1800])
    set(gca,'FontSize',13);
    scatter(scatterV, scatterN, '*', 'r')
	hold off;
    xlim([0, 2100])
    ylim([0, 1600])
    caxis([-0.06, 0.06])
    daspect([1, 1, 1]);
    
	% save figures
	saveas(gcf, strcat('/Volumes/Shield/CMAQ-DDM/Figure/', siteName{siteNum}, '_compare_log_quadratic.png'))
end