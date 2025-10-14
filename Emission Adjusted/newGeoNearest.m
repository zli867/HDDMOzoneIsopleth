T = readtable('/Users/zongrunli/Desktop/CMAQ-DDM/data/Observations.xlsx');
load('/Users/zongrunli/Desktop/CMAQ-DDM/data/QCoeff.mat');
siteName = {'Azusa', 'Glendora', 'LA North Main', 'Reseda', 'Pico Rivera', 'Pomona', 'Pasadena', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Rubidoux', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernardino'};
siteTableName = {'Azusa', 'Glendora', 'LA North Main', 'Reseda', 'Pico Rivera', 'Pomona', 'Pasadena', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Riverside', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernadino'};
baseNOx = 395.748;
baseVOC = 401.412;
ODV = [];
E_N = [];
E_V = [];
t = [];
a0 = [];
a1 = [];
a2 = [];
a3 = [];
a4 = [];
a5 = [];
for siteNum = 1: length(siteTableName)
    % CHANGE SELECT THE ODV DATA
    for i = 1: size(T, 1)
        temp = T(i, :);
        if temp.Year ~= 2222 & temp.Year ~= 3333 & temp.Year ~= 4444 & strcmp(temp.SiteName, siteTableName{siteNum})
            NOx = temp.NOx/baseNOx;
            VOC = temp.VOC/baseVOC;
            ODVTemp = temp.ODV;
            t = [t; temp.Year-1975];
            ODV = [ODV;ODVTemp];
			E_N = [E_N; NOx];
			E_V = [E_V; VOC];
            a0 = [a0; QCoeff(siteNum, 1)];
            a1 = [a1; QCoeff(siteNum, 2)];
            a2 = [a2; QCoeff(siteNum, 3)];
            a3 = [a3; QCoeff(siteNum, 4)];
            a4 = [a4; QCoeff(siteNum, 5)];
            a5 = [a5; QCoeff(siteNum, 6)];
        end
    end
    % adjusted_N = (E_N .* (beta(1) + beta(2) .* t))
    % adjusted_V = (E_V .* (beta(3) + beta(4) .* t))
end
allYear = unique(t);
allYear = sort(allYear, 'ascend');
N_adj = [];
V_adj = [];
E_N_ori = [];
E_V_ori = [];
for i = 1: length(allYear)
    thisYear = allYear(i);
    index = find(t == thisYear);
    NOx_tmp = E_N(index);
    VOC_tmp = E_V(index);
    ODV_tmp = ODV(index);
    a0_tmp = a0(index);
    a1_tmp = a1(index);
    a2_tmp = a2(index);
    a3_tmp = a3(index);
    a4_tmp = a4(index);
    a5_tmp = a5(index);
    fun = @(beta)ODV_tmp - (a0_tmp + a1_tmp .* beta(1) + a2_tmp .* beta(2) + a3_tmp .* beta(1) .^2 + a4_tmp .* beta(1) .* beta(2) + a5_tmp .* beta(2) .^2);
    x0 = [mean(NOx_tmp), mean(VOC_tmp)];
    result = lsqnonlin(fun, x0);
%     result = lsqnonlin(fun, x0, zeros(1, 2));
    N_adj = [N_adj; result(1)];
    V_adj = [V_adj; result(2)];
    E_N_ori = [E_N_ori; mean(NOx_tmp)];
    E_V_ori = [E_V_ori; mean(VOC_tmp)];
end
figure;
plot(allYear + 1975, N_adj, 'r*--'), hold on
plot(allYear + 1975, E_N_ori, 'b+--')
title('E_{NO_x}'), legend('adjusted', 'original')
figure;
plot(allYear+1975, V_adj, 'r*--'), hold on
plot(allYear+1975, E_V_ori, 'b+--')
title('E_{VOC}'), legend('adjusted', 'original')