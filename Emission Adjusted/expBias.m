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
    % adjusted_N = E_N .* exp(beta(1) + beta(2) .* t
    % adjusted_V = E_V .* exp(beta(3) + beta(4) .* t
end

fun = @(beta)((ODV - (a0 + a1 .* (E_N .* exp(beta(1) + beta(2) .* t)) + a2 .* (E_V .* exp(beta(3) + beta(4) .* t)) + ...
             a3 .* (E_N .* exp(beta(1) + beta(2) .* t)) .^2 + a4 .* (E_N .* exp(beta(1) + beta(2) .* t)) .* (E_V .* exp(beta(3) + beta(4) .* t)) + ...
             a5 .* exp(E_V .* (beta(3) + beta(4) .* t)) .^2))./ODV).^2;
%              (exp(beta(1) + beta(2) .* t - 1)).^2 + (exp(beta(3) + beta(4) .* t - 1)) .^2;
x0 = [0, 0, 0, 0];
for i = 1:100
    result = lsqnonlin(fun, x0);
    x0 = result;
end
beta = result;
ODV_adjusted = (a0 + a1 .* (E_N .* exp(beta(1) + beta(2) .* t)) + a2 .* (E_V .* exp(beta(3) + beta(4) .* t)) + a3 .* (E_N .* exp(beta(1) + beta(2) .* t)) .^2 + a4 .* (E_N .* exp(beta(1) + beta(2) .* t)) .* (E_V .* exp(beta(3) + beta(4) .* t)) + a5 .* (E_V .* exp(beta(3) + beta(4) .* t)) .^2);

figure(1)
maxValue = max(max(ODV), max(ODV_adjusted)) + 1;
scatter(ODV, ODV_adjusted), hold on
plot([0:1:maxValue], [0:1:maxValue], '--')
xlabel("ODV"), 
ylabel("Adjusted ODV")
hold off

adjusted_N = (E_N .* exp(beta(1) + beta(2) .* t));
adjusted_V = (E_V .* exp(beta(3) + beta(4) .* t));

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

figure(4)
time = [min(t):1:max(t)];
ratio_N = exp(beta(1) + beta(2) .* time);
ratio_V = exp(beta(3) + beta(4) .* time);
scatter(time, ratio_N)
xlabel('Year from 1975')
ylabel('ratio')
title("NO_x ratio = adjusted NO_x/ NO_x")

figure(5)
scatter(time, ratio_V)
xlabel('Year from 1975')
ylabel('ratio')
title("VOC ratio = adjusted VOC/ VOC")

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
    % adjusted_N = E_N .* exp(beta(1) + beta(2) .* t
    % adjusted_V = E_V .* exp(beta(3) + beta(4) .* t
end

fun = @(beta)((ODV - (a0 + a1 .* (E_N .* exp(beta(1) + beta(2) .* t)) + a2 .* (E_V .* exp(beta(3) + beta(4) .* t)) + ...
             a3 .* (E_N .* exp(beta(1) + beta(2) .* t)) .^2 + a4 .* (E_N .* exp(beta(1) + beta(2) .* t)) .* (E_V .* exp(beta(3) + beta(4) .* t)) + ...
             a5 .* exp(E_V .* (beta(3) + beta(4) .* t)) .^2))./ODV).^2;
             ((exp(beta(1) + beta(2) .* t) - 1)./exp(beta(1) + beta(2) .* t)).^2 + ((exp(beta(3) + beta(4) .* t) - 1)./exp(beta(3) + beta(4) .* t)) .^2;
x0 = [0, 0, 0, 0];
for i = 1:100
    result = lsqnonlin(fun, x0);
    x0 = result;
end
beta = result;
ODV_adjusted = (a0 + a1 .* (E_N .* exp(beta(1) + beta(2) .* t)) + a2 .* (E_V .* exp(beta(3) + beta(4) .* t)) + a3 .* (E_N .* exp(beta(1) + beta(2) .* t)) .^2 + a4 .* (E_N .* exp(beta(1) + beta(2) .* t)) .* (E_V .* exp(beta(3) + beta(4) .* t)) + a5 .* (E_V .* exp(beta(3) + beta(4) .* t)) .^2);

figure(1)
maxValue = max(max(ODV), max(ODV_adjusted)) + 1;
scatter(ODV, ODV_adjusted), hold on
plot([0:1:maxValue], [0:1:maxValue], '--')
xlabel("ODV"), 
ylabel("Adjusted ODV")
hold off

adjusted_N = (E_N .* exp(beta(1) + beta(2) .* t));
adjusted_V = (E_V .* exp(beta(3) + beta(4) .* t));

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

figure(4)
time = [min(t):1:max(t)];
ratio_N = exp(beta(1) + beta(2) .* time);
ratio_V = exp(beta(3) + beta(4) .* time);
scatter(time, ratio_N)
xlabel('Year from 1975')
ylabel('ratio')
title("NO_x ratio = adjusted NO_x/ NO_x")

figure(5)
scatter(time, ratio_V)
xlabel('Year from 1975')
ylabel('ratio')
title("VOC ratio = adjusted VOC/ VOC")

allYear = unique(t);
E_N_ori = [];
E_V_ori = [];
for i = 1: length(allYear)
    index = find(t == allYear(i));
    E_N_ori = [E_N_ori; mean(E_N(index))];
    E_V_ori = [E_V_ori; mean(E_V(index))];
end
adj_N = (E_N_ori .* exp(beta(1) + beta(2) .* allYear));
adj_V = (E_V_ori .* exp(beta(3) + beta(4) .* allYear));
figure;
plot(allYear + 1975, adj_N, 'r*--'), hold on
plot(allYear + 1975, E_N_ori, 'b+--')
title('E_{NO_x}'), legend('adjusted', 'original')
figure;
plot(allYear+1975, adj_V, 'r*--'), hold on
plot(allYear+1975, E_V_ori, 'b+--')
title('E_{VOC}'), legend('adjusted', 'original')