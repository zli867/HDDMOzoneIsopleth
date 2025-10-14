T = readtable('/Users/zongrunli/Desktop/CMAQ-DDM/data/Observations.xlsx');
load('/Users/zongrunli/Desktop/CMAQ-DDM/data/QCoeff.mat');
siteName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Rubidoux', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernardino'};
siteTableName = {'Azusa', 'Glendora', 'West LA', 'LA North Main', 'Reseda', 'Burbank', 'Pico Rivera', 'Pomona', 'Pasadena', 'Long Beach', 'LAX', 'Santa Clarita', 'Anaheim', 'Mission Viejo', 'La Habra', 'Banning', 'Perris', 'Riverside', 'Lake Elsinore', 'Crestline', 'Upland', 'Fontana', 'Redlands', 'San Bernadino'};
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
for siteNum = 1: 24
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

adjusted_N = [];
adjusted_V = [];
E_N_test = [];
E_V_test = [];
t_test = [];
for i = 1: size(ODV, 1)
    % Find the extremum
    syms x y;
    e1 = ODV(i) - (a0(i) + a1(i) * y + a2(i) * x + a3(i) * y^2 + a4(i) * x * y + a5(i) * x^2);
    e2 = (a2(i) + a4(i)*y + 2*a5(i)*x) * (E_N(i)-y) - (E_V(i) - x) * (a1(i) + 2*a3(i) * y + a4(i)*x);
    [x0, y0] = solve(e1, e2, x, y);
    res = double([x0, y0]);
    % find the nearest and real points
    flag = isreal(res(1, :)) || isreal(res(2, :)) || isreal(res(3, :)) || isreal(res(4, :));
    if flag
        dis = intmax;
        index = 1;
        for j = 1:4
            if(isreal(res(j, :)))
                disTmp = (res(j, 1) - E_V(i))^2 + (res(j, 2) - E_N(i))^2;
                if disTmp < dis
                    dis = disTmp;
                    index = j;
                end
            end
        end
        t_test = [t_test; t(i)];
        adjusted_N = [adjusted_N; res(index, 2)];
        adjusted_V = [adjusted_V; res(index, 1)];
        E_N_test = [E_N_test; E_N(i)];
        E_V_test = [E_V_test; E_V(i)];
    else
        continue
    end
end
% plot the original figures
figure;
scatter(E_N_test, adjusted_N)
xlabel('E_{NO_x}'), ylabel('Adjusted E_{NO_x}')

figure;
scatter(E_V_test, adjusted_V)
xlabel('E_{VOC}'), ylabel('Adjusted E_{VOC}')


% Average the emissions by years
N_changed = [];
V_changed = [];
N_original = [];
V_original = [];
time = unique(t_test);
time = sort(time);
for i = 1: size(time, 1)
    index = find(t_test == time(i));
    N_changed_tmp = mean(adjusted_N(index));
    V_changed_tmp = mean(adjusted_V(index));
    N_original = [N_original; mean(E_N_test(index))];
    V_original = [V_original; mean(E_V_test(index))];
    N_changed = [N_changed; N_changed_tmp];
    V_changed = [V_changed; V_changed_tmp];
end

ratio_N = N_changed ./ N_original;
ratio_V = V_changed ./ V_original;
% do not have solution in 1976 since 1976's maximum is lower than
% observation
figure
maxvalue = max(max(V_original(:)), max(V_changed(:)));
scatter(V_original, V_changed), hold on
plot([0: 0.01: maxvalue], [0: 0.01: maxvalue], '--')
xlabel('VOC Emission'), ylabel('Adjusted VOC Emission')

figure
maxvalue = max(max(N_original(:)), max(N_changed(:)));
scatter(N_original, N_changed), hold on
plot([0: 0.01: maxvalue], [0: 0.01: maxvalue], '--')
xlabel('NO_x Emission'), ylabel('Adjusted NO_x Emission')

figure
plot(time + 1975, N_original, '--*'), hold on
plot(time + 1975, N_changed, '--*')
xlabel('year'), ylabel('Normalized E_{NO_x}'), legend('Original Emission', 'Adjusted Emission')

figure
plot(time + 1975, V_original, '--*'), hold on
plot(time + 1975, V_changed, '--*')
xlabel('year'), ylabel('Normalized E_{VOC}'), legend('Original Emission', 'Adjusted Emission')

figure
plot(time + 1975, ratio_N, '--*')
xlabel('year'), ylabel('^{Adjusted E_{NO_x}}/_{Original E_{NO_x}}')

figure
plot(time + 1975, ratio_V, '--*')
xlabel('year'), ylabel('^{Adjusted E_{VOC}}/_{Original E_{VOC}}')

