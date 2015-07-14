%% Testing the score with example from Son et al.
format compact

time_pts = [1,2,3,4,5]'
x1 = [2,3,6,4,7]'
x2 = [1,2,3,5,3]'
y1 = [4,3,6,2,7]'
y2 = [5,2,3,1,3]'

figure()
subplot(1,2,1)
plot(time_pts, x1, 'ro-'), hold on;
plot(time_pts, x2, 'bo-'), hold off;

subplot(1,2,2)
plot(time_pts, y1, 'ro-'), hold on;
plot(time_pts, y2, 'bo-'), hold off;

%% simple correlation values
disp('Single correlation coefficients')
corr(x1, x2, 'type', 'spearman') % 0.667 Son
corr(x1, x2, 'type', 'pearson')  % 0.439 Son

%% ys1 and yr1 scores
disp('ys1 and yr1 coefficients')
ys1(x1, x2, time_pts) % 0.583 Son
ys1(y1, y2, time_pts) % 0.833 Son

yr1(x1, x2, time_pts) % 0.555 Son
yr1(y1, y2, time_pts) % 0.805 Son