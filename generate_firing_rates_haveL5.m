clc;close ALL;
% === 加载数据 ===
matpath = fullfile('Simulation_results','Test_sim_Svoboda_data','Test_sim_Svoboda_Simcolumn_dynthreshold_set_pertype_simulation_1.mat');
load(matpath, 'modelspt','cellinfo_all','inputspikes','cellinfo_input','simdata','simLen','V','modelspt_L5');

dt    = simdata.timestep;   % ms
Tmax  = simLen;             % ms
bin_ms = 1;
edges = 0:bin_ms:Tmax;           % 1 ms bin
visual_window = [100,200];

% === 只选第2个 barrel ✅ ===
barrel_id = 2;
is_in_barrel = cellinfo_all(:,5) == barrel_id;

% === 定义层级与兴奋/抑制分类 ===
is_L4   = cellinfo_all(:,4) == 2 | cellinfo_all(:,4) == 3;
is_L23  = cellinfo_all(:,4) == 5 | cellinfo_all(:,4) == 6;
is_exc  = cellinfo_all(:,6) == 1;
is_inh  = cellinfo_all(:,6) == -1;

N_L5 = size(modelspt_L5, 1);
N_L5_E = round(0.85 * N_L5);
L5Exc_idx = 1:N_L5_E;
L5Inh_idx = N_L5_E+1:N_L5;

groups = struct( ...
  'L4Exc',  find(is_L4  & is_exc  & is_in_barrel), ...
  'L4Inh',  find(is_L4  & is_inh  & is_in_barrel), ...
  'L23Exc', find(is_L23 & is_exc  & is_in_barrel), ...
  'L23Inh', find(is_L23 & is_inh  & is_in_barrel));

% 定义颜色变量
grey_color = [210, 210, 210] / 255;
red_color =[227,90,94] / 255;
blue_color = [61, 85, 160] / 255;
orange_color = [239, 164, 142] / 255;

% 预计算各组信息（调用局部函数 compute_group_info，定义在文件末尾）
info.L4Exc  = compute_group_info(groups.L4Exc,  modelspt, edges, bin_ms);
info.L4Inh  = compute_group_info(groups.L4Inh,  modelspt, edges, bin_ms);
info.L23Exc = compute_group_info(groups.L23Exc, modelspt, edges, bin_ms);
info.L23Inh = compute_group_info(groups.L23Inh, modelspt, edges, bin_ms);
info.L5Exc = compute_group_info(L5Exc_idx, modelspt_L5, edges, bin_ms);
info.L5Inh = compute_group_info(L5Inh_idx, modelspt_L5, edges, bin_ms);
% === 作图 ===
% 阴影颜色与透明度（t1 = onset 差，t2 = peak 差）
t1_color = [1.00, 0.90, 0.45]; % 浅黄
t2_color = [0.50, 0.95, 0.60];  % 浅绿
t3_color = [221, 161, 190]/255;  % 浅绿
alpha_sh = 0.35;


% --- 丘脑输入（保持不变） ---
subplot(4,1,1); hold on;
if exist('inputspikes','var') && ~isempty(inputspikes)
    % 只选择感兴趣的神经元行（这里 201:400）
    rows = 201:400;

    % 取出这些行并展平为时间向量
    th_spk = inputspikes(rows, :);
    th_spk = th_spk(:);
    th_spk = th_spk(th_spk>0 & th_spk<=edges(end));

    % 直方计数（用已有的 edges，确保 edges 对应于你想展示的时间范围）
    th_counts = histcounts(th_spk, edges);

    % 归一化：按所选神经元数归一到每个神经元的 firing rate (spikes / s)
    nCells = numel(rows);
    th_rate = th_counts / nCells * (1000/bin_ms);

    % 平滑并绘图（保留颜色、线宽、线型、图例等）
    th_rate = smoothdata(th_rate, 'gaussian', 1);
    plot(edges(1:end-1), th_rate, 'Color', grey_color, 'LineWidth', 2, 'LineStyle', '-', 'DisplayName', 'Thalamic Input');
    legend('Location', 'best');
end
xlabel('Time (ms)'); ylabel('Firing rate (Hz)');
title('Thalamic Input'); grid on; box on; xlim(visual_window);

% --- L4 子图：先绘制阴影（如果存在），再绘曲线以覆盖阴影 ---
ax = subplot(4,1,2); hold on;
% 获取 y 范围以填充阴影
yL4 = [0, max([max(info.L4Exc.rate_sm,[],'omitnan'), max(info.L4Inh.rate_sm,[],'omitnan'), 1]) * 1.05];
ylim(yL4);

% t1：onset_inh -> onset_exc
% if ~isnan(info.L4Inh.onset_time) && ~isnan(info.L4Exc.onset_time) && (info.L4Exc.onset_time ~= info.L4Inh.onset_time)
%     xs = sort([info.L4Inh.onset_time, info.L4Exc.onset_time]);
%     patch(ax, [xs(1) xs(2) xs(2) xs(1)], [yL4(1) yL4(1) yL4(2) yL4(2)], t1_color, ...
%         'FaceAlpha', alpha_sh, 'EdgeColor', 'none', 'DisplayName', 't1 (onset)');
% end

% t2：peak_inh -> peak_exc
% if ~isnan(info.L4Inh.first_peak_time) && ~isnan(info.L4Exc.first_peak_time) && (info.L4Exc.first_peak_time ~= info.L4Inh.first_peak_time)
%     xs = sort([info.L4Inh.first_peak_time, info.L4Exc.first_peak_time]);
%     patch(ax, [xs(1) xs(2) xs(2) xs(1)], [yL4(1) yL4(1) yL4(2) yL4(2)], t2_color, ...
%         'FaceAlpha', alpha_sh, 'EdgeColor', 'none', 'DisplayName', 't2 (first peak)');
% end

% 绘制 L4 曲线（画在阴影之上）
if ~isempty(groups.L4Exc)
    plot(edges(1:end-1), info.L4Exc.rate_sm, 'Color', red_color, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName','L4Exc');
end
if ~isempty(groups.L4Inh)
    plot(edges(1:end-1), info.L4Inh.rate_sm, 'Color', red_color, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName','L4Inh');
end

% % 额外标注 onset/peak 的竖线与点（可选）
% if ~isnan(info.L4Exc.onset_time)
%     plot(info.L4Exc.onset_time, interp1(edges(1:end-1), info.L4Exc.rate_sm, info.L4Exc.onset_time,'linear','extrap'), 'o', 'Color', red_color, 'MarkerFaceColor', red_color, 'HandleVisibility','off');
% end
% if ~isnan(info.L4Inh.onset_time)
%     plot(info.L4Inh.onset_time, interp1(edges(1:end-1), info.L4Inh.rate_sm, info.L4Inh.onset_time,'linear','extrap'), 'o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'HandleVisibility','off');
% end
% if ~isnan(info.L4Exc.first_peak_time)
%     xline(info.L4Exc.first_peak_time, '-','Color', red_color, 'LineWidth', 1, 'HandleVisibility','off');
% end
% if ~isnan(info.L4Inh.first_peak_time)
%     xline(info.L4Inh.first_peak_time, '--','Color', [0 0 0], 'LineWidth', 1, 'HandleVisibility','off');
% end

xlabel('Time (ms)'); ylabel('Firing rate (Hz)');
title('L4'); grid on; box on; xlim(visual_window);

% --- L2/3 子图：同样的阴影/曲线顺序 ---
ax2 = subplot(4,1,3); hold on;
yL23 = [0, max([max(info.L23Exc.rate_sm,[],'omitnan'), max(info.L23Inh.rate_sm,[],'omitnan'), 1]) * 1.05];
ylim(yL23);

% % t1（onset）
% if ~isnan(info.L23Inh.onset_time) && ~isnan(info.L23Exc.onset_time) && (info.L23Exc.onset_time ~= info.L23Inh.onset_time)
%     xs = sort([info.L23Inh.onset_time, info.L23Exc.onset_time]);
%     patch(ax2, [xs(1) xs(2) xs(2) xs(1)], [yL23(1) yL23(1) yL23(2) yL23(2)], t1_color, ...
%         'FaceAlpha', alpha_sh, 'EdgeColor', 'none', 'HandleVisibility','on', 'DisplayName', 't1 (onset)');
% end

% t2（first peak）
% if ~isnan(info.L23Inh.first_peak_time) && ~isnan(info.L23Exc.first_peak_time) && (info.L23Exc.first_peak_time ~= info.L23Inh.first_peak_time)
%     xs = sort([info.L23Inh.first_peak_time, info.L23Exc.first_peak_time]);
%     patch(ax2, [xs(1) xs(2) xs(2) xs(1)], [yL23(1) yL23(1) yL23(2) yL23(2)], t2_color, ...
%         'FaceAlpha', alpha_sh, 'EdgeColor', 'none', 'HandleVisibility','on', 'DisplayName', 't2 (first peak)');
% end

% 绘制 L2/3 曲线
if ~isempty(groups.L23Exc)
    plot(edges(1:end-1), info.L23Exc.rate_sm, 'Color', blue_color, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName','L23Exc');
end
if ~isempty(groups.L23Inh)
    plot(edges(1:end-1), info.L23Inh.rate_sm, 'Color', blue_color, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName','L23Inh');
end

% % 标注 onset/peak
% if ~isnan(info.L23Exc.onset_time)
%     plot(info.L23Exc.onset_time, interp1(edges(1:end-1), info.L23Exc.rate_sm, info.L23Exc.onset_time,'linear','extrap'), 'o', 'Color', blue_color, 'MarkerFaceColor', blue_color, 'HandleVisibility','off');
% end
% if ~isnan(info.L23Inh.onset_time)
%     plot(info.L23Inh.onset_time, interp1(edges(1:end-1), info.L23Inh.rate_sm, info.L23Inh.onset_time,'linear','extrap'), 'o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'HandleVisibility','off');
% end
% if ~isnan(info.L23Exc.first_peak_time)
%     xline(info.L23Exc.first_peak_time, '-','Color', blue_color, 'LineWidth', 1, 'HandleVisibility','off');
% end
% if ~isnan(info.L23Inh.first_peak_time)
%     xline(info.L23Inh.first_peak_time, '--','Color', [0 0 0], 'LineWidth', 1, 'HandleVisibility','off');
% end

xlabel('Time (ms)'); ylabel('Firing rate (Hz)');
title('L2/3'); grid on; box on; xlim(visual_window);


% L5
% N_L5 = size(modelspt_L5, 1);
% N_L5_E = round(0.85 * N_L5);
% L5Exc_idx = 1:N_L5_E;
% L5Inh_idx = N_L5_E+1:N_L5;
ax3 = subplot(4,1,4); hold on;
yL5 = [0, max([max(info.L5Exc.rate,[],'omitnan'), max(info.L5Inh.rate,[],'omitnan'), 1]) * 1.05];
ylim(yL23);
plot(edges(1:end-1), info.L5Exc.rate_sm, 'Color', green_color, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', 'L5Exc');
plot(edges(1:end-1), info.L5Inh.rate_sm, 'Color', green_color, 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'L5Inh');
xlabel('Time (ms)'); ylabel('Firing rate (Hz)');
title('L5'); box on; xlim(visual_window);

% --- 合并图例：取每个轴的可见句柄，避免重复 --- 
axes(ax); % 回到 L4 轴
hL = findobj(ax, '-property', 'DisplayName');
% 为了让图例显示 t1/t2，一般在绘制阴影时已设置 DisplayName
legend(ax, 'show', 'Location', 'best');
axes(ax2);
legend(ax2, 'show', 'Location', 'best');

% === 计算并打印 t1 与 t2（L4 与 L2/3） ===
t1_L4 = nan; t2_L4 = nan;
t1_L23 = nan; t2_L23 = nan;

if ~isnan(info.L4Exc.onset_time) && ~isnan(info.L4Inh.onset_time)
    t1_L4 = info.L4Exc.onset_time - info.L4Inh.onset_time;
end
if ~isnan(info.L4Exc.first_peak_time) && ~isnan(info.L4Inh.first_peak_time)
    t2_L4 = info.L4Exc.first_peak_time - info.L4Inh.first_peak_time;
end

if ~isnan(info.L23Exc.onset_time) && ~isnan(info.L23Inh.onset_time)
    t1_L23 = info.L23Exc.onset_time - info.L23Inh.onset_time;
end
if ~isnan(info.L23Exc.first_peak_time) && ~isnan(info.L23Inh.first_peak_time)
    t2_L23 = info.L23Exc.first_peak_time - info.L23Inh.first_peak_time;
end

% fprintf('L4: t1 (onset_exc - onset_inh) = %.2f ms, t2 (first-peak_exc - inh) = %.2f ms\n', t1_L4, t2_L4);
% fprintf('L2/3: t1 (onset_exc - onset_inh) = %.2f ms, t2 (first-peak_exc - inh) = %.2f ms\n', t1_L23, t2_L23);

if isnan(t1_L4) || isnan(t2_L4) || isnan(t1_L23) || isnan(t2_L23)
    warning('部分 onset 或 peak 无法确定（NaN）。通常是该群体无脉冲或无法检测到明确的峰值）。');
end

%% === 将 50 个 100ms 窗口映射到 0–100 ms，并跨实验平均 ===
win_len   = 100;                   % 100 ms 窗口长度
win_starts = 100:100:5000;         % 50 个窗口起点：100, 200, ..., 5000
local_t   = 0:win_len-1;           % 局部时间轴 [0,99] ms
t = edges(1:end-1);                % 与 rate_sm 对齐的全局时间轴

% 1) 切片到矩阵 (50 x 100)，每行是一个窗口映射到 0–100ms
M_L4E  = slice_to_matrix(info.L4Exc.rate_sm,  t, win_starts, win_len);
M_L4I  = slice_to_matrix(info.L4Inh.rate_sm,  t, win_starts, win_len);
M_L23E = slice_to_matrix(info.L23Exc.rate_sm, t, win_starts, win_len);
M_L23I = slice_to_matrix(info.L23Inh.rate_sm, t, win_starts, win_len);
M_L5E  = slice_to_matrix(info.L5Exc.rate_sm,  t, win_starts, win_len);
M_L5I  = slice_to_matrix(info.L5Inh.rate_sm,  t, win_starts, win_len);

% 2) 跨 50 次实验逐点取平均，得到新的平均 firing rate (0–100ms)
mean_L4E  = mean(M_L4E,  1, 'omitnan');
mean_L4I  = mean(M_L4I,  1, 'omitnan');
mean_L23E = mean(M_L23E, 1, 'omitnan');
mean_L23I = mean(M_L23I, 1, 'omitnan');
mean_L5E  = mean(M_L5E,  1, 'omitnan');
mean_L5I  = mean(M_L5I,  1, 'omitnan');

% 3) 每个实验（窗口）内的主峰时间差（Exc - Inh），基于“最大值所在时刻”
dL4  = win_peak_time_diff(M_L4E,  M_L4I,  local_t);
dL23 = win_peak_time_diff(M_L23E, M_L23I, local_t);
dL5 = win_peak_time_diff(M_L5E, M_L5I, local_t);

% 4) 计算峰值时间差
mean_dL4   = mean(dL4,  'omitnan');
std_dL4    = std(dL4,   0, 'omitnan');

mean_dL23  = mean(dL23, 'omitnan');
std_dL23   = std(dL23,  0, 'omitnan');

mean_dL5  = mean(dL5, 'omitnan');
std_dL5   = std(dL5,  0, 'omitnan');

% 4) 打印 50 个差值并给出均值
fprintf('\n=== 每个100ms窗口的主峰时间差（Exc - Inh，单位 ms）===\n');
for i = 1:numel(win_starts)
    a = win_starts(i); b = a + win_len;
    fprintf('窗口 %2d [%4d,%4d):  L4 Δ=%.2f,  L2/3 Δ=%.2f,  L5 Δ=%.2f\n', i, a, b, dL4(i), dL23(i), dL5(i));
end

fprintf('—— Δt 统计（忽略 NaN）——\n');
fprintf(['L4: mean = %.2f ms, std = %.2f ms; ' ...
        'L2/3: mean = %.2f ms, std = %.2f ms; '...
        'L5: mean = %.2f ms, std = %.2f ms\n'],...
        mean_dL4, std_dL4, mean_dL23, std_dL23,mean_dL5, std_dL5);

% ===== 计算平均曲线的峰值时间（0-100 ms 内） =====
%% ===== 1) 先算峰值时间 =====
% L4 Exc
[pks_L4E, locs_L4E] = findpeaks(mean_L4E, 'MinPeakProminence', 0.8*max(mean_L4E));
if isempty(locs_L4E)
    [pks_L4E, locs_L4E] = findpeaks(mean_L4E);
end
if ~isempty(locs_L4E)
    t_peak_L4E = local_t(locs_L4E(1));
else
    t_peak_L4E = NaN;
end

% L4 Inh
[pks_L4I, locs_L4I] = findpeaks(mean_L4I, 'MinPeakProminence', 0.8*max(mean_L4I));
if isempty(locs_L4I)
    [pks_L4I, locs_L4I] = findpeaks(mean_L4I);
end
if ~isempty(locs_L4I)
    t_peak_L4I = local_t(locs_L4I(1));
else
    t_peak_L4I = NaN;
end

% L2/3 Exc
[pks_L23E, locs_L23E] = findpeaks(mean_L23E, 'MinPeakProminence', 0.8*max(mean_L23E));
if isempty(locs_L23E)
    [pks_L23E, locs_L23E] = findpeaks(mean_L23E);
end
if ~isempty(locs_L23E)
    t_peak_L23E = local_t(locs_L23E(1));
else
    t_peak_L23E = NaN;
end

% L2/3 Inh
[pks_L23I, locs_L23I] = findpeaks(mean_L23I, 'MinPeakProminence', 0.8*max(mean_L23I));
if isempty(locs_L23I)
    [pks_L23I, locs_L23I] = findpeaks(mean_L23I);
end
if ~isempty(locs_L23I)
    t_peak_L23I = local_t(locs_L23I(1));
else
    t_peak_L23I = NaN;
end

% L5 Exc
[pks_L5E, locs_L5E] = findpeaks(mean_L5E, 'MinPeakProminence', 0.8*max(mean_L5E));
if isempty(locs_L5E)
    [pks_L5E, locs_L5E] = findpeaks(mean_L5E);
end
if ~isempty(locs_L5E)
    t_peak_L5E = local_t(locs_L5E(1));
else
    t_peak_L5E = NaN;
end

% L5 Inh
[pks_L5I, locs_L5I] = findpeaks(mean_L5I, 'MinPeakProminence', 0.8*max(mean_L5I));
if isempty(locs_L5I)
    [pks_L5I, locs_L5I] = findpeaks(mean_L5I);
end
if ~isempty(locs_L5I)
    t_peak_L5I = local_t(locs_L5I(1));
else
    t_peak_L5I = NaN;
end

% === 可选：是否添加“输入”的平均子图 ===
SHOW_INPUT = true;   % <--- 这里改 true/false 控制是否显示输入子图

% 如果展示输入子图，则也做窗口映射与跨实验平均
if SHOW_INPUT
    M_TH   = slice_to_matrix(th_rate, t, win_starts, win_len);
    mean_TH = mean(M_TH, 1, 'omitnan');
end

% === 分层展示平均 firing rate（按开关决定 2 或 3 个子图）===
nrows = 3 + double(SHOW_INPUT);  % true=>3, false=>2
figure('Name', sprintf('Aligned mean firing rate by layer (0–100 ms)'), 'Color','w');

row_idx = 1;

% (可选) 输入子图（放在最上方）
if SHOW_INPUT
    subplot(nrows,1,row_idx); hold on;
    plot(local_t, mean_TH, '-', 'Color', grey_color, 'LineWidth', 2, 'DisplayName','Input');
    xlabel('Time (ms)');
    ylabel('Firing rate (Hz)');
    title('INPUT');
    grid on; box on; xlim([0, 80]); legend('Location','best');
    set(gca, 'FontSize', 14, 'FontName', 'Arial');
    row_idx = row_idx + 1;
end

% --- L4 子图 ---
subplot(nrows,1,row_idx); hold on;
plot(local_t, mean_L4E, '-', 'Color', red_color, 'LineWidth', 2, 'DisplayName','Exc');
plot(local_t, mean_L4I, '--', 'Color', red_color, 'LineWidth', 2, 'DisplayName','PV+');

% 获取 L4 当前 y 轴范围
ylims_L4 = ylim;

% 计算填充区间并 patch
if ~isnan(t_peak_L4E) && ~isnan(t_peak_L4I)
    x_fill_L4 = sort([t_peak_L4E, t_peak_L4I]);
    xv = [x_fill_L4(1), x_fill_L4(2), x_fill_L4(2), x_fill_L4(1)];
    yv = [ylims_L4(1), ylims_L4(1), ylims_L4(2), ylims_L4(2)];
    patch(xv, yv, t1_color, 'FaceAlpha', alpha_sh, 'EdgeColor','none', ...
          'DisplayName', '\Deltat');
    % 若想让曲线显示更清晰，可再画一次曲线
    plot(local_t, mean_L4E, '-', 'Color', red_color, 'LineWidth', 2, 'HandleVisibility','off');
    plot(local_t, mean_L4I, '--', 'Color', red_color, 'LineWidth', 2, 'HandleVisibility','off');
end

xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
title('L4');
grid on; box on; xlim([0, 80]); legend('Location','best');
set(gca, 'FontSize', 14, 'FontName', 'Arial');
row_idx = row_idx + 1;

% --- L2/3 子图 ---
subplot(nrows,1,row_idx); hold on;
plot(local_t, mean_L23E, '-', 'Color', blue_color, 'LineWidth', 2, 'DisplayName','Exc');
plot(local_t, mean_L23I, '--', 'Color', blue_color, 'LineWidth', 2, 'DisplayName','PV+');

% 获取 L2/3 当前 y 轴范围
ylims_L23 = ylim;

% 计算填充区间并 patch
if ~isnan(t_peak_L23E) && ~isnan(t_peak_L23I)
    x_fill_L23 = sort([t_peak_L23E, t_peak_L23I]);
    xv = [x_fill_L23(1), x_fill_L23(2), x_fill_L23(2), x_fill_L23(1)];
    yv = [ylims_L23(1), ylims_L23(1), ylims_L23(2), ylims_L23(2)];
    patch(xv, yv, t2_color, 'FaceAlpha', alpha_sh, 'EdgeColor','none', ...
          'DisplayName', '\Deltat');
    plot(local_t, mean_L23E, '-', 'Color', blue_color, 'LineWidth', 2, 'HandleVisibility','off');
    plot(local_t, mean_L23I, '--', 'Color', blue_color, 'LineWidth', 2, 'HandleVisibility','off');
end

xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
title('L2/3');
grid on; box on; xlim([0, 80]); legend('Location','best');
set(gca, 'FontSize', 14, 'FontName', 'Arial');
row_idx = row_idx + 1;

% --- L5 子图 ---
subplot(nrows,1,row_idx); hold on;
plot(local_t, mean_L5E, '-', 'Color', green_color, 'LineWidth', 2, 'DisplayName','Exc');
plot(local_t, mean_L5I, '--', 'Color', green_color, 'LineWidth', 2, 'DisplayName','PV+');

% 获取 L5 当前 y 轴范围
ylims_L5 = ylim;

% 计算填充区间并 patch
if ~isnan(t_peak_L5E) && ~isnan(t_peak_L5I)
    x_fill_L5 = sort([t_peak_L5E, t_peak_L5I]);
    xv = [x_fill_L5(1), x_fill_L5(2), x_fill_L5(2), x_fill_L5(1)];
    yv = [ylims_L5(1), ylims_L5(1), ylims_L5(2), ylims_L5(2)];
    patch(xv, yv, t3_color, 'FaceAlpha', alpha_sh, 'EdgeColor','none', ...
          'DisplayName', '\Deltat');
    plot(local_t, mean_L5E, '-', 'Color', green_color, 'LineWidth', 2, 'HandleVisibility','off');
    plot(local_t, mean_L5I, '--', 'Color', green_color, 'LineWidth', 2, 'HandleVisibility','off');
end

xlabel('Time (ms)');
ylabel('Firing rate (Hz)');
title('L5');
grid on; box on; xlim([0, 80]); legend('Location','best');
set(gca, 'FontSize', 14, 'FontName', 'Arial');

%% === 四子图：同层上下（列=层），Inh 单独编号且用自己的 y 轴范围 ===
twin   = [100 180];   % 时间窗口(ms)
xtick_values = linspace(100, 180, 5);  % Original x-axis tick positions
xtick_labels = linspace(0, 80, 5);   % New x-axis labels from 0 to 100 ms
dotsz  = 5;
L4_inh_color = [212, 97, 83] / 255;
L23_inh_color =[78,102,145] / 255;
L5_inh_color =[232,173,137] / 255;

figure('Name','列=层，行=Exc/Inh');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
ax = gobjects(6,1);

% ---------- 1列：L4 ----------
% (1,1) L4 Exc —— y 轴范围只按 Exc 数量
ax(1) = nexttile(1); hold on;
rowsL4_exc = make_rows_single(groups.L4Exc);
plot_raster_group(modelspt, groups.L4Exc, rowsL4_exc, twin, red_color, dotsz);
xlim(twin); ylim([0.5, max(1,numel(groups.L4Exc))+0.5]);
ylabel('Neuron number'); title('L4 Exc');
set(gca, 'XTick', xtick_values, 'XTickLabel', xtick_labels); 
set(gca, 'FontSize', 14, 'FontName', 'Arial');
grid on; box on;

% (2,1) L4 Inh —— y 轴范围只按 Inh 数量
ax(4) = nexttile(4); hold on;
rowsL4_inh = make_rows_single(groups.L4Inh);
plot_raster_group(modelspt, groups.L4Inh, rowsL4_inh, twin, L4_inh_color, dotsz);
xlim(twin); ylim([0.5, max(1,numel(groups.L4Inh))+0.5]);
xlabel('Time (ms)'); ylabel('Neuron number'); title('L4 PV+');
set(gca, 'XTick', xtick_values, 'XTickLabel', xtick_labels);
set(gca, 'FontSize', 14, 'FontName', 'Arial');
grid on; box on;

% ---------- 2列：L2/3 ----------
% (1,2) L2/3 Exc —— y 轴范围只按 Exc 数量
ax(2) = nexttile(2); hold on;
rowsL23_exc = make_rows_single(groups.L23Exc);
plot_raster_group(modelspt, groups.L23Exc, rowsL23_exc, twin, blue_color, dotsz);
xlim(twin); ylim([0.5, max(1,numel(groups.L23Exc))+0.5]);
title('L2/3 Exc');ylabel('Neuron number');
set(gca, 'XTick', xtick_values, 'XTickLabel', xtick_labels);
set(gca, 'FontSize', 14, 'FontName', 'Arial');
grid on; box on;

% (2,2) L2/3 Inh —— y 轴范围只按 Inh 数量
ax(5) = nexttile(5); hold on;
rowsL23_inh = make_rows_single(groups.L23Inh);
plot_raster_group(modelspt, groups.L23Inh, rowsL23_inh, twin, L23_inh_color, dotsz);
xlim(twin); ylim([0.5, max(1,numel(groups.L23Inh))+0.5]);
xlabel('Time (ms)'); title('L2/3 PV+');ylabel('Neuron number');
set(gca, 'XTick', xtick_values, 'XTickLabel', xtick_labels);
set(gca, 'FontSize', 14, 'FontName', 'Arial');
grid on; box on;

% ---------- 3列：L5 ----------
% (1,3) L5 Exc
ax(3) = nexttile(3); hold on;
rowsL5_exc = make_rows_single(L5Exc_idx);  % 这里是 1:N_L5_E，相当于局部行号
plot_raster_group(modelspt_L5, L5Exc_idx, rowsL5_exc, twin, orange_color, dotsz);
xlim(twin); ylim([0.5, max(1,numel(L5Exc_idx))+0.5]);
title('L5 Exc'); ylabel('Neuron number');
set(gca, 'XTick', xtick_values, 'XTickLabel', xtick_labels);
set(gca, 'FontSize', 14, 'FontName', 'Arial');
grid on; box on;

% (2,3) L5 Inh
ax(6) = nexttile(6); hold on;
rowsL5_inh = make_rows_single(L5Inh_idx);
plot_raster_group(modelspt_L5, L5Inh_idx, rowsL5_inh, twin, L5_inh_color, dotsz);
xlim(twin); ylim([0.5, max(1,numel(L5Inh_idx))+0.5]);
xlabel('Time (ms)'); title('L5 PV+'); ylabel('Neuron number');
set(gca, 'XTick', xtick_values, 'XTickLabel', xtick_labels);
set(gca, 'FontSize', 14, 'FontName', 'Arial');
grid on; box on;
% 只联动 x 轴（时间对齐）；y 轴彼此独立
linkaxes(ax,'x');


%% === 局部函数：计算组信息 ===
function info = compute_group_info(idx, modelspt, edges, bin_ms)
    info = struct('counts',[],'rate',[],'rate_sm',[],'onset_time',nan,'first_peak_time',nan);
    if isempty(idx)
        return;
    end
    spk = modelspt(idx,:);
    spk = spk(spk>0 & spk<=edges(end));
    counts = histcounts(spk, edges); % 每 ms 的脉冲数（长度为 length(edges)-1）
    rate = counts / numel(idx) * (1000/bin_ms); % 转换为 Hz（群体平均）
    rate_sm = smoothdata(rate, 'gaussian', 1);
    info.counts = counts;
    info.rate = rate;
    info.rate_sm = rate_sm;
    % onset：第一个 counts>0 的 bin（未平滑）
    onset_bin = find(counts>0, 1, 'first');
    if ~isempty(onset_bin)
        info.onset_time = edges(onset_bin); % edges 为 bin 起点
    end
    % first peak：对平滑速率使用 findpeaks，优先使用一定的 prominence，若未找到再退化到无 prominence
    if max(rate_sm) > 0
        prom = 1 * max(rate_sm);
        [pks, locs] = findpeaks(rate_sm, 'MinPeakProminence', prom);
        if isempty(locs)
            [pks, locs] = findpeaks(rate_sm);
        end
        if ~isempty(locs)
            info.first_peak_time = edges(locs(1));
        end
    end
end
% ====== 辅助函数 ======
function [rows_map, order_ids] = make_layer_rows(exc_ids, inh_ids)
% 为“同一层”建立统一的 y 行号映射，使 Exc/Inh 在同一层对齐
    if nargin<1 || isempty(exc_ids), exc_ids = []; end
    if nargin<2 || isempty(inh_ids), inh_ids = []; end
    order_ids = unique([exc_ids(:); inh_ids(:)], 'stable');  % 层内统一顺序
    if isempty(order_ids)
        rows_map = containers.Map('KeyType','int32','ValueType','int32');
        return;
    end
    rows_map = containers.Map(num2cell(int32(order_ids)), num2cell(int32(1:numel(order_ids))));
end

function plot_raster_group(modelspt, idx_list, rows_map, twin, color_, dotsz)
% 把给定细胞集合按 rows_map 映射到统一 y 行号，并在 twin 时间窗内打点
    if isempty(idx_list), return; end
    X = []; Y = [];
    tmin = twin(1); tmax = twin(2);
    idx_list = idx_list(:)';
    for id = idx_list
        if ~isKey(rows_map, int32(id)), continue; end
        row = rows_map(int32(id));
        spk = modelspt(id,:);                      % 第 id 个细胞的所有脉冲时刻
        spk = spk(spk>0 & spk>=tmin & spk<=tmax);  % 时间窗内
        if isempty(spk), continue; end
        X = [X, spk(:)'];                          %#ok<AGROW>
        Y = [Y, repmat(row,1,numel(spk))];         %#ok<AGROW>
    end
    if ~isempty(X)
        plot(X, Y, '.', 'Color', color_, 'MarkerSize', dotsz, 'HandleVisibility','off');
    end
end

% function plot_raster_group(modelspt, idx_list, rows_map, twin, color_, dotsz)
% % 把给定细胞集合按 rows_map 映射到统一 y 行号，并在 twin 时间窗内打点
%     if isempty(idx_list), return; end
%     X = []; Y = [];
%     tmin = twin(1); tmax = twin(2);
%     idx_list = idx_list(:)';
%     
%     for id = idx_list
%         if ~isKey(rows_map, int32(id)), continue; end
%         row = rows_map(int32(id));
%         spk = modelspt(id,:);                      % 第 id 个细胞的所有脉冲时刻
%         spk = spk(spk>0 & spk>=tmin & spk<=tmax);  % 时间窗内
%         
%         if isempty(spk), continue; end
%         
%         % 将脉冲时间四舍五入到最近的3ms倍数
%         spk_3ms = round(spk / 3) * 3;
%         
%         % 只保留3ms内的脉冲时间
%         spk_3ms = spk_3ms(spk_3ms >= tmin & spk_3ms <= tmax);
%         
%         if isempty(spk_3ms), continue; end
%         
%         X = [X, spk_3ms(:)'];                          %#ok<AGROW>
%         Y = [Y, repmat(row, 1, numel(spk_3ms))];         %#ok<AGROW>
%     end
%     
%     if ~isempty(X)
%         plot(X, Y, '.', 'Color', color_, 'MarkerSize', dotsz, 'HandleVisibility','off');
%     end
% end


function rows_map = make_rows_single(idx_list)
    idx_list = idx_list(:)';                 % 保持顺序即可（也可按你想要的规则排序）
    if isempty(idx_list)
        rows_map = containers.Map('KeyType','int32','ValueType','int32'); return;
    end
    rows_map = containers.Map(num2cell(int32(idx_list)), ...
                              num2cell(int32(1:numel(idx_list))));
end

% FOR 50TIMTS
function M = slice_to_matrix(rate_sm, t, win_starts, win_len)
% 将每个窗口 [s, s+win_len) 的速率切出并映射到 0..win_len-1
% 输出 M 大小为 (nwin x win_len)，不足部分以 NaN 填充
    nwin = numel(win_starts);
    M = nan(nwin, win_len);
    if isempty(rate_sm) || isempty(t), return; end
    for i = 1:nwin
        a = win_starts(i); b = a + win_len;
        idx = (t >= a) & (t < b);
        seg = rate_sm(idx);
        L = min(numel(seg), win_len);
        if L > 0
            M(i, 1:L) = seg(1:L); % 对齐到局部 0..L-1 ms
        end
    end
end

function deltas = win_peak_time_diff(M_exc, M_inh, local_t)
% 对每个窗口，取 Exc 与 Inh 在该窗口内的“主峰”（最大值）时间差（Exc - Inh）
% 若某窗口缺数据或全为 NaN，则返回 NaN
    nwin = size(M_exc,1);
    deltas = nan(nwin,1);
    for i = 1:nwin
        re = M_exc(i,:); ri = M_inh(i,:);
        if all(~isfinite(re)) || all(~isfinite(ri))
            continue;
        end
        [~, ie] = max(re); [~, ii] = max(ri);
        if isfinite(re(ie)) && isfinite(ri(ii))
            deltas(i) = local_t(ie) - local_t(ii);
        end
    end
end
