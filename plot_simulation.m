function plot_simulation(Nsim, initdata, ConData, savename_input, varargin)
% Loads the results of simulations, and plots
% INPUT:
% * initdata: structure with initial settings, in '-filename--_initialsettings.mat';
% * ConData: structure with connectivity data, in '-filename--_ConData.mat'
% * Nsim: # simulations to plot
% * 'simvec' (optional): vector with which simulations to use

% exaple use: 
% * plot_simulation(10, initdata, ConData) plots the first 10 simulations
% * plot_simulation(10, initdata, ConData, )

f = filesep;

%% Plot example
if ~exist('Nsim','var')
    Nsim = input('Which simulation(s) do you want to plot? (give array)');
end

simvec = 1:Nsim;
% Look for 'varargin' inputs
len = length(varargin);
% check "len" for even number
if mod(len,2) > 0
    error('Wrong arguments: must be name-value pairs.');
end
for i = 1:2:len
    switch lower(varargin{i})
        case 'simvec'
            simvec=varargin{i+1};
            Nsim = length(simvec);
        otherwise
            % neglect invalid option
            disp(['Ignoring invalid input ' varargin{i}])
    end
end


whiskers = ones(1,Nsim);
for nt = 1:Nsim
    disp(['Loading simulation ' num2str(simvec(nt)) ' for plotting'])
    try
        if initdata.Vthresdyn
            thresholdname = 'dynthreshold_set';
        else
            thresholdname = 'fixthreshold_set';
        end
        load([ConData.savefolder ConData.FnametoSave '_Simcolumn_' thresholdname '_' initdata.setVthres.type '_simulation_' num2str(simvec(nt))  ]);
        load([ConData.savefolder savename_input '_Thalamic_Spike_Trains.mat'  ]);
    catch
        if nt>1
            keyboard
        end
        savefolder = input('What is the folder where the simulations were saved?', 's');
        if ~strcmp(savefolder(end), f)
            savefolder(end+1) = f ;
        end
        animal = input('What was the animal name?','s');
        savename = input('What were the file names?', 's');
        thresholdtype = input('What kind of threshold was used? (d/f)', 's');
        if strcmp(thresholdtype, 'd')
            thresholdname = 'dynthreshold_set';
        else
            thresholdname = 'fixthreshold_set';
        end
        thresholdset = input('How were the thresholds set? (d/t/i)','s');
        if strcmp(thresholdset, 'd')
            thresholdsetname = 'distribution';
        elseif strcmp(thresholdset, 't')
            thresholdsetname = 'pertype';
        else
            thresholdsetname = 'individual';
        end
        try
            if ~isempty(animal)
                ConData.savefolder = savefolder;
                ConData.FnametoSave = [savename '_' animal];
                initdata.setVthres.type = thresholdsetname;
                load([ConData.savefolder ConData.FnametoSave '_Simcolumn_' thresholdname '_' initdata.setVthres.type '_simulation_' num2str(simvec(nt))  ]);
                load([ConData.savefolder ConData.FnametoSave '_Thalamic_Spike_Trains.mat'  ], 'WhiskerTrace');
            else
                ConData.savefolder = savefolder;
                ConData.FnametoSave = savename;
                initdata.setVthres.type = thresholdsetname;
                load([ConData.savefolder ConData.FnametoSave '_Simcolumn_' thresholdname '_' initdata.setVthres.type '_simulation_' num2str(simvec(nt))  ]);
                load([ConData.savefolder ConData.FnametoSave '_Thalamic_Spike_Trains.mat'  ], 'WhiskerTrace');
            end
        catch
            load([savefolder savename num2str(simvec(nt)) '.mat'  ]);
            savenameth = input('What were the file names for the thalamic spike trains?', 's');
            load([savefolder savenameth], 'WhiskerTrace');
        end
    end
    
    [barrelnr, ind_barrel] = sort(cellinfo_all(:,5)); % sort by barrel
    celltypevec = cellinfo_all(ind_barrel,4);
    Nbarrel = max(barrelnr);
    Nall = length(modelsc);
    
    % === 仅绘制 L2/3 pyramidal (cell type = 5) 的膜电位 ===
    % 在“按 barrel 排序”的索引空间里先筛选出 L2/3 pyramidal
    is_L23_pyr_sorted = (cellinfo_all(ind_barrel,4) == 6) & (cellinfo_all(ind_barrel,5) == 2) ; 
    l23pyr_idx_sorted = find(is_L23_pyr_sorted);

    if isempty(l23pyr_idx_sorted)
        warning('未找到 L2/3 pyramidal 细胞，跳过该仿真绘制。');
    else
        % 你想画多少个细胞：默认取前 4 个；若不足 4 个则按实际数量
        Nplot = min(4, numel(l23pyr_idx_sorted));
        pick_sorted = l23pyr_idx_sorted(1:Nplot);   % 也可以改成 randperm(...) 随机抽样

        % 还原到“原全局索引”（与 V/U/VT 的行索引一致）
        pick_global = ind_barrel(pick_sorted);

        % 时间轴（注意：用 V 的列数来稳妥获取长度）
        time = (1:size(V,2)) * simdata.timestep;

        disp(['Plotting simulation ' num2str(simvec(nt)) ', L2/3 pyramidal n=' num2str(Nplot)])
        figure

        nrow = ceil(Nplot/2); ncol = 2; % 2x2 布局，少于4个也适配
        for ii = 1:Nplot
            gid = pick_global(ii);  % 全局 cell id（对应 V/U/VT 的行）
            % 取出该细胞的元数据（便于标题显示）
            pos_rc   = cellinfo_all(gid,1);  % rostral-caudal
            pos_dv   = cellinfo_all(gid,2);  % dorsal-ventral
            depth    = cellinfo_all(gid,3);  % 距离 pia
            ctype    = cellinfo_all(gid,4);  % 应为 5
            barrelID = cellinfo_all(gid,5);  % 桶号
            ei_flag  = cellinfo_all(gid,6);  % 1=兴奋

            subplot(nrow, ncol, ii)
            plot(time, V(gid,:), 'LineWidth', 1.0); hold on
            try
                plot(time, U(gid,:), 'LineWidth', 0.8);
            catch
                % U 未保存则忽略
            end
            if initdata.Vthresdyn
                try
                    plot(time, VT(gid,:), 'LineWidth', 0.8);
                catch
                    % VT 未保存则忽略
                end
            end
            ylim([-80, 20]); grid on; box on

            % 图例
            if initdata.Vthresdyn
                legend({'V','u','V_T'}, 'Location','best'); 
            else
                legend({'V','u'}, 'Location','best');
            end

            % 标题含关键信息：细胞ID / barrel / 坐标 / 深度 PV+ fast spike pyramidal
            title(sprintf('L2/3 pyramidal | id=%d | ',gid), 'Interpreter','none')
            xlabel('time (ms)'); ylabel('mV')
        end
        % 1) 构造群体信号（用你挑出的 L2/3 pyramidal：pick_global）
        Fs  = 1000 / simdata.timestep;       % 采样率 Hz
        lfp = mean(V(pick_global, :), 1);
        lfp = detrend(lfp);

        % 2) 全频段功率谱（1–150 Hz）
        win = round(Fs * 1.0);               % 1 s 窗，可调
        [pxx,f] = pwelch(lfp, hamming(win), [], [], Fs);
        band_all = (f >= 1) & (f <= 150);
        f_all = f(band_all); p_all = pxx(band_all);

        % 找全局主峰
        [pk_all, idx_all] = max(p_all);
        f_peak = f_all(idx_all);

        % 计算峰的“prominence”和半高宽（粗略）
        % 用邻域包络估计背景
        bk = movmedian(p_all, max(11, round(Fs/10)));   % 慢变背景
        prom = 10*log10(pk_all) - 10*log10(bk(idx_all));

        % FWHM：左右找功率=峰值一半处的频率
        half = pk_all/2;
        iL = find(p_all(1:idx_all) <= half, 1, 'last');
        iR = idx_all-1 + find(p_all(idx_all:end) <= half, 1, 'first');
        if isempty(iL), iL = 1; end
        if isempty(iR), iR = numel(p_all); end
        fwhm = f_all(iR) - f_all(iL);
        Q = f_peak / max(fwhm, eps);          % 品质因子

        % γ 带功率占比
        gamma_band = (f_all >= 30 & f_all <= 80);
        pow_gamma = trapz(f_all(gamma_band), p_all(gamma_band));
        pow_total = trapz(f_all, p_all);
        ratio_gamma = pow_gamma / pow_total;

        % 3) 判定规则（你可以调阈值）
        is_gamma = (f_peak >= 30 && f_peak <= 80) ...
                && (prom >= 3) ...            % ≥3 dB 高于局部背景
                && (Q >= 2) ...               % 峰要相对“尖”
                && (ratio_gamma >= 0.1);      % γ 功率占比≥10%

        fprintf('f_peak = %.1f Hz, prom=%.1f dB, Q=%.2f, γ-ratio=%.2f --> %s\n', ...
                f_peak, prom, Q, ratio_gamma, string(is_gamma));

        % 可视化
        figure;
        plot(f_all,10*log10(p_all),'LineWidth',1.2); hold on; grid on; box on
        xlabel('Frequency (Hz)'); ylabel('Power (dB)'); xlim([0 120])
        y = 10*log10(pk_all);
        plot(f_peak, y, 'ro', 'MarkerFaceColor','r');
        xline(30,'--k'); xline(80,'--k');
        title(sprintf('Peak = %.1f Hz | prom=%.1f dB | Q=%.2f | γ? %s', ...
              f_peak, prom, Q, string(is_gamma)));

    end

    
    
%     figure
%     if whiskers(nt)==1
%         for na=1:2
%             subplot(6,1,na)
%             plot((1:length(WhiskerTrace.Recording{na,simvec(nt)}))*WhiskerTrace.binsize, WhiskerTrace.Recording{na,simvec(nt)}, 'k', 'LineWidth',2)
%             xlim([1, length(WhiskerTrace.Recording{na,simvec(nt)})*WhiskerTrace.binsize])
%             if na==1
%                 title('Whisker Angle')
%             else
%                 title('Whisker Curvature')
%             end
%             set(gca, 'XGrid','on')
%         end
%     end
% 
%     [Nneuronth, ~] = size(inputspikes);
%     subplot(3,1,2)
%     mint = 1;
%     maxt = 0;
%     for nn = 1:Nneuronth
%         hold all
%         plot(inputspikes(nn,:), nn*ones(size(inputspikes(nn,:))), '.k')
%         u = unique(inputspikes(nn,:)); % remove 0
%         if length(u)>1
%             if u(2)<mint
%                 mint = u(2);
%             end
%         end
%         if max(inputspikes(nn,:))>maxt
%             maxt = max(inputspikes(nn,:));
%         end
%     end
%     NCB = 0;
%     [barrelnrth, ~] = sort(cellinfo_input(:,5)); % sort by barrel
%     for nb = 1:Nbarrel
%         % solid lines between barrels
%         NCB = NCB+sum(barrelnrth==nb);
%         try
%             plot([mint maxt],[NCB,NCB], '-b')
%         catch
%             plot([0 maxt],[NCB,NCB], '-b')
%         end
%     end
%     ylim([1 Nneuronth])
%     title('Thalamic input spikes')
%     box on
%     
%     subplot(3,1,3)
%     [Nneuron, ~] = size(modelspt);
%     for nn = 1:Nneuron
%         hold all
%         plot(modelspt(ind_barrel(nn),:), nn*ones(size(modelspt(ind_barrel(nn),:))), '.k')
%         u = unique(modelspt(ind_barrel(nn),:)); % remove 0
%         if length(u)>1
%             if u(2)<mint
%                 mint = u(2);
%             end
%         end
%         if max(modelspt(ind_barrel(nn),:))>maxt
%             maxt = max(modelspt(ind_barrel(nn),:));
%         end
%     end
%     
%     NCB = 0;
%     NCB_old = 0;
%     for nb = 1:Nbarrel
%         % solid lines between barrels
%         NCB = NCB+sum(barrelnr==nb);
%         try
%             plot([mint maxt],[NCB,NCB], '-b')
%         catch
%             plot([0 maxt],[NCB,NCB], '-b')
%         end
%         % dotted lines between cell-types
%         NCT = NCB_old;
%         for celltype = 1:15
%             NCT = NCT+sum(celltypevec(NCB_old+1:NCB)==celltype);
%             try
%                 plot([mint maxt],[NCT,NCT], '--r')
%             catch
%                 plot([0 maxt],[NCT,NCT], '--r')
%             end
%         end
%         NCB_old = NCB;
%     end
%     ylim([1 Nneuron])
%     title('Model spikes')
%     
%     for na=2:3
%         subplot(3,1,na)
%         xlim([mint, maxt])
%         xlabel('time (ms)')
%         ylabel('Neuron number')
%         set(gca, 'XGrid','on')
%         box on
%     end
end

end