function SpikeTrainStruct = make_thalamic_spike_trains_svoboda_recordings(savefolder, savename_input, SvobodaStruct, Barrelstruct, make_new_thalamic_kernels)
% Make Thalamic spike trains from Svoboda recording.

% Input:
% * savefolder & savename: folder where to store file with spike trains
% * SvobodaStruct with fields
%       * loadfolder (string): folder with Svoboda data files
%       * dataname (string): file name suffix
%       * volume (integer > 2): trials from which volume to use
%       * animal (string): animal name
%       * sessionvec (string array): array of sessions to load
%       * window (struct) with fields
%           *.start: 'pole in reach' or 'first touch' or 'first'
%           *.window (2x1 array with window (ms) around start, first touch or pole in reach trial
%       * trialvec (optional): vector of which trials to use
% * BarrelStruct{barrelidX, barrelIDY} with fields
%       * .mainbarrel (1,2,3): main, secondary or tertiary barrel
% * make_new_thalamic_kernels: 0 (use existing file) or 1 (make new kernels)

disp('Making new Thalamic spike trains')
    
f = filesep;
    %% Load Svoboda recordings
    addpath(genpath(SvobodaStruct.loadfolder));

    plotcheck = 0;
    [whiskermat, ~, ~, binsize_whisker] = load_data_across_sessions(SvobodaStruct.loadfolder, SvobodaStruct.animal, SvobodaStruct.sessionvec, SvobodaStruct.dataname, SvobodaStruct.volume, SvobodaStruct.window, plotcheck);
%     fprintf('binsize_whisker=%d\n',binsize_whisker);
    [~, ~, Ntrialtot] = size(whiskermat);
    if isfield(SvobodaStruct, 'trialvec')
        if ~isempty(SvobodaStruct.trialvec)
            trialvec = SvobodaStruct.trialvec;
        else
            trialvec = 1:Ntrialtot;
        end
    else
      trialvec = 1:Ntrialtot;
    end
    Ntrial = length(trialvec);
    for nt = 1:Ntrial
%       WhiskerTrace.Recording{1,nt} = pi*squeeze(whiskermat(1,:,trialvec(nt)))/180;    % from degree to radian (base angle)
%       WhiskerTrace.Recording{2,nt} = squeeze(whiskermat(2,:,trialvec(nt))); 
%         WhiskerTrace.Recording{1,nt} = zeros(size(squeeze(whiskermat(1,:,trialvec(nt)))));    % set angle=0
        
        % #Tag:输入的修改:1 针对弯曲度原始数据 
%         % 方法1:
%         % 设置弯曲的参数
%         stimulus_start = 1000;  % 弯曲开始时间(单位:ms)
%         stimulus_duration = 10;  % 弯曲持续时间(单位:ms)
%         stimulus_amplitude = 8e-3;  % 弯曲程度(10^-3)
%         tau_fast = 2;  % 快速过程的时间常数 (ms) - 同时用于上升和下降
%         % 获取当前试次的弯曲度数据
%         curvature_data = squeeze(whiskermat(2, :, trialvec(nt)));
%         
%         % 创建一个新的弯曲度数据序列，将其设置为零
%         new_curvature_data = zeros(1, length(curvature_data));
%         
%         % 设置时间序列
%         t = 1:length(curvature_data);  % 时间向量
%         
%         % 创建刺激信号，使用相同的时间常数模拟快速上升和快速下降过程
%         stimulus_signal = zeros(1, length(curvature_data));
%         
%         % 上升部分：快速上升
%         for i = stimulus_start:stimulus_start + stimulus_duration
%             if i <= length(stimulus_signal)
%                 time_from_start = i - stimulus_start;
%                 stimulus_signal(i) = stimulus_amplitude * (1 - exp(-time_from_start / tau_fast));
%             end
%         end
%         
%         % 下降部分：快速下降
%         for i = stimulus_start + stimulus_duration + 1:length(stimulus_signal)
%             time_from_end = i - (stimulus_start + stimulus_duration);
%             stimulus_signal(i) = stimulus_amplitude * exp(-time_from_end / tau_fast);
%         end
%         
%         % 将修改后的弯曲度数据存储回 WhiskerTrace 结构体
%         WhiskerTrace.Recording{2, nt} = stimulus_signal;

        % 方法2:
        % 设置弯曲的参数
%         stimulus_start = 1000;  % 弯曲开始时间(单位:ms)
%         stimulus_total_duration = 4;  % 刺激总持续时间(单位:ms) - 包括上升和下降
%         stimulus_amplitude = 8e-3;  % 弯曲程度(10^-3)
%         tau_fast = 2;  % 快速过程的时间常数 (ms) - 同时用于上升和下降
% 
%         % 获取当前试次的弯曲度数据
%         curvature_data = squeeze(whiskermat(2, :, trialvec(nt)));
% 
%         % 创建一个新的弯曲度数据序列，将其设置为零
%         new_curvature_data = zeros(1, length(curvature_data));
% 
%         % 设置时间序列
%         t = 1:length(curvature_data);  % 时间向量
% 
%         % 创建刺激信号
%         stimulus_signal = zeros(1, length(curvature_data));
% 
%         % 计算上升和下降的时间分配
%         % 这里我们让上升和下降各占一半时间
%         rise_time = stimulus_total_duration / 2;
%         fall_time = stimulus_total_duration / 2;
% 
%         % 上升部分：从0到峰值
%         for i = stimulus_start:stimulus_start + rise_time - 1
%             if i <= length(stimulus_signal)
%                 time_from_start = i - stimulus_start;
%                 stimulus_signal(i) = stimulus_amplitude * (1 - exp(-time_from_start / tau_fast));
%             end
%         end
%         % 这里我们不保持峰值，直接开始下降
%         % 下降部分：从峰值到0
%         for i = stimulus_start + rise_time:stimulus_start + rise_time + fall_time - 1
%             if i <= length(stimulus_signal)
%                 time_from_peak = i - (stimulus_start + rise_time);
%                 stimulus_signal(i) = stimulus_amplitude * exp(-time_from_peak / tau_fast);
%             end
%         end
% 
%         % 将修改后的弯曲度数据存储回 WhiskerTrace 结构体
%         WhiskerTrace.Recording{2,nt} = stimulus_signal;
        
        % 方法三
        % single
%         curvature_data = squeeze(whiskermat(2, :, trialvec(nt)));
%         stimulus_start = 50;      % ms
%         tau = 2;                    % ms
%         stimulus_amplitude = 6e-3;
% 
%         N = length(curvature_data);
%         t = 1:N;
%         stimulus_signal = zeros(1,N);
% 
%         x = (t - stimulus_start);        % 从开始时刻计时
%         idx = x >= 0;
%         u = x(idx) / tau;                % 无量纲时间
%         alpha = u .* exp(1 - u);         % α(t) = (t/τ)*exp(1 - t/τ), t>=0
%         alpha(alpha < 0) = 0;            % 理论上不会<0，这里仅保险
%         stimulus_signal(idx) = stimulus_amplitude * alpha;
%         stimulus_signal(t >(stimulus_start + tau)) = 0;
%         WhiskerTrace.Recording{2,nt} = stimulus_signal;
        
        % multi
        % 1连续型脉冲刺激
%         curvature_data = squeeze(whiskermat(2, :, trialvec(nt)));
%         tau = 6;
%         stimulus_amplitude = 6e-3;
%         N = length(curvature_data);
%         t = 1:N;
%         stimulus_signal = zeros(1,N);
%         num_pulses = 50;
%         for k = 1:num_pulses
%             stimulus_start = k * 50;          % 本段起始：50,100,150,...
%         
%             % 仍用你的写法：x 从本段起点计时，但只在 [0, tau] 内写入
%             x = (t - stimulus_start);
%             idx = (x >= 0) & (x <= tau);           % 仅当前时间段窗口
%         
%             if any(idx)
%                 u = x(idx) / tau;                  % 无量纲时间
%                 alpha = u .* exp(1 - u);           % α(t) = (t/τ)*exp(1 - t/τ)
%                 alpha(alpha < 0) = 0;              % 保险
%                 % 仅对当前窗口做“叠加写入”，不触碰其它时间点
%                 stimulus_signal(idx) = stimulus_signal(idx) + stimulus_amplitude * alpha;
%             end
%         end
%         WhiskerTrace.Recording{2,nt} = stimulus_signal;
        % 2矩形脉冲参数
%         curvature_data = squeeze(whiskermat(2, :, trialvec(nt)));
%         disp('length(curvature_data)');
%         disp(length(curvature_data)); %
        start_time = 100;        % 脉冲起始时间点
        duration = 1;           % 脉冲持续时间
        amplitude_A = 3e-2;     % 脉冲幅值

%         N = length(curvature_data);
        N = SvobodaStruct.window.window(2)-SvobodaStruct.window.window(1)
        t = 1:N;
        stimulus_signal = zeros(1,N);
        num_pulses = 20;

        for k = 1:num_pulses
            % 计算当前脉冲的起始点
            pulse_start = k * start_time;
            
            % 计算当前脉冲的结束点
            pulse_end = pulse_start + duration;
            
            % 找到在当前脉冲持续时间内的索引
            idx = (t >= pulse_start) & (t < pulse_end);
            
            % 在脉冲持续时间内设置幅值为A
            if any(idx)
                stimulus_signal(idx) = stimulus_signal(idx) + amplitude_A;
            end
        end
        WhiskerTrace.Recording{2,nt} = stimulus_signal;
        WhiskerTrace.Recording{1,nt} = zeros(size(stimulus_signal));    % set angle=0

%        WhiskerTrace.Recording{2,nt} = zeros(size(squeeze(whiskermat(2,:,trialvec(nt)))));    % set curvature=0
    end
    WhiskerTrace.binsize = 1; % ms binsize_whisker=2
    WhiskerTrace.quantity = {'base_angle','curvature'};                                 % See Whisker_Recording class for recognized quantities and units
    WhiskerTrace.unit = {'radian','mm-1'};
    
    %% Make / load kernels
    if make_new_thalamic_kernels
        disp('Making new Thalamic kernels')
        % Make new kernels
        thalamickernelfolder = ['..' f 'Make_New_Thinput' f];
        addpath(genpath(thalamickernelfolder));
        
        binsize = 1; % ms
        [Nbx, Nby] = size(Barrelstruct);
        Nbarrel = Nbx*Nby;
        KernelStruct = cell(Nbx, Nby);
        for nbx = 1:Nbx
            for nby = 1:Nby
                KernelStruct{nbx,nby} = make_kernels_angle_curvature_Svobodaexp(SvobodaStruct.Nkernel_ba, SvobodaStruct.Nkernel_c, SvobodaStruct.Nkernel_m, binsize, [savefolder SvobodaStruct.savename '_Thalamic_Kernels_barrel_' num2str(nbx) '_' num2str(nby)]);
            end
        end
        save([savefolder SvobodaStruct.savename '_Thalamic_Kernels'], 'KernelStruct')
        % delete temp files
        for nbx = 1:Nbx
            for nby = 1:Nby 
                delete([savefolder SvobodaStruct.savename '_Thalamic_Kernels_barrel_' num2str(nbx) '_' num2str(nby) '.mat'])
            end
        end   
    else
        % Load kernels from file
         load([savefolder SvobodaStruct.savename '_Thalamic_Kernels']);
         [Nbx, Nby] = size(KernelStruct);
         Nbarrel = Nbx*Nby;
    end
    
    %% Make spike trains for each barreloid
    SpikeTrainStruct = cell(Nbx,Nby);
    SpikeGenStruct = cell(Nbx,Nby);
    nb = 0;
    for nbx = 1:Nbx
        for nby = 1:Nby    
            nb = nb+1;
            disp(['Making Thalamic spike trains for barrel ' num2str(nb) '/' num2str(Nbarrel)])
            SpikeGenStruct{nbx,nby}.refra             = 3; % refractory period (ms) # Tag 不应期,表示每次生成尖峰后,神经元需要等待3毫秒才能再次发放尖峰,是否可以通过修改让响应变得干净？
            SpikeGenStruct{nbx,nby}.Ntrial_pertrace   = 1; % # trials for each Deflection trace 每个胡须动作对应 1 个试次
            SpikeGenStruct{nbx,nby}.binsize           = 1; % binsize spike trains (ms) 每 1 毫秒会记录一次尖峰的发生情况
            if Barrelstruct{nbx,nby}.mainbarrel == 1
                % principal barreloid
                SpikeGenStruct{nbx,nby}.delay             = 0; % (ms)
                SpikeGenStruct{nbx,nby}.scaling           = 1; % Scaling of PSTH
            elseif Barrelstruct{nbx,nby}.mainbarrel == 2
                % secondary barreloid
                SpikeGenStruct{nbx,nby}.delay             = 2.5; % (ms)
                SpikeGenStruct{nbx,nby}.scaling           = .3;
            elseif Barrelstruct{nbx,nby}.mainbarrel == 3
                % tertiary barreloid: no spikes
                SpikeGenStruct{nbx,nby}.delay             = 0; % (ms)
                SpikeGenStruct{nbx,nby}.scaling           = 0;
            end
            % Generate spike trains
            plotyn = 1;
            SpikeTrainStruct{nbx,nby} = kernel_recording_to_spiketrain(WhiskerTrace, KernelStruct{nbx,nby}, SpikeGenStruct{nbx,nby}, [savefolder savename_input '_Thalamic_Spike_Trains_barrel_' num2str(nbx) '_' num2str(nby)], plotyn);
            
        end
    end
    save([savefolder savename_input '_Thalamic_Spike_Trains'], 'SpikeTrainStruct', 'WhiskerTrace', 'SvobodaStruct')
    
    %% delete temp files
    for nbx = 1:Nbx
        for nby = 1:Nby 
            delete([savefolder savename_input '_Thalamic_Spike_Trains_barrel_' num2str(nbx) '_' num2str(nby) '.mat'])
        end
    end    
end