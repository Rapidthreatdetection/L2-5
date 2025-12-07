function Kernels = make_kernels_curvature(Nkernel, binsize, params)
% Example use
% params_ckernels.mf = 10;
% params_ckernels.spreadf = 40;
% params_ckernels.mv = 2; % normal
% params_ckernels.spreadv = 0.2;
% params_ckernels.mk = 80; % normal
% params_ckernels.spreadk = 40;
% params_ckernels.mb = 40; % uniform
% params_ckernels.spreadb = 1;
% params_ckernels.mq = 0; % uniform
% params_ckernels.spreadq = 2;
% params_ckernels.noiseamp = 0;                   % activation functions    
% params_ckernels.kerneltime = (0:100)*binsize;
% Kernels = make_kernels_curvature(Nkernel, binsize, params_ckernels); 

Kernels.ActivationFunction.function = @activation_Petersen_2008;%activation_Petersen_2008
%% 包含随机性
% for nk = 1:Nkernel;
%                 
%     % Make sin kernel #Tag:输入的修改:3 针对卷积核形状
% %     fkernel = (params.mf+params.spreadf*rand)/1000; 
% %     taukernel = (2*rand)*(1/(1*fkernel)); % 原始时间常数范围
% %     if rand<0.5
% %         signkernel = -1;
% %     else
% %         signkernel = 1;
% %     end
% %     kernel = signkernel*(heaviside(params.kerneltime).*(1-heaviside(params.kerneltime-1/fkernel))).*exp(-params.kerneltime/taukernel).*sin(2*pi*fkernel*params.kerneltime);
% %     kernel = kernel./sqrt(sum(kernel.^2)*binsize);
%     % 修改为单峰核 存在随机
%     fkernel = (params.mf+params.spreadf*rand)/1000;
%     taukernel = (0.04 + 0.06*rand) * (1 / fkernel);   % 缩短时间常数 
%     kernel = (params.kerneltime / taukernel) .* exp(-params.kerneltime / taukernel);
%     kernel = kernel / sqrt(sum(kernel.^2)*binsize); % 归一化
% 
%     Kernels.Kernels{nk,1} = kernel;
% 
%     Kernels.ActivationFunction.Params{nk}.v = [params.mv+params.spreadv*randn,0];
%     Kernels.ActivationFunction.Params{nk}.k = [params.mk+params.spreadk*randn,0];
%     Kernels.ActivationFunction.Params{nk}.b = [params.mb+params.spreadb*rand,0];
%     Kernels.ActivationFunction.Params{nk}.q = [params.mq+params.spreadq*rand,0];
%     Kernels.ActivationFunction.Params{nk}.noiseamp = params.noiseamp;
%                 
% end
% 
% Kernels.kerneltime = params.kerneltime;
% 
% % Reorganize
% new_position_vec = randperm(Nkernel);
% Kernels.ActivationFunction.Params = Kernels.ActivationFunction.Params(new_position_vec);
% Kernels.Kernels = Kernels.Kernels(new_position_vec);

%% 固定卷积核和激活函数参数版本
% for nk = 1:Nkernel
%     % 固定频率
%     fkernel = params.mf / 1000;
%     % 固定时间常数
%     taukernel = 0.05 / fkernel;  % 可以根据需要调整常量
%     % 固定卷积核形状
%     kernel = (params.kerneltime / taukernel) .* exp(-params.kerneltime / taukernel);
%     kernel = kernel / sqrt(sum(kernel.^2) * binsize);  % 归一化
%     % 存储卷积核
%     Kernels.Kernels{nk,1} = kernel;
%     % 固定激活函数参数
%     Kernels.ActivationFunction.Params{nk}.v = [params.mv, 0];
%     Kernels.ActivationFunction.Params{nk}.k = [params.mk, 0];
%     Kernels.ActivationFunction.Params{nk}.b = [params.mb, 0];
%     Kernels.ActivationFunction.Params{nk}.q = [params.mq, 0];
%     Kernels.ActivationFunction.Params{nk}.noiseamp = params.noiseamp;
% end
% 
% % 不再随机排序
% Kernels.kerneltime = params.kerneltime;

%% 高斯卷积核版本
t = (0:20)*binsize;          % 时间轴，比如 (0:100)*binsize
% 设定高斯中心和宽度（与 tau 相关，方便你控制形状）
mu    = 10;        % 高斯峰位置
sigma = 1;        % 高斯宽度
k     = 1.5;     % 缩放系数

for nk = 1:Nkernel
    % 固定卷积核形状
    kernel =k * 1/(sigma*sqrt(2*pi)) * exp(-(t - mu).^2 ./ (2 * sigma.^2));
    Kernels.Kernels{nk,1} = kernel;
    % 固定激活函数参数
    Kernels.ActivationFunction.Params{nk}.v = [3, 0];
    Kernels.ActivationFunction.Params{nk}.k = [200, 0];
    Kernels.ActivationFunction.Params{nk}.b = [500, 0];
    Kernels.ActivationFunction.Params{nk}.q = [4, 0];
    Kernels.ActivationFunction.Params{nk}.noiseamp = 0;
end

% 不再随机排序
Kernels.kerneltime = t;

end