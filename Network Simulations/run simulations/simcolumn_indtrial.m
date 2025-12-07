function ConData = simcolumn_indtrial(ConData, initdata, simdata, inputspikes, seeds, WhPara, varargin)
% Simcolumn_indtrial 
% * runs single trial simulation (needed: 
% * sets for this trial initial Vr, Vt dynamic threshold, V and u
% * saves each trial into a mat file

% INPUTS:
% * Data: structure with connectivity data and parameter values from reorganize_conmat
% * initdata: structure with initial conditions: Vrest(mean resting
% membrane potential; can be list); stdVrest; Vthres (mean threshold
% potential; stdVthres; V0 (mean starting membrane potential); stdV0; u0;
% (mean starting u); stdu0; Vthresdyn (dynamic threshold (1) or not (0));
% setVthres.type (how to set the threshold: 'distribution', 'pertype' or
% 'individual') and optional setindcell (.Vrest; .Vthres, .V0, .u0) if one wants to set
% them for each cell individually. In this case, set the desired value in
% setindcell, in a matrix of (Number of neurons NAll x total # simulations) ,  and (optional) ExternalInput = NAll x SimLen
% array with random inputs to each cell 
% * simdata: structure with simulation data: Tsim (total length simulation)
% and timestep, both in ms
% * inputspikes: cell structure with input spikes inputspikes{# iterations per condition}(Ncell x Nspikes_max) 
% * seeds (optional): structure with seeds to re-do simulations: initseed = integer (seed for initial conditions); runseed = integer (seed for running: synaptic failures and
% amplitude)
% * WhPara (optional): structure with base-angles and amplitudes whiskers for whisker
% modulation l2/3 (Crochet 2011)

% GPU version of simulation
% the strategy is to put everything on GPU and run from there
%run the simulation with full matrix operations

% seeds are optional: seeds.initseed for initial values membrane potential
% and threshold, seeds.runseed for running the simulations (synaptic
% failures and random amplitude synapse)




%% neuronal model parameters
a=ConData.Neuron_Para(:,1);
b=ConData.Neuron_Para(:,2);
% b(Data.Cellinfo_All(:,4)>5) = b(Data.Cellinfo_All(:,4)>5) + 0.3;
b(ConData.Cellinfo_All(:,4)==5) = b(ConData.Cellinfo_All(:,4)==5) - 0.15;
c=ConData.Neuron_Para(:,3);
d=ConData.Neuron_Para(:,4);
d(ConData.Cellinfo_All(:,4)>5) = d(ConData.Cellinfo_All(:,4)>5) + 2;

Tau_plas = 120; % time constant for short term synaptic dynamics in ms; 

% parameters for dynamic threshold model
DynThresMat = nan*ones(ConData.NAll, 6);
if initdata.Vthresdyn == 1
    for nt = 1:length(ConData.NtAll)
        cells = find(ConData.Cellinfo_All(:,4)==nt);
        ncells = length(cells);
        DynThresMat(cells,1) = ConData.Neuron_Vt.alpha(nt) *ones(ncells,1);
        DynThresMat(cells,2) = ConData.Neuron_Vt.Vi(nt)    *ones(ncells,1);
        DynThresMat(cells,3) = ConData.Neuron_Vt.Vmin(nt)  *ones(ncells,1);
        DynThresMat(cells,4) = ConData.Neuron_Vt.Ka(nt)    *ones(ncells,1);
        DynThresMat(cells,5) = ConData.Neuron_Vt.Ki(nt)    *ones(ncells,1);
        DynThresMat(cells,6) = ConData.Neuron_Vt.tau(nt)   *ones(ncells,1);
    end
end

%% simulation time
step = simdata.timestep;
simLen = simdata.TSim;
Nsim = simLen/step;

%% Crochet Whisker modulation
% the membrane potential of L2/3 neurons is directly modulated by whisking
% (Crochet et al, Neuron 2011)
% * model: Data.WhiskModel made in 'reorganize_conmat'
% * data: structure WhPara with base-angles (degrees), amplitudes and phase of the whiskers
% Assumptions:
% * Whisker angle correlated inputs to L2/3 population; 
% * not all neurons are modulated by whisking 
% * the distribution of whisking modulation is depth dependent (90% in L3 and 20% in L2); 
% * whisker modulated Neurons have evenly distributed phase preference; 
% * modulation strength is normally distributed with 0.08 +- 0.02 mV/deg; 
% * additional white noise inputs to reduce the r2 (in petersen paper the r2 is low; not
% implemented now)
% * only neurons in the principal barrel are modulated
MIn = zeros(ConData.NAll, Nsim);
if simdata.whiskermodulation
    if (exist('WhPara', 'var') && (~isempty(WhPara) && isfield(ConData, 'WhiskModel')))
        disp('Calculating parameters whisker modulation')
        % needed: whisker base angle trace for this trial in degrees,
        % Amplitude and Phase, it is calculated with
        % WhiskerPara_direct_modulation
        WhPara.MeanAm = nanmean(WhPara.Baseangles_degrees);
        WhPara.MinAm = nanmin(WhPara.Baseangles_degrees);
        WhPara.MeanProtraAm = nanmean(WhPara.Amplitude);
        if ~(length(WhPara.Phase) == Nsim)
            error('Please use a whisker phase-trace that has the same length as the curvature trace')
        else
            MIn = 0.5*WhPara.MeanProtraAm*repmat(ConData.WhiskModel.k,[1,Nsim]).*(sin(repmat(WhPara.Phase,[ConData.NAll,1]) - repmat(ConData.WhiskModel.Pha0, [1,Nsim])));
        end
        if ~(size(MIn,1) == ConData.NAll)
            error('Size matrix direct whisker modulation does not have correct number of neurons')
        end
        if ~(size(MIn,2) == Nsim)
            error('Size matrix direct whisker modulation does not have correct number of time steps')
        end
    elseif ~exist('WhPara', 'var') &&  isfield(ConData, 'WhiskModel')
        disp('Whisker modulation model (Crochet 2011) defined, but no data given; ignoring model.')

    elseif  exist('WhPara', 'var') && ~isfield(ConData, 'WhiskModel')
        disp('Whisker modulation data given, but no model (Crochet 2011) defined; ignoring data.')
    else
        disp('No direct whisker modulation')
    end
end

%% Random inputs to network
if isfield(initdata, 'ExternalInput')
    disp('Using given external input')
    RIn = initdata.ExternalInput;
else
    disp('No external input given')
    RIn = zeros(ConData.NAll,Nsim);
end

%% assign initial parameters for vr, vt, v, u
if (exist('seeds','var') && (~isempty(seeds) && isfield(seeds, 'initseed')))
    rng(seeds.initseed)
else
    rng('shuffle')
    scurr = rng;
    seeds.initseed = scurr.Seed;
end


vr=initdata.Vrest*ones(ConData.NAll,1) + initdata.stdVrest*randn(ConData.NAll,1);     % resting membrane potential

% threshold for spiking
if strcmp(initdata.setVthres.type, 'distribution')
    vt=initdata.Vthres*ones(ConData.NAll,1) + initdata.stdVthres*randn(ConData.NAll,1);   
else
    if strcmp(initdata.setVthres.type, 'pertype')
        % NB Neurons are in order of type, so one can loop over types
        Nt_temp = [0, cumsum(ConData.NtAll)];
        for i = 1:length(Nt_temp) - 1
            vt(Nt_temp(i)+1:Nt_temp(i+1),1) = ConData.Neuron_Vt.avg(i) + initdata.Vthresvar*ConData.Neuron_Vt.std(i)*randn(ConData.NtAll(i), 1);
            % The next parameters are only for the dynamic threshold model (fixed
            % threshold does not use them)
            % currently spike threshold of excitatory neurons are dynamic, while those
            % of inhibitory neurons are static
            VtModel.sl(Nt_temp(i)+1:Nt_temp(i+1), 1) = ConData.Neuron_Vt.sl(i);
            VtModel.a(Nt_temp(i)+1:Nt_temp(i+1), 1) = ConData.Neuron_Vt.a(i);
        end
    elseif strcmp(initdata.setVthres.type, 'individual')
        vt = initdata.setindcell.Vthres;
    else
        error('Choose method to set threshold: initdata.setVthres.type should be distribution, pertype or individual')
    end
end
VT0 = vt;
    
% initial values for neuron model
v=initdata.V0.*ones(ConData.NAll,1)+initdata.stdV0.*randn(ConData.NAll,1); % Initial values of v
u=initdata.u0.*ones(ConData.NAll,1)+initdata.stdu0.*randn(ConData.NAll,1); % Initial values of u
% required initial value for computation

if isfield(initdata, 'setindcell')
    namevec = fieldnames(initdata.setindcell);
    for ii = 1:length(namevec)
        if strcmp(namevec{ii}, 'Vrest')
            vr = initdata.setindcell.Vrest;
        elseif strcmp(namevec{ii}, 'Vthres')
            % do nothing; already done
        elseif strcmp(namevec{ii}, 'V0')
            v = initdata.setindcell.V0;
        elseif strcmp(namevec{ii}, 'u0')
            u = initdata.setindcell.u0;
        else
            error('Non-existing field for setindcell')
        end
    end
end

V0 = v;
U0 = u;

%% Preallocate
inputsc = single(zeros(ConData.NIn,1));
modelsc = single(zeros(ConData.NAll,1));
inputin = single(zeros(ConData.NAll,1));
modelin = single(zeros(ConData.NAll,1));
Finput = [];
Fmodel = [];


Vminus = v;
V = nan*ones(ConData.NAll, Nsim);
U = nan*ones(ConData.NAll, Nsim);
if initdata.Vthresdyn == 1
    VT = nan*ones(ConData.NAll, Nsim);
end

% Preallocate L5 time series storage (ensure N_L5 exists)
if ~exist('N_L5','var') || isempty(N_L5)
    N_L5 = 1000; % default fallback (will be overwritten later if defined)
end
V_L5 = nan*ones(N_L5, Nsim);
U_L5 = nan*ones(N_L5, Nsim);
if initdata.Vthresdyn == 1
    VT_L5 = nan*ones(N_L5, Nsim);
end
dS = zeros(ConData.NAll, 1);
vt_sp=zeros(ConData.NAll,1);
vtct=zeros(ConData.NAll,3);
modelvt=zeros(ConData.NAll, 1);
modelspt=zeros(ConData.NAll, 1);

inputspikes(size(inputspikes,1)+1:ConData.NIn,:) = 0;
inputspikes(ConData.NIn+1:end, :) = [];

%% ------------------ L5 subnetwork (minimal invasive implementation) ------------------
% User requested L5 with 85% pyramidal (E, bursting) and 15% PV (I, fast-spiking).
% Default size set by user confirmation
N_L5 = 1000; % total L5 cells (user-specified)
fracE_L5 = 0.85;
N_L5_E = round(N_L5 * fracE_L5);
N_L5_I = N_L5 - N_L5_E;
fired_L5_prev = [];   % buffer: 上一步的 L5 spikes
L5_refrac_timer = zeros(N_L5,1);   % refractory counter in ms
L5_refrac_period = 3;              % absolute refractory = 3 ms


% Pyramidal (bursting) excitatory cell parameters (L5E)
aE = 0.01;  % Self-adaptation
bE = 0.5;   % Recovery variable
cE = -65;   % Membrane potential after spike
dE = 8;     % Spike-triggered adaptation

% PV (FS) inhibitory cell parameters (L5I)
aI = 0.08;   % Self-adaptation parameter
bI = 0.85;   % Recovery variable gain
cI = -65;   % Membrane potential after spike
dI = 10;     % Spike-triggered adaptation

% build parameter vectors for L5
a_L5 = [aE*ones(N_L5_E,1); aI*ones(N_L5_I,1)];
b_L5 = [bE*ones(N_L5_E,1); bI*ones(N_L5_I,1)];
c_L5 = [cE*ones(N_L5_E,1); cI*ones(N_L5_I,1)];
d_L5 = [dE*ones(N_L5_E,1); dI*ones(N_L5_I,1)];

% initial membrane variables for L5 (sampled similarly to main network)
% ----- L5 excitatory (pyramidal, bursting) and inhibitory (PV) - unified defaults -----
% 为提高 L5_exc 的放电性，统一 E/I 的起始电位/阈值等（a/b/c/d 保持不同）
Vrest_common = -72;   % common resting potential for L5
V0_common    = -70;   % common initial membrane potential
stdV0_common = 2;

Vrest_E = Vrest_common;
V0_E    = V0_common;
stdV0_E = stdV0_common;

Vrest_I = Vrest_common;
V0_I    = V0_common;
stdV0_I = stdV0_common;

% Build vectors
vr_L5 = [Vrest_E*ones(N_L5_E,1); Vrest_I*ones(N_L5_I,1)];
v_L5  = [V0_E*ones(N_L5_E,1) + stdV0_E*randn(N_L5_E,1);
         V0_I*ones(N_L5_I,1) + stdV0_I*randn(N_L5_I,1)];
u_L5  = zeros(N_L5,1);   % FS 和 bursting 初始 u=0 是对的


% threshold for L5 
% 将 E/I 起始阈值统一以避免抑制侧偏（较高的阈值会抑制 E 的放电）
vt_L5_common = -40; % common spike threshold (mV)
vt_L5_E = vt_L5_common * ones(N_L5_E, 1);
vt_L5_I = vt_L5_common * ones(N_L5_I, 1);
vt_L5 = [vt_L5_E; vt_L5_I];
VT0_L5 = vt_L5;

% spike tracking for L5
modelsc_L5 = zeros(N_L5,1,'single');
modelspt_L5 = zeros(N_L5, 6000, 'single'); % prealloc small matrix; will grow if needed

% synaptic state for simplified L5 synapses (per-post conductance/current)
% Using simplified single-exponential synapse for speed/clarity
tau_syn_E = 3; % ms
tau_syn_I = 6; % ms
s_L5 = zeros(N_L5,1); % synaptic drive (sum of contributions)
% per-post tau for L5
tau_post_L5 = [tau_syn_E*ones(N_L5_E,1); tau_syn_I*ones(N_L5_I,1)];

% Connection probabilities (defaults)
p_L5E_to_L23E = 0.12; % L5 excitatory -> L2/3 excitatory
p_L5I_to_L23E = 0.08; % L5 inhibitory -> L2/3 excitatory


% amplitude defaults (increased to give stronger L5 drive and feedback)
% 注意：这些值是经验性调整，单位与模型中 I/g 的标度相关，可在需要时微调
amp_EE = 0.12; amp_EI = 0.20; amp_IE = -0.40; amp_II = -0.15;

% construct mapping from existing network to L5
post_L5_E_idx = 1:N_L5_E; post_L5_I_idx = N_L5_E+1:N_L5;


% ---------------------------------------------------------------------
% 统一化 L5 连接构建：使用 L5ConnParams（在上面已定义）来抽样
% - p_to_L5: 行 = presyn 类型编号（用户表 1..15），列 = [to L5_E, to L5_I]
% - p_L5_to_L23 and p_internal 在 L5ConnParams 中提供 L5->其他和 L5 内部概率
% 如果 L5ConnParams 中未定义某些域，则回退到原有的 p_* 变量
% 确保 L5ConnParams 在此文件内已定义（仅本文件内使用，不写回 ConData）
if ~exist('L5ConnParams','var') || isempty(L5ConnParams)
    maxType_for_table = 15;
    p_to_L5 = zeros(maxType_for_table, 2);
    % 默认建议值（行 = 用户给定类型序号 1..15，列 = [to L5_E, to L5_I])
    % p_to_L5(1,:)  = [0.20, 0]; % L4 Exc 1 (increased drive to L5_E)
    p_to_L5(2,:)  = [0.03, 0.05]; % L4 Exc 2 (increased drive to L5_E)
    p_to_L5(3,:)  = [0.05, 0.10]; % L4 PV+ FS
    % p_to_L5(4,:)  = [0.15, 0.20]; % L4 Non-FS
    p_to_L5(5,:)  = [0.08, 0.10]; % L2/3 Exc (increased drive to L5_E)
    p_to_L5(6,:)  = [0.05, 0.10]; % L2/3 PV+
    % p_to_L5(7,:)  = [0.50, 0.20]; % L2/3 PV+
    % p_to_L5(8,:)  = [0.45, 0.35]; % PV+ bursting
    % p_to_L5(9,:)  = [0.25, 0.20]; % Martinotti / SOM
    % p_to_L5(11,:) = [0.10, 0.30]; % VIP+/CR-
%     p_to_L5(12,:) = [0.08, 0.25]; % CR+ bipolar
%     p_to_L5(14,:) = [0.08, 0.25]; % CR+ multipolar
    % p_to_L5(15,:) = [0.20, 0.20]; % Neurogliaform

    L5ConnParams.p_to_L5 = p_to_L5;
    L5ConnParams.maxType = maxType_for_table;
    % L5 -> L2/3 probabilities (E row; I row)
    L5ConnParams.p_L5_to_L23 = [p_L5E_to_L23E; p_L5I_to_L23E];
    % internal L5 connectivity: rows pre E/I, cols post E/I
    % set more balanced internal probabilities to avoid runaway inhibition
    L5ConnParams.p_internal = [0.05, 0.25; 0.15, 0.25]; % [E->E, E->I; I->E, I->I]
    L5ConnParams.labels = { 'type1','type2','type3','type4','type5','type6','type7','type8','type9','type10','type11','type12','type13','type14','type15' };
end

if ~isfield(L5ConnParams, 'p_L5_to_L23')
    L5ConnParams.p_L5_to_L23 = [p_L5E_to_L23E; p_L5I_to_L23E]; % [Erow; Irow]
end
if ~isfield(L5ConnParams, 'p_internal')
    L5ConnParams.p_internal = [0.08, 0.40; 0.40, 0.20];
    % rows: pre E(1)/I(2), cols: post E(1)/I(2)
end

% create cell arrays: for each existing pre neuron, list of L5 posts it connects to
pre2L5 = cell(ConData.NAll,1);
for pi = 1:ConData.NAll
    typ = ConData.Cellinfo_All(pi,4);
    if typ < 1 || typ > L5ConnParams.maxType
        % unknown type -> no connections to L5 by default
        pre2L5{pi} = [];
        continue
    end
    % probability to L5_E and to L5_I based on presyn type
    pE = L5ConnParams.p_to_L5(typ,1);
    pI = L5ConnParams.p_to_L5(typ,2);

    % sample targets independently for E and I subpopulations
    maskE = rand(N_L5_E,1) < pE;
    maskI = rand(N_L5_I,1) < pI;
    pre2L5{pi} = [post_L5_E_idx(maskE)'; post_L5_I_idx(maskI)'];
end
% end of pre2L5 construction (now using L5ConnParams)

% create L5 -> L2/3 mapping (only target L2/3 excitatory cells, type 5)
post_idx_L23E = find(ConData.Cellinfo_All(:,4) == 5);
L5toL23 = cell(N_L5,1);
% get probabilities for E->L23E and I->L23E from L5ConnParams
p_L5E_to_L23E_use = L5ConnParams.p_L5_to_L23(1);
p_L5I_to_L23E_use = L5ConnParams.p_L5_to_L23(2);
for j = 1:N_L5
    if j <= N_L5_E
        mask = rand(length(post_idx_L23E),1) < p_L5E_to_L23E_use;
        L5toL23{j} = post_idx_L23E(mask);
    else
        mask = rand(length(post_idx_L23E),1) < p_L5I_to_L23E_use;
        L5toL23{j} = post_idx_L23E(mask);
    end
end

% internal L5 connectivity (cell -> list)
L5internal = cell(N_L5,1);
for j = 1:N_L5
    if j <= N_L5_E
        pEE = L5ConnParams.p_internal(1,1);
        pEI = L5ConnParams.p_internal(1,2);
        maskE = rand(N_L5_E,1) < pEE;
        maskI = rand(N_L5_I,1) < pEI;
        L5internal{j} = [post_L5_E_idx(maskE)'; post_L5_I_idx(maskI)'];
    else
        pIE = L5ConnParams.p_internal(2,1);
        pII = L5ConnParams.p_internal(2,2);
        maskE = rand(N_L5_E,1) < pIE;
        maskI = rand(N_L5_I,1) < pII;
        L5internal{j} = [post_L5_E_idx(maskE)'; post_L5_I_idx(maskI)'];
    end
end

% 在 L5internal 构建完成后运行（一次）
num_E_to_I = 0;
for j = 1:N_L5_E
    posts = L5internal{j};
    num_E_to_I = num_E_to_I + sum(posts > N_L5_E); % posts > N_L5_E 属于 I
end
fprintf('L5 E->I edges: %d (per E ~ %.2f)\n', num_E_to_I, num_E_to_I / N_L5_E);

% storage for final connectivity (simple adjacency lists)
final_connectivity_L5.pre2L5 = pre2L5;
final_connectivity_L5.L5toL23 = L5toL23;
final_connectivity_L5.L5internal = L5internal;

% (removed nested helper function to avoid MATLAB nested/local function conflicts)

%% ---------------------------------------------------------------------------

%% matrix for running variables
% these varibles run in GPU
STD_input = (ConData.PMat_IntoAll.STD);
STD_model = (ConData.PMat_AlltoAll.STD);

g_input = (ConData.PMat_IntoAll.g);
s_input = (ConData.PMat_IntoAll.s);
g_model = (ConData.PMat_AlltoAll.g);
s_model = (ConData.PMat_AlltoAll.s);

CN_input = (ConData.PMat_IntoAll.CN);
CN_model = (ConData.PMat_AlltoAll.CN);

Trise_input = (ConData.PMat_IntoAll.Trise);
Tfall_input = (ConData.PMat_IntoAll.Tfall);
Trise_model = (ConData.PMat_AlltoAll.Trise);
Tfall_model = (ConData.PMat_AlltoAll.Tfall);

inputspt_temp = inputspikes;

%% variables needed to run STDP
spikepairs_IntoAll = zeros(length(ConData.PMat_IntoAll.Am), 4);
spikepairs_AlltoAll = zeros(length(ConData.PMat_AlltoAll.Am), 4);
Ampost_IntoAll = ConData.PMat_IntoAll.Am;
Ampost_AlltoAll = ConData.PMat_AlltoAll.Am;
initial_connectivity = Ampost_AlltoAll;
% Ampre_IntoAll = Data.PMat_IntoAll.Am;
% Ampre_AlltoAll = Data.PMat_AlltoAll.Am;
% celltype is organized as postType, postEI, preType, preEI
Celltype_IntoAll = [ConData.Cellinfo_All(ConData.PMat_IntoAll.preCell(:,1), 4), ConData.Cellinfo_All(ConData.PMat_IntoAll.preCell(:,1), 6), ...
    ConData.Cellinfo_In(ConData.PMat_IntoAll.preCell(:,2), 4), ConData.Cellinfo_In(ConData.PMat_IntoAll.preCell(:,2), 6)];
Celltype_AlltoAll = [ConData.Cellinfo_All(ConData.PMat_AlltoAll.preCell(:,1), 4), ConData.Cellinfo_All(ConData.PMat_AlltoAll.preCell(:,1), 6), ...
    ConData.Cellinfo_All(ConData.PMat_AlltoAll.preCell(:,2), 4), ConData.Cellinfo_All(ConData.PMat_AlltoAll.preCell(:,2), 6)];

%% simulation starts
if exist('seeds','var') && isfield(seeds, 'runseed')
    rng(seeds.runseed)
else
    rng('shuffle')
    scurr = rng;
    seeds.runseed = scurr.Seed;
end
% seeds for synaptic failures and Am
timestepseed_input = randi(2^32, 1, Nsim);
timestepseed_model = randi(2^32, 1, Nsim);

% feedback buffer from L5 to main network (per-post synaptic drive)
s_L23_feedback = zeros(ConData.NAll,1);

for t=1:Nsim % simulation of Ni ms

    tm = t*step;

    % Optional: reset neuron states at specified intervals (ms)
    % If simdata.reset_interval is set (>0), then at each multiple of that
    % interval (excluding t=0) we restore membrane variables to their
    % initial values V0/U0 (as computed before the loop) and reset spike
    % threshold trackers. This preserves RNG and connectivity behavior.
    if isfield(simdata, 'reset_interval') && ~isempty(simdata.reset_interval) && simdata.reset_interval>0
        ri = simdata.reset_interval;
        tol = max(1e-9, step/2);
        if tm>0 && abs(mod(tm, ri)) < tol
            % 恢复神经元膜电位、恢复突触状态、输入状态等
            v = V0;
            u = U0;
            Vminus = v;
            vt_sp = zeros(ConData.NAll,1);
            vtct = zeros(ConData.NAll,3);
            inputin = single(zeros(ConData.NAll,1));
            modelin = single(zeros(ConData.NAll,1));
            Finput = [];
            Fmodel = [];
            s_input = (ConData.PMat_IntoAll.s);
            s_model = (ConData.PMat_AlltoAll.s);
            g_input = (ConData.PMat_IntoAll.g);
            g_model = (ConData.PMat_AlltoAll.g);
            STD_input = (ConData.PMat_IntoAll.STD);
            STD_model = (ConData.PMat_AlltoAll.STD);
            % 如果有动态阈值
            if initdata.Vthresdyn == 1
                vt = VT0;
            end
            % 可选：如果有外部输入、白噪声等，也可在此处重置
            % RIn(:,t) = ...  % 通常外部输入不变
            disp(['Reset neuron & synapse states at t = ' num2str(tm) ' ms (reset_interval = ' num2str(ri) ' ms)']);
        end
    end


    [ind_input,ind_sc]=find(inputspt_temp <= tm & inputspt_temp > 0);
    inputspt_temp(sub2ind(size(inputspikes), ind_input, ind_sc)) = 0;
    
    inputpre = [];
    
    if isempty(ind_input) == 0
        for i = 1:length(ind_input)
            inputpre = [inputpre; ConData.PMat_IntoAll.postIdx{ind_input(i)}];
        end
        inputsc(ind_input) = inputsc(ind_input)+1;
        clear ind_input ind_sc
        
        %now we can directly get all the parameters (for now just F4)
        %we need Am, Trise, Tfall, Plas and synaptic delay
        %Am will be generated from CV and Am matrix
        
        %  check for input spikes from input layer 
        [Finput, preSpikes_input, timestepseed_input(t)] = simcolumn_synapse_onestep(inputpre, Ampost_IntoAll, STD_input, ConData.PMat_IntoAll.CV, ...
            ConData.PMat_IntoAll.Fail, ConData.PMat_IntoAll.Delay, ConData.PMat_IntoAll.preID, ConData.NAll, Finput, tm, timestepseed_input(t));

        
    end
    
    % update running matrix for the synaptic current model
    if isempty(Finput) == 0
        [s_input, STD_input, Finput] = simcolumn_updatesynap_ver2(s_input, ...
            STD_input, ConData.PMat_IntoAll.Plas, Finput, tm);
    end
    
    idx = find(v>vt);%neuron reaches spike threshold
    idx_spingen = idx(vtct(idx, 3) == modelsc(idx));
    vt_sp(idx_spingen) = vt(idx_spingen);
    vtct(idx_spingen, 1) = 1; % keep track of neurons with spike generation in progress
    vtct(idx_spingen, 2) = tm; % time of first spike threshold cross
    vtct(idx_spingen, 3) = vtct(idx_spingen, 3) + 1;
    
    modelfired = find(v>0);    % find cells those are firing in layer2/3
    
    % check failed spikes
    idx = tm - vtct(:, 2) > 15 & vtct(:, 1) == 1;
    vtct(idx, 1) = 0;
    vtct(idx, 3) = vtct(idx, 3) - 1;
    
    modelpre= [];
    modelpost = [];
    
    if isempty(modelfired) == 0
        for i = 1:length(modelfired)
            % modelpre uses pre-post direction, need to use postIdx (which
            % contains all postsynaptic partner to one presynaptic neuron)
            % modelpost uses post-pre direction to calculate STDP, so use
            %  preIdx instead 
            modelpre  = [modelpre; ConData.PMat_AlltoAll.postIdx{modelfired(i)}];
            modelpost = [modelpost; ConData.PMat_AlltoAll.preIdx{modelfired(i)}];
            modelsc(modelfired(i)) = modelsc(modelfired(i)) + 1;
            modelspt(modelfired(i), modelsc(modelfired(i))) = tm - step;
            modelvt(modelfired(i), modelsc(modelfired(i))) = vt_sp(modelfired(i));
        end
        % reset spike generation tracker
        vtct(modelfired, 1) = 0;
        
        % check for input spikes from simulated layer
        [Fmodel, preSpikes_model, timestepseed_model(t)] = simcolumn_synapse_onestep(modelpre, Ampost_AlltoAll, STD_model, ...
        ConData.PMat_AlltoAll.CV, ConData.PMat_AlltoAll.Fail, ConData.PMat_AlltoAll.Delay, ...
        ConData.PMat_AlltoAll.preID, ConData.NAll, Fmodel, tm, timestepseed_model(t));
        
                

        % update spike pairs, which is used to run STDP 
        spikepairs_AlltoAll(preSpikes_model(:,1), 1) = preSpikes_model(:,2);
        spikepairs_AlltoAll(preSpikes_model(:,1), 3) = 1;
        % post synaptic spike timing
        spikepairs_AlltoAll(modelpost, 2) = tm;
        spikepairs_AlltoAll(modelpost, 4) = 1;

        % ------------------ Route model spikes to L5 subnetwork ------------------
        % For each model neuron that fired, deliver a simple synaptic event to L5 posts
        if exist('pre2L5','var')
            for pf = 1:length(modelfired)
                preid = modelfired(pf);
                posts = pre2L5{preid};
                if ~isempty(posts)
                    % only excitatory pres (ConData.Cellinfo_All(:,6)==1) should drive L5 excitatory per spec
                    preEIflag = ConData.Cellinfo_All(preid, 6);
                    if preEIflag == 1
                        % add excitatory amplitude
                        s_L5(posts) = s_L5(posts) + amp_EE;
                    else
                        % inhibitory 上层输入应该给 L5E 更小 或 不给
                        % 为避免抑制 E 过强，需要降低
                        s_L5(posts) = s_L5(posts) + 0.0; % 不给 I input
                    end
                end
            end
        end
        % --- Debug logging (插入点：路由完外部 spikes -> s_L5 后，L5 更新前) ---
%         if (mod(t,10) == 1) % 每 10 个时间步打印一次，避免太多输出
%             % 总体看 L5 s 输入被加到 E 和 I 上的瞬时值
%             sE_total = sum(s_L5(1:N_L5_E));
%             sI_total = sum(s_L5(N_L5_E+1:end));
%             fprintf('t=%.2f ms | sE_total=%.4g | sI_total=%.4g | mean_vE=%.4g mean_vI=%.4g\n', ...
%                 tm, sE_total, sI_total, mean(v_L5(1:N_L5_E)), mean(v_L5(N_L5_E+1:end)));
%         end
        % ------------------ end routing ------------------
    end
    
    %% STDP rule
    % right now only run it for L4-L23 network; i.e. assume no plasticity
    % in thalamo-cortical synapses
    if simdata.STDP
        spikeidx_AlltoAll = [modelpre; modelpost];
        if ~isempty(spikeidx_AlltoAll)
            [spikepairs_AlltoAll, Ampost_AlltoAll] = simcolumn_Plasticity_STDP...
                (spikepairs_AlltoAll, spikeidx_AlltoAll, Celltype_AlltoAll, Ampost_AlltoAll, 3);
        end
    end
    
    
    %% update running matrix for the synaptic current model
    if isempty(Fmodel) == 0
        [s_model, STD_model, Fmodel] = simcolumn_updatesynap_ver2(s_model, ...
            STD_model, ConData.PMat_AlltoAll.Plas, Fmodel, tm);
    end
    
    
    
    % reset the neuron model after firing
    inter=find(v>0);
    v(inter)=c(inter);
    u(inter)=u(inter)+d(inter);

%% ------------------ L5 subnetwork dynamics (per time step) ------------------

    % (1) decay synaptic drives
    if exist('tau_post_L5','var')
        s_L5 = s_L5 + step .* ( - s_L5 ./ tau_post_L5 );
    else
        s_L5 = s_L5 + step .* ( - s_L5 ./ tau_syn_E );
    end

    % (2) deliver internal E->I and I->E synaptic input from LAST timestep
    if ~isempty(fired_L5_prev)
        for jf = 1:length(fired_L5_prev)
            jid = fired_L5_prev(jf);

            posts = L5internal{jid};
            if ~isempty(posts)
                isEpost = posts <= N_L5_E;
                postsE = posts(isEpost);
                postsI = posts(~isEpost);

                if jid <= N_L5_E
                    % E pre -> E, I
                    if ~isempty(postsE)
                        s_L5(postsE) = s_L5(postsE) + amp_EE;
                    end
                    if ~isempty(postsI)
                        s_L5(postsI) = s_L5(postsI) + amp_EI;
                    end
                else
                    % I pre -> E, I
                    if ~isempty(postsE)
                        s_L5(postsE) = s_L5(postsE) + amp_IE;
                    end
                    if ~isempty(postsI)
                        s_L5(postsI) = s_L5(postsI) + amp_II;
                    end
                end
            end
        end
    end

    % --- Apply refractory: neurons in refractory do NOT integrate membrane ---
    mask_refrac = (L5_refrac_timer > 0);
    v_L5(mask_refrac) = c_L5(mask_refrac);     % clamp voltage at reset value
    u_L5(mask_refrac) = u_L5(mask_refrac);     % adaptation allowed to decay normally

    % (3) integrate L5 membrane voltage
    vtemp_L5 = v_L5;
    v_L5 = v_L5 + step .* ( 0.04 .* (v_L5 - vr_L5) .* (v_L5 - vt_L5) - u_L5 + s_L5 );
    u_L5 = u_L5 + step .* a_L5 .* ( b_L5 .* (vtemp_L5 - vr_L5) - u_L5 );

    % (4) detect NEW spikes
    % block spikes for neurons still in refractory period
    fired_L5 = find(v_L5 > 0 & L5_refrac_timer == 0);


    % -------- 正确统计 L5 spike 的位置 --------
    %spikes_E_thisstep = sum(fired_L5 <= N_L5_E);
    %spikes_I_thisstep = sum(fired_L5 >  N_L5_E);
    %fprintf("    L5 spikes: E=%d I=%d\n", spikes_E_thisstep, spikes_I_thisstep);

    if ~isempty(fired_L5)
        for jf = 1:length(fired_L5)
            jid = fired_L5(jf);

            % log spike
            modelsc_L5(jid) = modelsc_L5(jid) + 1;
            if modelsc_L5(jid) > size(modelspt_L5,2)
                modelspt_L5(:, size(modelspt_L5,2)+1:modelsc_L5(jid)) = 0;
            end
            modelspt_L5(jid, modelsc_L5(jid)) = tm - step;

            % feedback to L2/3
            posts_main = L5toL23{jid};
            if ~isempty(posts_main)
                if jid <= N_L5_E
                    s_L23_feedback(posts_main) = s_L23_feedback(posts_main) + amp_EE;
                else
                    s_L23_feedback(posts_main) = s_L23_feedback(posts_main) + amp_IE;
                end
            end
        end

        % reset L5 membrane states
        v_L5(fired_L5) = c_L5(fired_L5);
        u_L5(fired_L5) = u_L5(fired_L5) + d_L5(fired_L5);
        L5_refrac_timer(fired_L5) = L5_refrac_period;
    end

    % (5) NOW update fired_L5_prev (this must be last)
    fired_L5_prev = fired_L5;

    % (6) decay feedback buffer to L2/3
    tau_feedback = 6;
    s_L23_feedback = s_L23_feedback + step .* ( - s_L23_feedback ./ tau_feedback );

    % (7) record L5 states
    V_L5(:, t) = v_L5;
    U_L5(:, t) = u_L5;
    if exist('VT_L5','var')
        VT_L5(:, t) = vt_L5;
    end
    L5_refrac_timer = max(L5_refrac_timer - step, 0);


    % running synaptic current model with Euler method
    % the differential function is
    %dg/dt = ((tau2 / tau1) ** (tau1 / (tau2 - tau1))*s-g)/tau1 
    %ds/dt = -s/tau2 
    g_input = g_input + step*((CN_input.*s_input-g_input)./Tfall_input);
    s_input = s_input + step*( -s_input ./ Trise_input);
    
    g_model = g_model + step*((CN_model.*s_model-g_model)./Tfall_model);
    s_model = s_model + step*( -s_model ./ Trise_model);
    
    % running STD model with Euler method
    % the differential function is
    % dW/dt = -(W-1)/tau_plas
    STD_input = STD_input + step * (-(STD_input-1)./Tau_plas);
    STD_model = STD_model + step * (-(STD_model-1)./Tau_plas);
    
    % calculate synaptic current
    inputin = simcolumn_calcualteI_Mat_ver2(g_input, ConData.PMat_IntoAll.preEI, v);
    modelin = simcolumn_calcualteI_Mat_ver2(g_model, ConData.PMat_AlltoAll.preEI, v);
    

    % include feedback from L5 to main network
    I = inputin + modelin + s_L23_feedback;
    
    % running neuron model with Euler method
    v_temp = v;
    v=v+step*(0.04*(v-vr).*(v-vt)-u + I + RIn(:,t) + MIn(:,t));
    u=u+step*a.*(b.*(v_temp-vr)-u); 
    
    V(:,t) = v;
    U(:,t) = u;
    
    %% spike threshold model from Fontaine et al Plos Comp Bio 2014
    if initdata.Vthresdyn == 1    
        Vt_thetainf = DynThresMat(:,1).*(v - DynThresMat(:,2))+DynThresMat(:,3)+DynThresMat(:,4).*log(1+exp((v - DynThresMat(:,2))./DynThresMat(:,5)));
        vt = vt + (Vt_thetainf - vt).*step./DynThresMat(:,6);
        
        VT(:,t) = vt;
    end
    % fix spike threshold for neurons under spike generation: 'freeze' the vt value once v crosses the vt. Otherwise the vt keeps changing and sometimes results in false negative spikes 
    vt(vtct(:, 1) == 1) = vt_sp(vtct(:,1) == 1);
    
end


%% Save
final_connectivity = Ampost_AlltoAll;

V = single(V);
U = single(U);
if initdata.Vthresdyn == 1
    VT = single(VT);
end




cellinfo_all = ConData.Cellinfo_All;
cellinfo_input = ConData.Cellinfo_In;
Neuron_Para = ConData.Neuron_Para;

a = exist(ConData.savefolder, 'dir');
if ~(a==7) 
    mkdir(ConData.savefolder)
end

if isfield(simdata, 'whattosave')
    whattosave = simdata.whattosave;
else
    whattosave = {'Neuron_Para', 'Tau_plas', 'cellinfo_all', 'cellinfo_input', ...
        'modelsc', 'modelspt', 'inputsc', 'inputspikes', 'simLen', 'step', ...
        'V0','U0','VT0', 'vr', 'seeds', 'timestepseed_input', ...
        'timestepseed_model', 'initdata','simdata', 'final_connectivity', 'initial_connectivity', ...
        'V_L5', 'U_L5', 'VT_L5', 'modelsc_L5', 'modelspt_L5', 'final_connectivity_L5'};
end

% If user provided simdata.whattosave, ensure L5 variables are included so
% that generated files always contain L5 results when present in workspace.
if isfield(simdata, 'whattosave') && ~isempty(simdata.whattosave)
    l5vars = {'V_L5', 'U_L5', 'VT_L5', 'modelsc_L5', 'modelspt_L5', 'final_connectivity_L5'};
    for k = 1:length(l5vars)
        if ~ismember(l5vars{k}, whattosave)
            whattosave{end+1} = l5vars{k};
        end
    end
end

if initdata.Vthresdyn == 1
    thresholdname = ['dynthreshold_set'];
else
    thresholdname = ['fixthreshold_set'];
end


savename = [ConData.savefolder,ConData.FnametoSave '_Simcolumn_' thresholdname '_' initdata.setVthres.type '_simulation_' num2str(initdata.nsim)];
savename = strrep(savename, '\', '/');
save(savename, whattosave{:}, '-v7.3')

%% Overwrite Am in ConData (this changes in case of STDP, otherwise it stays the same)
if simdata.STDP
    ConData.PMat_AlltoAll.Am = final_connectivity;
end

