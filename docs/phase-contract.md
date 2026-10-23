# Experimental PF-FLIP phase contract

2026-10-23. This is a uniform-grid experiment for the unresolved NumPy/Bifrost
2D stability gate, not a new production default. The old graph (`phase_model=0`)
keeps its weighted phase-tag path, and the old `Sim` keeps its original path.

The two explicit experimental modes share the following mathematical definition
between an independent `PhaseContractSim` and the actual fused C++ graph:

1. Deposit **raw sums** of primary-particle mass with PF-FLIP Eq.6,
   `w=max(1-distance²/r²,0)³`, with `r=dx`. Auxiliary spray has zero mass and
   contributes to no phase or velocity channel. Momentum uses the same kernel.
2. Use Eq.7, `phi=min(sqrt(max(raw-log(rho_l/rho_g)*rho_g*rho0,0) /
   (alpha_phi*rho0*rho_l)),1)`. `rho0` is the calibrated kernel sum of a uniform
   reference seeding, not a per-face normalization by the instantaneous particle
   count. `alpha_phi=1` is the paper's default.
3. Mode 1 stops here. Mode 2 also deposits mass at cell centres with the same
   kernel and evaluates Eq.7 there. Cells with phi >= 0.5 are liquid. Interior
   faces between two liquid cells are set to 1; between two air cells to 0;
   differing neighbours retain their face-sampled Eq.7 value. Closed boundary
   faces retain their raw phase value but have zero pressure coupling/velocity.

The final pure-face operation is required by PF-FLIP §3.3. **The separate
cell-centred splat is our explicit sampling choice**: the paper states the cell
tagging rule but does not specify how to derive cell samples from face samples.
We do not claim this choice is the authors' unpublished implementation. The
tested 2D seeding is the same regular 2x2 subcell lattice used for face calibration.

The momentum and mean-mass float32 ports retain their existing semantics.
Experimental phase ports instead carry the compressed/snap-adjusted phase,
computed from double raw accumulators before a float32 port conversion. The
independent Python implementation rounds particle momentum and the same field
boundaries to float32 explicitly when `port_rounding=True`. Disabling that switch
is a separate precision experiment, never an unrecorded change to a comparison.
Momentum is retained as its own loop state, rounded from the unrounded final
particle velocity, just as the C++ output ports are. Recomputing it from the
already-rounded velocity creates a second rounding and is not the same boundary.
The first 6-second diagnostic arms predate this harness correction and remain
retained as failed comparisons; they must not be presented as corrected runs.

ST weights, adaptive particle volumes and kernel support, interior solid support,
moving boundaries, forcing and adaptive pressure grids require separate contracts
and gates. This first implementation **refuses** those combinations rather than
applying the uniform calibration silently. Their intended raw deposit is a sum of
`kernel * rho_phase * represented_volume * temporal_weight`, but the calibration,
kernel support near solid/refinement boundaries and pure-cell tagging must be
verified before these modes are enabled. Existing legacy features remain available
with phase_model=0.

Acceptance starts with Eq.7 analytic points, the pure-face truth table, actual
P2G field comparison, full graph single steps, and unsupported-mode refusal. Then
compare mode 1 and mode 2 on identical still-water and plunge seeds, without
emission first, measuring water/air speeds, particle energy and distribution.
Tolerance and timestep changes are separate arms. Agreement is necessary but is
not sufficient: both implementations must also satisfy the physical controls.

The run ledger and retained control outputs live under
`out/runs/phase_contract_20261023/`. Existing 40-second simulation caches and
render presets are not inputs to this job and are not changed.

## 2026-10-23 验证记录与当前结论

Xue 允许继续核查后，先实现上述统一相场的显式实验模式，再执行静水、
撞水、求解容差和时间步对照。**相场构造与短单步对应关系已通过；撞水的
完整独立序列仍未通过本项目预设的相似性门槛，不能进入生产。** 本轮没有
打开额外喷雾来改变画面，也没有进入粗 3D、正式海浪、渲染或删除旧缓存。

这一轮刻意关闭 escape、droplet、emission、ST 与 adapt，先隔离基础两相
的相场、转移和推进。因此它不是完整 PF-FLIP 场景的验收，也不能证明
开启逃逸处理后仍会有相同的能量变化。phase_model=1/2 只供研究；默认
phase_model=0 保留，逃逸分类与喷雾机制没有因本次相场模式而自动获得论文对应性。

### 已通过的局部与静水检查

- Eq.7 的独立解析点、alpha 缩放、纯相邻接真值表通过。
- 纯水、纯空气、混合界面，两个模式的实际 fused P2G 共 6 组相场与 NumPy
  float32 输出逐位一致。三步真实 graph 的位置逐位一致，速度最大差
  `3.35e-8` 格/秒，独立动量端口最大差 `2.98e-6`。
- ST、粒子自适应和非法 alpha 的实验配置被实际 graph 拒绝；这些组合
  尚未建立新相场校准，不能直接将这个 uniform 模式用于生产自适应场景。
- 与修改前保存控制比较，165 个旧 P2G 输出数组逐位一致。额外执行旧喷雾
  配置的 41 步真实 graph，保存的第 0/40 步位置、速度、相位、半径、喷雾、
  逃逸标记共 12 个数组，与修改前的 20 秒实验对应输出逐位一致。
- 首轮静水：48×40，初始 7,680 粒子，5 秒。只有 Eq.7 时最终水 RMS 速度约
  `0.0665` 格/秒；加入纯相面后约 `0.00942` 格/秒。两边一致，水重心漂移
  小于 `0.011` 格，采样总能量峰值增长约 `0.172%`。初始边界粒子被推进器
  夹到 0.51 格的位置会改变初始势能，不能把这个小变化全部归给压力投影。
  静水运行早于下述动量端口修正，原始文件及这个事实一并保留。

### 对照程序中发现的独立动量端口差异

C++ 由未舍入的最终速度分别导出 float32 速度和 float32 动量。最初的
Python 对照却先舍入速度，再重算动量；这会产生第二次舍入。一个保存状态
有 3,053 个动量分量不同，最大绝对差为 `0.0009765625`。现在实验参考把动量
作为独立循环状态保存，在调用方舍入速度之前计算。它仍独立求解全部物理，
没有把 C++ 的速度或压力灌给 Python。

原始 `plunge_raw`、`plunge_pure` 的失败记录保留；修正后另写
`plunge_pure_momentum`。C++ 的最终位置、速度和两个动量数组在修正前后
逐位相同，证实该轮只修改了 Python 边界对应方式。这个修正没有让完整
撞水对照通过，不能把全部后期差异归咎于它。

早期分析把稀疏保存的最后一帧当作最终状态；后续分析明确使用另存的实际
final 状态。6 秒、每 0.1 秒保存的旧末帧约为 5.9025 秒，而非 6 秒。后续
门槛同时覆盖全部稀疏帧与实际最终状态，没有因此放宽门槛。

### 6 秒撞水与数值敏感性

64×80，初始 20,480 个粒子，水柱撞入水池。预先固定的独立实现相似性门槛为：
4 格分箱水分布总变差 <0.05、水重心差 <0.2 格、水高度第 99 百分位差 <1 格。
这是项目的诊断标准，不是论文提供的定理或任意湍流场景的普适标准。

| 配置 | 水分布差最大值 | 重心差最大值，格 | 高度 P99 差最大值，格 | 采样粒子总能量峰值／初始，Python / C++ |
| --- | ---: | ---: | ---: | ---: |
| 纯相面、修正动量端口、dt=0.0025、tol=1e-6 | 0.02214 | 0.05436 | 2.6187 | 1.0272 / 1.0302 |
| 只收紧 tol 至 1e-10 | 0.02031 | 0.05941 | 2.4644 | 1.0434 / 1.0338 |
| dt=0.00125、tol=1e-10、按 Eq.13 缩放混合量 | 0.02896 | 0.06880 | 2.9596 | 1.0475 / 1.0771 |

三个配置均因高度门槛失败。整体水分布与重心较接近，并不意味着飞出的
高处水粒子已一致。能量按固定粒子质量计算动能加重力势能，不含新的能量
来源；瞬态超调值得调查，但不是仅凭一个数就证明代码不稳定或已经 blow-up。

时间步减半时，依据论文 Eq.13 保持 `(1-alpha)/dt` 近似不变：水 alpha 从
0.97 改为 0.985，空气从 0.9 改为 0.95。否则固定 alpha 会约加倍人工黏性，
容易把额外阻尼误当成数值收敛。这里仅做一个层级的敏感性检查；混合量公式
本身是近似，且两套稀疏输出时间相差 0.00125 秒，不能据此估计收敛阶。
收紧容差和缩小步长都没有消除当前差异与能量超调；不能以此证明所有误差
均与压力无关，也不能把继续拧紧 tol 当作已找到的修复方案。

### 在相同撞击状态下检查能量发生了什么

取修正边界后 C++ 缓存中的约 2.8、3.0、3.2、3.4、3.6 秒状态，每个分别
用 phase_model=1/2 执行一次真实 fused graph 和一次独立 Python 单步。
这些稀疏状态没有保存独立动量端口，因此本审计明确从保存的 float32 速度
重建一次共同动量，交给两边；这是共同状态实验，不是原序列的精确续跑。

10 个单步全部满足真实压力残差 <`8.54e-11`；按压力系数计算的网格动能在
投影时均减少。两边最终粒子总能量相对差 <`1.34e-9`，速度最大差
`1.31e-5` 格/秒。因此这个共同状态结果不是只观察 Python 后推测 C++。

约 3.2 秒的同一输入，在纯相模式中，网格投影降低动能，但粒子 G2P 后的
动能增加约 2.09e6，随后位置推进的势能增加约 4.18e5，总机械能增加约
2.51e6（该输入总机械能的约 0.226%）。只用 raw Eq.7 的该单步总机械能反而
减少约 1.56e6。数值使用粒子 Vp=1 的一致单位；网格与粒子的能量归一化不同，
这里比较各自变化的正负，不能把两种绝对能量直接相减。

这把重点进一步移到**相场密度、转移和粒子推进组合的能量行为**，而不是
“压力矩阵已经解错”的单一解释。它尚未识别某一行根因，也没有证明纯相处理
本身应删除。不能靠关闭论文步骤或加大阻尼让图好看来代替对应关系核查。

### 论文对应性仍有明确缺口

重新核对 PF-FLIP §3.3–3.7 与实际调用路径：

- 本次 cell-centre raw splat 是公开说明的实现选择，作者的具体 cell tagging
  采样方式仍需确认；不能把这一选择自动等同于作者实现。
- 当前主流体推进是二阶中点法；论文 §3.7 描述第三阶 Runge–Kutta 和局部
  自适应时间步。当前 uniform/non-ST 路径没有这个三阶局部推进。
- 主求解路径没有把水相速度向空气侧外推数格后再用于水粒子推进；论文
  §3.7 明确提到它。`reference2d/validation.py` 中单相对照的 `_extrapolate`
  不是主 PF-FLIP 路径拥有该步骤的证据。
- 论文 §3.5 的体积漂移校正仍未启用／实现于这里的主路径；原参考文件亦
  明写 no divergence correction。需测量粒子覆盖／体积漂移，不能直接把
  当前能量现象归咎于缺少这一项。
- 当前逃逸机制使用对侧相粒子局部密度，尚未与论文 §3.4 的 phi 阈值分类
  完成新的统一契约；本轮故意关闭该机制，并未验证新相场与喷雾联动。

下一步应围绕共同状态和已暴露的能量峰值，逐项核对水相推进速度、外推、
分类与局部积分，并继续独立 NumPy／真实 C++ 验证。先完成这一层，再恢复
逃逸、额外喷雾、ST 与自适应的分项联动，之后才做 20 秒完整对照和粗 3D。
不能声称本轮已经修复原 20 秒视频后半段的混乱喷雾。

完整证据为 [phase_contract_20261023.json](../results/phase_contract_20261023.json)。
复现入口为 `gates/validate_phase_contract.py`、`experiments/phase_contract_2d.py`、
`gates/audit_phase_energy.py` 与 `gates/regress_phase_legacy.py`。
只读分析脚本 `experiments/analyze_phase_contract.py` 生成曲线、对照图和 6 秒
诊断视频；蓝色为水，灰色为空气，逃逸和额外喷雾均关闭。输出和所有失败
分支保留在本 worktree 的 `out/runs/phase_contract_20261023/`。

## 2026-10-23 后续：逐步定位分叉和增能来源

Xue 要求找出确定原因，不能用“流体混沌”代替检查。本轮将问题拆成三件事：
原始两套实现的物理定义差异、统一定义后的数值分叉、两边共有的能量行为。
以下结论属于实验 phase_model=2、uniform/non-ST、逃逸和额外发射关闭的
64×80 撞水案例，不自动推广到原始 20 秒完整喷雾或生产 3D。

### 相同投影速度能否消除完整序列的差异

新增仅在 `PFFLIP_FLAT_STAGE_TRACE` 环境变量设置时使用的诊断导出。
实际 C++ graph 逐步保存 P2G 后、压力投影后和相场的内部 double 面数组。
文件按步骤追加；调用者必须选择新文件，不能并发共用路径。它只作用于
实验相场的 flat/non-ST/non-adaptive 路径，默认运行不写文件，也不改变物理。

先重现此前 `plunge_pure_strict` 的完整 2,400 步：导出打开后最终 13 个状态
数组与修改前控制逐位相同。随后 Python 从相同初始粒子开始，仍自己做
P2G、G2P、位置推进和独立动量端口舍入，仅在每一步使用 C++ 的投影后速度。

**全部 2,400 步的位置和速度逐位相同；逐步 P2G 速度、相场及最终动量也
全部逐位相同。** 再共享 P2G 的冗余对照得到相同结果。因此在此范围内，
不必假设另一个隐藏的 P2G/G2P/推进差异来解释独立运行时的分叉；替换压力
投影这一阶段已经足以消除它。该实验故意取消压力求解的独立性，是因果
隔离，不是独立 solver 验收，不能拿它把之前的失败改成通过。

默认 phase_model=0 的旧喷雾短序列也重新检查：41 步中的 12 个保留位置、
速度、相位、半径、喷雾和逃逸数组与修改前逐位相同。

### 首个差异及一次明确的放大机制

再让 Python 独立求解压力，逐步对照同一 C++ trace；它在所有此前保留采样
时刻及实际最终状态，均逐位重现原独立 Python 轨迹。两边的不同迭代算法、
预条件器和浮点运算顺序，不会产生逐位相同的近似压力解，即使都满足容差。

| 时刻／步骤（从 0 计） | 测量到的事件 |
| --- | --- |
| 第 0 步 | 投影前速度及相场完全相同；独立投影后 double 网格速度最大差 `4.18645e-12` 格/秒。导出 float32 粒子速度最大差 `9.31323e-10`。 |
| 第 30 步结束，约 0.0775 秒 | 第一次出现位置差：一个坐标分量相差 `1.90735e-6` 格。此时相场仍完全相同。 |
| 第 1010 步开始，约 2.525 秒 | 第一次 cell 水/空气标签不同：格子 `[55,26]` 的 phi 分别为 `0.49942720465` 和 `0.50043018569`，分处 0.5 阈值两侧。 |
| 固定此刻同一输入速度的干预 | 只更换这一个 cell 标签，保留自己的 raw face mass 和其他所有输入，单步投影速度最大改变 `0.75021773` 格/秒。 |

cell phi 从各自保存的位置按同一独立中心 splat 重建；面 phi 则直接核对
C++ 实际 trace。第一处标签改变使一个面 phi 最大跳变约 `0.79442`。
若替换所有 C++ 面 phi，速度最大差为 `0.75060342`；只替换标签已复现这个
量级，二者速度场最大差 `0.00342642`。这将**分类阈值的离散跳变**与
普通连续的面质量变化分开了。不能将该干预当作修复，也不能据此删除论文
的纯相面处理：cell tagging 的具体采样及整套界面处理仍需共同验证。

这条分叉链现在有可复现证据：压力近似解的微小差异 → float32 循环状态和
粒子位置逐渐不同 → 跨过 cell 分类阈值 → 面密度／压力系数跳变 → 下一步
速度差显著放大。它不是“压力没有收敛”的证明，也不是所有增能现象的解释。

### 单独定位粒子压力动能来自哪些面

此项使用前述 `energy_replay/step_1280_phase_2` 的**共同输入**单步，来源为
另一个容差 1e-6 的已存轨迹，再以 1e-10 重算；它不是上面 strict 序列中
同一时刻的状态，不能把两个实验拼成一次连续观测。

将压力插值增量分为纯空气面、混合面、纯水面三部分。若总速度增量为
`dv=sum(dv_i)`，将动能功按 `m*(v_before+dv/2) dot dv_i` 分配；这样包含交叉项，
各部分严格加回总动能变化，数值重构最大误差约 `4.66e-10`。

6,008 个水粒子的压力更新动能净增长为 `2,662,463.525`（Vp=1 单位），其中：

| 插值来源 | 对水粒子的压力动能净贡献 |
| --- | ---: |
| 纯空气网格面 | +2,633,455.136（净总量的约 98.9%） |
| 混合面 | +36,329.836 |
| 纯水面 | −7,321.447 |

大部分增长集中在粒子采样 phi 为 0.5–0.9 的界面水粒子，且空间上集中在
两侧壁：距侧壁 3 格内的水粒子净增 `2,882,138.026`，其余水粒子净减
`219,674.501`。最大几个增长粒子在左侧约 x=1.2、y=33–35 格的位置。
因此这不是整个水体均匀被加热，也不是主要由 phi<0.3 的孤立水粒子造成。

当前 G2P 使用同一速度场的双线性插值。压力更新按面密度计算，纯空气面
密度低，可能有较大的速度变化；邻近水粒子的插值会读到这些面，并按水粒子
自身质量形成动能。网格的密度加权动能下降，不保证这种非伴随的粒子转移
也降低动能。上面的分解实测了这个路径；它**没有**证明所有空气到水的动量
传递都是错误的，也没有证明直接删掉空气贡献是物理正确的修复。

接下来修复验收应针对水相界面速度／推进、侧壁附近的支持和分类、以及
论文要求的液相速度外推做共同状态检查，明确哪些速度用于 G2P、哪些用于
推进。论文 §3.7 明确说外推用于推进，不能擅自声称它规定了相同的 G2P
外推。应重新检查粒子能量、界面覆盖与相似性，然后恢复逃逸和发射，并做
完整 20 秒对照。当前没有改动生产默认值，没有将 6 秒失败门槛放宽，亦未
声称原视频后半段的混乱喷雾已修复。

复现入口：`gates/trace_phase_causality.py`、`gates/trace_phase_sequence.py`、
`gates/trace_phase_independent.py`、`gates/intervene_phase_tag.py`、
`gates/analyze_phase_causality.py`。小型证据为
[phase_causality_20261023.json](../results/phase_causality_20261023.json)，原始数据
在 `out/runs/phase_causality_20261023/`。首轮只读分析已写完数字结果，但绘图
因 Matplotlib 参数位置报错；修正后另写 r2，完整侧壁范围另写 r3，旧文件保留。

### 2026-10-23 关于分类容差和导数的讨论

Xue 提出合理的工程问题：允许近似数值解，为何分类仍是硬阈值；能否使用
容差带或导数来判断更倾向水还是空气？这是一条待验证的算法方向，本轮
没有修改分类规则。随后只完成独立的[压力失败防护](pressure-replay-audit.md)，
原严格撞水 2,400 步最终物理状态逐位不变，因此不能将它说成分类放大已修好。

相场本来就是连续指示量，但不是概率，也不能直接当作精确体积分数。
PF-FLIP §3.3 明确按 0.5 分类 cell，然后把同相邻接面归为纯相；这是此次
离散跳变的位置，不代表所有面相场都从来只有 0/1。
来源：[PF-FLIP 论文](https://animation.rwth-aachen.de/media/papers/94/2025-Siggraph-Adaptive_Phase_Field_FLIP.pdf)。

工程比较器常用滞回：高于上阈值切成水，低于下阈值切成空气，中间保留
上一标签。例如 0.49/0.51 仅是讲解用数字，未经校准，不能直接用于本项目。
它减少阈值附近来回切换，但引入历史依赖、延迟及初始化/细化/restart 状态；
两条轨迹还可能在新的边界再次分叉。压力残差的容差不能直接转换成相场
容差，因为它们是不同量，且 phi 还受粒子采样、网格与历史误差影响。
参考：[TI 对比较器噪声与滞回的解释](https://www.ti.com/document-viewer/lit/html/SSZT122/GUID-E7C418E2-6753-4FD9-85EC-9D09D3464D30)，
这是工程类比，并非适用 PF-FLIP 的正确性证明。

连续过渡带也可消除该处跳跃，但会改变密度和压力耦合；应核查论文原有
界面厚度、alpha_phi、采样和边界支持，不能随意抹平密度以让视频更相似。
必须用静水、界面平移/撞击、近壁/离壁和网格/时间步对照检查偏差，不能只看
同一 demo 更顺滑。当前优先验证论文对应性，避免新增经验设计。

空间一阶导数提供界面方向与陡峭程度，时间导数提供变化趋势；纯水内部和
纯空气内部都可能梯度为零，所以导数不能单独替代相位判断。它可以帮助
定义不确定带、界面距离或一致的邻域重建。局部差分通常是一遍 O(N) 的
邻域计算，不等于额外压力求解；实际成本需测量，且求导也会放大采样噪声。
若试验此类方法，应保留原 0.5 实现作控制，不改变当前失败门槛，并明确
区分论文复现修补与另加的算法改进。
