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
