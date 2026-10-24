# 重新核查喷雾分辨率与实现边界

2026-10-22，针对已完成的 40 秒海浪缓存。

**PF-FLIP 第七节要求更高分辨率的喷雾模拟和渲染，不能只理解为把已有粒子栅格化到更细的体积。** 但原文没有要求再做一次同等高分辨率的压力求解，也没有完整给出细喷雾的出生采样、粒子权重和光学重建规则。我们的 3D 自适应并非完全没有做：本次运行使用了空气粒子合并与拆分，已有的 3D 自适应压力路径则没有启用，也不能直接支持这个实体海滩。现有代码具备液滴动力学和额外喷雾发射，但这不等于已经验证了论文整套多分辨率喷雾流程。已有缓存可以支持一个近似的离线喷雾方案；它缺少逐步投影速度场，不能承诺精确重放。应把这些不同层次分别验收，再决定需要改实现、改配置，还是补充缓存，而不是继续把外观差距全部交给材质参数。([PF-FLIP，第七节](https://ge.in.tum.de/download/Adaptive_Phase_Field_FLIP_preprint.pdf))

## 2–4 倍指喷雾的空间尺度，不是第二套压力求解

第七节把逃逸液体转为拉格朗日液滴，计算空气阻力、重力和蒸发，并在界面速度梯度大的地方增加喷雾粒子。其细密度场采用 MSBG 稀疏存储，粒子在相应尺度上模拟，随后形成渲染体密度。海浪表格中的压力尺寸为 3072×1536×1024，喷雾场为 12288×6144×4096：**4 倍是每轴倍率，完整包围盒的潜在体素数变为 64 倍，不是粒子数必须变成 64 倍。** 重建只在出图时执行，也不意味着液滴运动只在出图时积分。([PF-FLIP，第七节 p.15、表 7 p.20](https://ge.in.tum.de/download/Adaptive_Phase_Field_FLIP_preprint.pdf))

对本次基础格长 0.35 m，2–4 倍对应 0.175–0.0875 m 的细场格长。这是尺度换算，不是已经通过验证的推荐参数。拉格朗日粒子可以在比压力格更细的尺度上采样位置、粒径和轨迹，并从基础速度场取得空气运动；**原文没有说明另建细压力或细速度求解器**。因此，“我们没有细速度场”不能直接当作未实现论文的证据。真正未能建立对应关系的是细粒子的出生分布、权重、积分约束与密度重建。插值和更细体素都不能创造基础流场里没有解析出的涡旋。

主流体空间自适应、液滴局部时间步、喷雾细场和 ST-FLIP 时间采样是不同机制。ST-FLIP 附录 B 还另有两倍分辨率的水面重建，它改善的是出图几何，不是喷雾轨迹。应避免把这些“细化”合并成一个开关。([PF-FLIP，第五至七节](https://ge.in.tum.de/download/Adaptive_Phase_Field_FLIP_preprint.pdf); [ST-FLIP，p.9、附录 B](https://ge.in.tum.de/download/ST-FLIP.pdf))

## 3D 粒子自适应已在运行，自适应压力尚未用于这片海

代码审计对象是主工作区 `1c659316acbcebe79c282ea08b41dfb329a805e9`。本文所在文档分支从 `21b28c78e64d44837e7474439f810d80be02c3b4` 建立；下列已链接代码、图定义及历史记录在两个提交间无差异。生产输入与缓存只在本机 `out/runs/force_rear_3d_035_20261020/` 中，不属于公开仓库。已读保存配置、代码、既有 gate 结果和缓存头；没有重新运行 gate，也没有把当前源码等同于每一历史 chunk 的已验证二进制哈希。

| 层次 | 实现与本次运行证据 | 对当前瓶颈意味着什么 |
| --- | --- | --- |
| 3D 粒子自适应 | `adapt=1`；远离水面的空气粒子八合一，接近界面再拆为八个。末帧有 18,646,826 个 scale 1、5,386,927 个 scale 2 粒子 | 确实在节省深空气区粒子；没有把水或喷雾细分到基础尺度以下 |
| 3D 自适应压力 | 已有 octree、naive/cleanup 路径及 3D gate；本次明确设置 `pressure_solver=0` | 本次使用均匀压力，不是请求自适应后被静默忽略 |
| 实体碰撞与自适应压力 | 当前自适应路径遇到切割/封闭网格面的实体会回退 uniform | 不能直接切换开关就让当前海滩使用自适应压力 |
| MSBG 整体架构 | 当前仍有细均匀场、稠密 owner map 和每步 hash/tree 构建 | 不是论文 MSBG 的完整内存与性能架构 |

八合一和拆分来自实际 3D 分支，末帧计数来自 `dumps/c7/scl.0236.npy`，运行开关来自 `inputs/plan.json`。自适应压力已有均匀网格退化、线性压力、3D 旋流、z 不变退化及静水 gate；它们证明受限条件下的实现，不证明实体海滩支持或整套 MSBG 已完成。([粒子 split/merge](https://github.com/cuix-github/Anisim/blob/21b28c78e64d44837e7474439f810d80be02c3b4/nodes/src/PFFlipSolve3.cpp#L1519); [压力选择与实体限制](https://github.com/cuix-github/Anisim/blob/21b28c78e64d44837e7474439f810d80be02c3b4/nodes/src/PFFlipSolve3.cpp#L847); [3D gate](https://github.com/cuix-github/Anisim/blob/21b28c78e64d44837e7474439f810d80be02c3b4/gates/validate_adaptive_3d.py); [稠密 owner map](https://github.com/cuix-github/Anisim/blob/21b28c78e64d44837e7474439f810d80be02c3b4/nodes/src/PFFlipAdaptive.cpp#L124))

这对硬件预算有实际影响：论文通过整个管线的空间自适应节约远场成本，我们目前只取得其中一部分。不过在不改变最细格长时，启用自适应压力主要改变资源分配，不会自动增加浪尖喷雾细节。历史无障碍测试中自适应路径还比 uniform 慢，不能未经同场测量就承诺它会加速这片海。([历史成本记录及限制](https://github.com/cuix-github/Anisim/blob/21b28c78e64d44837e7474439f810d80be02c3b4/docs/history.md#L324); [成本数据](https://github.com/cuix-github/Anisim/blob/21b28c78e64d44837e7474439f810d80be02c3b4/results/3d/adaptive/cost_indicative.json))

## 已有喷雾动力学，但细采样和光学权重仍需建立对应关系

当前代码已实现额外喷雾、粒径相关阻力、重力、蒸发及局部子步。额外喷雾是无质量、单向耦合的液体标记粒子，不参加主流场 splat，回到水里或蒸发后删除。生产配置启用了 `escape`、`droplet_model` 和 `emit_model`；已有 gate 证明发射开关不改变非喷雾粒子的结果，这属于机制验证，不能替代 reference 外观验证。([标记和删除](https://github.com/cuix-github/Anisim/blob/21b28c78e64d44837e7474439f810d80be02c3b4/nodes/src/PFFlipSolve3.cpp#L1108); [液滴积分](https://github.com/cuix-github/Anisim/blob/21b28c78e64d44837e7474439f810d80be02c3b4/nodes/src/PFFlipSolve3.cpp#L1286); [喷雾 gate 记录](https://github.com/cuix-github/Anisim/blob/21b28c78e64d44837e7474439f810d80be02c3b4/results/3d/spray/receipt.json))

当前发射规则是本项目的明确选择：合格界面水粒子按 `min(1, rate × (G − G0) × dt)` 的概率出生，每来源粒子每步最多一个，并设每基础格上限。新粒子复制来源的位置与速度；该循环没有单独的两倍或四倍出生位置采样。局部喷雾步长按基础格 CFL 限制，并采样主步投影后的 MAC 速度。**增加 rate 会增加样本，但不能保证增加不同的轨迹结构；这仍是待验证的瓶颈假设。**([梯度与出生循环](https://github.com/cuix-github/Anisim/blob/21b28c78e64d44837e7474439f810d80be02c3b4/nodes/src/PFFlipSolve3.cpp#L1422))

Ihmsen 的被引方法包含对出生位置和速度的采样，并允许从预先计算的流体进行单向后处理；PF 的短段落没有说照搬其全部规则。该差别值得用于检查我们的出生分布，但不能据此直接添加一套论文未要求的模型。当前 `SURF_MIST_DETAIL`、云团半径、密度倍率和 VDB 实验属于已有粒子的表示与光学控制，不生成新轨迹。物理液滴半径、一个数值粒子代表的水量和栅格化核半径必须分别记录，才能判断“白点”“胖云团”“细而连贯的喷溅”到底在哪一步形成。([Ihmsen，§3](https://cg.informatik.uni-freiburg.de/publications/2012_CGI_sprayFoamBubbles.pdf); [当前 mist 构建](https://github.com/cuix-github/Anisim/blob/21b28c78e64d44837e7474439f810d80be02c3b4/maya/build_surface_scene.py#L130))

`esc_phi` **有连接、有作用，0 也是有效阈值**。正式计划未覆盖它，而当前图默认 0；escape gate 使用 0.3。这是配置与调用链推断，尚不是所有历史 chunk 的运行端口记录。代码以另一相的独立粒子计数和 `rho0_face` 校准后判断逃逸；气液计数没有在局部除以两者之和，因此不能保证互补。论文液体判据默认得到 0.3 cutoff；只有在我们的空气计数确实等于 `1−φ` 时，`esc_phi=0.3` 才能直接对应它，不能把论文的 `δφ=0.2` 原样填入同名直觉参数。阈值与相场对应关系是本项目应验证的事情，不应请作者调试我们的接线。([阈值使用](https://github.com/cuix-github/Anisim/blob/21b28c78e64d44837e7474439f810d80be02c3b4/nodes/src/PFFlipSolve3.cpp#L1099); [图默认值](https://github.com/cuix-github/Anisim/blob/21b28c78e64d44837e7474439f810d80be02c3b4/graphs/sim_3d_fused.json#L183); [gate 覆盖](https://github.com/cuix-github/Anisim/blob/21b28c78e64d44837e7474439f810d80be02c3b4/gates/validate_escape.py#L32); [PF-FLIP，§3.4](https://ge.in.tum.de/download/Adaptive_Phase_Field_FLIP_preprint.pdf))

## 现有缓存能支持近似离线喷雾，两个缓存可在同一场景渲染

8 个五秒 chunk 各有 60 份常规样本，合计 480 份，保存位置、速度、相别、尺度、逃逸、额外喷雾和半径。它包含空气和水，不只有水面粒子；每个 chunk 末还有 13 类完整重启状态数组。但常规缓存频率是 **12 Hz，主步是 48 Hz**，没有逐步保存投影后的 `u/v/w`、压力或 Eulerian 相场。出图位置经过 ST 时间同步，不是完整抖动积分状态；常规样本缺少 tau/weights 和稳定粒子 ID，不能把相邻数组同一行当作同一个粒子。([缓存与重启契约](https://github.com/cuix-github/Anisim/blob/21b28c78e64d44837e7474439f810d80be02c3b4/graphs/loop3d.py#L35); [ST 出图同步](https://github.com/cuix-github/Anisim/blob/21b28c78e64d44837e7474439f810d80be02c3b4/nodes/src/PFFlipSolve3.cpp#L1737); 本机完成记录 `receipts/cache_check.json`)

工程上可以在每个已存时刻从两相粒子重建占据场、空气速度及发射候选，然后在时间上插值，单独推进新喷雾。这保留主水体轨迹，属于需要实现和验证的**近似单向离线模拟**。更小的喷雾子步不会恢复 12 Hz 样本之间缺失的流体事件。更接近原求解的路径是从已有五秒重启点回放一个短区间，逐主步导出投影速度、相场和必要的边界信息；这仍会调用 solver，必须另行约定，不能当作免费后处理。论文没有给出离线喷雾 pass 的完整等价性保证，尤其需要确认液滴回落并入主流体时的处理。

渲染组织本身并不复杂：同一个物理时刻，**主缓存生成水面，喷雾缓存生成稀疏密度体积，二者一起进入同一个 Arnold 场景**，共享镜头和光照，参与遮挡、透射和散射。不需要先分别渲染两张 beauty 再简单相加。联合时必须统一坐标、单位、时间和 shutter；本次原始位置/速度是格单位，乘 0.35 转为米/米每秒，半径已经是米。还要明确粒子归属：若新 pass 替换旧额外喷雾，就不能同时再画旧 `spr=1` 集合；原始逃逸液滴、额外发射粒子和表面泡沫不能混作同一计数。这是基于现有数据接口的渲染工程方案，不是声称论文使用 Arnold。

## 文献库把下载、阅读和实现证据分开记录

新[文献索引](../../references/README.md)按核心方法、喷雾与液滴、水面重建归组，机器可读的[目录](../../references/catalog.json)记录来源、版本、SHA-256、页数及实际阅读范围。当前有 **5 份完整论文 PDF、1 份明确标注的书籍预览、2 项未取得 PDF 的来源记录**。PDF 只作为本地阅读副本，不进入公开 Git；原有模拟和渲染结果保持原位。

| 资料 | 对应问题 | 本次阅读范围与取得情况 |
| --- | --- | --- |
| PF-FLIP / ST-FLIP | 两相方法、自适应、时间采样、喷雾及重建 | 完整作者 PDF；定向重读相关章节，未宣称逐页精读 |
| Gissler 2017 | 相对空气速度和阻力 | 完整作者 PDF；摘要、§3–4.2，重点 pp.3–5 |
| Song 2005 | 拉格朗日液滴先例 | 作者上传网页 §5.1–5.4；PDF 下载被拒，未绕过 |
| Ihmsen 2012 | 额外粒子出生及单向后处理 | 完整作者 PDF；摘要、引言、§3.1–3.2 |
| Babinsky–Sojka 2002 | 液滴粒径分布 | 出版方摘要和引言；未找到合法公开全文 PDF |
| Sirignano 2010 | 蒸发和半径平方规律 | 授权 49 页 PDF 预览；只读相关 printed pp.21–24，非全书 |
| Bhatacharya 2011 | 粒子水面的 level-set 重建 | 正确的 SCA 2011 作者 PDF；摘要、§3.1 和 §3.2 开头，未以 2015 扩展版替代 |

代码注释与历史记录证明此前参考过第七节及 Ihmsen，但没有找到逐篇精读的完整证据，不能补写成“这些文章以前都读完并完整实现了”。被引用论文包含的所有功能，也不自动变成 PF 的实现要求。例如，Song 的压力梯度力与反馈，或 Gissler 的变形/遮挡模型，必须先确认 PF 是否采用。下一轮应围绕出生采样、阻力单位、蒸发生命、粒子权重及密度核建立逐项对应，而不是只依据论文标题添功能。([Gissler](https://cg.informatik.uni-freiburg.de/publications/2017_CAG_generalizedDragForce_v2.pdf); [Song 作者文本](https://www.researchgate.net/publication/220183677_Stable_but_nondissipative_water); [本项目已有喷雾记录](https://github.com/cuix-github/Anisim/blob/21b28c78e64d44837e7474439f810d80be02c3b4/docs/open-questions.md#L116))

## 下一步先确定可验证的对应关系

给作者的关键问题应集中在论文没有公开说明的四件事：细喷雾“模拟分辨率”实际改变哪些采样；一个粒子代表一个液滴还是带权水量；物理半径怎样映射成密度核和消光；离线重放需要什么场与时间频率。我们的图默认值、筛选、单位和实体支持则由本项目自行核查。

后续最有价值的验收对象是一个短的浪尖或碰撞区间：同时记录发射位置与数量、粒径及寿命、轨迹、密度积分和透射率，用它区分出生不足、运动缺失和光学表现。这样才能决定是否只需重建、是否需要一个新喷雾 pass，或是否要回放主模拟的一小段。本次只整理和审计文献、代码与已有缓存，没有运行模拟、喷雾重放、密度构建或渲染，也没有保证现有主流场一定能复现 reference 的全部结构。
