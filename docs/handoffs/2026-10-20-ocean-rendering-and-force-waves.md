# 海浪场景与渲染交接 2026年10月20日

这次工作的目标是用已实现的 PF-FLIP / ST-FLIP Bifrost solver，尽可能接近论文 Breaking Ocean Waves demo 的观感。渲染已经有明显进展，薄层体积泡沫的几个版本值得保留；当前工作转向重新设计造浪方式。**下一步首先要交付一个真正能在 Maya 里看清并播放的 2D 粒子预演。正式 3D 场景、生产分辨率、启动器和启动时间都没有确定，也没有已排期的新正式模拟。**

本文件同时记录成功结果、失败尝试和未完成事项。不要把生成了 cache、通过了 headless bounds 检查，理解成 Maya 显示已经验收。上一轮正是在这里交付失误：用户打开场景只能看见零星几个点，无法观察水体。最新显示修复仍未完成 GUI 验收。

## 接手时先看这里

技术仓库在 `C:/Users/5d149/source/Anisim`，远端为 https://github.com/cuix-github/Anisim 。不要被会话的起始工作目录或旧项目名误导。先在这个仓库 `git pull --ff-only`，读取 `CLAUDE.md`、`README.md` 和 `docs/open-questions.md`，再看本文件与 `docs/force-wave-preview.md`。如果 checkout 有未提交修改，先检查内容，不能直接 reset 或覆盖。

交接前的技术分支是 `codex/force-wave-previews`，基于 `origin/main` 的 `d186ed8`。截至整理前的最近提交是 `7f73e65`，其后只有 `maya/build_wave_2d_scene.py` 的显示修复尚未提交；本次交接会保存它，明确作为未验收的候选。以最终 Git 状态和 PR 为准，不要假设新 agent 切到 main 就能看到这些代码。

本次交接保存于 **[Anisim draft PR 43](https://github.com/cuix-github/Anisim/pull/43)**，其中也包含此前的新施力实现和未验收的 Maya 显示候选。它保持 draft，不能当成已验收的生产功能。

此前仍开放的 PR：

| PR | 内容 | 注意事项 |
| --- | --- | --- |
| [Anisim PR 41](https://github.com/cuix-github/Anisim/pull/41) | 旧场景喷雾续跑工具、2D 诊断和 Maya builder | 旧 18–30 秒生产续跑已经取消；代码存在不等于任务仍有效 |
| [Anisim PR 42](https://github.com/cuix-github/Anisim/pull/42) | 新施力场景的时间与磁盘预算 | 是规划，不是可直接启动的正式配置；其中“尚未实现造浪器”的历史状态已被本文件更新 |

当前分支复用了 PR 41 中部分 helper 和 loop runner 改动。合并时要核对重叠，不能覆盖新接口。仓库采用 squash merge；每次 push 前确认对应 PR 仍 OPEN，避免往已合并的分支继续推。新工作走 PR，用户审阅合并。本次 handoff 不授权合并其他开放 PR。

## 用户的工作方式和当前授权

用户愿意讨论物理和方法，希望有依据地设计实验，先做控制和验证，再花时间运行。不要无目的地一个参数一个参数试。可复现的证据、真实可见的结果，比乐观地宣告完成更重要。

当前顺序是：新的施力造浪器 → 2D 粒子预演 → 粗 3D 粒子预演 → 用户确认浪形 → 正式模拟的分辨率和时间预算 → 再决定何时启动。2D 不能复现横向交叉、斜向碰撞和沿浪脊变化；它用于看纵向传播、破碎和冲岸，粗 3D 才检查方向性。

最近的具体要求是先不管粗 3D，把 2D 做到像之前能正常观看的 demo。用户允许为此花几分钟到十几分钟模拟，不要求极端低精度。当前新一轮 2D 已算完，不要因为显示不好就继续盲目重算。

此前一小时是实施第一批预演的预算，不是未来所有工作的统一授权。用户反复要求：昂贵的模拟或渲染启动前，要把范围、张数和预估时间说清楚并达成一致；预计超时，提前告知。修复和检查现有文件可以继续，但本次 handoff 请求只授权整理，不授权新运行。

旧数据保留。虽然用户曾说磁盘紧张时可以考虑删除旧 cache，之后明确改为“旧数据先保留”。不要用较早的条件授权开始清理。不要删除失败实验，它们用于追查矛盾。

用户重视先还原两篇论文的演示，再考虑 incremental improvements。自定义 hybrid 粒子分类器、解析液滴替代、新的物理模型等不能悄悄加入。现有 clustered foam 是得到认可的外观实验，但不是论文原方法的等价实现。

电脑控制曾两次因用户按 Escape 被中止。中止后没有继续操作；不能据此推断画面已验收。接手 agent 需要重新获得当前窗口状态，不能复用旧窗口句柄或坐标。

## 目标参考与仍未解决的观感

参考视频是 [Breaking Ocean Waves 所在视频](https://youtu.be/nt9BohngvoE?t=145)，主要看 2:26 起的场景。用户提供的截图标注 domain 1 kilometre、2 billion particles、3072×1536×1024。其目标不是逐项追平这个算量，而是海面、浪峰破碎、向后飘散的水滴和细碎白水给人的规模感。

论文入口：

- [PF-FLIP preprint，第七节渲染讨论](https://ge.in.tum.de/download/Adaptive_Phase_Field_FLIP_preprint.pdf)
- [ST-FLIP，第四节与附录 B](https://ge.in.tum.de/download/ST-FLIP.pdf)

用户的核心感受是：当前画面像一米鱼缸里的水，晶莹、金属或玻璃感过强，大片水体和孤立碎片显得厚重；参考更像广阔海域，细节层次和白水的连续/碎裂关系更丰富。不要将此简单解释为“粒子不够”，也不要因为本项目单元的实际米数较小，就承诺必然能捕获参考里的全部细节。域范围、强浪形成过程、可解析尺度、重建、白水表现、视角与画面尺度都有关。

渲染改善后，用户怀疑更根本的问题在场景：原 30 秒模拟里，浪产生得晚、强度不足，只有很短时间能观察卷浪和飞溅。于是提出换掉旧造浪设计。**渲染成果应保留；新场景是否有效仍需要观察。**

## 当前最值得保留的渲染配置

“最终认定”应准确理解为目前经过验证、用户认可方向的候选集合，而不是已经选定一套最终海洋预设。共同基础已明确：**corrected-baseline 水体与灯光 + 独立薄层体积泡沫，mist 关闭。**

完整机器可读注册表为 [`maya/presets/breaker_validated_20261019.json`](../../maya/presets/breaker_validated_20261019.json)，说明为 [`maya/presets/README.md`](../../maya/presets/README.md)。注册表保存实际场景设置、图像/场景/foam ASS 的 hash、相机、收据和图像检查。它是 frozen-scene registry，**不是新增的 `-Preset mix_u60_c40` 命令行选项**。

用户此前偏向 clustered，觉得 uniform 有整体连续感但局部太白、太不透明。25% uniform + 75% clustered 的图得到明确的 “This is WOW!”。之后授权 20/80 的六视图，再要求同镜头的 60/40 和远距离透视。没有明确宣布 60/40 永久取代其他比例。

agent 的独立判断是：纯 clustered 比纯 uniform 更能保留颗粒、孔洞和边缘变化；25/75 是很有价值的折中；20/80 的 24.9 秒近景是六视图中较好的展示候选。60/40 更饱满、连续，也更容易形成 matte 白片。应该保留这个比较，不把用户一次查看某张图解读成最终选择。

所有下列目录相对于 `C:/Users/5d149/source/Anisim/`。展示图优先打开 `beauty_view.png`，不要拿未经显示转换的线性输出进行外观判断。

| ID | uniform / clustered | 代表镜头 | 完整结果目录 |
| --- | --- | --- | --- |
| `uniform` | 100 / 0 | 24.9 s，18 m oblique | `out/runs/breaker_foam_views_20261019/oblique/A/` |
| `clustered` | 0 / 100 | 同一 18 m 镜头 | `out/runs/breaker_foam_views_20261019/oblique/B/` |
| `mix_u25_c75` | 25 / 75 | 同一 18 m 镜头，用户强烈认可 | `out/runs/breaker_foam_mix_20261019/C/` |
| `mix_u20_c80` | 20 / 80 | 24.9 s，22 m oblique | `out/runs/breaker_six_views_20261019/shots/break_detail/C/` |
| `mix_u60_c40` | 60 / 40 | 同一 22 m 镜头 | `out/runs/breaker_foam_mix60_20261019/C/` |
| `mix_u60_c40` 远距透视 | 60 / 40 | 24.9 s，1920×1080，距离目标 100 m | `out/runs/breaker_mix60_perspective_20261020/` |

前面三项是相同相机/几何的一组，后面两个 22 m 镜头是另一组。18 m 和 22 m 之间，镜头与水体重建范围都变了；不能把全部差异归因于泡沫比例。远距透视图与 22 m 60/40 的几何、材质、泡沫已核对相同，改变的是摄影机和输出分辨率。

每个可复现目录应一起保留 `candidate.mb`、`foam.ass`、`state.json`、`receipt.json`、`image_check.json` 和 EXR/PNG。仅复制 `.mb` 不保证外部 stand-in 还能加载。修改渲染实验时复用冻结网格，不要无意中重新计算水面让比较变量失控。

### 水材质与光照的实际值

这些值来自验证过的保存场景，不能用旧 launcher 的默认值替换。

| 项目 | 数值 |
| --- | --- |
| 环境 | Maya 2027，MtoA 5.6.1.1，Arnold 7.5.1.1，CPU |
| 常规预览 | 960×540，AA 4，OIDN |
| 工作/显示色彩 | ACEScg；ACES 1.0 SDR-video (sRGB) |
| Water base / metalness | 0 / 0 |
| Water specular / transmission | 1 / 1 |
| IOR / specular roughness | 1.33 / 0.03 |
| Transmission color | (0.3, 0.55, 0.6) |
| Transmission depth | 40 个保存场景坐标单位，对应 10 个模拟米 |
| Transmission scatter | (0, 0, 0) |
| Thin-walled | false |
| Dome color / intensity | (0.03, 0.035, 0.045) / 12 |
| Dome camera visibility | off，反射/照明仍存在 |
| Sun color / intensity | (1, 0.95, 0.86) / 6 |
| Sun rotation | (-40, 160, 0) degrees |
| Specular / transmission / total ray depth | 3 / 8 / 10 |
| Volume depth / volume samples | 3 / 2 |
| Mist | off |

Specular=1 不是“把水变成镜子”的证明；IOR 的 Fresnel 响应仍然参与。曾试 roughness 0.12，反射变柔，但没有因此解决大水域尺度感，也没有取代 0.03 baseline。

### 水面重建与泡沫的实际值

| 项目 | 值与意义 |
| --- | --- |
| 原始 hero simulation cell | 0.25 m |
| 粒子 sphere radius / surface voxel | 各 0.5 cell，即 0.125 m |
| Gaussian smoothing | **0 次**，不要再误读 generator 默认值 |
| Masked mean-curvature flow | 10 次 |
| Mesher smoothing | 0.1 |
| Mask sigma / theta / zeta / floor | 2 cells / 2 / 5 / 0 |
| Foam 表现 | 在曲率/高度选中的既有水面上采样，以 Arnold volume points 形成独立薄层 |
| Foam 球支持半径 | 0.0225 m，即 2.25 cm |
| Foam 中心向外偏移 | 0.0075 m，即 0.75 cm |
| 名义支持厚度 | 水面内 1.5 cm、外 3 cm；不是物理气泡直径 |
| Optical-depth 参数 | 0.8 |
| Surface quadrature area | 0.014² m² |
| Volume step | 0.005 m |
| Scatter / color / anisotropy | 0.98 / (1,1,1) / 0.2 |
| Cluster 空间长度 | 0.24 m 与 0.065 m |
| Cluster seeds | 774、991 |
| Density export | 64 个常量密度分组，零密度不输出 |

混合是在空间 extinction density 上做：

`density_C = uniform_weight * density_A + clustered_weight * density_B`

它不是两张 PNG 叠加，不是“60% 水面被泡沫覆盖”，也不是每个粒子有 60% 概率被选中。相同几何上 A/B/C 保持总 primitive extinction 积分预算一致，但透射是非线性的，画面亮度仍可不同。补少量 uniform 分量，会在 clustered 的空隙间加入稀薄连接；这就是此前讨论的“保留颗粒，同时有整体连续感”。

这层泡沫不需要重新跑主 solver。它也没有气泡生成/合并/寿命/输运历史，因而好看的单帧不能证明动画连续或能留下冲岸后的泡沫尾迹。它是静态外观近似，不是新的 whitewater simulation。

### 必须保留的单位约定

旧 hero 渲染采用 **4 个场景坐标单位 / 模拟米**。父 Maya ASS 与子 foam ASS 都必须声明 `meters_per_unit=0.01`。这两个概念不能混为一谈：后者是 Arnold procedural 导入的 DCC 元数据约定，前者是这套冻结场景实际沿用的模拟比例。

曾经子 ASS 声明 1，而父 ASS 声明 0.01，造成导入后泡沫放大/偏移 100 倍，整个 layer 没在水面上。原生 ASS 回读成功、transform 看起来 identity，均没有发现这个问题。后来用 **AiProceduralExpand 检查实际父场景展开结果**才定位。

成功版本先验证每个导入点、半径、变换链、距离水面 anchor 的偏移，再看渲染 volume AOV。不要只看导出节点存在就说“散射生效”。新 2D 粒子场景的米/厘米问题是另一条 Maya viewport 路径，不能直接套用旧 hero 的 4 units/m。

### 摄影机和图像尺度

18 m 组：rotation `(-10,145,0)`，target `(145.5,13.5,70)` 模拟米。22 m 组：rotation `(-15,140,0)`，target `(145.5,13.5,58)` 模拟米。这两组是正交近景；“18/22 m”指取景宽度，不是相机离水面的距离。

已经完成的远距透视 60/40 图：

| 项目 | 值 |
| --- | --- |
| Camera perspective | true，DOF off |
| Eye，模拟米 | (205.0885, 38.3819, -25.9942) |
| Target，模拟米 | (143, 12.5, 48) |
| 到 target 的直线距离 | 100 m |
| 高于静水面 | 26.3819 m |
| Target 平面取景宽度 | 56 m |
| Rotation / focal length | (-15,140,0) / 64.2857 mm |
| Sensor aperture | 36 × 20.25 mm |
| 输出 | 1920×1080，AA 4 |
| 实际水体重建范围 | 约 50×68 m 的现有局部重建，**不是整个 172×80 m 水域** |
| 计时 | Arnold call 370.00 s；process 380.25 s |

这张图已通过图像检查，AOV closure relative L1 约 8.83e-8。它用 594,229 vertices、1,188,078 triangles 的既有网格。不能用它的成本直接声称整片水域渲染成本相同。

六视图在 `out/runs/breaker_six_views_20261019/`：23.025 s 的 pre、24.9 s 的 break、27.025 s 的 shore，各有 wide/detail；`comparison_full.png` 保留原分辨率。wide 56 m，detail 22–24 m。第一张 65° 高视角把浪压平；之后两张 wide 改为 40°。总共六次渲染，已完成，不要重新执行整批脚本来“查看结果”。

## 渲染尝试的来龙去脉

按以下顺序阅读 [`docs/rendering/README.md`](../rendering/README.md) 链接的 study，可以避免重复踩坑。

| 阶段 | 结论与后续 |
| --- | --- |
| 原水体、matte 几何检查、roughness 比较 | 单纯改材质无法解释所有厚块与棱角。`breaker_inspect` 适合诊断，不是最终真实海水 preset |
| ACES 显示修正 | 同一 render 的显示转换明显改善颜色；用户明确觉得更舒服。检查 `beauty_view.png` |
| Deep Water + Physical Sky | 有依据的候选，但太蓝，用户更喜欢 corrected baseline；没有被选为替代方案 |
| Neutral water + finer mist | 颜色与已有 mist 的尺度实验，没有建立新的默认 baseline；孤立白点仍存在 |
| 高曲率界面 SSS | 第一版数值节点把 mask 变成常量 0，修正后散射才生效；外观仍像给大块水体涂白 |
| Thin foam volume A/B | 前两次 A 因 ASS 单位问题失效；`breaker_foam_layer_20261019_v6` 才是成功的单位对齐版本 |
| Uniform A / clustered B | 同样支持位置/半径、相同 extinction 预算，仅密度分布不同；A 连续，B 细碎，均保留 |
| 25/75、20/80、60/40 | 经检查的中间版本；25/75 获用户明确好评，20/80 做六视图，60/40 做近景和远距透视 |
| 角状碎片审计 | 一个被圈出的 fragment 在 raw/current mesh 都只有 6 vertices、8 triangles；轮廓问题已在重建几何中，不能靠灯光或法线圆化 |
| Surface correspondence | 单纯细分 mesher 没修复该 6/8 碎片；局部 volume resampling 增加面数但也明显改变形状，未替换 baseline |

重要的失败解释：早期体积层没贡献时曾怀疑 per-point density 不受支持，后来真实证据是单位错配。不能把这段历史写成“Arnold 不支持 per-point volume data”。现行 64 bins 是验证成功的工程方案，不是该推论的证明。

旧 reconstruction 的内存峰值也不意味着只能重跑模拟。一次 public-node bounds 审计发现，约 8.4 m 高粒子范围的 `volume_bounds` 给出了 108 m 高的范围，随后 dense mask field_to_volume 分配过大。只把 **mask field** 限制到粒子 bounds + 8 m padding，保持 0.5-cell 网格对齐，复现控制网格所有点/三角形不变。六视图成功使用的是这条 bounded 路径，不是之前试过的 tiled 拼接路径。

论文的 implicit density 0.5 isocontour、特征保护曲率流和本实现的 signed-distance 演化/meshing 并非完全等价。论文常见 30 次与目前 10 次迭代不能仅按数字判断平滑强度一致。此前 proposal 的 hybrid classifier / analytic droplet 路径仍被搁置。

### Whitewater 和 spray 要分清

Spray 是空中水滴；bubbles 是水中空气；foam 是水面泡沫。主两相 solver 能表示逃逸水粒子和被困空气，不等于当前渲染已经完整表现这三类。

PF-FLIP 第七节还描述额外的梯度触发喷雾与液滴空气阻力、重力、蒸发；不能把参考里所有飞溅都归功于主两相求解。项目已实现 droplet / extra emission 路径，原 hero run 没开启这两条新模型；新 force previews 已打开。表面泡沫目前是渲染近似，水下 bubbles 也没有完整可见的渲染实现。

“mist 关闭”是渲染可见性控制，与 solver 的 `emit_model` 开关不同。已有 cache 中逃逸粒子、兼容路径的 radius 标记、新生成 spray 不是同一计数。此前局部 24.9 秒 cache 中提到的 1,192 个 radius-tagged 粒子不能当成所有逃逸粒子的总数。更细 mist voxel 只能更细地表示已有粒子云，不会自动新增液滴，也不会自动缩小 cloud radius。

## 数据和工具在哪里

本机 gallery：`C:/Users/5d149/source/Anisim/out/gallery/index.html`。它是离线索引，不是远程网站。用户非常喜欢，保持原图与链接可用。

| 位置 | 用途 |
| --- | --- |
| `maya/presets/README.md` | 参数解释和已验证候选 |
| `maya/presets/breaker_validated_20261019.json` | 五个 frozen setups 的精确状态、hash、receipts |
| `docs/rendering/studies/` | 各轮经过、失败与结论 |
| `experiments/rendering/breaker/` | 分组后的渲染实验脚本 |
| `experiments/rendering/breaker/README.md` | 脚本导航 |
| `results/rendering/` | 可提交的小型数字证据 |
| `out/runs/` | 原 cache、渲染、场景、图像、失败实验，未删除 |
| `experiments/rendering/build_gallery.py` | 只扫描已有媒体生成 gallery，不运行 solver 或 Arnold |

旧的 `experiments/breaker_*.py` 入口作为兼容 shim 保留。整理时核对了原媒体大小/时间戳、保存的 Maya/ASS hash 和引用，没有迁移历史大文件。不要为了 handoff 或美化目录再移动 cache。

`breaker_foam_mix.py`、`breaker_foam_mix60.py`、`breaker_six_views.py` 等是当时具体实验的仪器，带固定输出路径、先 commit 和防重跑断言；不是随便再执行一遍的通用 preset。需要新 render 时读清 mode、输出目录、receipt 和已经消耗的授权张数。

## 为什么放弃旧 wave zone

旧 `wavesea4` 是 688×96×320、0.25 m cell，即约 172×24×80 m，dt=0.025 s，30 s 共 1,200 steps。五段运行的 receipts 共约 16.893 h；不是尚未算完的模拟。

旧 coast-engineering wave zone 在源区域内追随/改写速度。用户不喜欢这个控制方式，尤其此前不完整的区域处理造成过堆水与异常流动。现在要求源只施力，由 solver 的投影、输运和边界处理决定后续流动。**施力本身不保证不会有边界事故。**

用户同意约一千个不同振幅、周期、相位、方向的施力分量，借鉴 Tessendorf 的频谱组织思想，但不把 FFT 海面直接当流体结果。两三个主要 shoreward wave groups 携带主要能量，背景频段提供层次。1000 是分量数量，不是 1000 个随机猛烈冲击点。

原先次日 7 点从 18 秒继续旧场景到 30 秒的任务已经取消。`anisim-18-30` automation 已删除，准备好的旧 job 标记 cancelled，没有执行。不要因 PR 41 还开放就重新启动它。

## 新施力实现与验证

实现和入口：

- `nodes/src/PFFlipSolve.h` / `nodes/src/PFFlipSolve3.cpp`：通用力与阻尼端口。
- `graphs/generate_sim_graph.py`、`graphs/sim_3d.json`、`graphs/sim_3d_fused.json`、`graphs/loop3d.py`：Bifrost 编排与 absolute-step 读取。
- `reference2d/external_forces.py`：独立参考。
- `gates/validate_external_forces.py`、`validate_external_forces_3d.py`、`validate_force_resume.py`：验证。
- `experiments/force_wave_preview.py`：频谱配方、离线 acceleration bake、run 和粒子 cache。
- `experiments/force_wave_audit.py` / `force_wave_health.py`：粒子诊断与 checkpoint 检查。
- `maya/build_wave_2d_scene.py`：Maya cache 播放及旁置的实际 solver graph。

新增 `external_forces`、`liquid_acceleration_u/v/w` 和 `linear_drag_u/v/w`。加速度以 cells/s² 表示，阻尼以 1/s 表示。支持 scalar broadcast 或对应 MAC face 数组；尺寸错误、非有限值、负阻尼会报错。公式为：

`u_new = (u_old + dt * actual_liquid_fraction * acceleration) / (1 + dt * damping)`

它位于重力之后、压力投影之前，保留 FLIP 所需的旧 P2G 速度。加速度按实际液体相分数加权，阻尼作用于两相；没有投影后的目标速度覆盖。新场景 `zone_length=0`。

当前频谱源是 **Python 离线 bake → Bifrost 每步读取 → 自定义 solver**。它没有被实现成 Maya 里可实时编辑频谱的独立 compound。保存场景里带着真正的 `sim_3d_fused` graph，但其输出未接入播放路径，避免按 Play 就开始重新模拟。不要把它宣传成 live interactive solver。

已经完成的验证：

| 检查 | 结果 |
| --- | --- |
| 旧 DLL control → 新 DLL disabled path | 比较的 grid/particle 六类输出 bit-identical，原控制保留 |
| 独立 NumPy 2D force gate | grid max error 4.27e-7，position 3.82e-6，velocity 3.93e-7；NaN 输入拒绝 |
| 3D discrete-curl anchor | 16×12×8，散度/法向边界为零的 force，投影保持 dtF；最大速度误差 4.66e-10 |
| Resume clock | 16 steps 与 8+8 steps 的 13 个 state channels bit-identical；错误重置时钟的 negative control 有差异 |
| 单波校准 | H=1 m、T=6 s，gain 从 1 调到 0.7；35 s 运行，x=75/95 m 拟合高度 1.127/1.007 m |
| 单波波长 | 测得 52.65 m，线性理论 50.69 m，偏差 3.85% |

单波用 5 m 带内粒子高度的 97.5 percentile 拟合，不是重建水面。机制 gates、单波校准、混合海况数值稳定性、生产场景成功是不同层面的结论。不要把小 gate 当成所有 nonlinear breaking 已验证。

## 新场景配方与三个已完成的 cache

共同物理配方：180 m 长、32 m 高；静水深 12 m；beach toe x=70 m，坡度 1:5，初始 shoreline x=130 m。施力源中心 x=38 m，支持 x=14–62 m，4 s smooth ramp。seed=20，3 个 lead components 加 997 个背景分量。

Lead 的目标高度/周期为 4 m / 8 s、3 m / 9.5 s、1.5 m / 6.5 s；3D heading +15°、−18°、+5°。背景目标 Hs=1.5 m，周期 3–7.5 s，3D heading 限制 ±35°。gain=0.7。相位按线性估计聚焦 x=100 m、z=36 m、t=18 s。目标高度不是实际浪高保证。

2D 设所有 heading 为零并重新计算 lead phase，保留振幅；它不是 3D 中央切片的精确预报。源会向两边发波，x<20 m 的阻尼与 3D 侧面阻尼还没有建立可靠吸收边界。

全部新 previews 启用 ST、escape、adaptive air particles、droplet model、gradient emission：`emit_grad_threshold=10/s`、`emit_rate=0.05`、`emit_max_per_cell=8`。**有实体海滩时 pressure 仍走 uniform 路径，不能说 adaptive pressure 也打开了。**

| Run tag | 网格与步长 | 数据 | 真实 workflow 用时 |
| --- | --- | --- | --- |
| `force_sea_2d_20261020` | 180×32×1，dx=1 m，dt=1/24 s | 30 s，360 帧 / 12 fps，主体水约 9,600 | 42.7 s |
| `force_sea_3d_20261020` | 180×32×72，dx=1 m，dt=1/24 s | 30 s，360 帧 / 12 fps，初始水约 691,203，总粒子约 110 万 | 594.5 s，32 threads |
| **`force_sea_2d_half_20261020`** | **360×64×1，dx=0.5 m，dt=1/48 s** | **30 s，720 帧 / 24 fps，主体水约 38,400** | **103.3 s，8 threads** |

计时包含 bake/export，graph-call 时间还包含启动和 I/O；不是纯 solver kernel benchmark。最后一项只加密了 2D，粗 3D 没有一起重算，因此现在不能再把最新 2D 与旧粗 3D 称为相同离散精度的 matched pair。

采样帧中 x>68 m 区域同时存在的 escaped water / extra spray 峰值，分别为 95 / 1、2,747 / 7、306 / 3。这些不是累计发射量，也不等于画面所有飞溅。额外喷雾仍稀少。旧 diagnostic 曾把 emission rate 提到 1.0；新的匹配配方没有继承这种 20 倍加强，不能悄悄再加回去制造“成功”。

所有 step 的压力 residual 检查通过 1e-6 限制，finite 检查通过。每 5 s checkpoint 的液体位移上界没有超出约 0.72 cell/step，但没有由此证明所有中间步的 CFL 或收敛。最新 2D 在 20 s、x<10 m 的最高液体粒子达到 **23.73 m**，初始静水仅 12 m，后墙反射/堆积问题仍然明显。最高粒子位置不是 crest height 的可靠定义。

## 目前最紧急的 Maya 显示问题

最新应该接着检查的文件是：

```text
C:/Users/5d149/source/Anisim/out/runs/force_sea_2d_half_20261020/open_in_maya.bat
C:/Users/5d149/source/Anisim/out/runs/force_sea_2d_half_20261020/scenes/wave_2d_points.ma
```

启动器已指向 `wave_2d_points.ma`，不是同目录中较早、仍显示失败的 `wave_2d.ma`。时间轴 1–720、24 fps，初始选帧 470，约 19.5625 s。前一版 1 m preview 的启动器/场景不是最新交付。

最新 builder 保存的结构：

- cache/full 的位置是模拟米；额外 `cache/maya_cm/` 是乘 100 的显示副本，DAG parent identity。原 cache 没有被改写。
- 四个 layer：water、droplets、spray、air；air 默认隐藏。
- `particle_play` 直接输出 points，使用 Maya wireframe；尝试设置每层 RGB override。
- scene 有 `wave_viewport_setup` Python scriptNode，仅在加载后设相机、关闭 grid、设置 wireframe；不运行 solver。
- actual solver graph 在 `wave_solver_OPEN_IN_BIFROST` 下，输出未连接；timeline 只读 cache。
- headless 六帧 bounds 检查通过；**颜色、粒子屏幕尺寸、完整 GUI 播放速率尚未验收。**

诊断经过及其证据边界：

1. 原来约 9,600 水粒子的 cache 不是只有几个点；headless 可以读到完整位置并得到正确世界 bounds。
2. 加 `point_scope` 后的 shaded display 几乎全黑。`point_size`、shape、颜色、单位转换是调查方向，没有锁定一个已验证的根因。
3. Maya VNN float3 默认值需要 `"{r,g,b}"` 格式；逗号字符串没有按期望写入。现在 builder 对颜色做 query 回读。这个修正**没有单独解决可见性**。
4. `c.rename(parent, ...)` 可能同时改 shape 名；不能硬写 `waterParticlesShape`。实际节点应从 `c.listRelatives(parent, shapes=True)` 获取。
5. 默认摄影机 far clip 太近已修正至 1000 m。这个修正也不等于全问题解决。
6. GUI 中观察到 direct-output + wireframe 能显示一整片白色水点；切换 smoothShaded 又显著丢失可见点。此时的控制是 1 m preview，不是最后的 0.5 m 版。
7. 将位置显式转厘米、取消 parent scale、改 scope sphere 的版本，headless 通过但 GUI 仍失败，不能称为根因修复。
8. 最后改回 direct-output + wireframe、厘米显示副本、layer override，生成了新的 `wave_2d_points.ma`。在最后打开/播放验收前，电脑控制被 Escape 中止。**没有实际看到最终这个场景正常播放。**

不能再重复“保存成功，所以可播放”的交付。下一位先读取启动器和 scene 内容，确认当前打开的是哪一份；在 GUI 看清主体水、完整浪形和时间轴变化，再向用户交付。需要保留一份独立于 shading 的粒子分布检查，但它不能替代 Maya 实测。

`receipts/maya_scene.json` 当前 status=passed 指 headless evaluation/bounds，不含 GUI acceptance。`results/3d/force_wave_preview/maya_2d.json` 和 maya_3d 是早期收据，不是最新半米版的可视验收。对这些“passed”的含义必须保持明确。

粗 3D 的旧 `out/runs/force_sea_3d_20261020/scenes/wave_3d.ma` 同样没有可用显示验收，且用户明确暂时先不管它。

## 正式模拟的硬件预算与尚未准备的东西

硬件是 Ryzen 9 9950X3D，16 cores / 32 logical threads，约 62 GiB usable RAM，RTX 5090。现有 custom CPU solver 的主要预算看 CPU/RAM/I/O；不要把显卡余量换算成未经实现的 solver 加速。

用户提出的正式目标是约 15 小时内得到 **30 s** 模拟和可查看 cache；更早的 60 s 想法被这个阶段的 30 s 目标取代。原 hero run 已使用 32 threads，不能再给它乘一次所谓“启用全线程加速”。

PR 42 的线性容量筛选：180×32×72 m，dt=1/48 s，0.35 m cell 约 9.76 M cells，计入 50% 假设开销和 1 h cache 工作约 15.04 h；1/3 m cell 约 11.20 M cells，在同样假设下约 17.11 h。**这不是新场景实测，不是最后承诺。** 0.35 m 比旧 0.25 m 粗，是运行时与细节的取舍。

正式决定需要在接受浪形之后，串行测 startup、传播、破碎与 spray 活跃阶段，包括 I/O 和粒子增长，不能仅用安静前 5 秒外推。30 s、dt=1/48 共 1,440 steps；留 1 h 给 cache 时，14 h 算力预算约要求平均 35 s/step。不要通过放松 pressure tolerance 或临时关功能偷偷满足时限。

本次 handoff 只读检查 C 盘约剩 **1.226 TB，约 1,142 GiB**。之前较完整盘点中 Anisim out 约 635.86 GiB，`wavesea4/dumps` 439.31 GiB，而 `wavesea4/cache` 仅 9.47 GiB。原 dumps 是证据与再导出来源，不是随手可删的缩略缓存。

规划曾提出 650 GiB 新 run cap、保留 250 GiB free reserve；这些是尚待落实的预算，不是已经实现的所有保护。输出 cadence、粒子增长、checkpoint 和显示副本都要计入。没有已准备好的最终 production launcher，也没有新的 watchdog/故障自动修复任务在替它工作。旧 7 点 automation 已取消。

## 接手后的具体顺序

1. 保留当前代码、cache 和所有失败文件，确认 branch/PR 状态；不要先清目录或切换旧模型。
2. 只处理最新半米 2D 的 Maya 显示，实测完整主体水和起浪到冲岸的 playback。验证 layer、单位、相机、frame mapping，区分 UI 问题和数据问题。
3. 与用户看浪形：何时破、破多久、是否足够有力、后墙是否污染主观察区、逃逸粒子是否足够可见。用实际数据，不用“1000 分量”本身证明复杂性。
4. 如果要改物理配方，先解释改什么、预期和预算；分开测源强、吸收边界与分辨率，保留控制。不要改写速度重新伪装成施力。
5. 用户愿意后再恢复粗 3D，检查斜向碰撞、沿脊破碎和边界。先验收可见性，再谈大模拟。
6. 最后做正式 30 s 的速度/内存/磁盘抽样，确定分辨率与启动时间，再准备 launcher、checkpoint 和必要的低频故障处理。

本交接期间没有启动模拟、重建、Maya 或 Arnold，没有删除任何数据。渲染候选与新场景是两条不同进度线：前者已经有值得保留的结果，后者的可视预演仍是眼下的阻塞。
