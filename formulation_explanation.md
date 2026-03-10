# Section 2. Formulation 公式详细解释

> 本文档对 aba0c8.md 第 2 节（Formulation）中的所有公式进行逐一解释，给出必要的推导过程，并阐明公式之间的内在联系。

---

## 2.1 运动学（Kinematics）

### 公式 (1)：线元映射

\[ d\mathbf{x} = \mathbf{F}(\mathbf{X},t)\,d\mathbf{X} \]

**含义：** 参考构型 \(\Omega_0\) 中的无穷小线元 \(d\mathbf{X}\) 经过变形后映射为当前构型 \(\Omega\) 中的线元 \(d\mathbf{x}\)，线性映射的媒介正是**变形梯度张量** \(\mathbf{F}\)。位移场 \(\mathbf{u} = \mathbf{x} - \mathbf{X}\)，即当前位置与参考位置之差，刻画了质点的平动；但要描述局部变形（拉伸、剪切、转动），需要对 \(\mathbf{x}\) 关于 \(\mathbf{X}\) 求梯度。

---

### 公式 (2)：变形梯度

\[ \mathbf{F}(\mathbf{X},t) = \frac{\partial \mathbf{x}}{\partial \mathbf{X}} \]

**推导：** 由 \(\mathbf{x} = \mathbf{X} + \mathbf{u}(\mathbf{X},t)\)，对参考坐标 \(\mathbf{X}\) 求偏导：

\[ \mathbf{F} = \frac{\partial \mathbf{x}}{\partial \mathbf{X}} = \mathbf{I} + \frac{\partial \mathbf{u}}{\partial \mathbf{X}} \]

其中 \(\mathbf{I}\) 是二阶单位张量，\(\frac{\partial \mathbf{u}}{\partial \mathbf{X}}\) 是位移梯度。公式 (1) 即由此定义直接得到：\(d\mathbf{x} = \frac{\partial \mathbf{x}}{\partial \mathbf{X}} d\mathbf{X} = \mathbf{F}\,d\mathbf{X}\)。

**Jacobian 行列式：** \(J = \det\mathbf{F} > 0\) 表示体积比（当前体积微元与参考体积微元之比），\(J = 1\) 对应不可压缩情形。

---

### 公式 (3)：右 Cauchy-Green 变形张量

\[ \mathbf{C} = \mathbf{F}^T\mathbf{F} \]

**来源：** 变形后线元长度平方为：

\[ |d\mathbf{x}|^2 = d\mathbf{x}\cdot d\mathbf{x} = (\mathbf{F}\,d\mathbf{X})\cdot(\mathbf{F}\,d\mathbf{X}) = d\mathbf{X}\cdot\mathbf{F}^T\mathbf{F}\,d\mathbf{X} = d\mathbf{X}\cdot\mathbf{C}\,d\mathbf{X} \]

因此 \(\mathbf{C} = \mathbf{F}^T\mathbf{F}\) 直接描述变形后线元长度的改变，是纯粹的应变度量（不含刚体转动）。\(\mathbf{C}\) 为对称正定张量。

---

### 公式 (4)：\(\mathbf{C}\) 的主不变量

\[ I_1 = \mathrm{tr}(\mathbf{C}),\quad I_2 = \frac{1}{2}\left[(\mathrm{tr}(\mathbf{C}))^2 - \mathrm{tr}(\mathbf{C}^2)\right],\quad I_3 = \det(\mathbf{C}) \]

**说明：** 三个主不变量在坐标系旋转下保持不变（客观性），因此本构模型常以 \(I_1, I_2, I_3\) 为自变量写出。三者的物理意义：
- \(I_1 = \lambda_1^2+\lambda_2^2+\lambda_3^2\)（主拉伸比的平方和，刻画总体拉伸）；
- \(I_2 = \lambda_1^2\lambda_2^2+\lambda_2^2\lambda_3^2+\lambda_3^2\lambda_1^2\)（刻画面积变化）；
- \(I_3 = (\lambda_1\lambda_2\lambda_3)^2 = J^2\)（体积变化的平方）。

其中 \(\lambda_1,\lambda_2,\lambda_3\) 是 \(\mathbf{C}\) 的三个特征值（主伸长率的平方）。

---

### 公式 (5)：链拉伸（chain stretch）

\[ \lambda_{ch} = \sqrt{\frac{I_1}{3}} \]

**来源（8-链模型）：** Arruda-Boyce 8-链模型将弹性体聚合物网络理想化为一个正立方体，8 条链从中心伸向 8 个角点。当立方体受到任意变形时，利用对称性可以证明每条链的端到端长度 \(r\) 满足：

\[ r^2 = \frac{r_0^2}{3}(I_1) \]

因此平均链拉伸为：

\[ \lambda_{ch} = \frac{r}{r_0} = \sqrt{\frac{I_1}{3}} \]

\(\lambda_{ch}\) 是连接宏观变形度量 \(\mathbf{C}\)（通过 \(I_1\)）与分子尺度链行为的桥梁，为后续定义损伤驱动力提供了物理基础。

---

## 2.2 热力学考量（Thermodynamic Considerations）

### 公式 (6)：内部机械功率

\[ P_{int} = \int_{\Omega_0}\left[\mathbf{P}:\nabla \mathbf{u} + f_p\dot{p} + f_{\tilde{\lambda}}\dot{\tilde{\lambda}} + \boldsymbol{\xi}_{\tilde{\lambda}}\cdot\nabla\dot{\tilde{\lambda}} + f_d\dot{d}\right]dV \]

**说明：** 采用虚功原理框架，内部功率由五项之和构成：
- \(\mathbf{P}:\nabla\mathbf{u}\)：第一 Piola-Kirchhoff 应力张量 \(\mathbf{P}\) 与速度梯度的双点积，对应经典弹性能输入率；
- \(f_p\dot{p}\)：压力场 \(p\) 的功率共轭 \(f_p\) 与压力变化率的乘积，处理近不可压约束；
- \(f_{\tilde{\lambda}}\dot{\tilde{\lambda}}\)：非局部链拉伸 \(\tilde{\lambda}\) 的标量共轭力与其变化率之积；
- \(\boldsymbol{\xi}_{\tilde{\lambda}}\cdot\nabla\dot{\tilde{\lambda}}\)：非局部链拉伸梯度的矢量共轭力（微形态项），引入后使模型具有内禀长度尺度；
- \(f_d\dot{d}\)：损伤场 \(d\) 的热动力学共轭力与损伤速率之积。

---

### 公式 (7)：外部机械功

\[ W_{ext} = \int_{\partial\Omega_0}\mathbf{T}\cdot\mathbf{u}\,dA + \int_{\partial\Omega_0}\tilde{\iota}\,\tilde{\lambda}\,dA \]

**说明：** 外部功由宏观面力 \(\mathbf{T}\) 在宏观位移 \(\mathbf{u}\) 上做的功，以及微形态牵引力 \(\tilde{\iota}\) 在非局部链拉伸 \(\tilde{\lambda}\) 上做的功两部分构成。第二项是由于将 \(\tilde{\lambda}\) 作为独立状态变量所需引入的额外边界功。

---

### 公式 (8)：外部机械功率

\[ P_{ext} = \dot{W}_{ext} = \int_{\partial\Omega_0}\mathbf{T}\cdot\dot{\mathbf{u}}\,dA + \int_{\partial\Omega_0}\tilde{\iota}\,\dot{\tilde{\lambda}}\,dA \]

**推导：** 直接对公式 (7) 关于时间求导即得。

---

### 公式 (9)：力学平衡方程与边界条件

\[ \nabla\cdot\mathbf{P} = \mathbf{0}\quad\text{in}\;\Omega_0,\qquad \mathbf{u} = \bar{\mathbf{u}}\quad\text{on}\;\partial_D\Omega_0,\qquad \mathbf{P}\cdot\mathbf{n}_0 = \mathbf{T}\quad\text{on}\;\partial_N\Omega_0 \]

**推导：** 对虚功原理 \(\delta P_{int} = \delta P_{ext}\)，取位移虚变量 \(\delta\mathbf{u}\)，对公式 (6) 中的 \(\mathbf{P}:\nabla\delta\mathbf{u}\) 用散度定理积分：

\[ \int_{\Omega_0}\mathbf{P}:\nabla\delta\mathbf{u}\,dV = \int_{\partial\Omega_0}(\mathbf{P}\cdot\mathbf{n}_0)\cdot\delta\mathbf{u}\,dA - \int_{\Omega_0}(\nabla\cdot\mathbf{P})\cdot\delta\mathbf{u}\,dV \]

与公式 (7) 右侧对应项对比，由 \(\delta\mathbf{u}\) 的任意性得到体积分项 \(\nabla\cdot\mathbf{P} = \mathbf{0}\) 和面积分项 \(\mathbf{P}\cdot\mathbf{n}_0 = \mathbf{T}\)。

---

### 公式 (10)：压力共轭力

\[ f_p = 0 \]

**含义：** 压力场 \(p\) 是一个 Lagrange 乘子，仅用于约束不可压条件，其热动力学共轭力恒为零。这也是后续本构关系推导的一个出发点。

---

### 公式 (11)：微力平衡方程与边界条件

\[ \nabla\cdot\boldsymbol{\xi}_{\tilde{\lambda}} - f_{\tilde{\lambda}} = 0\quad\text{in}\;\Omega_0,\qquad \tilde{\lambda} = \bar{\tilde{\lambda}}\quad\text{on}\;\partial_D\Omega_0,\qquad \boldsymbol{\xi}_{\tilde{\lambda}}\cdot\mathbf{n}_0 = \tilde{\iota}\quad\text{on}\;\partial_N\Omega_0 \]

**推导：** 类比公式 (9) 的推导，对非局部链拉伸取虚变量 \(\delta\tilde{\lambda}\)，对公式 (6) 中的 \(\boldsymbol{\xi}_{\tilde{\lambda}}\cdot\nabla\delta\tilde{\lambda}\) 施用散度定理：

\[ \int_{\Omega_0}\boldsymbol{\xi}_{\tilde{\lambda}}\cdot\nabla\delta\tilde{\lambda}\,dV = \int_{\partial\Omega_0}(\boldsymbol{\xi}_{\tilde{\lambda}}\cdot\mathbf{n}_0)\delta\tilde{\lambda}\,dA - \int_{\Omega_0}(\nabla\cdot\boldsymbol{\xi}_{\tilde{\lambda}})\delta\tilde{\lambda}\,dV \]

再与 \(f_{\tilde{\lambda}}\dot{\tilde{\lambda}}\) 合并，由任意性得微力平衡方程。

---

### 公式 (12)：损伤共轭力

\[ f_d = 0 \]

**含义：** 损伤场 \(d\) 的虚功原理贡献为零（即不存在独立的损伤"惯性力"），损伤演化完全由耗散不等式驱动（见公式 (16)）。

---

### 公式 (13)：Helmholtz 自由能变化率

\[ \frac{d\Psi}{dt} = \frac{\partial\Psi}{\partial\mathbf{F}}:\dot{\mathbf{F}} + \frac{\partial\Psi}{\partial p}\dot{p} + \frac{\partial\Psi}{\partial\tilde{\lambda}}\dot{\tilde{\lambda}} + \frac{\partial\Psi}{\partial\nabla\tilde{\lambda}}:\nabla\dot{\tilde{\lambda}} + \frac{\partial\Psi}{\partial d}\dot{d} \]

**来源：** 设 \(\Psi = \Psi(\mathbf{F},p,\tilde{\lambda},\nabla\tilde{\lambda},d)\)，对时间直接应用链式法则展开全微分即得。每一项代表对应状态变量变化率对自由能变化率的贡献。

---

### 公式 (14)：Clausius-Planck 耗散不等式

\[ D = \left(\mathbf{P} - \frac{\partial\Psi}{\partial\mathbf{F}}\right):\dot{\mathbf{F}} + \left(f_p - \frac{\partial\Psi}{\partial p}\right)\dot{p} + \left(f_{\tilde{\lambda}} - \frac{\partial\Psi}{\partial\tilde{\lambda}}\right)\dot{\tilde{\lambda}} + \left(\boldsymbol{\xi}_{\tilde{\lambda}} - \frac{\partial\Psi}{\partial\nabla\tilde{\lambda}}\right)\cdot\nabla\dot{\tilde{\lambda}} + \left(f_d - \frac{\partial\Psi}{\partial d}\right)\dot{d} \geq 0 \]

**推导：** 热力学第二定律要求总内部耗散率 \(D = P_{int} - \dot{\Psi} \geq 0\)。将公式 (6) 减去公式 (13) 即得上式，其中 \(P_{int}\) 中的 \(\mathbf{P}:\nabla\mathbf{u} = \mathbf{P}:\dot{\mathbf{F}}\)（因为 \(\dot{\mathbf{F}} = \nabla\dot{\mathbf{u}}\)）。

---

### 公式 (15)：本构关系（Coleman-Noll 程序）

\[ \mathbf{P} = \frac{\partial\Psi}{\partial\mathbf{F}},\quad f_p = \frac{\partial\Psi}{\partial p},\quad \boldsymbol{\xi}_{\tilde{\lambda}} = \frac{\partial\Psi}{\partial\nabla\tilde{\lambda}},\quad f_{\tilde{\lambda}} = \frac{\partial\Psi}{\partial\tilde{\lambda}} \]

**推导：** 在 Coleman-Noll 程序中，耗散不等式 (14) 对任意过程均须满足。对于前四项，由于其对应的变化率（\(\dot{\mathbf{F}},\dot{p},\dot{\tilde{\lambda}},\nabla\dot{\tilde{\lambda}}\)）可以独立任意选取，各括号内的因子必须单独为零，从而得到公式 (15) 中的四个本构关系。

---

### 公式 (16)：损伤耗散不等式

\[ D = -\frac{\partial\Psi}{\partial d}\,\dot{d} \geq 0 \]

**推导：** 将公式 (15) 代入公式 (14)，前四项贡献消失。再利用公式 (12)（\(f_d = 0\)），第五项化简为：

\[ D = \left(0 - \frac{\partial\Psi}{\partial d}\right)\dot{d} = -\frac{\partial\Psi}{\partial d}\,\dot{d} \geq 0 \]

---

### 公式 (17)：损伤耗散总能量

\[ I = -\int_0^T\int_{\Omega}\frac{\partial\Psi}{\partial d}\,\dot{d}\,dV\,dt \]

**含义：** 对公式 (16) 中的局部耗散率在空间 \(\Omega\) 和时间 \([0,T]\) 上积分，得到损伤累积与裂纹扩展所消耗的总能量，是断裂能的宏观度量。

---

### 公式 (18)：损伤演化条件

\[ \frac{\partial\Psi}{\partial d} \leq 0\quad\text{in}\;\Omega_0 \]

**推导：** 由公式 (16) 和损伤不可逆性（\(\dot{d} \geq 0\)），直接得 \(-\frac{\partial\Psi}{\partial d}\dot{d} \geq 0\)，从而 \(\frac{\partial\Psi}{\partial d} \leq 0\)。这要求自由能关于损伤变量单调递减，即材料随损伤积累而软化，具有明确的物理意义。

---

## 2.3 模型专化（Model Specialization）

### 2.3.1 Helmholtz 自由能

#### 公式 (19)：含压力场的自由能（近不可压）

\[ \hat{\psi}(\mathbf{F},p) = \psi(\mathbf{F}) - p(J-1) - \frac{p^2}{2\kappa} \]

**推导（扰动 Lagrange 乘子法）：**
- 纯 Lagrange 乘子法严格约束 \(J=1\)（不可压），在数值计算中容易导致锁死（locking）问题。
- 扰动法引入压力场 \(p\) 作为独立变量，并附加罚项 \(-\frac{p^2}{2\kappa}\)（相当于允许轻微压缩，\(\kappa\) 控制体积模量）。
- 对 \(p\) 取变分（Euler-Lagrange 方程）：

\[ \frac{\partial\hat{\psi}}{\partial p} = -(J-1) - \frac{p}{\kappa} = 0 \quad\Rightarrow\quad p = -\kappa(J-1) \]

当 \(\kappa \to \infty\) 时，\(J \to 1\)，退化为完全不可压条件。

---

#### 公式 (20)：含损伤的局部自由能

\[ \Psi^{\mathrm{loc}}(\mathbf{F},p,d) = a(d)\psi(\mathbf{F}) - b(d)p(J-1) - \frac{p^2}{2\kappa} \]

**与公式 (19) 的联系：** 在公式 (19) 的基础上，用**退化函数** \(a(d)\) 和 \(b(d)\) 分别作用于弹性项和体积约束项，以描述损伤对材料力学响应的劣化：
- \(a(d)\) 降低剪切弹性模量；
- \(b(d)\) 降低不可压约束的严格性（裂纹张开时允许体积增大）；
- \(b\) 的阶次高于 \(a\)（即 \(b\) 随损伤衰减更快），保证在损伤区域裂纹可以物理地张开而不受不可压约束限制。

---

#### 公式 (21)：Neo-Hookean 自由能

\[ \psi(\mathbf{F}) = \frac{\mu}{2}\left(I_1(\mathbf{C}) - 3 - 2\ln J\right) \]

**来源：** Neo-Hookean 模型是最简单的超弹性本构，其中：
- \(I_1 - 3\) 是主拉伸比的无量纲化偏差度量（参考构型时 \(I_1 = 3\)）；
- \(-2\ln J\) 项保证在纯体积变形（\(J \neq 1\)）下也有能量贡献，且当 \(J=1\)（不可压）时该项为零；
- \(\mu\) 为剪切模量。

公式 (21) 即为 \(\psi(\mathbf{F})\) 的具体形式，代入公式 (20) 即完成 \(\Psi^{\mathrm{loc}}\) 的具体化。

---

#### 公式 (22)：微-宏耦合非局部自由能

\[ \Psi^{\mathrm{nloc}}_{\mathrm{micmac}}(\lambda_{ch},\tilde{\lambda}) = \frac{h_{nl}}{2}\left[\lambda_{ch} - \tilde{\lambda}\right]^2 \]

**含义：** 以平方惩罚项将局部链拉伸 \(\lambda_{ch}\)（来自宏观变形，公式 (5)）与非局部链拉伸 \(\tilde{\lambda}\)（由网络非均一性和非局部效应决定）耦合在一起。惩罚系数 \(h_{nl} > 0\) 控制两者之差的"允许量"。当 \(h_{nl} \to \infty\)，\(\tilde{\lambda} \to \lambda_{ch}\)，退化为纯局部模型。该项是引入内禀长度尺度的核心机制之一。

---

#### 公式 (23)：梯度非局部自由能

\[ \Psi^{\mathrm{nloc}}_{\mathrm{grd}}(\nabla\tilde{\lambda},d) = \frac{h_{nl}g(d)\ell^2}{2}\nabla\tilde{\lambda}\cdot\nabla\tilde{\lambda} \]

**含义：** 引入非局部链拉伸的梯度 \(\nabla\tilde{\lambda}\) 作为自由能的自变量，赋予模型空间非局部性（梯度增强效应）。长度尺度参数 \(\ell\) 控制非局部交互的影响域，具有真实的物理含义（对应聚合物网络链的平均链长或相关长度）。松弛函数 \(g(d)\) 使该梯度惩罚在损伤区域减弱（见公式 (25)）。

---

### 2.3.2 退化函数

#### 公式 (24)：退化函数定义

\[ a(d) = (1-k_\ell)(1-d)^2 + k_\ell,\qquad b(d) = (1-k_\ell)(1-d)^3 + k_\ell \]

**说明：**
- \(k_\ell \ll 1\) 是数值条件化参数（防止刚度矩阵奇异）；
- \(d=0\)（无损）时，\(a(0) \approx 1\), \(b(0) \approx 1\)，材料保持完整；
- \(d=1\)（完全损伤）时，\(a(1) = k_\ell \approx 0\), \(b(1) = k_\ell \approx 0\)，材料几乎失去承载能力；
- \(b(d)\) 的多项式阶次（3次）高于 \(a(d)\)（2次），确保不可压约束随损伤更快消失，允许裂纹张开。

---

### 2.3.3 松弛函数

#### 公式 (25)：幂次松弛函数

\[ g(d) = (1-d)^m,\quad m>0 \]

**含义：** 随着损伤 \(d\) 增大，\(g(d)\) 从 1 单调减小到 0，逐渐削弱公式 (23) 中的梯度非局部项对损伤区域的影响。参数 \(m\) 控制衰减速率：
- \(m\) 越小，衰减越陡（在损伤较大时才显著衰减）；
- 当 \(m \to 0\)，\(g(d) \to \mathbf{1}_{d<1}\)（不连续函数），对应 Poh-Sun（2017）的极限情形；
- 连续的 \(g(d)\) 形式在高度非线性数值计算中收敛性优于不连续情形。

#### 公式 (26)：指数型松弛函数（Poh-Sun 2017）

\[ g(d) = (1-R)\exp(-md) + R - \frac{\exp(-m)}{1-\exp(-m)} \]

**含义：** 这是文献中另一种松弛函数形式，参数 \(m\) 和 \(R\) 控制函数形态。对比公式 (25)，两者在不同参数下具有相似的衰减曲线，但公式 (26) 更加光滑，在 \(d=1\) 时不严格为零（有一个残差项），适合某些对数值稳定性要求更高的场景。

---

### 2.3.4 本构关系

#### 公式 (27)：总 Helmholtz 自由能

\[ \Psi(\mathbf{F},p,\tilde{\lambda},\nabla\tilde{\lambda},d) = a(d)\psi(\mathbf{F}) - b(d)p(J-1) - \frac{p^2}{2\kappa} + \frac{h_{nl}}{2}\left[\lambda_{ch}-\tilde{\lambda}\right]^2 + \frac{h_{nl}g(d)\ell^2}{2}\nabla\tilde{\lambda}\cdot\nabla\tilde{\lambda} \]

**推导：** 将 \(\Psi = \Psi^{\mathrm{loc}} + \Psi^{\mathrm{nloc}}_{\mathrm{micmac}} + \Psi^{\mathrm{nloc}}_{\mathrm{grd}}\) 展开，即将公式 (20)、(22)、(23) 相加得到。该式统一了局部弹性、不可压约束、微-宏耦合和梯度非局部四部分的能量贡献。

**验证公式 (18)：** 对损伤求偏导：

\[ \frac{\partial\Psi}{\partial d} = a'(d)\psi(\mathbf{F}) - b'(d)p(J-1) + \frac{h_{nl}g'(d)\ell^2}{2}|\nabla\tilde{\lambda}|^2 \]

由于 \(a'(d) < 0\)（\(a\) 随 \(d\) 增大而减小），\(-b'(d)p(J-1)\) 在通常情况下也满足符号要求，\(g'(d) < 0\)，三项均为非正，保证 \(\frac{\partial\Psi}{\partial d} \leq 0\)，满足耗散不等式。

---

#### 公式 (28)：具体本构关系

\[ \mathbf{P} = a(d)\frac{\partial\psi}{\partial\mathbf{F}} - b(d)pJ\mathbf{F}^{-T}\quad\text{in}\;\Omega_0 \]

\[ f_p = -b(d)(J-1) - \frac{p}{\kappa}\quad\text{in}\;\Omega_0 \]

\[ \boldsymbol{\xi}_{\tilde{\lambda}} = h_{nl}g(d)\ell^2\nabla\tilde{\lambda}\quad\text{in}\;\Omega_0 \]

\[ f_{\tilde{\lambda}} = -h_{nl}\left[\lambda_{ch}-\tilde{\lambda}\right]\quad\text{in}\;\Omega_0 \]

**推导：** 将公式 (27) 分别对 \(\mathbf{F}\)、\(p\)、\(\nabla\tilde{\lambda}\)、\(\tilde{\lambda}\) 求偏导，代入公式 (15)：

- **应力张量 \(\mathbf{P}\)：** \(\frac{\partial}{\partial\mathbf{F}}\left[a(d)\psi - b(d)p(J-1)\right] = a(d)\frac{\partial\psi}{\partial\mathbf{F}} - b(d)p\frac{\partial J}{\partial\mathbf{F}}\)。利用 \(\frac{\partial J}{\partial\mathbf{F}} = J\mathbf{F}^{-T}\) 得到公式 (28)\(_1\)。

- **压力方程：** \(f_p = \frac{\partial\Psi}{\partial p} = -b(d)(J-1) - \frac{p}{\kappa}\)，此即公式 (28)\(_2\)，同时这也是确定压力场 \(p\) 的方程（令 \(f_p = 0\) 即可求解 \(p\)）。

- **非局部矢量共轭力：** \(\boldsymbol{\xi}_{\tilde{\lambda}} = \frac{\partial\Psi}{\partial\nabla\tilde{\lambda}} = h_{nl}g(d)\ell^2\nabla\tilde{\lambda}\)，即公式 (28)\(_3\)。

- **非局部标量共轭力：** \(f_{\tilde{\lambda}} = \frac{\partial\Psi}{\partial\tilde{\lambda}} = \frac{h_{nl}}{2}\cdot 2(\lambda_{ch}-\tilde{\lambda})\cdot(-1) = -h_{nl}(\lambda_{ch}-\tilde{\lambda})\)，即公式 (28)\(_4\)。

---

## 2.4 损伤函数（Damage Function）

### 公式 (29)：损伤共轭力的本构定义

\[ f_d(d,\tilde{\lambda},\mathcal{H}) = \begin{cases} d & \text{if } \tilde{\lambda} < \lambda_{cr} \\ d - 1 + \dfrac{\lambda_{cr}-1}{\mathcal{H}(\tilde{\lambda})-1}\left(1 - c + c\,e^{-\gamma(\mathcal{H}(\tilde{\lambda})-\lambda_{cr})}\right) & \text{if } \tilde{\lambda} \geq \lambda_{cr} \end{cases} \]

**说明：** 这是对损伤共轭力 \(f_d\) 的本构规定（而非推导结果），类似于强度准则。参数含义：
- \(\lambda_{cr}\)：损伤起始阈值（临界链拉伸），\(\tilde{\lambda} < \lambda_{cr}\) 时无损伤；
- \(0 \leq c \leq 1\)：调节最大损伤；
- \(\gamma \geq 0\)：控制从无损到全损的过渡陡度；
- \(\mathcal{H}(\tilde{\lambda})\)：历史函数（公式 (30)），保证损伤的不可逆性。

---

### 公式 (30)：历史函数

\[ \mathcal{H}(\tilde{\lambda};\mathbf{X},t) = \max_{s\in[0,t]}\tilde{\lambda}(\mathbf{X},s) \]

**含义：** 历史函数记录了非局部链拉伸 \(\tilde{\lambda}\) 在每一质点 \(\mathbf{X}\) 处历史上达到的最大值。这直接保证了损伤的单调不可逆性：一旦损伤发生，即使后续卸载也不会恢复（因为 \(\mathcal{H}\) 不会减小）。

---

### 公式 (31)：损伤函数（强形式解）

\[ d = \begin{cases} 0 & \text{if } \tilde{\lambda} < \lambda_{cr} \\ 1 - \dfrac{\lambda_{cr}-1}{\mathcal{H}(\tilde{\lambda})-1}\left(1 - c + c\,e^{-\gamma(\mathcal{H}(\tilde{\lambda})-\lambda_{cr})}\right) & \text{if } \tilde{\lambda} \geq \lambda_{cr} \end{cases} \]

**推导：** 由公式 (12)（\(f_d = 0\)），将公式 (29) 中的 \(f_d\) 令为零，直接求解 \(d\)：
- 当 \(\tilde{\lambda} < \lambda_{cr}\)：\(f_d = d = 0\)，故 \(d = 0\)；
- 当 \(\tilde{\lambda} \geq \lambda_{cr}\)：令 \(f_d = 0\)，则：

\[ d = 1 - \frac{\lambda_{cr}-1}{\mathcal{H}(\tilde{\lambda})-1}\left(1 - c + c\,e^{-\gamma(\mathcal{H}(\tilde{\lambda})-\lambda_{cr})}\right) \]

**性质验证：**
- \(\tilde{\lambda} = \lambda_{cr}\)（刚达阈值）时，\(\mathcal{H} = \lambda_{cr}\)，指数项 \(e^0 = 1\)，括号内为 1，分式为 1，\(d = 0\)——连续；
- \(\mathcal{H} \to \infty\) 时，分式趋近于 0，\(d \to 1\)——材料完全损伤；
- 由于 \(\mathcal{H}\) 单调不减，\(d\) 也单调不减，满足不可逆性。

---

## 2.5 强形式与弱形式（Strong and Weak Forms）

### 公式 (32)：微力平衡方程（代入本构关系后）

\[ h_{nl}\left[\lambda_{ch} - \tilde{\lambda} + \nabla\cdot(g(d)\ell^2\nabla\tilde{\lambda})\right] = 0\quad\text{in}\;\Omega_0 \]

\[ h_{nl}g(d)\ell^2\nabla\tilde{\lambda}\cdot\mathbf{n}_0 = \hat{\iota}\quad\text{on}\;\partial_N\Omega_0 \]

**推导：** 将公式 (28)\(_3\) 和 (28)\(_4\) 代入公式 (11) 的微力平衡方程：

\[ \nabla\cdot\boldsymbol{\xi}_{\tilde{\lambda}} - f_{\tilde{\lambda}} = \nabla\cdot\left(h_{nl}g(d)\ell^2\nabla\tilde{\lambda}\right) - \left(-h_{nl}(\lambda_{ch}-\tilde{\lambda})\right) = h_{nl}\left[\lambda_{ch} - \tilde{\lambda} + \nabla\cdot(g(d)\ell^2\nabla\tilde{\lambda})\right] = 0 \]

边界条件类似代入 \(\boldsymbol{\xi}_{\tilde{\lambda}}\cdot\mathbf{n}_0 = \hat{\iota}\) 得到。

---

### 公式 (33)：约化后的微力平衡方程

\[ \lambda_{ch} - \tilde{\lambda} + \nabla\cdot(g(d)\ell^2\nabla\tilde{\lambda}) = 0\quad\text{in}\;\Omega_0 \]

\[ \nabla\tilde{\lambda}\cdot\mathbf{n}_0 = 0\quad\text{on}\;\partial_N\Omega_0 \]

**推导：** 由于 \(h_{nl} > 0\)，可以在公式 (32) 两侧除以 \(h_{nl}\)；自然边界条件中令 \(\hat{\iota} = 0\)（微形态牵引力在外边界可忽略），得到公式 (33)。

---

### 公式 (34)：完整强形式

**力学平衡（位移方程）：**

\[ \nabla\cdot\left[a(d)\frac{\partial\psi}{\partial\mathbf{F}} - b(d)pJ\mathbf{F}^{-T}\right] = \mathbf{0}\quad\text{in}\;\Omega_0 \]

**压力方程（近不可压约束）：**

\[ -b(d)(J-1) - \frac{p}{\kappa} = 0\quad\text{in}\;\Omega_0 \]

**非局部链拉伸方程：**

\[ \tilde{\lambda} - \lambda_{ch} - \nabla\cdot(g(d)\ell^2\nabla\tilde{\lambda}) = 0\quad\text{in}\;\Omega_0 \]

**推导：** 
- 第一式由公式 (9)\(_1\) 代入公式 (28)\(_1\) 得到；  
- 第二式由公式 (15) 中 \(f_p = \frac{\partial\Psi}{\partial p}\) 并令 \(f_p = 0\)（公式 (10)）得到公式 (28)\(_2\)，即 \(-b(d)(J-1)-\frac{p}{\kappa} = 0\)；  
- 第三式来自公式 (33)，整理符号后得。

---

### 公式 (35)：边界条件（与强形式配套）

\[ \mathbf{u} = \bar{\mathbf{u}}\quad\text{on}\;\partial_D\Omega_0,\qquad \mathbf{P}\cdot\mathbf{n}_0 = \mathbf{T}\quad\text{on}\;\partial_N\Omega_0 \]

\[ \tilde{\lambda} = \bar{\tilde{\lambda}}\quad\text{on}\;\partial_D\Omega_0,\qquad \nabla\tilde{\lambda}\cdot\mathbf{n}_0 = 0\quad\text{on}\;\partial_N\Omega_0 \]

---

### 公式 (36)：引入链拉伸上界的修正非局部方程

\[ \tilde{\lambda} - \lambda_{ch}^{\mathrm{max}} + \ll\lambda_{ch}^{\mathrm{max}} - \lambda_{ch}\gg - \nabla\cdot(g(d)\ell^2\nabla\tilde{\lambda}) = 0\quad\text{in}\;\Omega_0 \]

**推导动机：** 公式 (34)\(_3\) 中的驱动力 \(\lambda_{ch}\) 在裂纹张开时会持续增大，导致非物理的损伤区扩展（damage zone broadening）。统计力学告诉我们，当链拉伸超过临界值 \(\lambda_{ch}^{\mathrm{max}}\) 后，所有链已断裂，\(\lambda_{ch}\) 不再有效增大。

**具体修改：** 在公式 (22) 中，用 \(\lambda_{ch}^{\mathrm{max}} - \ll\lambda_{ch}^{\mathrm{max}} - \lambda_{ch}\gg\) 替代 \(\lambda_{ch}\)。利用公式 (37) 中 Macaulay 括号的定义：

\[ \lambda_{ch}^{\mathrm{max}} - \ll\lambda_{ch}^{\mathrm{max}} - \lambda_{ch}\gg = \begin{cases} \lambda_{ch}^{\mathrm{max}} - (\lambda_{ch}^{\mathrm{max}} - \lambda_{ch}) = \lambda_{ch}, & \text{if } \lambda_{ch} < \lambda_{ch}^{\mathrm{max}} \\ \lambda_{ch}^{\mathrm{max}}, & \text{if } \lambda_{ch} \geq \lambda_{ch}^{\mathrm{max}} \end{cases} \]

即当 \(\lambda_{ch} < \lambda_{ch}^{\mathrm{max}}\) 时，与原始公式 (34)\(_3\) 一致；当 \(\lambda_{ch} \geq \lambda_{ch}^{\mathrm{max}}\) 时，驱动力被截断为 \(\lambda_{ch}^{\mathrm{max}}\)，防止非物理增大。这一修改在代码中对应：`Lmbda_ch = conditional(sqrt(Ic/3) > Lmbda_ch_max, Lmbda_ch_max, sqrt(Ic/3))`。

---

### 公式 (37)：Macaulay 括号

\[ \ll x \gg = \begin{cases} x, & \text{if } x > 0 \\ 0, & \text{if } x \leq 0 \end{cases} \]

**含义：** Macaulay 括号（正值截断算子）等价于 \(\ll x \gg = \max(x,0)\)。在本文中用于实现链拉伸的物理上界截断（公式 (36)）。

---

### 公式 (38)：试验函数空间

\[ \mathbb{U} = \{\mathbf{u}\in H^1(\Omega_0);\;\mathbf{u} = \bar{\mathbf{u}}\;\text{on}\;\partial_D\Omega_0\} \]

\[ \mathbb{P} = \{p \in L^2(\Omega_0)\} \]

\[ \mathbb{L} = \{\tilde{\lambda}\in H^1(\Omega_0);\;\tilde{\lambda} = \bar{\tilde{\lambda}}\;\text{on}\;\partial_D\Omega_0\} \]

**含义：** 
- \(H^1(\Omega_0)\)（一阶 Sobolev 空间）：函数本身及其一阶弱导数均属 \(L^2\)，保证应变能有界；
- \(L^2(\Omega_0)\)（Lebesgue 平方可积空间）：压力场仅需平方可积，不要求导数有界（与混合有限元中压力的插值阶次一致）；
- Dirichlet 条件在试验函数空间中强制施加。

---

### 公式 (39)：检验函数空间

\[ \mathbb{V} = \{\mathbf{v}\in H_0^1(\Omega_0);\;\mathbf{v} = 0\;\text{on}\;\partial_D\Omega_0\} \]

\[ \mathbb{Q} = \{q \in L^2(\Omega_0)\} \]

\[ \mathbb{B} = \{\beta\in H_0^1(\Omega_0);\;\beta = 0\;\text{on}\;\partial_D\Omega_0\} \]

**含义：** 检验函数在 Dirichlet 边界上为零（齐次边界条件），使弱形式推导时边界积分自然消失。

---

### 公式 (40)：弱形式（Galerkin 变分方程）

**位移（力学平衡）弱形式：**

\[ \int_{\Omega_0}\left[a(d)\frac{\partial\psi}{\partial\mathbf{F}} - b(d)pJ\mathbf{F}^{-T}\right]:\nabla\mathbf{v}\,dV - \int_{\partial_N\Omega_0}\mathbf{T}\cdot\mathbf{v}\,dA = 0 \]

**压力（不可压约束）弱形式：**

\[ -\int_{\Omega_0}\left(b(d)(J-1) + \frac{p}{\kappa}\right)q\,dV = 0 \]

**非局部链拉伸弱形式：**

\[ \int_{\Omega_0}\tilde{\lambda}\beta\,dV - \int_{\Omega_0}\left[\lambda_{ch}^{\mathrm{max}} - \ll\lambda_{ch}^{\mathrm{max}} - \lambda_{ch}\gg\right]\beta\,dV + \int_{\Omega_0}g(d)\ell^2\nabla\tilde{\lambda}\cdot\nabla\beta\,dV = 0 \]

**推导：** 以第三个方程为例，从公式 (36) 出发，乘以检验函数 \(\beta\) 并在 \(\Omega_0\) 上积分：

\[ \int_{\Omega_0}\left[\tilde{\lambda} - \lambda_{ch}^{\mathrm{max}} + \ll\lambda_{ch}^{\mathrm{max}} - \lambda_{ch}\gg\right]\beta\,dV - \int_{\Omega_0}\nabla\cdot(g(d)\ell^2\nabla\tilde{\lambda})\,\beta\,dV = 0 \]

对最后一项用散度定理（Green 公式）分部积分：

\[ -\int_{\Omega_0}\nabla\cdot(g(d)\ell^2\nabla\tilde{\lambda})\beta\,dV = \int_{\Omega_0}g(d)\ell^2\nabla\tilde{\lambda}\cdot\nabla\beta\,dV - \int_{\partial\Omega_0}g(d)\ell^2(\nabla\tilde{\lambda}\cdot\mathbf{n}_0)\beta\,dA \]

自然边界条件 \(\nabla\tilde{\lambda}\cdot\mathbf{n}_0 = 0\)（公式 (33)）使边界积分项消失，最终得到公式 (40)\(_3\)。

前两个弱形式的推导类似：将强形式 (34)\(_1\) 和 (34)\(_2\) 乘以对应检验函数，散度定理 + Neumann 边界条件即可。

---

## 公式整体联系总结

下图给出本节各公式之间的逻辑脉络：

```
运动学描述
  (1)-(2) 变形梯度 F
       ↓
  (3) 右Cauchy-Green张量 C
       ↓
  (4) 主不变量 I1, I2, I3
       ↓
  (5) 链拉伸 λ_ch = √(I1/3)          ← 连接分子与连续介质尺度

热力学框架
  (6)-(8) 虚功原理（内/外功率）
       ↓
  (9)-(12) 平衡方程（散度定理）
       ↓
  (13) 自由能变化率（链式法则）
       ↓
  (14) Clausius-Planck耗散不等式
       ↓
  (15) 本构关系（Coleman-Noll）
       ↓
  (16)-(18) 损伤耗散不等式

模型具体化
  (19) 近不可压自由能（扰动Lagrange法）
       ↓
  (20)-(23) 含损伤的总自由能（局部+非局部）
       ↓
  (24) 退化函数 a(d), b(d)
  (25)-(26) 松弛函数 g(d)
       ↓
  (27) 总Helmholtz自由能 Ψ
       ↓
  (28) 代入(15)得到具体本构关系

损伤模型
  (29)-(30) 损伤共轭力本构 + 历史函数
       ↓
  (31) 损伤函数 d(𝒩(λ̃))  ← 来自 f_d = 0 (公式12)，即 d(ℋ(λ̃))

强/弱形式
  (28) 代入 (9)(11)(10) →  (32)-(35) 强形式
       ↓
  (36)-(37) 引入上界截断（物理改进）
       ↓
  (38)-(39) 函数空间
       ↓
  (40) 弱形式（Galerkin分部积分）→ 有限元离散化
```

---

*本文档所有公式均采用 LaTeX 渲染，行内变量使用 \(\cdot\) 格式，显示公式使用 \[\cdot\] 格式。*
