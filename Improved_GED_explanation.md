# Improved_GED.py 代码详细解释

> **说明**：本文档按主函数运行顺序逐段解释 `Improved_GED.py` 的每一行代码。该程序基于 FEniCS 框架实现了一种**改进梯度增强损伤（Modified/Improved Gradient-Enhanced Damage, GED）**模型，用于模拟超弹性材料的平面应变断裂力学问题。涉及计算的代码段均附有对应的数学公式。

---

## 目录

1. [导入库](#1-导入库-第-1319-行)
2. [FEniCS 与求解器参数配置](#2-fenics-与求解器参数配置-第-2256-行)
3. [局部投影辅助函数](#3-局部投影辅助函数-第-6073-行)
4. [边界子域类定义](#4-边界子域类定义-第-7792-行)
5. [损伤函数与材料参数](#5-损伤函数与材料参数-第-95118-行)
6. [几何与数值参数](#6-几何与数值参数-第-120136-行)
7. [网格生成与子域设置](#7-网格生成与子域设置-第-139168-行)
8. [边界标记与积分测度](#8-边界标记与积分测度-第-170186-行)
9. [函数空间与场变量定义](#9-函数空间与场变量定义-第-191230-行)
10. [位移边界条件](#10-位移边界条件-第-233243-行)
11. [运动学描述](#11-运动学描述-第-248257-行)
12. [本构模型与弹性能量泛函](#12-本构模型与弹性能量泛函-第-260290-行)
13. [损伤演化方程设置](#13-损伤演化方程设置-第-293319-行)
14. [主求解循环](#14-主求解循环-第-322453-行)
15. [后处理与结果输出](#15-后处理与结果输出-第-456488-行)

---

## 1. 导入库（第 13–19 行）

```python
13: from __future__ import division
14: from dolfin import *
15: from mshr import *
16: import os
17: import sympy
18: import numpy as np
19: import matplotlib.pyplot as plt
```

**解释：**

- **第 13 行**：`from __future__ import division` — 确保 Python 2 中整数相除返回浮点数（Python 3 中已为默认行为），保证数值计算精度。
- **第 14 行**：`from dolfin import *` — 导入 FEniCS 核心库 **DOLFIN**，提供有限元网格、函数空间、变分形式、求解器等全部功能。
- **第 15 行**：`from mshr import *` — 导入 FEniCS 的几何建模库 **mshr**，用于构造圆形、矩形等几何体并生成网格。
- **第 16 行**：`import os` — 标准库，用于创建目录、处理文件路径。
- **第 17 行**：`import sympy` — 符号计算库（本脚本中已导入但未显式调用，预留用于符号推导）。
- **第 18 行**：`import numpy as np` — 数值计算库，用于数组操作、向量运算、文件输出等。
- **第 19 行**：`import matplotlib.pyplot as plt` — 绘图库，用于可视化网格和后处理结果。

---

## 2. FEniCS 与求解器参数配置（第 22–56 行）

```python
22: set_log_level(LogLevel.WARNING)
26: parameters["form_compiler"]["representation"]="uflacs"
27: parameters["form_compiler"]["optimize"]=True
28: parameters["form_compiler"]["cpp_optimize"]=True
29: parameters["form_compiler"]["quadrature_degree"]=2
30: info(parameters,True)
```

**解释：**

- **第 22 行**：将 DOLFIN 的日志级别设置为 `WARNING`，只输出警告及以上级别的信息，减少终端噪音。
- **第 26–29 行**：配置 FFC（FEniCS 编译器）参数：
  - `"uflacs"` — 使用 UFL 代数计算系统（Unified Form Language Algebraic Compiler System）作为编译后端，适合复杂超弹性本构；
  - `optimize=True` / `cpp_optimize=True` — 开启生成 C++ 代码及运行时 JIT 编译优化；
  - `quadrature_degree=2` — 数值积分使用 2 阶 Gauss 求积（对二次位移场为精确积分）。
- **第 30 行**：打印所有参数树（调试用，实际运行时可注释）。

```python
34: solver_u_parameters = {"nonlinear_solver": "snes", ...}
48: solver_Lmbda_parameters = {"nonlinear_solver": "snes", ...}
```

**解释：**

- **第 34–45 行**：位移问题的非线性求解器参数字典：
  - 使用 PETSc **SNES**（Scalable Nonlinear Equations Solver）框架；
  - 线性求解器选用 **MUMPS**（直接法，适合小到中等规模问题）；
  - 非线性迭代方法：`newtontr`（Trust-Region Newton 方法），对大变形问题的病态矩阵有更好的鲁棒性；
  - 线搜索：`cp`（临界点线搜索）；
  - 各类容差均设为 $[ 10^{-6} ]$；最大迭代 300 次。（`error_on_nonconvergence: False` 表示即使未收敛也不中止程序。）
- **第 48–56 行**：损伤变量 $[ \Lambda ]$ 的求解器参数：
  - 方法为 `vinewtonssls`（**Variable-Inequality Newton Semi-Smooth Line Search**），支持**带界约束（Box Constraint）**的非线性问题，确保损伤变量满足下界 $[ \Lambda_{lb} ]$ 和上界 $[ \Lambda_{ub} ]$。

---

## 3. 局部投影辅助函数（第 60–73 行）

```python
60: def local_project(v, V, u=None):
61:     dv = TrialFunction(V)
62:     v_ = TestFunction(V)
63:     a_proj = inner(dv, v_)*dx
64:     b_proj = inner(v, v_)*dx
65:     solver = LocalSolver(a_proj, b_proj)
66:     solver.factorize()
67:     if u is None:
68:         u = Function(V)
69:         solver.solve_local_rhs(u)
70:         return u
71:     else:
72:         solver.solve_local_rhs(u)
73:         return
```

**解释：**

该函数执行**逐单元（element-wise）局部 $[ L^2 ]$ 投影**，将任意 UFL 表达式 $[ v ]$ 投影到函数空间 $[ V ]$ 中。

对应的变分投影问题为：在每个单元 $[ K ]$ 上求 $[ u_h \in V_h ]$，使得：

$$[ \int_K u_h \cdot \delta v \, \mathrm{d}x = \int_K v \cdot \delta v \, \mathrm{d}x \quad \forall \, \delta v \in V_h ]$$

由于使用 `DG`（Discontinuous Galerkin）空间（见后文），每个单元的质量矩阵是解耦的，可用 `LocalSolver` 高效逐单元求解，无需全局组装。此函数用于后处理时投影应力张量 $[ \mathbf{P} ]$、变形梯度 $[ \mathbf{F} ]$ 和体积比 $[ J ]$。

---

## 4. 边界子域类定义（第 77–92 行）

```python
77: class bot_boundary(SubDomain):
78:     def inside(self, x, on_boundary):
79:         return on_boundary and near(x[1], 0, 0.01*hsize)
81: class top_boundary(SubDomain):
82:     def inside(self, x, on_boundary):
83:         return on_boundary and near(x[1], H, 0.01*hsize)
85: class pin_point(SubDomain):
86:     def inside(self, x, on_boundary):
87:         return near(x[0], L/2, hsize) and near(x[1], 0, 0.01*hsize)
90: bot_boundary = bot_boundary()
91: top_boundary = top_boundary()
92: pin_point    = pin_point()
```

**解释：**

三个子域类继承自 FEniCS 的 `SubDomain`，通过重载 `inside` 方法定义空间位置判断逻辑：

- **`bot_boundary`（第 77–79 行）**：下边界，即满足 $[ x_2 = 0 ]$ 且在几何边界上（容差 $[ 0.01 h_{\text{size}} ]$）的点集。
- **`top_boundary`（第 81–83 行）**：上边界，即满足 $[ x_2 = H ]$ 且在几何边界上的点集。
- **`pin_point`（第 85–87 行）**：销钉点，位于 $[ x_1 = L/2, \, x_2 = 0 ]$ 附近，容差为 $[ h_{\text{size}} ]$，用于防止刚体运动（pin the rigid body mode）。
- **第 90–92 行**：将类实例化为对象，供后续标记使用。

$$[ \text{bot\_boundary} = \{ \mathbf{x} \in \partial\Omega \mid x_2 \approx 0 \}, \quad \text{top\_boundary} = \{ \mathbf{x} \in \partial\Omega \mid x_2 \approx H \} ]$$

---

## 5. 损伤函数与材料参数（第 95–118 行）

```python
95:  lambda_i = 1.2
96:  Lmbda_ch_max = 2.7
97:  alpha_damage = 1
98:  beta_damage = 20
101: mu    = 1
102: kappa = 1000
103: k_ell = 1e-5
104: eta = 0
107: body_force = Constant((0, 0))
108: load_min = 0
109: load_max = 0.6
110: load_steps = 100
```

**解释：**

**损伤参数（第 95–98 行）：**

- $[ \lambda_i = 1.2 ]$：损伤阈值，链段拉伸比 $[ \Lambda_{\mathrm{ch}} ]$ 超过此值时损伤开始演化。
- $[ \Lambda_{\mathrm{ch,max}} = 2.7 ]$：链段拉伸比的上限截断值，防止数值奇异。
- $[ \alpha_d = 1 ]$：损伤幅值系数（最大损伤量）。
- $[ \beta_d = 20 ]$：损伤演化曲线的锐度系数，值越大损伤跃变越剧烈。

损伤变量 $[ d ]$ 由链段拉伸比 $[ \Lambda ]$ 按如下公式计算（见第 369 行）：

$$[ d(\Lambda) = 1 - \frac{\lambda_i - 1}{\Lambda - 1} \left[ 1 - \alpha_d + \alpha_d \exp\!\left(-\beta_d (\Lambda - \lambda_i)\right) \right], \quad \Lambda \geq \lambda_i ]$$

当 $[ \Lambda < \lambda_i ]$ 时，$[ d = 0 ]$（无损伤）。

**材料参数（第 101–104 行）：**

- $[ \mu = 1 ]$：剪切模量（无量纲归一化）。
- $[ \kappa = 1000 ]$：体积模量（大值 $[ \kappa \gg \mu ]$ 使材料接近不可压缩）。
- $[ k_\ell = 10^{-5} ]$：残余刚度系数，防止损伤区域刚度矩阵奇异（数值稳定用）。
- $[ \eta = 0 ]$：粘性系数（本算例中关闭粘性正则化）。

**加载参数（第 107–110 行）：**

- `body_force` — 体力 $[ \mathbf{b} = \mathbf{0} ]$（本算例无体力）。
- 位移荷载从 $[ t = 0 ]$ 线性增加至 $[ t = 0.6 ]$，共 100 个荷载步。

---

## 6. 几何与数值参数（第 113–136 行）

```python
113: L, H = 1, 1
114: resolution = 40
115: segment = 20
116: hsize = 0.02
117: ell_multi = 2
118: ell = Constant(ell_multi*hsize)
121: maxiteration = 100
122: AM_tolerance = 2e-3
125: r2 = 0.47
126: r1 = 0.45
127: cirle_center = 0.5
129: crack_thickness = 0.04
130: crack_tip = 0.2
132: g_power = 0.28
135: simulation_params = "eta_%0.1f_ell_%0.2f_..."
136: savedir = simulation_params + "/"
```

**解释：**

- **第 113 行**：计算域为单位正方形 $[ \Omega = [0,1] \times [0,1] ]$，$[ L = H = 1 ]$。
- **第 114 行**：`resolution = 40` — mshr 网格生成的分辨率参数（控制最大网格单元尺寸）。
- **第 116 行**：$[ h_{\mathrm{size}} = 0.02 ]$ — 特征网格尺寸，用于边界识别容差与长度参数。
- **第 117–118 行**：相场正则化长度参数 $[ \ell = 2 \times h_{\mathrm{size}} = 0.04 ]$，定义损伤弥散区域宽度。在梯度损伤模型中，$[ \ell ]$ 控制损伤的空间梯度惩罚强度：

$$[ \ell = \ell_{\mathrm{multi}} \times h_{\mathrm{size}} ]$$

- **第 121–122 行**：交替最小化（Alternate Minimization）算法的最大迭代次数 100 次，收敛容差 $[ \varepsilon_{\mathrm{AM}} = 2 \times 10^{-3} ]$。
- **第 125–127 行**：J 积分计算所用的权函数环形区域：内圆半径 $[ r_1 = 0.45 ]$，外圆半径 $[ r_2 = 0.47 ]$，圆心 $[ (0.5, 0.5H) ]$。
- **第 129–130 行**：初始裂缝厚度 $[ 0.04 ]$，初始裂尖位置 $[ x_{\mathrm{tip}} = 0.2 ]$。
- **第 132 行**：松弛函数指数 $[ g_p = 0.28 ]$（见后文松弛函数 $[ \mathcal{R}(d) ]$）。
- **第 135–136 行**：根据关键参数生成输出目录名（便于参数研究时区分不同算例）。

---

## 7. 网格生成与子域设置（第 139–168 行）

```python
139: domain = Rectangle(Point(0.0, 0.0), Point(L, H))
140: domain.set_subdomain(1, Circle(...r2) - Circle(...r1))
142: domain.set_subdomain(2,  Rectangle(Point(0, 0.47*H),  Point(L, 0.475*H)))
143: domain.set_subdomain(3,  Rectangle(Point(0, 0.475*H), Point(L, 0.48*H)))
...
153: domain.set_subdomain(13, Rectangle(Point(0, 0.525*H), Point(L, 0.53*H)))
155: mesh = generate_mesh(domain, resolution)
156: boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim()-1, mesh.domains())
157: top_boundary.mark(boundary_markers, 1)
158: ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)
159: mf = MeshFunction("size_t", mesh, 2, mesh.domains())
160: dx = Measure("dx", domain=mesh, subdomain_data=mf)
163: if not os.path.exists(savedir):
164:     os.makedirs(savedir)
165: plt.figure()
166: plot(mesh)
167: plt.savefig(savedir + "mesh.pdf")
168: plt.close()
```

**解释：**

- **第 139 行**：创建矩形计算域 $[ \Omega = [0,L] \times [0,H] ]$。
- **第 140 行**：设置子域 1 为环形区域（外圆减内圆），用于 J 积分的域积分计算：

$$[ \Omega_J = \{ \mathbf{x} : r_1 \leq \|\mathbf{x} - \mathbf{x}_c\| \leq r_2 \} ]$$

  其中 $[ \mathbf{x}_c = (0.5, 0.5H) ]$ 为圆心。

- **第 142–153 行**：设置子域 2–13 为沿裂缝路径（$[ y \approx 0.5H ]$ 附近）的细长矩形带，用于裂缝区域的网格加密控制（mshr 对指定子域自动加密）。这 12 个水平带覆盖区间 $[ y \in [0.47H, 0.53H] ]$，共宽 $[ 0.06H ]$，对应裂缝可能扩展的区域。

- **第 155 行**：调用 mshr 在指定分辨率下生成非结构三角形网格。

- **第 156–160 行**：
  - 创建**边界网格函数** `boundary_markers`（维度 = 拓扑维度 - 1 = 1），继承子域标记，再将顶部边界标记为 1；
  - 创建**面网格函数** `mf`（维度 = 2），继承子域标记，用于区分不同积分子域；
  - 定义带子域数据的面积分测度 $[ \mathrm{d}x ]$ 和边界积分测度 $[ \mathrm{d}s ]$：

$$[ \int_{\Omega_i} (\cdot)\, \mathrm{d}x(i), \quad \int_{\Gamma_j} (\cdot)\, \mathrm{d}s(j) ]$$

- **第 163–168 行**：如输出目录不存在则创建，绘制并保存网格图像（PDF 格式）。

---

## 8. 边界标记与积分测度（第 172–186 行）

```python
172: lines  = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
173: points = MeshFunction("size_t", mesh, mesh.topology().dim()-2)
176: lines.set_all(0)
177: bot_boundary.mark(lines, 1)
178: top_boundary.mark(lines, 1)
179: file_results = XDMFFile(savedir + "/" + "lines.xdmf")
180: file_results.write(lines)
183: points.set_all(0)
184: pin_point.mark(points, 1)
185: file_results = XDMFFile(savedir + "/" + "points.xdmf")
186: file_results.write(points)
```

**解释：**

- **第 172 行**：创建边（线段）级别的网格函数，用于可视化边界线。
- **第 173 行**：创建点（顶点）级别的网格函数，用于可视化关键点。
- **第 176–180 行**：将上下边界线段标记为 1，写入 XDMF 文件供 ParaView 可视化验证。
- **第 183–186 行**：将销钉点标记为 1，写入 XDMF 文件。

这些操作主要用于**可视化调试**，确认边界条件施加位置正确。

---

## 9. 函数空间与场变量定义（第 191–230 行）

```python
191: T_DG0 = TensorFunctionSpace(mesh,'DG',0)
192: DG0   = FunctionSpace(mesh,'DG',0)
194: V_CG2 = VectorFunctionSpace(mesh, "Lagrange", 2)
196: CG1 = FunctionSpace(mesh, "Lagrange", 1)
197: V_CG2elem = V_CG2.ufl_element()
198: CG1elem = CG1.ufl_element()
200: MixElem = MixedElement([V_CG2elem, CG1elem])
202: V = FunctionSpace(mesh, MixElem)
205: w_p = Function(V)
206: u_p = TrialFunction(V)
207: v_q = TestFunction(V)
208: (u, p) = split(w_p)
209: (v, q) = split(v_q)
211: Lmbda  = Function(CG1)
212: Lmbda_trial = TrialFunction(CG1)
213: beta   = TestFunction(CG1)
214: damage = Function(CG1)
215: damage_previous = Function(CG1)
218: PTensor = Function(T_DG0, name="Nominal Stress")
219: FTensor = Function(T_DG0, name="Deformation Gradient")
220: JScalar = Function(CG1, name="Volume Ratio")
```

**解释：**

本程序采用**稳定混合有限元（Stabilized Mixed FEM）**处理近不可压缩超弹性问题，避免体积锁定（volumetric locking）。

- **第 191–192 行**：`T_DG0` — 0 阶间断 Galerkin（DG）张量函数空间，`DG0` — 标量版本，用于应力/变形梯度的逐单元后处理。
- **第 194 行**：`V_CG2` — 2 阶 Lagrange 向量空间，存储位移场 $[ \mathbf{u} ]$（二次插值）。
- **第 196 行**：`CG1` — 1 阶连续 Lagrange 标量空间，存储压力场 $[ p ]$ 和损伤场 $[ \Lambda ]$、$[ d ]$。
- **第 200–202 行**：定义混合元素 $[ \mathbf{V} = \mathbf{CG2} \times \mathrm{CG1} ]$，即 **Taylor-Hood 型**混合格式（$[ P_2/P_1 ]$），满足 inf-sup 稳定性条件：

$$[ \mathbf{V}_h = \{ (\mathbf{u}_h, p_h) \in [\mathcal{P}_2]^2 \times \mathcal{P}_1 \} ]$$

- **第 205–209 行**：`w_p` 是混合场函数（含位移和压力），`u_p`/`v_q` 分别是试验函数和测试函数；`split` 将混合场分离为位移分量 $[ \mathbf{u} ]$ 和压力分量 $[ p ]$。
- **第 211–215 行**：`Lmbda` — 网络链段拉伸比场 $[ \Lambda(\mathbf{x}) ]$；`damage`/`damage_previous` — 当前和前一时步损伤场 $[ d(\mathbf{x}) ]$（确保损伤不可逆性）。
- **第 218–220 行**：后处理用函数：名义应力张量 $[ \mathbf{P} ]$、变形梯度 $[ \mathbf{F} ]$、体积比 $[ J ]$。

**J 积分权函数（第 223–230 行）：**

```python
223: q_func_space = FunctionSpace(mesh, 'CG', 1)
224: q_expr = Expression('sqrt(...) < r1 ? 1.0 : (sqrt(...) > r2 ? 0.0 : (r2 - sqrt(...))/(r2-r1))', ...)
228: q_J_integral = interpolate(q_expr, q_func_space)
```

定义 J 积分的平滑权函数 $[ q(\mathbf{x}) ]$，在裂尖附近内圆内取 1，外圆外取 0，环形区域内线性插值：

$$[ q(\mathbf{x}) = \begin{cases} 1, & \|\mathbf{x} - \mathbf{x}_c\| \leq r_1 \\ \dfrac{r_2 - \|\mathbf{x} - \mathbf{x}_c\|}{r_2 - r_1}, & r_1 < \|\mathbf{x} - \mathbf{x}_c\| < r_2 \\ 0, & \|\mathbf{x} - \mathbf{x}_c\| \geq r_2 \end{cases} ]$$

---

## 10. 位移边界条件（第 233–243 行）

```python
233: u1 = Expression([0, "(-t/L)*x[0] + t"], t=0.0, L=L, degree=1)
234: u2 = Expression([0, "(t/L)*x[0] - t"], t=0.0, L=L, degree=1)
236: bc_u1 = DirichletBC(V.sub(0), u1, top_boundary)
237: bc_u2 = DirichletBC(V.sub(0), u2, bot_boundary)
238: bc_u  = [bc_u1, bc_u2]
241: bc_Lmbda_B = DirichletBC(CG1, 1, bot_boundary)
242: bc_Lmbda_T = DirichletBC(CG1, 1, top_boundary)
243: bc_Lmbda = [bc_Lmbda_B, bc_Lmbda_T]
```

**解释：**

- **第 233–238 行**：施加纯剪切荷载——上边界（$[ x_2 = H ]$）和下边界（$[ x_2 = 0 ]$）分别施加反对称的水平梯度位移：

上边界位移（$[ t ]$ 为荷载乘子，随时间步增加）：

$$[ \mathbf{u}\big|_{x_2=H} = \left(0,\; -\frac{t}{L} x_1 + t\right) = \left(0,\; t\left(1 - \frac{x_1}{L}\right)\right) ]$$

下边界位移：

$$[ \mathbf{u}\big|_{x_2=0} = \left(0,\; \frac{t}{L} x_1 - t\right) = \left(0,\; t\left(\frac{x_1}{L} - 1\right)\right) ]$$

两组边界条件形成**反对称剪切（pure mode-II shear）**加载模式，使试件中部产生剪切裂缝扩展。

- **第 241–243 行**：对链段拉伸比 $[ \Lambda ]$ 在上下边界施加 Dirichlet 条件 $[ \Lambda = 1 ]$（未变形参考状态），防止损伤在边界处演化。

---

## 11. 运动学描述（第 248–257 行）

```python
248: d = len(u)
249: I = Identity(d)
250: F = I + grad(u)
251: C = F.T*F
253: J  = det(F)
254: Ic = tr(C) + 1
255: Lmbda_ch = sqrt(Ic/3)
257: Lmbda_ch = conditional(sqrt(Ic/3) > Lmbda_ch_max, Lmbda_ch_max, sqrt(Ic/3))
```

**解释：**

以下为有限变形连续介质力学的基本运动学量：

- **变形梯度（第 250 行）**：

$$[ \mathbf{F} = \mathbf{I} + \nabla\mathbf{u} ]$$

  其中 $[ \mathbf{I} ]$ 为二阶单位张量，$[ \nabla\mathbf{u} ]$ 为位移梯度。

- **右 Cauchy-Green 变形张量（第 251 行）**：

$$[ \mathbf{C} = \mathbf{F}^T \mathbf{F} ]$$

- **体积比（第 253 行）**：

$$[ J = \det(\mathbf{F}) ]$$

  $[ J > 1 ]$ 表示体积膨胀，$[ J < 1 ]$ 表示体积压缩，近不可压缩材料约束 $[ J \approx 1 ]$。

- **第一不变量（修正，考虑平面应变第三方向，第 254 行）**：

$$[ I_c = \mathrm{tr}(\mathbf{C}) + 1 ]$$

  加 1 是因为平面应变时 $[ \mathbf{C} ]$ 是 $[ 2\times2 ]$ 矩阵，需补充第三主值（$[ \lambda_3 = 1/J ]$ 近似为 1）。

- **链段拉伸比（第 255–257 行）**：基于非高斯网络模型（如 8-chain 模型），链段拉伸比由不变量计算：

$$[ \Lambda_{\mathrm{ch}} = \sqrt{\frac{I_c}{3}} ]$$

  同时进行上限截断处理，防止数值发散：

$$[ \Lambda_{\mathrm{ch}} = \min\!\left(\sqrt{\frac{I_c}{3}},\; \Lambda_{\mathrm{ch,max}}\right) ]$$

---

## 12. 本构模型与弹性能量泛函（第 260–290 行）

```python
262: def w(damage):
263:     return damage
265: def a(damage):
266:     return (1.0-damage)**2
268: def b_sq(damage):
269:     return (1.0-damage)**3
271: def P(u, damage):
272:     return a(damage)*mu*(F - inv(F.T)) - b_sq(damage)*p*J*inv(F.T)
274: def energy_density_function(u, damage):
275:     return a(damage)*(mu/2.0)*(Ic-3.0-2.0*ln(J)) - b_sq(damage)*p*(J-1.0) - (1./(2.*kappa))*(p**2)
278: elastic_energy = ((1-k_ell)*a(damage)+k_ell)*(mu/2.0)*(Ic-3.0-2.0*ln(J))*dx \
279:                  - b_sq(damage)*p*(J-1.0)*dx - 1./(2.*kappa)*p**2*dx
281: external_work  = dot(body_force, u)*dx
282: elastic_potential = elastic_energy - external_work
285: F_u = derivative(elastic_potential, w_p, v_q)
286: J_u = derivative(F_u, w_p, u_p)
288: problem_u = NonlinearVariationalProblem(F_u, w_p, bc_u, J=J_u)
289: solver_u  = NonlinearVariationalSolver(problem_u)
290: solver_u.parameters.update(solver_u_parameters)
```

**解释：**

**退化函数（第 262–269 行）：**

损伤模型引入三个依赖于损伤变量 $[ d \in [0,1] ]$ 的退化函数：

$$[ w(d) = d \quad \text{（能量耗散函数）} ]$$

$$[ a(d) = (1-d)^2 \quad \text{（弹性能量退化函数）} ]$$

$$[ b^2(d) = (1-d)^3 \quad \text{（压力-体积耦合退化函数）} ]$$

注：`b_sq` 实际上返回的是 $[ b^3(d) = (1-d)^3 ]$，代码注释中 `b squared` 与实际幂次有出入，实质上是 $[ (1-d)^3 ]$，对应混合有限元压力项退化。

**名义应力张量（第 271–272 行）：**

基于 Helmholtz 自由能对变形梯度的偏导（第一 Piola-Kirchhoff/名义应力）：

$$[ \mathbf{P}(\mathbf{u}, d) = a(d)\, \mu \left(\mathbf{F} - \mathbf{F}^{-T}\right) - b^2(d)\, p\, J\, \mathbf{F}^{-T} ]$$

其中：
- $[ a(d)\, \mu \left(\mathbf{F} - \mathbf{F}^{-T}\right) ]$ — 损伤退化后的 Neo-Hookean 偏量应力贡献；
- $[ -b^2(d)\, p\, J\, \mathbf{F}^{-T} ]$ — 静水压力（Lagrange 乘子 $[ p ]$）对名义应力的贡献。

**弹性应变能密度（第 274–275 行）：**

基于 Neo-Hookean 模型的单位参考体积弹性应变能密度（含混合格式稳定项）：

$$[ \psi(\mathbf{u}, d, p) = a(d)\, \frac{\mu}{2}\left(I_c - 3 - 2\ln J\right) - b^2(d)\, p\,(J-1) - \frac{p^2}{2\kappa} ]$$

- 第一项：等容 Neo-Hookean 弹性能（$[ I_c - 3 - 2\ln J ]$ 为 Rivlin-Saunders 型等容部分）；
- 第二项：不可压缩约束 $[ J = 1 ]$ 的弱形式（Lagrange 乘子 $[ p ]$ 作为压力）；
- 第三项：压力稳定项（Augmented Lagrangian），当 $[ \kappa \to \infty ]$ 时退化为严格不可压缩约束。

**全弹性势能（第 278–282 行）：**

$$[ \Pi(\mathbf{u}, p, d) = \int_\Omega \left[\left((1-k_\ell) a(d) + k_\ell\right) \frac{\mu}{2}(I_c - 3 - 2\ln J) - b^2(d)\, p\,(J-1) - \frac{p^2}{2\kappa}\right] \mathrm{d}V - \int_\Omega \mathbf{b} \cdot \mathbf{u}\, \mathrm{d}V ]$$

其中 $[ (1-k_\ell) a(d) + k_\ell ]$ 保证完全损伤时仍有残余刚度 $[ k_\ell ]$，避免刚度矩阵奇异。

**变分方程（第 285–290 行）：**

对混合场 $[ (\mathbf{u}, p) ]$ 取方向导数得到弱形式（一阶最优性条件）：

$$[ F_u = \frac{\mathrm{d}\Pi}{\mathrm{d}(\mathbf{u},p)} \cdot (\mathbf{v}, q) = 0 \quad \forall \, (\mathbf{v}, q) \in \mathbf{V}_h ]$$

$$[ J_u = \frac{\mathrm{d}^2 \Pi}{\mathrm{d}(\mathbf{u},p)^2} \cdot [(\delta\mathbf{u}, \delta p), (\mathbf{v}, q)] \quad \text{（切线刚度矩阵）} ]$$

用 `NonlinearVariationalProblem` 封装后交给 SNES 非线性求解器求解。

---

## 13. 损伤演化方程设置（第 293–319 行）

```python
294: Lmbda_0 = interpolate(Expression("1", degree=0), CG1)
295: Lmbda_previous = interpolate(Expression("1", degree=0), CG1)
297: def Relaxation(damage):
298:     return (1-damage)**g_power
300: def Viscous_Relaxation(Lmbda, Lmbda_previous):
301:     return conditional(Lmbda > Lmbda_previous, 1, 0)
303: dt = 1
305: Lmbda_WF = dt*(Lmbda*beta)*dx \
306:          + Relaxation(damage)*dt*(ell**2)*inner(grad(Lmbda), grad(beta))*dx \
307:          - dt*Lmbda_ch*beta*dx \
308:          + Viscous_Relaxation(Lmbda, Lmbda_previous)*eta*(Lmbda - Lmbda_previous)*beta*dx
309: Lmbda_Jacobian = derivative(Lmbda_WF, Lmbda, Lmbda_trial)
311: Lmbda_lb = interpolate(Expression("x[0]>=0 & x[0]<=L/5 & near(x[1], H/2, 0.01*hsize) ? 1.5 : 1", ...), CG1)
313: Lmbda_ub = interpolate(Expression("100", degree=0), CG1)
316: problem_Lmbda = NonlinearVariationalProblem(Lmbda_WF, Lmbda, bc_Lmbda, J=Lmbda_Jacobian)
317: problem_Lmbda.set_bounds(Lmbda_lb, Lmbda_ub)
318: solver_Lmbda = NonlinearVariationalSolver(problem_Lmbda)
319: solver_Lmbda.parameters.update(solver_Lmbda_parameters)
```

**解释：**

- **第 294–295 行**：初始化 $[ \Lambda_0 = 1 ]$（参考状态，无拉伸）和 $[ \Lambda_{\mathrm{prev}} = 1 ]$（上一时步值，用于粘性正则化和不可逆性）。

**松弛函数（第 297–298 行）：**

梯度惩罚项的空间松弛函数，使损伤区域的梯度正则化效果减弱（模拟损伤区域的软化）：

$$[ \mathcal{R}(d) = (1 - d)^{g_p} ]$$

  其中 $[ g_p = 0.28 ]$。当 $[ d \to 1 ]$ 时，$[ \mathcal{R}(d) \to 0 ]$，梯度项消失。

**粘性松弛（第 300–301 行）：**

$$[ \mathcal{V}(\Lambda, \Lambda_{\mathrm{prev}}) = \begin{cases} 1 & \Lambda > \Lambda_{\mathrm{prev}} \\ 0 & \Lambda \leq \Lambda_{\mathrm{prev}} \end{cases} ]$$

仅在链段拉伸增加时施加粘性阻尼（当 $[ \eta = 0 ]$ 时该项为零，本算例中粘性项不起作用）。

**损伤演化弱形式（第 305–308 行）：**

链段拉伸比 $[ \Lambda ]$ 的非局部梯度增强演化方程（弱形式）为：

$$[ \int_\Omega \Lambda\, \beta\, \mathrm{d}V + \mathcal{R}(d)\, \ell^2 \int_\Omega \nabla\Lambda \cdot \nabla\beta\, \mathrm{d}V - \int_\Omega \Lambda_{\mathrm{ch}}\, \beta\, \mathrm{d}V + \mathcal{V} \cdot \eta \int_\Omega (\Lambda - \Lambda_{\mathrm{prev}})\, \beta\, \mathrm{d}V = 0 ]$$

  对应的强形式偏微分方程（Euler-Lagrange 方程）为：

$$[ \Lambda - \mathcal{R}(d)\, \ell^2\, \Delta\Lambda = \Lambda_{\mathrm{ch}} \quad \text{in } \Omega ]$$

  这是一个**Helmholtz 型非局部方程**，将局部链段拉伸 $[ \Lambda_{\mathrm{ch}} ]$ 平滑为非局部场 $[ \Lambda ]$，正则化长度 $[ \ell ]$ 控制平滑半径。

- **第 311–313 行**：设置界约束（Box Constraint）：
  - 下界：在初始裂缝线（$[ x_1 \in [0, L/5],\; x_2 \approx H/2 ]$）处设 $[ \Lambda_{lb} = 1.5 ]$（预置损伤），其余区域 $[ \Lambda_{lb} = 1 ]$（无损伤初始状态）；
  - 上界：$[ \Lambda_{ub} = 100 ]$（实际上不起限制作用）。

- **第 316–319 行**：封装带界约束的损伤变分问题，使用支持界约束的 SNES 求解器 `vinewtonssls`。

---

## 14. 主求解循环（第 322–453 行）

```python
322: load_multipliers = np.linspace(load_min, load_max, load_steps)
326: energies   = np.zeros((len(load_multipliers), 5))
327: iterations = np.zeros((len(load_multipliers), 2))
330: (u, p) = w_p.split()
332: file_tot = XDMFFile(MPI.comm_world, savedir + "/results.xdmf")
334: file_tot.parameters["rewrite_function_mesh"] = False
335: file_tot.parameters["functions_share_mesh"]  = True
336: file_tot.parameters["flush_output"]          = True
338: J_integral_list = []
339: crack_length = []
340: traction_x_list = []
341: traction_y_list = []
```

**解释（第 322–341 行）：**

- **第 322 行**：生成均匀分布的 100 个荷载步 $[ t_k = k \cdot t_{\max} / (N-1),\; k=0,1,\ldots,99 ]$，从 0 到 0.6。
- **第 326–327 行**：初始化存储各时步能量和迭代次数的数组。
- **第 330 行**：将混合场 `w_p` 的视图分离为位移 `u` 和压力 `p`（用于后处理输出）。
- **第 332–336 行**：创建 XDMF 输出文件，配置为追加模式（`flush_output`），多场共享同一网格（节省存储）。
- **第 338–341 行**：初始化 J 积分列表、裂缝长度列表、牵引力列表。

### 14.1 荷载循环（第 343–381 行）

```python
343: for (i_t, t) in enumerate(load_multipliers):
351:     iteration = 1
352:     err_Lmbda = 1
355:     while err_Lmbda > AM_tolerance and iteration < maxiteration:
357:         solver_u.solve()
359:         solver_Lmbda.solve()
361:         Lmbda_numeric = project(Lmbda, CG1)
362:         damage_values = np.zeros_like(...)
363:         for i, value in enumerate(Lmbda_numeric.vector().get_local()):
364:             if value < lambda_i:
366:                 new_damage = 0
367:                 damage_values[i] = max(new_damage, damage_previous....)
368:             else:
369:                 new_damage = 1 - ((lambda_i-1)/(value-1)) * (1 - alpha_damage + alpha_damage * np.exp(-beta_damage * (value - lambda_i)))
370:                 damage_values[i] = max(new_damage, damage_previous....)
372:         damage.vector().set_local(damage_values)
374:         Lmbda_error = Lmbda.vector() - Lmbda_0.vector()
375:         err_Lmbda = Lmbda_error.norm('linf')
377:         Lmbda_0.assign(Lmbda)
380:         volume_ratio = assemble(J/(L*H)*dx)
381:         iteration = iteration + 1
```

**解释：**

**交替最小化算法（Alternate Minimization，AM）：**

对每个荷载步 $[ t_k ]$，固定损伤场求解位移，再固定位移场求解损伤，交替迭代直到收敛。这是相场断裂力学中的标准求解策略：

$$[ \begin{cases} \text{Step 1: } \min_{(\mathbf{u},p)} \Pi(\mathbf{u}, p;\, d^{(n)}) & \Rightarrow \text{更新} (\mathbf{u},p) \\ \text{Step 2: } \min_{\Lambda} \mathcal{F}(\Lambda;\, \mathbf{u}^{(n+1)}) & \Rightarrow \text{更新} \Lambda \\ \end{cases} ]$$

- **第 357 行**：固定 $[ \Lambda ]$（即固定 $[ d ]$）求解位移问题（非线性弹性）。
- **第 359 行**：固定 $[ \mathbf{u} ]$ 求解链段拉伸比 $[ \Lambda ]$（带界约束的 Helmholtz 问题）。
- **第 361–372 行**：由计算所得的 $[ \Lambda ]$ 值逐节点更新损伤场 $[ d ]$：

  当 $[ \Lambda < \lambda_i ]$ 时（无损伤阶段）：

  $$[ d = 0 ]$$

  当 $[ \Lambda \geq \lambda_i ]$ 时（损伤演化阶段）：

  $$[ d(\Lambda) = 1 - \frac{\lambda_i - 1}{\Lambda - 1}\left[(1 - \alpha_d) + \alpha_d \exp\!\left(-\beta_d(\Lambda - \lambda_i)\right)\right] ]$$

  同时强制损伤不可逆性（第 367、370 行）：

  $$[ d^{(k)} = \max\!\left(d_{\mathrm{new}},\; d^{(k-1)}\right) ]$$

- **第 374–376 行**：用 $[ L^\infty ]$ 范数衡量 $[ \Lambda ]$ 的收敛误差：

  $$[ e_\Lambda = \|\Lambda^{(n+1)} - \Lambda^{(n)}\|_{L^\infty} = \max_{\mathbf{x} \in \Omega}\left|\Lambda^{(n+1)}(\mathbf{x}) - \Lambda^{(n)}(\mathbf{x})\right| ]$$

  当 $[ e_\Lambda < \varepsilon_{\mathrm{AM}} = 2\times10^{-3} ]$ 时认为该荷载步的 AM 迭代收敛。

- **第 380 行**：计算当前时步的平均体积比（监控近不可压缩约束满足情况）：

  $$[ \bar{J} = \frac{1}{|\Omega|}\int_\Omega J\, \mathrm{d}V ]$$

### 14.2 荷载步后处理（第 384–453 行）

```python
385: Lmbda_previous.assign(Lmbda)
386: damage_previous.assign(damage)
389: local_project(P(u, damage), T_DG0, PTensor)
390: local_project(F, T_DG0, FTensor)
391: local_project(J, CG1, JScalar)
```

**解释（第 384–391 行）：**

- **第 385–386 行**：将收敛后的 $[ \Lambda ]$ 和 $[ d ]$ 保存为"上一时步"值，用于下一荷载步的不可逆性约束和粘性正则化。
- **第 389–391 行**：将名义应力 $[ \mathbf{P} ]$、变形梯度 $[ \mathbf{F} ]$、体积比 $[ J ]$ 投影到 DG0/CG1 空间，便于后处理和可视化。

```python
394: u.rename("Displacement", "u")
395: p.rename("Pressure", "p")
396: Lmbda.rename("Stretch", "Lmbda")
397: damage.rename("Damage", "damage")
400: Lmbda_ch_projected = project(Lmbda_ch, CG1)
401: Lmbda_ch_projected.rename("Lambda_ch", "Lmbda_ch_projected")
402: file_tot.write(u, t)
...
409: file_tot.write(JScalar,t)
412: Vct_space = VectorFunctionSpace(mesh, "Lagrange", 1)
413: laplacian_Lmbda = grad(Lmbda)
414: grad_term = (ell**2) * laplacian_Lmbda
415: grad_term_proj = project(grad_term, Vct_space)
416: grad_term_proj.rename("Grad_Lambda_Term", "grad_term_proj")
417: file_tot.write(grad_term_proj,t)
420: u1.t = t
421: u2.t = t
```

**解释（第 394–421 行）：**

- 对各场变量重命名并写入 XDMF 文件（用于 ParaView 可视化）。
- **第 400 行**：将局部链段拉伸 $[ \Lambda_{\mathrm{ch}} = \min(\sqrt{I_c/3}, \Lambda_{\mathrm{ch,max}}) ]$ 投影到 CG1 空间输出，便于与非局部 $[ \Lambda ]$ 比较。
- **第 413–417 行**：计算并输出非局部梯度项 $[ \ell^2 \nabla\Lambda ]$，便于分析正则化效果：

$$[ \boldsymbol{g} = \ell^2 \nabla\Lambda ]$$

- **第 420–421 行**：更新 Dirichlet 边界条件表达式中的时间参数 $[ t ]$，使边界位移在下一荷载步正确取值。

**裂缝长度计算（第 425–436 行）：**

```python
425: damage_values = damage.compute_vertex_values(mesh)
429: for vertex in vertices(mesh):
430:     if damage_values[vertex.index()] >= 0.95:
431:         if vertex.point().x() > rightmost_x:
432:             rightmost_vertex = vertex
433:             rightmost_x = vertex.point().x()
435: crack_length.append(rightmost_x)
```

**解释：**

逐顶点扫描损伤场，寻找损伤值 $[ d \geq 0.95 ]$ 的最右侧节点的 $[ x ]$ 坐标，作为当前时步的**裂缝尖端位置**（即裂缝长度）：

$$[ a(t) = \max\!\left\{ x_1(\mathbf{x}) : d(\mathbf{x}) \geq 0.95 \right\} ]$$

**J 积分计算（第 439–445 行）：**

```python
439: F_1 = F[0, 0]; F_2 = F[1, 0]
440: F_1_vector = as_vector([F_1, F_2])
442: J_expression = -1*(energy_density_function(u, damage)*grad(q_J_integral)[0] \
443:                     - inner(P(u, damage), outer(F_1_vector, grad(q_J_integral))))
444: J_integral = assemble(J_expression*dx(1))
```

**解释：**

采用**域积分形式**计算裂缝驱动力 J 积分（等效于能量释放率 $[ G ]$）。在有限变形框架下，J 积分的域积分表达式为：

$$[ J = -\int_{\Omega_J} \left[ \psi\, \frac{\partial q}{\partial X_1} - \mathbf{P} : \left(\mathbf{F}_{,1} \otimes \nabla q\right) \right] \mathrm{d}V ]$$

其中：
- $[ \psi ]$ — 弹性应变能密度；
- $[ \mathbf{F}_{,1} = \partial\mathbf{F}/\partial X_1 ]$ — 变形梯度对参考坐标 $[ X_1 ]$ 的梯度（代码中用第一列 $[ (F_{11}, F_{21})^T ]$ 表示）；
- $[ q(\mathbf{x}) ]$ — 裂尖权函数（见第 9 节）；
- $[ \Omega_J ]$ — J 积分环形积分域（子域 1，由 `dx(1)` 指定）；
- $[ \mathbf{P} : (\mathbf{F}_{,1} \otimes \nabla q) = P_{iJ} F_{i1,\cdot} \nabla q_J ]$ — 双点缩并（`inner`）。

**牵引力计算（第 449–453 行）：**

```python
449: P_12 = P(u, damage)[0, 1]
450: P_22 = P(u, damage)[1, 1]
451: traction_y = assemble(P_22*ds(1))
452: traction_y_list.append(traction_y)
```

**解释：**

计算上边界（标记 1）的法向牵引力（总反力），用于绘制载荷-位移曲线：

$$[ F_y = \int_{\Gamma_{\mathrm{top}}} P_{22}\, \mathrm{d}S ]$$

其中 $[ P_{22} ]$ 为名义应力张量 $[ \mathbf{P} ]$ 的 (2,2) 分量（法向应力）。

---

## 15. 后处理与结果输出（第 456–488 行）

```python
457: crack_length_arr = np.array(crack_length)
458: J_integral_arr   = np.array(J_integral_list)
459: traction_y_arr   = np.array(traction_y_list)
460: np.savetxt(savedir + '/J_disp.txt',
461:            np.column_stack((crack_length_arr, J_integral_arr)),
462:            header='Crack Length | J Integral List', fmt='%f', delimiter=' | ')
463: np.savetxt(savedir + '/traction_disp.txt',
464:            np.column_stack((load_multipliers, traction_y_arr)),
465:            header='displacament | traction', fmt='%f', delimiter=' | ')
```

**解释（第 457–463 行）：**

将全过程计算结果保存为文本文件：
- `J_disp.txt`：裂缝长度 $[ a(t_k) ]$ 与 J 积分 $[ J(t_k) ]$ 的对应关系，反映裂缝扩展驱动力随裂缝长度的变化；
- `traction_disp.txt`：施加位移 $[ t_k ]$ 与顶部反力 $[ F_y(t_k) ]$ 的对应关系（载荷-位移曲线）。

```python
466: num_plot = load_steps
467: plt.figure(1)
468: plt.plot(crack_length[1:num_plot], J_integral_list[1:num_plot], label='total')
469: plt.xlabel('Crack length')
470: plt.ylabel('J')
471: plt.title('J Integral')
472: plt.legend()
473: plt.savefig(savedir + '/J_Integral.pdf', transparent=True)
474: plt.show()
477: plt.figure(2)
478: plt.plot(load_multipliers[1:num_plot], traction_y_list[1:num_plot])
479: plt.xlabel('Displacement')
480: plt.ylabel('Total force')
481: plt.title('Traction')
482: plt.savefig(savedir + '/traction_disp.pdf', transparent=True)
```

**解释（第 466–482 行）：**

绘制并保存两张曲线图：

1. **J 积分曲线（第 467–474 行）**：横轴为裂缝长度 $[ a ]$，纵轴为能量释放率 $[ J ]$，体现断裂韧性曲线（R-曲线）特征：

$$[ J = J(a) ]$$

2. **载荷-位移曲线（第 477–482 行）**：横轴为边界位移 $[ t ]$，纵轴为总反力 $[ F_y ]$，体现材料的宏观力学响应：

$$[ F_y = F_y(t) = \int_{\Gamma_{\mathrm{top}}} P_{22}\, \mathrm{d}S ]$$

两图均跳过第 0 步（`[1:num_plot]`），避免零荷载步的无意义数据影响图形。结果以透明背景 PDF 格式保存，便于嵌入论文。

---

## 附：主要符号与公式汇总

| 符号 | 含义 | 定义 |
|------|------|------|
| $[ \mathbf{u} ]$ | 位移场 | 未知量 |
| $[ p ]$ | 静水压力（Lagrange 乘子） | 未知量 |
| $[ \mathbf{F} ]$ | 变形梯度 | $[ \mathbf{F} = \mathbf{I} + \nabla\mathbf{u} ]$ |
| $[ J ]$ | 体积比 | $[ J = \det\mathbf{F} ]$ |
| $[ \mathbf{C} ]$ | 右 Cauchy-Green 张量 | $[ \mathbf{C} = \mathbf{F}^T\mathbf{F} ]$ |
| $[ I_c ]$ | 第一不变量（修正） | $[ I_c = \mathrm{tr}(\mathbf{C}) + 1 ]$ |
| $[ \Lambda_{\mathrm{ch}} ]$ | 局部链段拉伸比 | $[ \Lambda_{\mathrm{ch}} = \sqrt{I_c/3} ]$ |
| $[ \Lambda ]$ | 非局部链段拉伸比 | 由 Helmholtz 方程求解 |
| $[ d ]$ | 损伤变量 | $[ d \in [0,1] ]$ |
| $[ a(d) ]$ | 弹性退化函数 | $[ a(d) = (1-d)^2 ]$ |
| $[ b^2(d) ]$ | 压力退化函数 | $[ b^2(d) = (1-d)^3 ]$ |
| $[ \mathbf{P} ]$ | 名义应力张量 | $[ \mathbf{P} = a(d)\mu(\mathbf{F}-\mathbf{F}^{-T}) - b^2(d)\,pJ\mathbf{F}^{-T} ]$ |
| $[ \psi ]$ | 弹性应变能密度 | $[ \psi = a(d)\frac{\mu}{2}(I_c-3-2\ln J) - b^2(d)\,p(J-1) - \frac{p^2}{2\kappa} ]$ |
| $[ \ell ]$ | 正则化长度参数 | $[ \ell = 2 h_{\mathrm{size}} ]$ |
| $[ \mathcal{R}(d) ]$ | 松弛函数 | $[ \mathcal{R}(d) = (1-d)^{g_p} ]$ |

---

*文档作者：根据 Mohammad Mousavi 的原始代码（Cornell University, 06/01/2024）整理。*
