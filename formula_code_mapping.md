# 论文公式与实现代码对应关系

本文档将论文 *A chain stretch-based gradient-enhanced model for damage and fracture in elastomers*（aba0c8.md）中的每条公式与 `Improved_GED.py` 中的对应代码一一对应，公式在前，代码在后，并附有简要说明。

---

## 1. 运动学（Kinematics）

### 公式 (1)：变形映射

\[
d\mathbf{x} = \mathbf{F}(\mathbf{X},t)\,d\mathbf{X}
\]

变形梯度 \(\mathbf{F}\) 将参考构型中的线元 \(d\mathbf{X}\) 映射到当前构型中的线元 \(d\mathbf{x}\)。

### 公式 (2)：变形梯度定义

\[
\mathbf{F}(\mathbf{X},t) = \frac{\partial\mathbf{x}}{\partial\mathbf{X}}
\]

利用位移 \(\mathbf{u} = \mathbf{x} - \mathbf{X}\) 可得 \(\mathbf{F} = \mathbf{I} + \nabla\mathbf{u}\)。

**对应代码**（第 249–251 行）：

```python
d = len(u)
I = Identity(d)          # 单位张量
F = I + grad(u)          # 变形梯度 F = I + ∇u
```

---

### 公式 (3)：右 Cauchy-Green 变形张量

\[
\mathbf{C} = \mathbf{F}^T\mathbf{F}
\]

**对应代码**（第 252 行）：

```python
C = F.T*F                # 右 Cauchy-Green 张量
```

---

### 公式 (4)：不变量

\[
\begin{array}{c}
I_1 = \mathrm{tr}(\mathbf{C}),\\[4pt]
I_2 = \dfrac{1}{2}\!\left((\mathrm{tr}(\mathbf{C}))^2 - \mathrm{tr}(\mathbf{C}^2)\right),\\[4pt]
I_3 = \det(\mathbf{C}).
\end{array}
\]

雅可比行列式 \(J = \det\mathbf{F} > 0\) 描述体积比。代码采用平面应变假设, 令 \(I_1 = \mathrm{tr}(\mathbf{C}) + 1\)（附加平面外分量），即代码中的 `Ic`。

**对应代码**（第 253–254 行）：

```python
J  = det(F)              # 体积比 J = det(F)
Ic = tr(C) + 1           # 平面应变下的修正 I₁ = tr(C) + 1
```

---

### 公式 (5)：链拉伸（Arruda-Boyce 8 链模型）

\[
\lambda_{ch} = \sqrt{\frac{I_1}{3}}
\]

链拉伸 \(\lambda_{ch}\) 描述弹性体材料点处聚合物链的平均拉伸程度。代码中另加了 \(\lambda_{ch}^{\max}\) 上界截断：

\[
\lambda_{ch} = \min\!\left(\sqrt{\frac{I_1}{3}},\ \lambda_{ch}^{\max}\right)
\]

**对应代码**（第 256–257 行）：

```python
Lmbda_ch = sqrt(Ic/3)
Lmbda_ch = conditional(sqrt(Ic/3) > Lmbda_ch_max, Lmbda_ch_max, sqrt(Ic/3))
```

参数 `Lmbda_ch_max = 2.7`（第 96 行）对应论文中的 \(\lambda_{ch}^{\max}\)，用于防止局部链拉伸无界增大从而避免损伤区拓宽（见论文公式 (36)）。

---

## 2. 自由能（Helmholtz Free Energy）

### 公式 (21)：Neo-Hookean 自由能密度

\[
\psi(\mathbf{F}) = \frac{\mu}{2}\!\left(I_1(\mathbf{C}) - 3 - 2\ln J\right)
\]

其中 \(\mu\) 为剪切模量。

### 公式 (19)：近不可压缩约束

\[
\hat{\psi}(\mathbf{F},p) = \psi(\mathbf{F}) - p(J-1) - \frac{p^2}{2\kappa}
\]

其中 \(\kappa\) 为体积模量，\(p\) 为压力场（Lagrange 乘子）。

### 公式 (20)：含损伤的局部自由能密度

\[
\Psi^{\mathrm{loc}}(\mathbf{F},p,d) = a(d)\,\psi(\mathbf{F}) - b(d)\,p(J-1) - \frac{p^2}{2\kappa}
\]

**对应代码**（第 274–279 行）：

```python
def energy_density_function(u, damage):
    return a(damage)*(mu/2.0)*(Ic-3.0-2.0*ln(J)) \
           - b_sq(damage)*p*(J-1.0) \
           - (1./(2.*kappa))*(p**2)

# 弹性势能积分（含 k_ell 数值稳定项）
elastic_energy = ((1-k_ell)*a(damage)+k_ell)*(mu/2.0)*(Ic-3.0-2.0*ln(J))*dx \
                 - b_sq(damage)*p*(J-1.0)*dx \
                 - 1./(2.*kappa)*p**2*dx
```

`elastic_energy` 中第一项的系数 `((1-k_ell)*a(damage)+k_ell)` 与论文公式 (24) 的完整降解函数完全一致；`energy_density_function` 用于后处理（J 积分计算），内部 `a(damage)` 为简化形式 \((1-d)^2\)。

---

## 3. 降解函数（Degradation Functions）

### 公式 (24)：降解函数 \(a(d)\) 和 \(b(d)\)

\[
a(d) = (1-k_\ell)(1-d)^2 + k_\ell,\qquad b(d) = (1-k_\ell)(1-d)^3 + k_\ell
\]

其中 \(k_\ell\) 为数值调节参数，保证全损伤状态下有残余刚度。\(b(d)\) 比 \(a(d)\) 有更高次多项式，使体积模量随损伤比剪切模量下降更快，从而允许裂纹张开。

**对应代码**（第 265–272 行）：

```python
def a(damage):           # 剪切降解：≈ (1-d)²（含 k_ell 时见 elastic_energy）
    return (1.0-damage)**2

def b_sq(damage):        # 体积降解：≈ (1-d)³
    return (1.0-damage)**3

def P(u, damage):        # 名义应力张量
    return a(damage)*mu*(F - inv(F.T)) - b_sq(damage)*p*J*inv(F.T)
```

参数 `k_ell = 1e-5` （第 103 行）。函数名 `b_sq` 源自论文注释（"b squared"），实际计算时直接使用 \((1-d)^3\) 而非其平方根。

---

## 4. 松弛函数（Relaxation Function）

### 公式 (25)：松弛函数

\[
g(d) = (1-d)^m,\quad m > 0
\]

该函数随损伤积累逐渐削弱非局部相互作用，是防止损伤区拓宽的两个关键建模特征之一（另一个是 \(\lambda_{ch}^{\max}\) 上界）。

**对应代码**（第 297–298 行）：

```python
def Relaxation(damage):
    return (1-damage)**g_power
```

参数 `g_power = 0.28` （第 132 行）对应论文中的 \(m = 0.28\)。

---

## 5. 本构关系（Constitutive Relations）

### 公式 (28)：名义应力与微力

由自由能 (27) 对状态变量求偏导得：

\[
\begin{aligned}
\mathbf{P} &= a(d)\,\frac{\partial\psi}{\partial\mathbf{F}} - b(d)\,p\,J\,\mathbf{F}^{-T} & &\mathrm{in}\ \Omega_0,\\
f_p &= -b(d)(J-1) - \frac{p}{\kappa} & &\mathrm{in}\ \Omega_0.
\end{aligned}
\]

对于 Neo-Hookean 材料，\(\partial\psi/\partial\mathbf{F} = \mu(\mathbf{F} - \mathbf{F}^{-T})\)，代入得：

\[
\mathbf{P} = a(d)\,\mu\!\left(\mathbf{F} - \mathbf{F}^{-T}\right) - b(d)\,p\,J\,\mathbf{F}^{-T}
\]

**对应代码**（第 271–272 行）：

```python
def P(u, damage):        # 名义应力（第一 Piola-Kirchhoff 应力）
    return a(damage)*mu*(F - inv(F.T)) - b_sq(damage)*p*J*inv(F.T)
```

压力的本构方程 \(f_p = 0\)（公式 (10)）通过混合弱形式隐式求解（见公式 (40)₂）。

---

## 6. 损伤函数（Damage Function）

### 公式 (31)：损伤演化

\[
d =
\begin{cases}
0 & \text{若}\ \tilde{\lambda} < \lambda_{cr}, \\[4pt]
1 - \dfrac{\lambda_{cr}-1}{\mathcal{H}(\tilde{\lambda})-1}
\!\left(1 - c + c\,e^{-\gamma\!\left(\mathcal{H}(\tilde{\lambda})-\lambda_{cr}\right)}\right)
& \text{若}\ \tilde{\lambda} \geq \lambda_{cr}.
\end{cases}
\]

其中：
- \(\lambda_{cr}\) 为损伤起始临界拉伸；
- \(c \in [0,1]\) 控制最大损伤量（取 1 时允许完全损伤）；
- \(\gamma \geq 0\) 控制从无损到完全损伤的过渡陡度；
- \(\mathcal{H}(\tilde{\lambda};\mathbf{X},t) = \max_{s\in[0,t]}\tilde{\lambda}(\mathbf{X},s)\) 为历史最大值函数，保证损伤不可逆。

**对应代码**（第 362–372 行）：

```python
Lmbda_numeric = project(Lmbda, CG1)
damage_values = np.zeros_like(Lmbda_numeric.vector().get_local())
for i, value in enumerate(Lmbda_numeric.vector().get_local()):
    if value < lambda_i:
        new_damage = 0
        damage_values[i] = max(new_damage, damage_previous.vector().get_local()[i])
    else:
        new_damage = 1 - ((lambda_i-1)/(value-1)) * \
                     (1 - alpha_damage + alpha_damage * np.exp(-beta_damage * (value - lambda_i)))
        damage_values[i] = max(new_damage, damage_previous.vector().get_local()[i])
damage.vector().set_local(damage_values)
```

参数对应关系：

| 论文符号 | 代码变量 | 取值 |
|---------|---------|------|
| \(\lambda_{cr}\) | `lambda_i` | 1.2 |
| \(c\) | `alpha_damage` | 1 |
| \(\gamma\) | `beta_damage` | 20 |
| \(\mathcal{H}(\tilde{\lambda})\) | `value`（历史最大，由 SNES 单调约束保证） | — |

`max(new_damage, damage_previous...)` 强制执行损伤不可逆条件（对应历史函数 \(\mathcal{H}\)）。

---

## 7. 强形式与弱形式（Strong and Weak Forms）

### 公式 (36)：含上界约束的非局部链拉伸方程（强形式）

\[
\tilde{\lambda} - \lambda_{ch}^{\max} + \left\langle\!\left\langle\lambda_{ch}^{\max} - \lambda_{ch}\right\rangle\!\right\rangle - \nabla\cdot\!\left(g(d)\,\ell^2\,\nabla\tilde{\lambda}\right) = 0\quad\mathrm{in}\ \Omega_0.
\]

Macaulay 括号定义为 \(\langle\!\langle x\rangle\!\rangle = \max(x,0)\)，使得当 \(\lambda_{ch} < \lambda_{ch}^{\max}\) 时驱动力为 \(\lambda_{ch}\)，当 \(\lambda_{ch} \geq \lambda_{ch}^{\max}\) 时驱动力被限制为 \(\lambda_{ch}^{\max}\)。代码通过 `conditional` 在计算 `Lmbda_ch` 时直接截断（见第 5.1 节公式 (5) 对应代码），等效实现此操作。

### 公式 (40)：弱形式（含人工黏度，见附录 C 公式 C.3）

\[
\begin{aligned}
&\int_{\Omega_0}\!\!\left[a(d)\frac{\partial\psi}{\partial\mathbf{F}} - b(d)pJ\mathbf{F}^{-T}\right]:\nabla\boldsymbol{v}\,dV
 - \int_{\partial_N\Omega_0}\!\!\mathbf{T}\cdot\boldsymbol{v}\,dA = 0,\\[4pt]
&-\int_{\Omega_0}\!\!\left(b(d)(J-1)+\frac{p}{\kappa}\right)q\,dV = 0,\\[4pt]
&\int_{\Omega}\!\!\eta\,\dot{\tilde{\lambda}}\,\beta\,dV
+\int_{\Omega}\!\!\tilde{\lambda}\,\beta\,dV
-\int_{\Omega}\!\!\left[\lambda_{ch}^{\max} - \left\langle\!\left\langle\lambda_{ch}^{\max}-\lambda_{ch}\right\rangle\!\right\rangle\right]\beta\,dV
+\int_{\Omega}\!\!g(d)\,\ell^2\,\nabla\tilde{\lambda}\cdot\nabla\beta\,dV = 0.
\end{aligned}
\]

#### 弱形式 (40)₁₊₂：位移与压力问题

**对应代码**（第 278–290 行）：

```python
elastic_energy    = ((1-k_ell)*a(damage)+k_ell)*(mu/2.0)*(Ic-3.0-2.0*ln(J))*dx \
                    - b_sq(damage)*p*(J-1.0)*dx - 1./(2.*kappa)*p**2*dx
external_work     = dot(body_force, u)*dx
elastic_potential = elastic_energy - external_work

# 对 w_p 求方向导数（梯度）= 弱形式残差
F_u = derivative(elastic_potential, w_p, v_q)
J_u = derivative(F_u, w_p, u_p)
problem_u = NonlinearVariationalProblem(F_u, w_p, bc_u, J=J_u)
solver_u  = NonlinearVariationalSolver(problem_u)
solver_u.parameters.update(solver_u_parameters)
```

`derivative(elastic_potential, w_p, v_q)` 自动计算关于试探函数 \((\boldsymbol{u},p)\) 的方向导数，即弱形式 (40)₁₊₂ 的残差。

#### 弱形式 (40)₃：非局部链拉伸问题（含人工黏度，附录 C）

\[
\int_{\Omega}\eta\,\dot{\tilde{\lambda}}\,\beta\,dV
+ \int_{\Omega}\tilde{\lambda}\,\beta\,dV
- \int_{\Omega}\lambda_{ch}\,\beta\,dV
+ \int_{\Omega}g(d)\,\ell^2\,\nabla\tilde{\lambda}\cdot\nabla\beta\,dV = 0
\]

**对应代码**（第 303–319 行）：

```python
dt = 1

Lmbda_WF = dt*(Lmbda*beta)*dx \
         + Relaxation(damage)*dt*(ell**2)*inner(grad(Lmbda), grad(beta))*dx \
         - dt*Lmbda_ch*beta*dx \
         + Viscous_Relaxation(Lmbda, Lmbda_previous)*eta*(Lmbda - Lmbda_previous)*beta*dx

Lmbda_Jacobian = derivative(Lmbda_WF, Lmbda, Lmbda_trial)

# 下界约束：裂缝区域 λ̃ ≥ 1.5，其他区域 λ̃ ≥ 1
Lmbda_lb = interpolate(Expression(
    "x[0]>=0 & x[0]<=L/5 & near(x[1], H/2, 0.01*hsize) ? 1.5 : 1",
    hsize=hsize, L=L, H=H, degree=0), CG1)
Lmbda_ub = interpolate(Expression("100", degree=0), CG1)

problem_Lmbda = NonlinearVariationalProblem(Lmbda_WF, Lmbda, bc_Lmbda, J=Lmbda_Jacobian)
problem_Lmbda.set_bounds(Lmbda_lb, Lmbda_ub)
solver_Lmbda = NonlinearVariationalSolver(problem_Lmbda)
solver_Lmbda.parameters.update(solver_Lmbda_parameters)
```

各项对应关系：

| 弱形式项 | 代码表达式 |
|---------|-----------|
| \(\int\tilde{\lambda}\,\beta\,dV\) | `dt*(Lmbda*beta)*dx` |
| \(\int g(d)\ell^2\nabla\tilde{\lambda}\cdot\nabla\beta\,dV\) | `Relaxation(damage)*dt*(ell**2)*inner(grad(Lmbda), grad(beta))*dx` |
| \(-\int\lambda_{ch}\,\beta\,dV\) | `-dt*Lmbda_ch*beta*dx`（`Lmbda_ch` 已含截断） |
| \(\int\eta\,\dot{\tilde{\lambda}}\,\beta\,dV\) | `Viscous_Relaxation(...)*eta*(Lmbda - Lmbda_previous)*beta*dx` |

其中 `Viscous_Relaxation` 函数（第 300–301 行）为条件性人工黏度：

```python
def Viscous_Relaxation(Lmbda, Lmbda_previous):
    return conditional(Lmbda > Lmbda_previous, 1, 0)
```

仅在 \(\tilde{\lambda}\) 递增时施加黏度, 确保黏度项与损伤不可逆性一致。当 `eta = 0` （第 104 行）时, 该项不起作用 （对应无人工黏度的主要算例）。

---

## 8. 边界条件（Boundary Conditions）

### 公式 (34)/(35)：边界条件

位移 Dirichlet 条件（三角形位移，Section 4 第一边值问题）：

\[
\bar{u}_2^{\mathrm{top}} = t\!\left(1 - \frac{X_1}{L}\right),\qquad \bar{u}_2^{\mathrm{bot}} = -t\!\left(1 - \frac{X_1}{L}\right)
\]

**对应代码**（第 234–243 行）：

```python
# 施加位移边界条件（顶面向内，底面向外）
u1 = Expression([0, "(-t/L)*x[0] + t"], t=0.0, L=L, degree=1)
u2 = Expression([0, "(t/L)*x[0] - t"], t=0.0, L=L, degree=1)
bc_u1 = DirichletBC(V.sub(0), u1, top_boundary)
bc_u2 = DirichletBC(V.sub(0), u2, bot_boundary)
bc_u = [bc_u1, bc_u2]

# 非局部链拉伸边界条件：顶底边界处 λ̃ = 1（无损伤）
bc_Lmbda_B = DirichletBC(CG1, 1, bot_boundary)
bc_Lmbda_T = DirichletBC(CG1, 1, top_boundary)
bc_Lmbda = [bc_Lmbda_B, bc_Lmbda_T]
```

\(t\) 为载荷乘子 （在 `load_multipliers` 中从 0 线性增大到 0.6）。顶面 \(X_1\) 方向固定 （表达式 x 分量为 0）, \(X_2\) 分量按三角形施加。\(\tilde{\lambda} = 1\) 的 Dirichlet 条件 （公式 (35)₃）对应无损伤边界。

---

## 9. 交替最小化方案（Staggered Solution Scheme）

论文第 3 节描述了交替最小化迭代格式，外层循环检验收敛准则 \(|\tilde{\lambda}^j - \tilde{\lambda}^{j-1}|_\infty < 2\times10^{-3}\)。

**对应代码**（第 350–381 行）：

```python
iteration = 1
err_Lmbda = 1

while err_Lmbda > AM_tolerance and iteration < maxiteration:
    # 第一内层循环：固定 λ̃，求解 u 和 p
    solver_u.solve()
    # 第二内层循环：固定 u 和 p，求解 λ̃（含约束）
    solver_Lmbda.solve()

    # 根据 λ̃ 更新损伤场（公式 (31)）
    Lmbda_numeric = project(Lmbda, CG1)
    ...  # 损伤计算见第 6 节

    # 外层收敛检验：∞-范数误差
    Lmbda_error = Lmbda.vector() - Lmbda_0.vector()
    err_Lmbda = Lmbda_error.norm('linf')
    Lmbda_0.assign(Lmbda)
    iteration = iteration + 1

# 更新历史变量（不可逆性）
Lmbda_previous.assign(Lmbda)
damage_previous.assign(damage)
```

参数 `AM_tolerance = 2e-3` （第 121 行）, `maxiteration = 100` （第 120 行）。

---

## 10. J 积分（J Integral）

### 公式 (A.1)：路径 J 积分

\[
J = -\frac{d\Psi}{d\tilde{a}} = \int_{S}\!\!\left(\Psi\,N_1 - P_{IJ}\,N_J\,\frac{\partial x_i}{\partial X_1}\right)dS
\]

### 公式 (A.2)：域 J 积分

\[
J = \int_{A_1}\!\!\left(-\Psi\,\frac{\partial q}{\partial X_1} + P_{IJ}\,F_{i1}\,\frac{\partial q}{\partial X_J}\right)dA
\]

其中 \(A_1\) 为裂纹尖端周围的环形区域 （内外半径分别为 \(r_1\) 和 \(r_2\)）, 权函数 \(q\) 在内边界 \(C_1\) 处为 1、在外边界 \(C_3\) 处为 0。

**对应代码**（第 439–444 行）：

```python
# 取变形梯度第一列（对应 F_{i1}，i=1,2）
F_1 = F[0, 0]
F_2 = F[1, 0]
F_1_vector = as_vector([F_1, F_2])

# 域 J 积分：∫(-Ψ ∂q/∂X₁ + P:outer(F₁, ∇q)) dA
J_expression = -1*(energy_density_function(u, damage)*grad(q_J_integral)[0] \
                   - inner(P(u, damage), outer(F_1_vector, grad(q_J_integral))))
J_integral = assemble(J_expression*dx(1))   # 仅在圆环区域 dx(1) 积分
```

权函数 `q_J_integral` （第 224–228 行）为环形区域上的平滑函数: 内圆半径 `r1=0.45`, 外圆半径 `r2=0.47`, 在 \(r < r_1\) 处取 1, 在 \(r > r_2\) 处取 0, 中间线性插值。

---

## 11. 参数汇总

| 论文符号 | 代码变量 | 取值 | 含义 |
|---------|---------|------|------|
| \(\mu\) | `mu` | 1 | 剪切模量 |
| \(\kappa\) | `kappa` | 1000 | 体积模量 |
| \(k_\ell\) | `k_ell` | 1e-5 | 数值残余刚度 |
| \(\ell\) | `ell` | \(2h = 0.04\) | 非局部长度参数 |
| \(\lambda_{cr}\) | `lambda_i` | 1.2 | 损伤起始临界拉伸 |
| \(\lambda_{ch}^{\max}\) | `Lmbda_ch_max` | 2.7 | 链拉伸驱动力上界 |
| \(c\) | `alpha_damage` | 1 | 最大损伤系数 |
| \(\gamma\) | `beta_damage` | 20 | 损伤函数陡度参数 |
| \(m\) | `g_power` | 0.28 | 松弛函数幂次 |
| \(\eta\) | `eta` | 0 | 人工黏度系数（主算例关闭）|

---

*本文档对应代码版本: `Improved_GED.py`, 论文: aba0c8.md （A chain stretch-based gradient-enhanced model for damage and fracture in elastomers, Mousavi et al.）*
