### A Pluto.jl notebook ###
# v0.16.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 5867632c-fff5-11eb-3a19-2f309efd424a
begin
	using LinearAlgebra
	using PlutoUI, LaTeXStrings
	using HypertextLiteral
	using Plots,  PlotThemes
	using SymPy
	using FileIO
	using Images # , ImageIO, ImageMagick, ImageShow 
	# using CalculusWithJulia
end

# ╔═╡ 5928498a-cee5-4e32-a2c6-d46d3b9dcca4
TableOfContents(depth=4)

# ╔═╡ 77f286eb-1757-4044-86dd-6a550affea75
md"""
# MATH333: Methods of Applied Mathematics 
## [Syllabus](https://www.dropbox.com/s/zq5pn2yca17y0bu/Syllabus%20Math333-S211.pdf?dl=0)
"""

# ╔═╡ 79a8d881-0e07-44b1-a4cd-783b97db124e
md"""


## 9.1 Vector Functions

__Parametric curves:__ Recall that we say that ``C`` is a __parametric curve in space__ if the ``x-, y-`` and ``z-``coordinates of a point on the curve are defined by a set of ordered triples of functions 

```math
x=f(t),\quad  y=g(t), \quad z=h(t)
```
that are continuous on some interval ``a\leq t\leq b``.

__*Vector-Valued Functions*__ The functions ``r: \mathbb{R}\to \mathbb{R}^2 ``  
```math
{\bf{r}}(t) = \langle f(t), g(t)\rangle = f(t){\bf{i}} + g(t){\bf{j}}
```

and ``r: \mathbb{R}\to \mathbb{R}^3 ``  
```math
{\bf{r}}(t) = \langle f(t), g(t), h(t)\rangle = f(t){\bf{i}} + g(t){\bf{j}} +h(t){\bf k}
```
are called __vector-valued functions__ or simply __`vector function`__
"""

# ╔═╡ b691771d-5ce5-4937-a249-cb01da6243ac
begin
	g(t) = [2*cos(t), 2*sin(t), t]
	default(label=nothing)
	plt91_1 = plot(t->2*cos(t),t->2*sin(t),t->t,0,4π, 
			frame_style=:origin,
			aspect_ratio=1
		)
	# scatter!(t->1,t->0,t->0,0,1)
	# arrow!([0,0,0],g(π/2))
	# arrow!([0,0,0],g(π))
	# arrow!([0,0,0],g(2π))
	
	md"""
	__Example__ (Circular Helix)
	
	Graph the curve traced by the vector function
	```math
	\mathbf{r}(t) = 2\cos t \mathbf{i}+2\sin t \mathbf{j} + t \mathbf{k}, \quad t\geq 0.
	```
	
	$plt91_1
	"""
end

# ╔═╡ 30c6cfaa-abca-4099-8bee-c6cc582ef3b2
md"""
### Limits, Continuity, and Derivatives

__Definition (Limit of a Vector Function)__

If ``\lim_{t\to a} f(t),\lim_{t\to a} g(t),`` and ``\lim_{t\to a} h(t)`` exist, then
```math

\lim_{t\to a} \mathbf{r}(t)=\left\langle \lim_{t\to a} f(t),\lim_{t\to a} g(t),\lim_{t\to a} h(t)\right\rangle.
```

*Remark*
- ``t\to a`` can be replaced by ``t\to a^+``, ``t\to a^-``, ``t\to\infty`` or ``t\to-\infty``

__Theorem (Properties of Limits)__

If ``\lim_{t\to a} \mathbf{r}_1(t)=L_1`` and ``\lim_{t\to a} \mathbf{r}_2(t)=L_2``, then
```math
\begin{array}{lll}
\text{(i)} & \lim_{t\to a}c\mathbf{r}_1(t) = c L_1, \quad c \text{ a scalar}\\ \\
\text{(ii)} & \lim_{t\to a}[\mathbf{r}_1(t)+\mathbf{r}_2(t)] =  L_1+L_2\\ \\
\text{(iii)} & \lim_{t\to a}\mathbf{r}_1(t)\cdot \mathbf{r}_2(t) =  L_1\cdot L_2\\ \\
\end{array}
```

__Definition (Continuity of a Vector Function)__

A vector function ``\mathbf{r}`` is said to be continuous at ``t=a`` if

(i) ``\mathbf{r}(a)``,  

(ii) ``\lim_{t\to a}\mathbf{r}(t)`` exists, and 

(iii) ``\lim_{t\to a}\mathbf{r}(t)=\mathbf{r}(a)``.

*Remark*
* ``\mathbf{r}`` is continuous at ``t=a`` if and only if ``f,g`` and ``h`` are continuous there.

"""

# ╔═╡ 2d9f3647-327d-4b3f-b002-a4e1e9dbd59c
md"""

__Definition (Derivative of a Vector Function)__

The __derivative__ of a vector function ``\mathbf{r}`` is
```math
\mathbf{r}'(t) = \lim_{\Delta t\to 0}\frac{\mathbf{r}(t+\Delta t)-\mathbf{r}(t)}{\Delta t}
```
for all ``t`` for which the limit exists.


__Theorem (Differentiation of Components)__
If ``\mathbf{r}(t)=\langle f(t), g(t), h(t) \rangle``, where ``f, g`` and ``h`` are differentiable, then
```math
\mathbf{r}'(t) = \langle f'(t), g'(t), h'(t) \rangle.
```

_Remarks_
- _Smooth Curves_: A vector function ``\mathbf{r}`` have continuous first derivatives and ``\mathbf{r}'(t)\not=0`` for all ``t\in (a,b)``, then ``\mathbf{r}`` is said to be __smooth function__ and the curve ``C`` traced by ``\mathbf{r}`` is called a __smooth curve__.
- _Geometric Interpretation of ``\mathbf{r}'(t) `` at a point ``P``_ represents the direction of the tangent line at ``P``.
"""

# ╔═╡ 5437877a-5788-463c-97b3-3a203b878917
md"""
__Example__

Find parametric equations of the tangent line to the graph of the curve ``C`` whose parametric equations are 
```math
x=t^2, \quad y=t^2-t, \quad z=-7t
```
at ``t=3``.

__Solution__ in class
"""

# ╔═╡ 90db9f5e-908a-45ca-acb5-7fded8d28183
md"""
### Higher-Order Derivatives
Higher-order derivatives of a vector function are also obtained by differentiating its components.

__Theorem (Chain Rule)__

If ``\mathbf{r}`` is a differentiable vector function and ``s=u(t)`` is a differentiable scalar function, then the derivative of ``mathbf{r}(s)`` with respect to ``t`` is
```math
\frac{d\mathbf{r}}{dt}=\frac{d\mathbf{r}}{ds}\frac{ds}{dt}=\mathbf{r}'(s)u'(t).
```


__Theorem (Rules of Differentiation)__

Let ``\mathbf{r}_1`` and ``\mathbf{r}_2`` be differentiable vector functions and ``u(t)`` a differentiable scalar function

- ``\frac{d}{dt}\left[\mathbf{r}_1(t)+\mathbf{r}_2(t)\right]=\mathbf{r}_1'(t)+\mathbf{r}_2'(t)``.

- ``\frac{d}{dt}\left[u(t)\mathbf{r}_1(t)\right]=u(t)\mathbf{r}_1(t)'+u'(t)\mathbf{r}_1(t)``.


- ``\frac{d}{dt}\left[\mathbf{r}_1(t)\cdot\mathbf{r}_2(t)\right]=\mathbf{r}_1(t)\cdot\mathbf{r}_2'(t)+\mathbf{r}_1'(t)\cdot\mathbf{r}_2(t)``.


- ``\frac{d}{dt}\left[\mathbf{r}_1(t)\times\mathbf{r}_2(t)\right]=\mathbf{r}_1(t)\times\mathbf{r}_2'(t)+\mathbf{r}_1'(t)\times\mathbf{r}_2(t)``.

"""

# ╔═╡ cadb7973-85ae-431e-9a61-448a2243b72e
md"""
### Integral of a Vector Function
```math
\begin{array}{lcl}
\int \mathbf{r}(t) dt &=& \left[\int f(t) dt\right]\mathbf{i} 
+ \left[\int g(t) dt\right]\mathbf{j}
+ \left[\int h(t) dt\right]\mathbf{k}
\\ \\
\int_a^b \mathbf{r}(t) dt &=& \left[\int_a^b f(t) dt\right]\mathbf{i} 
+ \left[\int_a^b g(t) dt\right]\mathbf{j}
+ \left[\int_a^b h(t) dt\right]\mathbf{k}
\\ \\
\end{array}
```

__Length of a Space Curve__
```math
s=\int_a^b \sqrt{[f'(t)]^2+[g'(t)]^2+[h'(t)]^2} dt = \int_a^b\|\mathbf{r}'(t)\| dt
```

__Arc Length as a Parameter__ A curve is the plane or in space can be parameterized in terms of the arc length ``s``

__Example__
In (Circular Helix)
```math
	\mathbf{r}(t) = 2\cos t \mathbf{i}+2\sin t \mathbf{j} + t \mathbf{k}, \quad t\geq 0.
```

"""

# ╔═╡ d1013f55-85ea-4e0b-9499-830511054374
begin
	fig9_5_2=load(download("https://www.dropbox.com/s/kiidfzkegm7p25b/fig9_5_2.png?dl=0"))
md"""
## 9.5: Directional Derivative
	
__*The Gradient of a Function*__ 
```math
\nabla f(x,y,z) = \frac{\partial f}{\partial x} \mathbf{i} +\frac{\partial f}{\partial y} \mathbf{j} + \frac{\partial f}{\partial z} \mathbf{k} =\underset{\text{vector differential operator}}{\left(\frac{\partial }{\partial x} \mathbf{i}+\frac{\partial }{\partial y} \mathbf{j}+\frac{\partial }{\partial z} \mathbf{k}\right)}f(x,y,z)
```
__Example__ 
	If ``F(x, y, z) = xy^2 + 3x^2 - z^3`` , find ``\nabla F(x, y, z)`` at ``(2, 1, 4)``.
	
### A Generalization of Partial Differentiation: Directional Derivative
	
$fig9_5_2
	
"""
end

# ╔═╡ e8cb1b44-72df-4544-8c60-fa9413699391
md"""
__Definition__: The directional derivative of ``z = f(x, y)`` in the direction of a unit vector ``\mathbf{u}=\cos θ \mathbf{i} + \sin θ\mathbf{j}`` is
```math
D_{\mathbf{u}} f(x, y) = \lim_{h\to 0}\frac{f(x + h \cos θ, y + h \sin θ) - f(x, y)}{h}
```
provided the limit exists.

- __Functions of Three Variables__ For a function ``w = F(x, y, z)`` the __directional derivative__ is defined by
```math
D_{\mathbf{u}} F(x, y,z) = \lim_{h\to 0}\frac{F(x + h \cos α, y + h \cos β,z+h\cosγ) - F(x, y,z)}{h}
```
where ``α , β``, and ``γ`` are the direction angles of the unit vector ``\mathbf{u}`` measured relative to the positive ``x-``, ``y-``, and ``z-``axes, respectively.

__Computing a Directional Derivative__
If ``z = f (x, y)`` is a differentiable function of ``x`` and ``y`` and ``\mathbf{u}=\cos θ\mathbf{i} + \sin θ \mathbf{j}``, then
```math
 D_{\mathbf{u}}f(x, y) =\nabla f (x, y) \cdot \mathbf{u}
```

- For ``w=F(x,y,z)``,
```math
 D_{\mathbf{u}}F(x, y,x) =\nabla f (x, y,z) \cdot \mathbf{u}
```
"""

# ╔═╡ ca122aa2-cc1a-45c5-bf96-024eeca7b34e
md"""
__EXAMPLE__ (Directional Derivative)

Find the directional derivative of ``F(x, y, z) = xy^2 - 4x^2y + z^2`` at ``(1, -1, 2)`` in the direction of ``6 \mathbf{i} + 2 \mathbf{j} + 3 \mathbf{k}``

"""

# ╔═╡ 9d1bddd0-45bb-4789-a388-46719159f497
md"""
### Maximum Value of the Directional Derivative
Let ``f`` represent a function of either two or three variables. Then
```math
 D_{\mathbf{u}} f = ||\nabla f||||\mathbf{u}||\cos \phi = ||\nabla f||\cos \phi, \quad (||\mathbf{u}||=1)
```
where ``f`` is the angle between ``\nabla f`` and ``\mathbf{u}``. Because ``0 \leq \phi\leq π``, we have ``-1 \leq  \cos \phi \leq 1`` and, 
consequently, ``-||\nabla f|| \leq D_{\mathbf{u}} f \leq ||\nabla f ||``. In other words:

- The __maximum value__ of the directional derivative is ``||\nabla f ||`` and it occurs when ``\mathbf{u}`` has the same direction as ``\nabla f`` (when ``\cos\phi= 1``).

- The __minimum value__ of the directional derivative is ``-||\nabla f||`` and it occurs when ``\mathbf{u}`` and ``\nabla f`` have opposite directions (when ``\cos \phi =-1``)

__Remark__
- ``\nabla f`` points in the direction in which ``f`` __increases__ most rapidly, whereas 
- ``-\nabla f`` points in the direction of in which ``f`` has the most rapid __decrease__.
"""

# ╔═╡ fa1b0c79-0c40-4e88-ba65-e2565fc9b604
md"""
__Example__: The temperature in a rectangular box is approximated by
```math
T(x, y, z) = xyz(1 - x)(2 - y)(3 - z),\quad  0≤x≤1, 0 ≤ y ≤ 2, 0 ≤ z ≤ 3.
```

If a *mosquito* is located at ``({1\over 2}, 1, 1)``, in which direction should it fly to cool off as rapidly as possible?

"""

# ╔═╡ 0f489937-3266-4d83-880d-c4ec0ebd2d34
md"## 9.7 Curl and Divergence"

# ╔═╡ 392780ee-6c31-42fb-b311-a7e7df34b1dc
html"<iframe width=\"380\" height=\"220\" src=\"https://www.youtube.com/embed/rB83DpBJQsE\" title=\"YouTube video player
\" frameborder=\"0\" allow=\"accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture\" allowfullscreen></iframe>"


# ╔═╡ 751bff35-70ab-4d5a-909c-692157146226

md"""

__Vector Fields__ Vector functions of two and three variables,
```math
\begin{array}{lcl} 
F(x, y) &=& P(x, y) \mathbf{i} + Q(x, y) \mathbf{j}, \\ \\
F(x, y, z) &=& P(x, y, z) \mathbf{i} + Q(x, y, z) \mathbf{j} + R(x, y, z) \mathbf{k}
\end{array}
```
are also called __vector fields__.
"""

# ╔═╡ 509d1cbe-f592-4768-a29c-03881e34441c
begin
	xs = -1:0.5:1
	ys = -1:0.5:1
	
	
	df(x, y) = normalize([-y, x]) ./ 10
	
	xxs = [-y for x in xs for y in ys]
	yys = [x for x in xs for y in ys]
	
	plt=quiver(xxs, yys, quiver=df)
	plot!(;frame_style=:origin)
	md"""
	__Example__: Graph the two-dimensional vector field ``F(x, y) = -y \mathbf{i} + x \mathbf{j}``.
	
	$plt
	"""
end

# ╔═╡ 506abf7e-04b1-4220-9e10-bc9e01383ec3
md"""
The del ``\nabla`` operator
```math
\nabla = \frac{\partial}{\partial x}\mathbf{i}
+ \frac{\partial}{\partial y}\mathbf{j}
+ \frac{\partial}{\partial z}\mathbf{k}
```
combined with a scalar function ``\phi(x, y, z)`` gives the __gradient of ``ϕ``__
```math
F(x,y,z)=\nabla\phi  = \frac{\partial\phi}{\partial x}\mathbf{i}
+ \frac{\partial\phi}{\partial y}\mathbf{j}
+ \frac{\partial\phi}{\partial z}\mathbf{k}
```
combined with a vector function 
```math
F(x, y, z) = P(x,y,z)\mathbf{i} + Q(x,y,z)\mathbf{j}+R(x,y,z)\mathbf{k}
```
gives a vector field and a scalar function
* the __curl of ``F``__
```math
\text{curl } F= \nabla \times F
```
* the __divergence of ``F``__
```math
\text{div } F= \nabla \cdot F
```

 
"""

# ╔═╡ 3e2e4f04-b687-47e7-96d8-6774724a9207
md"""
__EXAMPLE__  If ``F=  (x^2y^3 - z^4)\mathbf{i} \;+\; 4x^5y^2z\mathbf{j} \;-\; y^4
z^6 \mathbf{k}``, find 

(a) ``\text{curl } F``, 

(b) ``\text{div } F``, 

(c) ``\text{div}(\text{curl } F)``
"""

# ╔═╡ 330410cd-40e8-4398-938d-ee5e4be6f650
md"""
__Properties__
- If f is a scalar function with continuous second partial derivatives, then
```math
\text{curl}(\text{grad } f) = \nabla \times \nabla f = 0
```
- If F is a vector field having continuous second partial derivatives, then
```math
\text{div}(\text{curl } F) = \nabla \cdot \left(\nabla \times F\right) = 0
```


"""

# ╔═╡ b7452e63-29ae-4c5e-87fa-5faecedf986c
md"""
## 9.8 Line Integrals
__*Terminology*__ Suppose ``C`` is a curve parameterized by ``x = f (t), y = g(t), a \leq t \leq b``, and ``A`` and ``B`` are the points ``( f (a), g(a))`` and ``( f (b), g(b))``, respectively. We say that

1. ``C`` is a __smooth curve__ if ``f'`` 	 and ``g'``	 are continuous on the closed interval ``[a, b]`` and not simultaneously zero on the open interval ``(a, b)``.
2. ``C`` is __piecewise smooth__ if it consists of a finite number of smooth curves ``C_1, C_2, \cdots , C_n`` joined end to end—that is, ``C = C_1 \cup C_2 \cup \cdots \cup C_n``.
3. ``C`` is a __closed curve__ if ``A = B``.
4. ``C`` is a __simple closed curve__ if ``A = B`` and the curve does not cross itself.
5. If ``C`` is not a closed curve, then the __positive direction__ on ``C`` is the direction corresponding to increasing values of ``t``.

"""

# ╔═╡ 0de4e7c3-216a-432a-b998-493a963901cf
md"""
__Definition__ (*Line Integrals in the Plane*)

Let ``G`` be a function of two variables ``x`` and ``y`` defined on a region of the plane containing a smooth curve ``C``.
* The __line integral of ``G`` along ``C`` from ``A`` to ``B`` with respect to ``x``__ is
```math
\int_C G(x, y) dx = \lim_{||P||\to 0}\sum_{k=1}^{n}{G(x_k^*,y_k^*)}\Delta x_k.
```
* The __line integral of ``G`` along ``C`` from ``A`` to ``B`` with respect to ``y``__ is
```math
\int_C G(x, y) dy = \lim_{||P||\to 0}\sum_{k=1}^{n}{G(x_k^*,y_k^*)}\Delta y_k.
```
* The __line integral of ``G`` along ``C`` from ``A`` to ``B`` with respect to arc length ``s``__ is
```math
\int_C G(x, y) ds = \lim_{||P||\to 0}\sum_{k=1}^{n}{G(x_k^*,y_k^*)}\Delta s_k.
```
"""

# ╔═╡ 8451fe29-cc80-426e-bc65-43ca5eafb5d3
Resource("https://www.dropbox.com/s/yjmo0w2fcyg6b5m/Line%20Integral.png?raw=1")

# ╔═╡ 4844492f-27b9-4350-bbd6-66d97b87e909
Resource("https://www.dropbox.com/s/2o6fus9653k39gd/Line%20Integral%202.png?raw=1",:width=>600)

# ╔═╡ a3f94e01-c503-40f7-a4b2-d40d25401fe3
@htl("""
<div style="color: red;font-weight:800;">Method of Evaluation—Curve Defined Parametrically</div>
<div style="margin: 20px 20px;">$(md"""
If ``C`` is a smooth curve parameterized by ``x= f (t),y = g(t), a \leq t \leq b``,
```math
{\large
\begin{array}{lcl}
\int_C G(x,y)dx &=& \int_a^b G(f(t),g(t))f'(t) dt\\ \\
	

\int_C G(x,y)dx &=& \int_a^b G(f(t),g(t))g'(t) dt \\ \\
\int_C G(x,y)dx &=& \int_a^b G(f(t),g(t))\sqrt{\left[f'(t)\right]^2+\left[g'(t)\right]^2} dt\\
	\end{array}
}
```
""")</div>
""")

# ╔═╡ 52f1f384-ba7f-4a1d-b1fb-311a10ea21a9
md"""
A line integral along a __piecewise-smooth curve__ ``C``. For example, if ``C`` is composed of smooth curves ``C_1`` and ``C_2``, then
```math
\int_C G(x, y) ds =\int_{C_1} G(x, y) ds +\int_{C_2} G(x, y) ds
```
__*Notation*__ In many applications, line integrals appear as a sum
```math
\int_C P(x, y) dx +\int_C Q(x, y) dy.
```
It is common practice to write this sum as one integral without parentheses as
```math
\int_C P(x, y) dx + Q(x, y) dy \quad \text{ or simply }\quad \int_C P dx + Q dy. 
```

A line integral along a __closed curve ``C``__ is very often denoted by
```math
\oint_C P dx + Q dy.
```
"""

# ╔═╡ f185e1e7-fc9e-43d2-9f41-43be18c444c3
md"""
__Example__ Evaluate ``\oint_C xy^2 dx``, where ``C`` is the quarter-circle ``x = 4\cos t, y = 4\sin t, 0 ≤ t ≤ {π\over 2}``

__Example__ Evaluate ``\oint_C x dx``, where ``C`` is the circle ``x = \cos t, y = \sin t, 0 ≤ t ≤ 2π``

"""

# ╔═╡ c6ffb2c9-3b92-4032-ab96-1804b11290f5
begin
	xEx984=0:0.1:2
	pltEx984=plot(xEx984,xEx984.^2,frame_style=:origin;c=:blue)
	plot!(pltEx984;xlims=(-1,3))
	plot!(pltEx984,repeat([2],11),0:.4:4,c=:blue)
	plot!(pltEx984,0:.2:2,repeat([0],11),c=:blue)
	annotate!(pltEx984,
			[(1.2,0,("^",20,-90.0,:top,:red)),
			 (2,1.2,("^",20,0.0,:bottom,:red)),
			(1.2,1.25,("^",20,135.0,:left,:red)),
			(1.1,2,(L"y=x^2",14,0.0,:bottom,:red))
			])
	md"""
	__Example__ Evaluate ``\oint_C  y^2 dx - x^2 dy`` on the closed curve ``C`` that is shown below
	
	$pltEx984
	"""
end

# ╔═╡ 62afd674-74c5-429d-9d21-1673fdb7b385
md"""
__Remark__: 
```math
\int_{-C} Pdx + Q dy = -\int_{C} Pdx + Q dy
```
"""

# ╔═╡ 4c4bd58e-04d5-4da0-8943-6d128a6b749f
md"""
__Line Integrals in Space__ 
If ``C`` is a smooth curve in 3-space defined by the parametric 
equations ``x = f (t), y = g(t), z = h(t), a \leq t \leq b``, then 
```math
\int_C G(x, y, z) dz = \int^b_a G( f(t), g(t), h(t)) h'(t) dt.
```
and
```math
\int_C G(x, y, z) ds = \int^b_a G( f(t), g(t), h(t)) 
\sqrt{\left[f'(t)\right]^2+\left[g'(t)\right]^2+\left[h'(t)\right]^2}
dt.
```
"""

# ╔═╡ db44a1a7-705d-4fc9-a16b-cd3032353b3a
md"""
__Example__ Evaluate 
```math
\int_C y dx + x dy + z dz,
``` 
where ``C`` is the helix ``x = 2\cos t, y = 2 \sin t, z = t, 0 \leq t \leq 2\pi``.
"""

# ╔═╡ 17ab5bf0-65b0-4079-99a0-ccd19721d5c6
md"""
__*Remark*__ Let ``F(x,y,z)=P(x,y,z) \mathbf{i}+Q(x,y,z)+\mathbf{j} + R(x,y,z)\mathbf{k}`` be defined along a curve ``C: x = f (t)``, ``y = g(t)``, ``z=h(t)``, ``a \leq t \leq b``, and suppose ``\mathbf{r}(t) = f (t) \mathbf{i} + g(t) \mathbf{j} +h(t) \mathbf{k}`` is the position vector of points on ``C``. We can write the line integral as
```math
\int_C P dx + Qdy+Rdz = \int_C F\cdot d\mathbf{r}

```



"""

# ╔═╡ 07f6f180-ac84-4952-a24f-234990c35191
md"""
### Work
We define the __work__ done by a force field ``F(x,y)=P(x,y)\mathbf{i}+Q(x,y)\mathbf{j}`` along a smooth curve ``C; x=f(t), y=g(t), a\leq t\leq b`` as the line integral
```math
W = \int_C P(x, y) dx + Q(x, y) dy\quad \text{or} \quad  W = \int_C F\cdot dr.
```
"""

# ╔═╡ b942c102-5698-4086-8056-a053b0ca9d43
md"""
__Example__ Find the work ``W = \int_C F · dr`` done by the force ``F = -y^2\mathbf{i} + xy\mathbf{j}`` acting along the curve ``C`` defined by ``x = 2t, y = t^3, 0 ≤ t ≤ 2``

"""

# ╔═╡ f84606c4-3de8-4156-a4ea-3df383974de5
md"""
## 9.9 Independence of the Path
In ``2-``space, if 
```math
\mathbf{F}(x,y) = P(x,y)\mathbf{i} +Q(x,y)\mathbf{j} \quad \text{is a vector field.}
```
and ``C`` is a __path__ defined as
```math
\mathbf{r}(t) = f(t)\mathbf{i} + g(t)\mathbf{j}, \quad a\leq t\leq b
```
then
```math
\int_C Pdx +Qdy = \int_C\mathbf{F}\cdot d\mathbf{r},
```
where
```math
d\mathbf{r} = dx\mathbf{i} + dy\mathbf{j}
```
__*Remark*__

The value of a line integral ``\int_C\mathbf{F}\cdot d\mathbf{r}`` depends on the path of integration. 
"""

# ╔═╡ f389ba80-79ea-45a8-a8cb-3fdbdcfb1971
Ex991Paths = @bind ex991path Radio(["c1"=>"C₁", "c2"=>"C₂","c3"=>"C₃","c4"=>"C₄" ], default="c1");html""

# ╔═╡ c3222e2a-cf12-42ba-bc9a-0f914dbfb529
md"__Example__ Evaluate ``\int_C  y dx + x dy`` on each path shown below

$Ex991Paths 
"

# ╔═╡ 6363ee6a-df6d-4150-b4c3-5e959bdf8c6c
begin
	pltToDisp = plot(;frame_style=:none)
	indx = parse(Int64,SubString(ex991path,2:2))
	ex991Solutions = [
			L"\int_{C_1}  y dx + x dy=\int_0^1 x^2dx+2x^2dx=\int_0^1 3x^2 dx=1" 
			L"\int_{C_2}  y dx + x dy=\int_0^1 xdx+xdx=\int_0^1 2x dx=1" 
			L"\int_{C_3}  y dx + x dy=\int_0^1 y\cdot 0+0\cdot 0+\int_0^1 1 dx+x\cdot 0 =\int_0^1 dx=1" 
			L"\int_{C_3}  y dx + x dy=\int_0^1 0\cdot dx+1\cdot 0+\int_0^1 y\cdot 0+1\cdot dy =\int_0^1 dy=1" 
		]
	annotate!(pltToDisp,
			[
			 (.5,.2,(ex991Solutions[indx]))
			])
	md"""
	
	
	"""
end

# ╔═╡ 60db1dcb-ef1e-4a4b-b4ad-e2d5415601d8
begin
	xEx991=0:0.1:1
	ticks = [-0.5,0,1,2]
	ticklabels = [ "$x" for x in ticks ]
	pltEx991=plot(frame_style=:origin;
			c=:blue, xlims=(-0.5,2), 
			ylims=(-0.5,1.2), 
			yticks=(0:1),
			xticks=(ticks,ticklabels))
	annotate!(pltEx991, 
				[ (0,-0.2,(L"(0,0)",8)),
				  (1.2,1,(L"(1,1)",8)),
				  (0,1.3,(L"y",12)),
				  (2,-0.3,(L"x",12))
				]
			 )
	scatter!(pltEx991,[0 1],[0 1], c=:black)
	
	pltEx991a=plot(pltEx991,xEx991,xEx991.^2,c=:blue)
	annotate!(pltEx991a,[(0.4,0.7,(L"y=x^2",12,0.0,:top,:black)),
			(0.5,0.25,("^",12,-45.0,:top,:red)),
			(0.8,0.8^2,("^",12,-45.0,:top,:red))] )
	
		
	pltEx991b=plot(pltEx991,xEx991,xEx991,c=:blue)
	annotate!(pltEx991b,[(0.6,0.3,(L"y=x",12,0.0,:top,:black)),
			(0.3,0.3,("^",12,-45.0,:top,:red)),
			(0.65,0.65,("^",12,-45.0,:top,:red))] )
	
	
	pltEx991c=plot(pltEx991,repeat([0],11),0:.1:1,c=:blue)
	plot!(pltEx991c,0:.1:1,repeat([1],11),c=:blue)
	annotate!(pltEx991c,[(0,0.5,("^",12,0,:top,:red)),
			(0.65,1,("^",12,-90.0,:top,:red))] )
	
	pltEx991d=plot(pltEx991,repeat([1],11),0:.1:1,c=:blue)
	plot!(pltEx991d,0:.1:1,repeat([0],11),c=:blue)
	annotate!(pltEx991d,[(0.5,0,("^",12,-90.0,:top,:red)),
			(1,0.65,("^",12,0.0,:top,:red))] )
	
	l = @layout([a b;c d; e])
	plt991=plot(pltEx991a,pltEx991b,pltEx991c,pltEx991d,pltToDisp, layout=l, 
		title=[L"C_1" L"C_2" L"C_3" L"C_4"])
	md"""
	$plt991
	"""
end

# ╔═╡ ca86c2c9-7416-4d45-bcec-a644c151a43f
md"""
### Conservative Vector Field
__Definition__: A vector function ``\mathbf{F}`` in ``2-`` or ``3-``space is said to be __conservative__ if ``\mathbf{F}`` can be written as the 
gradient of a scalar function ``\phi``. The function ``\phi`` is called a __potential function__ for ``\mathbf{F}``.

__*Example*__:  Is ``\mathbf{F}(x,y)=y\mathbf{i}+x\mathbf{j}`` conservative?
"""

# ╔═╡ efd8e89b-6f7d-44d2-b798-5e0c06fd3622
md"""
### Path Independence 
If the value of a line integral is the same for every path in a region 
connecting the initial point ``A`` and terminal point ``B``, then the integral is said to be __independent of the path__.

"""

# ╔═╡ fba21f02-4aa0-46a1-b900-386647c20450
md"""
__Theorem__ Fundamental Theorem for Line Integrals

Suppose ``C`` is a path in an open region ``R`` of the ``xy-``plane and is defined by ``r(t) = x(t)\mathbf{i} + y(t)\mathbf{j}``,``a \leq t \leq b``. If ``\mathbf{F}(x, y) = P(x, y)\mathbf{i} + Q(x, y)\mathbf{j}`` is a conservative vector field in ``R`` and ``\phi`` is a potential function for ``\mathbf{F}``, then
```math
\int_C \mathbf{F}\cdot d\mathbf{r} =\int_C f\cdot \nabla \phi\cdot d\mathbf{r} = \phi(B) - \phi(A), 
```
where ``A = (x(a), y(a))`` and ``B = (x(b), y(b))``.

__*Example*__
Evaluate ``\int_C ydx + xdy``, where ``C`` is a path with initial point ``(0, 0)`` and terminal point ``(1, 1)``.
"""

# ╔═╡ ed3590a1-b43a-4212-81a8-bed8b21f0858
md"""
__Terminology__

- We say that a region (in the plane or in space) is __connected__ if every pair of points ``A`` and ``B`` in the region can be joined by a piecewise-smooth curve that lies entirely in the region. 
- A region ``R`` in the plane is __simply connected__ if it is connected and every simple closed curve ``C`` lying entirely within the region can be shrunk, or contracted, to a point without leaving ``R``. A simply connected region has no holes in it. 
- A region that is not connected is called __disconnected__.
- A region that is connected with multiple holes is called __multiply connected__.
- A region ``R`` is said to be __open__ if it contains no boundary points.

__Theorem__ 

In an open connected region ``R``, ``\int_C \mathbf{F}\cdot d\mathbf{r}`` is independent of the path ``C`` if and only if the vector field ``\mathbf{F}`` is conservative in ``R``.
"""

# ╔═╡ 0615d26b-858d-4a03-8595-e85da72adca7
md"""
###  Integrals Around Closed Paths
__Theorem__

In an open connected region ``R``, ``\int_C \mathbf{F}\cdot d\mathbf{r}`` is independent of the path if and only if ``\int_C \mathbf{F}\cdot d\mathbf{r}=0``
for every closed path ``C`` in ``R``.
"""

# ╔═╡ 36380af4-f8de-4b7e-97c8-65673754ff76
md"""
__Summary__

```math
\mathbf{F} \text{ conservative} \Longleftrightarrow \text{ path independence} \Longleftrightarrow 
\int_C \mathbf{F}\cdot d\mathbf{r}=0.
```

"""

# ╔═╡ 5febf63a-7834-4248-8884-7be76bb08202
md"""
### Test for a Conservative Field
__Theorem (*Test for a Conservative Field*)__

Suppose ``\mathbf{F}(x, y) = P(x, y)\mathbf{i} + Q(x, y)\mathbf{j}`` is a conservative vector field in an open region ``R``, and that ``P`` and ``Q`` are continuous and have continuous first partial derivatives in ``R``. Then
```math
\frac{\partial P}{\partial y} = \frac{\partial Q}{\partial x}
```
for all ``(x, y)`` in ``R``. Conversely, if the equality holds for all ``(x, y)`` in a simply connected region ``R``, then ``\mathbf{F}`` is conservative in ``R``.

"""

# ╔═╡ d1b8b4eb-e856-4745-9b83-808e2c55b49f
md"""
__Example__

Determine which of the vector fields is conservative
```math
\begin{array}{lcl}
\mathbf{F}(x, y) &=& (x^2-2y^3)\mathbf{i} + (x + 5y)\mathbf{j} \\
\mathbf{F}(x, y) &=& -ye^{-xy}\mathbf{i} + -xe^{-xy}\mathbf{j} \\
\end{array}
```
"""

# ╔═╡ b0a55996-3e23-4696-9a5c-94284b982155
md"""
__Example__ 

(a) Show that ``\int_C \mathbf{F}\cdot d\mathbf{r}``, where ``\mathbf{F}(x, y) = (y^2 - 6xy + 6)\mathbf{i} + (2xy - 3x^2 - 2y)\mathbf{j}``, is independent of the path ``C`` between ``(-1, 0)`` and ``(3, 4)``.

(b) Find a potential function ``\phi`` for ``\mathbf{F}``.

(c) Evaluate ``\int_{(-1,0)}^{(3,4)}\mathbf{F}\cdot d\mathbf{r}``.

"""

# ╔═╡ 02561b84-2447-4462-91e9-4872c44e1ebe
md"""
### Conservative Vector Fields in 3-Space
```math
\mathbf{F}(x, y,z) = P(x, y,z)\mathbf{i} + Q(x, y,z)\mathbf{j} + R(x, y,z)\mathbf{k}
```
If ``\mathbf{F}`` is conservative and ``P, Q,`` and ``R`` are continuous and have continuous first partial derivatives in some open region of 3-space, then
```math
\frac{\partial P}{\partial y}
= \frac{\partial Q}{\partial x},
\quad \frac{\partial P}{\partial z}
= \frac{\partial R}{\partial x},
\quad  \frac{\partial Q}{\partial z}=\frac{\partial R}{\partial y}

```
__Remark__
```math
\text{curl}(\mathbf{F}) = 0
```
"""

# ╔═╡ a0f19e1e-2ee6-440c-8c39-9e1f0d07b2d5
x, y = symbols("x,y", real=true);html""

# ╔═╡ 9badddfe-d618-421e-89c9-e54639b94423
md"""
## 9.12 Green’s Theorem
"""

# ╔═╡ b0a38158-bfdb-4b6b-bf70-7d6405da0c2b
md"""
### Line Integrals Along Simple Closed Curves
__Definition__ 
- The __positive direction__ around a simple closed curve ``C`` is that direction a point on the curve must move in order to keep the region ``R`` bounded by ``C`` to the left.
- The __positive__ and *negative* directions correspond to the __counterclockwise__ and *clockwise* directions, respectively.

__Theorem__ (*Green’s Theorem in the Plane*)
Suppose that ``C`` is a piecewise-smooth simple closed curve bounding a simply connected region ``R``. If ``P``, ``Q``, ``\partial P/\partial y``, and ``\partial Q/\partial x`` are continuous on ``R``, then
```math
\oint_C P dx + Q dy = \iint_R \left(\frac{\partial Q}{\partial x}-\frac{\partial P}{\partial y}\right) dA. 
```
"""

# ╔═╡ 506a2816-5d1f-4792-93ed-cbf526888897
md"""
__EXAMPLE__ Evaluate 
```math
\oint_C (x^2 - y^2) dx + (2y - x) dy,
```
where ``C`` consists of the boundary of the region in the 
first quadrant that is bounded by the graphs of 
```math
y = x^2 \quad \text{ and } y = x^3.
```

"""

# ╔═╡ 6bc1fefb-be6c-4e38-9299-13a632d2ee34
md"""
__EXAMPLE__ 
Evaluate (+)
```math
\oint_C (x^5 + 3y) dx + (2x - e^{y^3}) dy, 
```
where ``C`` is the circle ``(x - 1)^2 + ( y - 5)^2 = 4``.
"""

# ╔═╡ de3153ea-e97d-4fda-b70f-bd23b49cd8b7
begin
	sec912ex3 = Resource("https://www.dropbox.com/s/b4uwhqw0ih5kneh/sec912ex3.png?raw=1")
	md"""
	__EXAMPLE__ 
	Find the work done by the force 
	```math
	\mathbf{F} =(-16y - \sin x^2) \mathbf{i} + (4e^y + 3x^2)\mathbf{j}
	```
	acting along the simple closed curve C shown
	
	$sec912ex3
	"""
end

# ╔═╡ 4958f085-527e-4ea7-b004-82d1317e1ed7
begin
	sec912ex4 = Resource("https://www.dropbox.com/s/khnwio2m5kiy1j8/sec912ex4.png?raw=1")
	md"""
	__Example 💣__
	
	Let ``C`` be the closed curve consisting of the four straight line segments ``C_1``, ``C_2``, ``C_3``, ``C_4`` shown below. 
	
	$sec912ex4

	Green’s theorem is not applicable to the line integral
	```math
	\oint_C \frac{-y}{x^2+y^2} dx +  \frac{x}{x^2+y^2} dy
	```
	since ``P, Q, \partial P/\partial y``, and ``\partial Q/\partial x`` are not continuous at the origin.

	"""
end

# ╔═╡ 5aa87065-e94c-425a-bf24-d197aaacafd8
begin
	region_with_holes = Resource("https://www.dropbox.com/s/yvmj37mezm05o1u/region_with_holes.png?raw=1")	
	md"""
	### Region with Holes
	
	$region_with_holes
	
	
	"""
end

# ╔═╡ ea5a1ed0-4494-4850-8943-18d8029cf356
begin
	sec912ex5=Resource("https://www.dropbox.com/s/rnrf7u7fhynyodp/sec912ex5.png?raw=1")
	md"""
	__EXAMPLE__
	Evaluate  
	```math
	\oint_C \frac{-y}{x^2+y^2} dx + \frac{x}{x^2+y^2} dy
	```
	where ``C = C_1 \cup C_2`` is the boundary of the shaded region ``R`` shown 
	
	$sec912ex5
	"""
end

# ╔═╡ 00f35bc7-5456-4a97-a474-166bde07dea8
begin 
	P(x,y)=-y/(x^2+y^2)
	Q(x,y)=x/(x^2+y^2)
	Py=simplify(diff(P(x,y),y))
	Qx=factor(simplify(diff(Q(x,y),x)))
	Py,Qx
end

# ╔═╡ a53e17e1-bb2d-4d41-b8e0-1df381601af5
begin
	sec912rem = Resource("https://www.dropbox.com/s/t2496eg6087q1zs/sec912rem.png?raw=1")
md"""
__*Remark*__
	In the following figure, suppose that ``C_1`` and ``C_2`` are two nonintersecting 
piecewise-smooth simple closed paths that have the same counterclockwise orientation. 
	
$sec912rem
	
Suppose further that ``P`` and ``Q`` have continuous first partial derivatives such that
```math
\frac{\partial P}{\partial y}=\frac{\partial Q}{\partial x}
```
in the region ``R`` bounded between ``C_1`` and ``C_2``. Then
```math
\oint_{C_1} P dx + Q dy =\oint_{C_2} P dx + Q dy.
```

"""
end


# ╔═╡ fb03d35d-9a94-4824-bb95-5c48a9a6ad7f
md"__Example__: Solve example  💣"

# ╔═╡ 7fbbe932-5b35-4fa5-a649-46cba8b5df47
md"""
## 9.13 Surface Integrals
__Definition__ (Surface Area)

Let ``f`` be a function for which the first partial derivatives ``f_x`` and ``f_y`` are continuous on a closed region ``R``. Then the area of the surface over ``R`` is given by

```math
A(S) = \iint_R \sqrt{1 + \left[f_x(x, y)\right]^2 + \left[f_y(x, y)\right]^2} dA.
```
__Differential of Surface Area The function__
```math
dS = \sqrt{1 + \left[f_x(x, y)\right]^2 + \left[f_y(x, y)\right]^2} dA
```
is called the __differential of the surface area__.
"""

# ╔═╡ 5ab19faa-f029-48d4-b9d3-7cef666acd1e
begin
	sec913si=Resource("https://www.dropbox.com/s/kvzr3futqo4pv66/sec913sui.png?raw=1",:width=>500)
md"""
$sec913si
"""
end

# ╔═╡ 4fdfad2e-e860-4358-b888-b6642fa0b6c3
md"""
__Remark__ The surface area integral is a generalization of the arc length integral.
"""

# ╔═╡ 3aa6fe06-ae82-48dc-9999-877533fa21ae
md"""
### Surface Integral
The generalization of the line integral ``\int_C G(x, y) ds`` is called a __surface integral__
"""

# ╔═╡ b4305759-a665-429a-a888-fadf5280b6ba
md"""
__Definition__ (Surface Integral)


Let ``G`` be a function of three variables defined over a region of 3-space containing the surface ``S``. Then the __surface integral of ``G`` over ``S``__ is given by
```math
\iint_S G(x, y, z) dS = \lim_{\|P\|\to 0} \sum_{k=1}^n G(x^*_k, y^*_k, z^*_k) \Delta S_k.
```
__Method of Evaluation__
```math
\iint_S G(x, y, z) dS = \iint_R G(x, y, f(x,y))\sqrt{1 + \left[f_x(x, y)\right]^2 + \left[f_y(x, y)\right]^2} dA 
```
- If ``y=g(x,z)`` and ``R`` on ``xz-``plane, then
```math
\iint_S G(x, y, z) dS = \iint_R G(x, g(x,z),z)\sqrt{1 + \left[g_x(x, z)\right]^2 + \left[g_z(x, z)\right]^2} dA 
```
- If ``x=h(y,z)`` and ``R`` on ``yz-``plane, then
```math
\iint_S G(x, y, z) dS = \iint_R G(h(y,z),y,z)\sqrt{1 + \left[h_y(y, z)\right]^2 + \left[h_z(y, z)\right]^2} dA 
```
"""

# ╔═╡ dd993ad7-758c-4cee-93ad-f9d62277e90b
md"""
__Remark__

- __Mass of a Surface__ Suppose ``\rho(x, y, z)`` represents the density of a surface at any point, or mass per unit surface area; then the mass ``m`` of the surface is
```math
m = \iint_S \rho(x, y, z) dS.
```
__Example__ Find the mass of the surface of the paraboloid ``z = 1 + x^2 + y^2`` in the first octant for ``1 \leq z \leq 5`` if the density at a point ``P`` on the surface is directly proportional to its distance from the ``xy``-plane.
"""

# ╔═╡ 1e20d2a6-bea6-4c16-862c-eaec3bb99d9a
md"""
__Example__
Evaluate ``\iint_S xz^2 dS``, where ``S`` is that portion of the cylinder ``y = 2x^2 + 1`` in the first octant bounded by ``x = 0, x = 2, z = 4``, and ``z = 8``.
"""

# ╔═╡ ffd52892-0427-49b2-8b1e-891d64c6c688
md"""
__Orientable Surfaces__
* Roughly, an __orientable surface ``S``__, with two sides that could be painted different colors. 
* The Möbius strip is not an orientable surface and is one-sided. A person who starts to paint the surface of a Möbius strip at a point will paint the entire surface and return to the starting point

"""

# ╔═╡ 71e53648-d66d-400a-9f5c-1ee9a9b8b9c5
@htl("<iframe width='560' height='315' src='https://www.youtube.com/embed/pgbHR290RW8' title='YouTube video playe' frameborder='0' allow='accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture' allowfullscreen></iframe>")

# ╔═╡ b9c77a2e-c384-4313-bea8-2694ab4dfc56
md"""
We say a smooth surface ``S`` is __orientable__ or is an __oriented surface__ if there exists a continuous unit normal vector function ``\mathbf{n}`` defined at each point ``(x, y, z)`` on the surface. The vector field ``\mathbf{n}(x, y, z)`` is called the __orientation__ of ``S``.


"""

# ╔═╡ e997b0f3-98de-4fd6-8ea5-e4b3440f0f5a
Resource("https://www.dropbox.com/s/6o9vhuyls0lp03g/sec913mobus.png?raw=1",:width=>500)

# ╔═╡ bbe13ec2-af6b-40fa-b4f7-74a68d59298d
md"""
A surface ``S`` defined by ``z = f (x, y)`` has an upward orientation when the unit 
normals are directed upward—that is, have positive ``k`` components, and it has a downward orientation when the unit normals are directed downward—that is, have 
negative ``k`` components. 

If a smooth surface ``S`` is defined by ``g(x, y, z) = 0``, then recall that a unit normal is
```math
\mathbf{n} =\frac{1}{\|\nabla g\|}\nabla g
```
"""


# ╔═╡ 45aac5bb-e392-4eb0-8a5e-976fa6f17938
begin
	sec913flux = Resource("https://www.dropbox.com/s/60txes56q202gij/sec913flux.png?raw=1",:width=>300)
md"""
### Integrals of Vector Fields
If ``\mathbf{F}(x, y, z) = P(x, y, z) \mathbf{i} + Q(x, y, z) \mathbf{j} + R(x, y, z) \mathbf{k}`` is the velocity field of a fluid, then, as shown here,

$sec913flux
the volume of the fluid flowing through an element of surface area ``S`` per unit time is approximated by
```math
\text{(height)(area of base)} = (\text{comp}_n\mathbf{F})\Delta S = (F \cdot n) \Delta S
```
"""
end

# ╔═╡ 904f1f79-4866-4837-b128-c8a30567fa76
md"""
The total volume of a fluid passing through ``S`` per unit time is called the __flux of ``\mathbf{F}`` through ``S``__ and is given by
```math 
\text{flux} =\iint_S (\mathbf{F} \cdot \mathbf{n}) dS. \tag{*}
```
In the case of a closed surface ``S``, if ``\mathbf{n}`` is the outer (inner) normal, then (*) gives the volume of fluid flowing out (in) through ``S`` per unit time.
"""

# ╔═╡ ce1a2204-adfd-4862-a5fc-8637dcc1ce9d
md"""
__EXAMPLE__ Flux Through a Surface
Let 
```math 
\mathbf{F}(x, y, z) = z \mathbf{j} + z \mathbf{k}
```
represent the flow of a liquid. Find the flux of ``F`` through the surface ``S``
given by that portion of the plane 
```math
z = 6 - 3x - 2y
```
in the first octant oriented upward.
"""

# ╔═╡ 81441e8f-253c-44a3-b160-c0ac350791ae
md"""
__Remark__ 
Then the flux of a vector field ``\mathbf{F}`` out of the surface ``S`` is
```math
\iint_S \mathbf{F}\cdot \mathbf{n} dS =\iint_{S_1} \mathbf{F}\cdot \mathbf{n} dS
+\iint_{S_2} \mathbf{F}\cdot \mathbf{n} dS
```
where we take ``S_1`` oriented upward and ``S_2`` oriented downward.
"""

# ╔═╡ ea174feb-97f8-497d-8a47-2fcdfeef29c8
md"""
## 9.14 Stokes’ Theorem
Green’s theorem can be written in vector notation as
```math
\oint_C \mathbf{F}\cdot d\mathbf{r}=\oint_C \mathbf{F}\cdot\mathbf{T} ds
=\iint_R (\text{curl }\mathbf{F})\cdot \mathbf{k} dA
```
that is, the line integral of the tangential component of ``\mathbf{F}`` is the double integral of the normal component of ``\text{curl }\mathbf{F}``.
"""

# ╔═╡ 47751cc9-f4e0-400e-8b6e-5c37c954ce96
md"""
__Stokes’ Theorem__

Let ``S`` be a piecewise-smooth orientable surface bounded by a piecewise-smooth simple closed curve ``C``. Let 
```math
\mathbf{F}(x, y, z)= P(x, y, z) \mathbf{i} + Q(x, y, z) \mathbf{j} + R(x, y, z) \mathbf{k}
```
be a vector field for which ``P,Q``, and ``R`` are continuous and have continuous first partial derivatives in a region of 3-space containing ``S``. If ``C`` is traversed in the positive direction, then 
```math
\oint_C \mathbf{F}\cdot d\mathbf{r}=\oint_C (\mathbf{F}\cdot \mathbf{T}) dS 
= \iint_S (\text{curl }\mathbf{F})\cdot \mathbf{n} dS,
```
where ``\mathbf{n}`` is a unit normal to ``S`` in the direction of the orientation of ``S``.
"""

# ╔═╡ a89ccf32-03b9-45c3-be77-7029ee8fdfee
begin
sec913intro = Resource("https://www.dropbox.com/s/gtbyo6j5yek2gsv/sec913intro.png?raw=1")
md""" 

$sec913intro
"""
end

# ╔═╡ a0d30cd9-2fbf-45be-b20d-da3619755922
md"""
__EXAMPLE__

Let ``S`` be the part of the cylinder ``z = 1 - x^2`` for ``0 \leq x \leq 1``, 
``-2\leq y \leq 2``. Verify Stokes’ theorem for the vector field 
```math
\mathbf{F} = xy \mathbf{i} + yz \mathbf{j} + xz \mathbf{k}.
```
Assume ``S`` is oriented upward.

"""

# ╔═╡ b3f1155f-fc04-4a35-ab00-d9f7851c5b12
md"""
__EXAMPLE__
Evaluate 
```math
\oint_C z dx + x dy + y dz,
```
where ``C`` is the trace of the cylinder ``x^2 + y^2 = 1`` in the plane ``y + z = 2``. Orient ``C`` counterclockwise as viewed from above.

"""

# ╔═╡ 2111a86d-051a-406f-9cf5-934f4f35e2c3
md"""
## 9.16 Divergence Theorem
Recall that Green’s theorem (Tangent Form) can be written in vector notation as
```math
\oint_C \mathbf{F}\cdot d\mathbf{r}=\oint_C \mathbf{F}\cdot\mathbf{T} ds
=\iint_R (\text{curl }\mathbf{F})\cdot \mathbf{k} dA
```
that is, the line integral of the tangential component of ``\mathbf{F}`` is the double integral of the normal component of ``\text{curl }\mathbf{F}``.

We give now Green's theorem (normal form)
```math
\begin{array}{lcl}
\oint_C \mathbf{F}\cdot d\mathbf{r}&=&\oint_C \mathbf{F}\cdot\mathbf{n} ds\\
&=&\oint_C Pdy - Qdx \\
&=&\iint_R \left[\frac{\partial P}{\partial x}-\left(
-\frac{\partial Q}{\partial y}
\right)\right] dA\\
&=&\iint_R \left[\frac{\partial P}{\partial x}+\frac{\partial Q}{\partial y}\right] dA\\
&=&\iint_R (\text{div }\mathbf{F}) dA
\end{array}
```
"""

# ╔═╡ 05359c98-79e7-4c57-be02-eeda4e63d171
md"""
__Theorem__ (*Divergence Theorem*)

Let ``D`` be a closed and bounded region in 3-space with a piecewise-smooth boundary `S` that is oriented outward. Let 
```math
\mathbf{F}(x, y, z) = P(x, y, z) \mathbf{i} + Q(x, y, z) \mathbf{j} + R(x, y, z) \mathbf{k}
```
be a vector field for which ``P``, ``Q``, and ``R`` are continuous and have continuous first partial derivatives in a region of 3-space containing ``D``. Then
```math
\iint_S (\mathbf{F}\cdot n) dS = \iiint_D \text{div } \mathbf{F} dV.
```
"""

# ╔═╡ cff2af83-26f1-4cc7-b1c6-d0d9617496ad
md"""
__EXAMPLE__ (Verifying Divergence Theorem)

Let ``D`` be the region bounded by the hemisphere 
```math
x^2 + y^{\small 2} + (z - 1)^2 = 9, \quad 1 \leq z \leq 4, 
```
and the plane ``z = 1``. Verify the divergence theorem if 
```math
\mathbf{F} = x \mathbf{i} + y \mathbf{j} + (z - 1) \mathbf{k}.
```
"""

# ╔═╡ 59058a3f-6571-4090-90bf-b2cbac2dc709
md"""
__EXAMPLE__ (Using Divergence Theorem)

If 
```math
\mathbf{F} = xy \mathbf{i} + y^2z \mathbf{j} + z^3 \mathbf{k},
``` 
evaluate 
```math
\iint_S (F ⋅ n) dS,
``` 
where ``S`` is the unit cube defined by 
```math
0 \leq x \leq 1,\quad 0 \leq y \leq 1,\quad  0 \leq z \leq 1.
```
"""

# ╔═╡ 0cf4ab34-6479-43ba-aed4-c300ad7cf8a4
md"""
## 4.1 Definition of the Laplace Transform

__Definition__ (*Laplace Transform*)

Let ``f`` be a function defined for ``t \ge 0``. Then the integral
```math
\mathcal{L}\{f(t)\}=\int_0^{\infty}e^{-st}f(t) dt 
```
is said to be the __Laplace transform__ of ``f``, provided the integral converges.

__Notation__

```math
\mathcal{L}\{f (t)\} = F (s),\quad \mathcal{L}\{g(t)\} = G(s),\quad \mathcal{L}\{y(t)\} = Y(s), \quad \text{and} \quad \mathcal{L}\{H(t)\} = h(s).
```
"""

# ╔═╡ 7fb106db-664f-4139-9d29-96cd855a0605
md"""
__Examples__

Evaluate

1. ``\mathcal{L}\{1\}``.
1. ``\mathcal{L}\{t\}``.
1. ``\mathcal{L}\{e^{-3t}\}``.
1. ``\mathcal{L}\{\sin(2t)\}``.



"""

# ╔═╡ bfdd5cc2-8da9-4530-b44a-b2b45ffa119b
md"""
__``\mathcal{L}`` is a Linear Transform__ For a sum of functions, we can write
```math
\int_0^{\infty} e^{-st}\left[\alpha f(t)+\beta g(t)\right] dt = \alpha \int_0^{\infty} e^{-st} f(t)dt +\beta \int_0^{\infty} e^{-st} g(t) dt
```
whenever both integrals converge for ``s > c``. Hence it follows that
```math
\mathcal{L}\{\alpha f(t) + \beta g(t)\} = \alpha \mathcal{L}\{ f (t)\} + \beta\mathcal{L}\{g(t)\} = \alpha F (s) + \beta G(s).
```
Because of this property, ``\mathcal{L}`` is said to be a __linear transform__. 

"""

# ╔═╡ 776b5a94-b578-4fbd-8645-c5af41371b0a
md"""
__Theorem__ (*Transforms of Some Basic Functions*)
```math
\begin{array}{llcl}
\textbf{(a)} & \mathcal{L}\{1\} &=& \frac{1}{s}\\
\textbf{(b)} & \mathcal{L}\{t^n\} &=& \frac{n!}{s^{n+1}},\quad n=1,2,3,\cdots\\
\textbf{(c)} & \mathcal{L}\{e^{at}\} &=& \frac{1}{s-a}\\
\textbf{(d)} & \mathcal{L}\{\sin kt\} &=& \frac{k}{s^2+k^2}\\
\textbf{(e)} & \mathcal{L}\{\cos kt\} &=& \frac{s}{s^2+k^2}\\
\textbf{(f)} & \mathcal{L}\{\sinh kt\} &=& \frac{k}{s^2-k^2}\\
\textbf{(g)} & \mathcal{L}\{\cosh kt\} &=& \frac{s}{s^2-k^2}\\
\end{array}
```
"""

# ╔═╡ e39cb5f0-7110-4e98-af12-ed35a7b80f81
md"""
### Sufficient Conditions for Existence of ``\mathcal{L}\{ f (t )\}`` 

__Definition__ (*Exponential Order*)

A function ``f`` is said to be of __exponential order__ if there exist constants ``c, M > 0``, and ``T > 0`` such that 
```math
| f (t)| \leq Me^{ct} \quad \text{for all}\quad t > T.
```

__Theorem__ (*Sufficient Conditions for Existence*)

If ``f(t)`` is piecewise continuous on the interval ``[0, \infty)`` and of exponential order, then ``\mathcal{L}\{ f (t)\}`` exists for ``s > c``.
"""

# ╔═╡ 364a5da1-fb7f-4878-98b3-318935195882
md"""
__EXAMPLE__

Evaluate ``\mathcal{L}\{ f (t)\}`` for 
```math
f(t) =\left\{\begin{array}{lcl} 
0, &\text{  }&0 \leq t < 3 \\
2, &\text{  }& t \ge 3
\end{array}
\right.
```
"""


# ╔═╡ 86bf7675-87f8-4b4a-ba4d-99c3ba88c02e
md"""
## 4.2 The Inverse Transform and Transforms of Derivatives
If ``F(s)`` represents the __Laplace transform__ of a function ``f(t)``, 
that is, ``\mathcal{L}\{ f (t)\} = F (s)``, we then say ``f(t)`` is the __inverse Laplace transform__ of ``F(s)`` and write
```math
f(t) = \mathcal{L}^{-1}\{F(s)\}.
```

__Some Inverse Transforms__

```math
\begin{array}{lcl}
1 &=& \mathcal{L}^{-1}\left\{\frac{1}{s}\right\} \\
t^n &=& \mathcal{L}^{-1}\left\{\frac{n!}{s^{n+1}}\right\}, \quad n=1,2,3,\cdots \\
e^{at} &=& \mathcal{L}^{-1}\left\{\frac{1}{s-a}\right\} \\
\sin kt &=& \mathcal{L}^{-1}\left\{\frac{k}{s^2+k^2}\right\} \\
\cos kt &=& \mathcal{L}^{-1}\left\{\frac{s}{s^2+k^2}\right\} \\
\sinh kt &=& \mathcal{L}^{-1}\left\{\frac{k}{s^2-k^2}\right\} \\
\cosh kt &=& \mathcal{L}^{-1}\left\{\frac{s}{s^2-k^2}\right\} \\
\end{array}
```
"""

# ╔═╡ b71ec033-596d-4b1a-bfb1-aeb6e5e21914
md"""
Examples:
Evaluate
- ``\mathcal{L}^{-1}\left\{\frac{1}{s^5}\right\}``
- ``\mathcal{L}^{-1}\left\{\frac{1}{s^2+7}\right\}``
"""

# ╔═╡ e80855f9-d666-46cc-a96d-49fada0341cb
md"""
### ``\mathcal{L}^{-1}`` is a Linear Transform

```math
\mathcal{L}^{-1}\left\{\alpha F(s) +\beta G(s)\right\} =
\alpha \mathcal{L}^{-1}\left\{F(s)\right\} + \beta \mathcal{L}^{-1}\left\{G(s)\right\} 
```

"""

# ╔═╡ 9c5f7821-d77f-40ae-bcc9-e558727c034f
md"""
__Example__

Evaluate 
```math
\begin{array}{ll}
\text{✏}\quad &  \mathcal{L}^{-1}\left\{\frac{-2s+6}{s^2+4}\right\} \\
\text{✏}\quad & \mathcal{L}^{-1}\left\{\frac{s^2+6s+9}{(s-1)(s-2)(s+4)}\right\} 
\end{array}
```

"""

# ╔═╡ c41fae51-0d1b-4470-a2d9-507609b485f6
md"""
### Transforms of Derivatives
__Theorem__ (*Transform of a Derivative*)

If ``f, f', \cdots , f^{(n-1)}`` are continuous on ``[0, \infty)`` and are of exponential order and if ``f^{(n)}(t)`` is piecewise continuous on ``[0, \infty)``, then
```math
\mathcal{L}\{f^{(n)}(t)\} = s^n F(s) - s^{n-1}f(0) - s^{n-2} f'(0) - \cdots - f^{(n-1)}(0),
```
where ``F(s) = \mathcal{L}\{f (t)\}``.
"""

# ╔═╡ 295399f1-14b3-4cf4-83ab-349a92600337
HTML("""<h3>Solving Linear ODEs</h3>
<div style="color:darkblue;">The Laplace transform of a linear differential equation with constant coefficients becomes an algebraic equation in """, "Y(s)", """</div>
""")

# ╔═╡ 7b5c964d-c364-46b2-ab08-c817fb0bb08e
md"""
__EXAMPLEs__

- Use the Laplace transform to solve the initial-value problem
```math
\frac{dy}{dx} + 3y = 13\sin 2t, \quad y(0)=6
```
- Solve 
```math 
y''- 3y' + 2y = e^{-4t},\quad  y(0) =1,\quad y'(0) = 5.
```
"""

# ╔═╡ 3307e144-c20f-41a6-bba4-25b899079780
md"""
__Theorem__ (*Behavior of ``F(s)`` as ``s \to \infty``*)

If ``f`` is piecewise continuous on ``[0, \infty)`` and of exponential order, then 
```math
\lim_{s\to \infty} \mathcal{L}\{ f (t)\} = 0.
```
"""

# ╔═╡ 0c20d396-1a35-4aef-8c1b-866c9496b4af
begin
	t, s = symbols("t, s", real=true)
	Lp(f) = integrate(exp(-s*t)*f,(t,0,oo))
	Lp(sin(t))
	html""
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
Images = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlotThemes = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[compat]
FileIO = "~1.11.0"
HypertextLiteral = "~0.9.0"
Images = "~0.24.1"
LaTeXStrings = "~1.2.1"
PlotThemes = "~2.0.1"
Plots = "~1.20.0"
PlutoUI = "~0.7.1"
SymPy = "~1.0.52"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "485ee0867925449198280d4af84bdb46a2a404d0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.0.1"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[ArrayInterface]]
deps = ["IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "b6dec2ed4f10840e2cf836508525656450d4d289"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.1.26"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "a4d07a1c313392a77042855df46c5f534076fab9"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.0"

[[AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "d127d5e4d86c7680b20c35d40b503c74b9a39b5e"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.4"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c3598e525718abcc440f69cc6d5f60dda0a1b61e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.6+5"

[[CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "e2f47f6d8337369411569fd45ae5753ca10394c6"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.0+6"

[[CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "bdc0937269321858ab2a4f288486cb258b9a0af7"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.3.0"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "9995eb3977fbf67b86d0a0a0508e83017ded03f2"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.14.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "32a2b8af383f11cbb65803883837a149d10dfe8a"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.10.12"

[[ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "42a9b08d3f2f951c9b283ea427d96ed9f1f30343"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.5"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[CommonEq]]
git-tree-sha1 = "d1beba82ceee6dc0fce8cb6b80bf600bbde66381"
uuid = "3709ef60-1bee-4518-9f2f-acd86f176c50"
version = "0.2.0"

[[CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "344f143fa0ec67e47917848795ab19c6a455f32c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.32.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[Conda]]
deps = ["JSON", "VersionParsing"]
git-tree-sha1 = "299304989a5e6473d985212c28928899c74e9421"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.5.2"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "6d1c23e740a586955645500bbec662476204a52c"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.1"

[[CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[DataAPI]]
git-tree-sha1 = "ee400abb2298bd13bfc3df1c412ed228061a2385"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.7.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "7d9d316f04214f7efdbb6398d545446e246eff02"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.10"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "abe4ad222b26af3337262b8afb28fab8d215e9f8"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.3"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "a32185f5428d3986f47c2ab78b1f216d5e6cc96f"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.5"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "92d8f9f208637e8d2d28c664051a00569c01493d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.1.5+1"

[[EllipsisNotation]]
deps = ["ArrayInterface"]
git-tree-sha1 = "8041575f021cba5a099a456b4163c9a08b566a02"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.1.0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "LibVPX_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "3cc57ad0a213808473eafef4845a74766242e05f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.3.1+4"

[[FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "70a0cfd9b1c86b0209e38fbfe6d8231fd606eeaf"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.1"

[[FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "f985af3b9f4e278b1d24434cbb546d6092fca661"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.3"

[[FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3676abafff7e4ff07bbd2c42b3d8201f31653dcc"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.9+8"

[[FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "937c29268e405b6808d958a9ac41bfe1a31b08e7"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.11.0"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "35895cf184ceaab11fd778b4590144034a167a2f"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.1+14"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "cbd58c9deb1d304f5a245a0b7eb841a2560cfec6"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.1+5"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "dba1e8614e98949abfa60480b13653813d8f0157"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "182da592436e287758ded5be6e32c406de3a2e47"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.58.1"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d59e8320c2747553788e4fc42231489cc602fa50"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.58.1+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7bf67e9a481712b3dbe9cb3dac852dc4b1162e02"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+0"

[[Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "2c1cf4df419938ece72de17f368a021ee162762e"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "44e3b40da000eab4ccb1aecdc4801c040026aeb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.13"

[[HypertextLiteral]]
git-tree-sha1 = "72053798e1be56026b81d4e2682dbe58922e5ec9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.0"

[[IdentityRanges]]
deps = ["OffsetArrays"]
git-tree-sha1 = "be8fcd695c4da16a1d6d0cd213cb88090a150e3b"
uuid = "bbac6d45-d8f3-5730-bfe4-7a449cd117ca"
version = "0.3.1"

[[IfElse]]
git-tree-sha1 = "28e837ff3e7a6c3cdb252ce49fb412c8eb3caeef"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.0"

[[ImageAxes]]
deps = ["AxisArrays", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "794ad1d922c432082bc1aaa9fa8ffbd1fe74e621"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.9"

[[ImageContrastAdjustment]]
deps = ["ColorVectorSpace", "ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "2e6084db6cccab11fe0bc3e4130bd3d117092ed9"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.7"

[[ImageCore]]
deps = ["AbstractFFTs", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "db645f20b59f060d8cfae696bc9538d13fd86416"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.8.22"

[[ImageDistances]]
deps = ["ColorVectorSpace", "Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "6378c34a3c3a216235210d19b9f495ecfff2f85f"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.13"

[[ImageFiltering]]
deps = ["CatIndices", "ColorVectorSpace", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageCore", "LinearAlgebra", "OffsetArrays", "Requires", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "bf96839133212d3eff4a1c3a80c57abc7cfbf0ce"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.6.21"

[[ImageIO]]
deps = ["FileIO", "Netpbm", "OpenEXR", "PNGFiles", "TiffImages", "UUIDs"]
git-tree-sha1 = "13c826abd23931d909e4c5538643d9691f62a617"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.5.8"

[[ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils", "Libdl", "Pkg", "Random"]
git-tree-sha1 = "5bc1cb62e0c5f1005868358db0692c994c3a13c6"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.2.1"

[[ImageMagick_jll]]
deps = ["JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "1c0a2295cca535fabaf2029062912591e9b61987"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.10-12+3"

[[ImageMetadata]]
deps = ["AxisArrays", "ColorVectorSpace", "ImageAxes", "ImageCore", "IndirectArrays"]
git-tree-sha1 = "ae76038347dc4edcdb06b541595268fca65b6a42"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.5"

[[ImageMorphology]]
deps = ["ColorVectorSpace", "ImageCore", "LinearAlgebra", "TiledIteration"]
git-tree-sha1 = "68e7cbcd7dfaa3c2f74b0a8ab3066f5de8f2b71d"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.2.11"

[[ImageQualityIndexes]]
deps = ["ColorVectorSpace", "ImageCore", "ImageDistances", "ImageFiltering", "OffsetArrays", "Statistics"]
git-tree-sha1 = "1198f85fa2481a3bb94bf937495ba1916f12b533"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.2.2"

[[ImageShow]]
deps = ["Base64", "FileIO", "ImageCore", "OffsetArrays", "Requires", "StackViews"]
git-tree-sha1 = "832abfd709fa436a562db47fd8e81377f72b01f9"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.1"

[[ImageTransformations]]
deps = ["AxisAlgorithms", "ColorVectorSpace", "CoordinateTransformations", "IdentityRanges", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "e4cc551e4295a5c96545bb3083058c24b78d4cf0"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.8.13"

[[Images]]
deps = ["AxisArrays", "Base64", "ColorVectorSpace", "FileIO", "Graphics", "ImageAxes", "ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageShow", "ImageTransformations", "IndirectArrays", "OffsetArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "8b714d5e11c91a0d945717430ec20f9251af4bd2"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.24.1"

[[Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[IndirectArrays]]
git-tree-sha1 = "c2a145a145dc03a7620af1444e0264ef907bd44f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "0.5.1"

[[Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "61aa005707ea2cebf47c8d780da8dc9bc4e0c512"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.4"

[[IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[IterTools]]
git-tree-sha1 = "05110a2ab1fc5f932622ffea2a003221f4782c18"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.3.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "c7f1c695e06c01b95a67f0cd1d34994f3e7db104"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.2.1"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a4b12a1bd2ebade87891ab7e36fdbce582301a92"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.6"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[LibVPX_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "12ee7e23fa4d18361e7c2cde8f8337d4c3101bc7"
uuid = "dd192d2f-8180-539f-9fb4-cc70b1dcf69a"
version = "1.10.0+0"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "761a393aeccd6aa92ec3515e428c26bf99575b3b"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+0"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "5455aef09b40e5020e1520f551fa3135040d4ed0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2021.1.1+2"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "0fb723cd8c45858c22169b2e42269e53271a6df7"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.7"

[[MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "4ea90bd5d3985ae1f9a908bd4500ae88921c5ce7"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.0"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[Netpbm]]
deps = ["ColorVectorSpace", "FileIO", "ImageCore"]
git-tree-sha1 = "09589171688f0039f13ebe0fdcc7288f50228b52"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.1"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "c0f4a4836e5f3e0763243b8324200af6d0e0f90c"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.5"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "520e28d4026d16dcf7b8c8140a3041f0e20a9ca8"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.7"

[[PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "646eed6f6a5d8df6708f15ea7e02a7a2c4fe4800"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.10"

[[Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "2276ac65f1e236e0a6ea70baff3f62ad4c625345"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.2"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "477bf42b4d1496b454c10cce46645bb5b8a0cf2c"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.2"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "501c20a63a34ac1d015d5304da0e645f42d91c9f"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.11"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs"]
git-tree-sha1 = "e39bea10478c6aff5495ab522517fae5134b40e3"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.20.0"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "Logging", "Markdown", "Random", "Suppressor"]
git-tree-sha1 = "45ce174d36d3931cd4e37a47f93e07d1455f038d"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.1"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "afadeba63d90ff223a6a48d2009434ecee2ec9e8"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.1"

[[PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "169bb8ea6b1b143c5cf57df6d34d022a7b60c6db"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.92.3"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[Ratios]]
deps = ["Requires"]
git-tree-sha1 = "7dff99fbc740e2f8228c6878e2aad6d7c2678098"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.1"

[[RecipesBase]]
git-tree-sha1 = "b3fb709f3c97bfc6e948be68beeecb55a0b340ae"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.1"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "2a7a2469ed5d94a98dea0e85c46fa653d76be0cd"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.3.4"

[[Reexport]]
deps = ["Pkg"]
git-tree-sha1 = "7b1d07f411bc8ddb7977ec7f377b97b158514fe0"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "0.2.0"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[Rotations]]
deps = ["LinearAlgebra", "StaticArrays", "Statistics"]
git-tree-sha1 = "2ed8d8a16d703f900168822d83699b8c3c1a5cd8"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.0.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["OpenSpecFun_jll"]
git-tree-sha1 = "d8d8b8a9f4119829410ecd706da4cc8594a1e020"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "0.10.3"

[[StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[Static]]
deps = ["IfElse"]
git-tree-sha1 = "62701892d172a2fa41a1f829f66d2b0db94a9a63"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.3.0"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3240808c6d463ac46f1c1cd7638375cd22abbccb"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.12"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "1958272568dc176a1d881acb797beb909c785510"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.0.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "fed1ec1e65749c4d96fc20dd13bea72b55457e62"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.9"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "000e168f5cc9aded17b6999a560b7c11dda69095"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.0"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[SymPy]]
deps = ["CommonEq", "CommonSolve", "LinearAlgebra", "Markdown", "PyCall", "RecipesBase", "SpecialFunctions"]
git-tree-sha1 = "1ef257ecbcab8058595a68ca36a6844b41babcbd"
uuid = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
version = "1.0.52"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "d0c690d37c73aeb5ca063056283fde5585a41710"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.5.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TiffImages]]
deps = ["ColorTypes", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "OffsetArrays", "OrderedCollections", "PkgVersion", "ProgressMeter"]
git-tree-sha1 = "03fb246ac6e6b7cb7abac3b3302447d55b43270e"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.4.1"

[[TiledIteration]]
deps = ["OffsetArrays"]
git-tree-sha1 = "52c5f816857bfb3291c7d25420b1f4aca0a74d18"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.3.0"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[VersionParsing]]
git-tree-sha1 = "80229be1f670524750d905f8fc8148e5a8c4537f"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.2.0"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll"]
git-tree-sha1 = "2839f1c1296940218e35df0bbb220f2a79686670"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.18.0+4"

[[WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "59e2ad8fd1591ea019a5259bd012d7aee15f995c"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.3"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "acc685bcf777b2202a904cdcb49ad34c2fa1880c"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.14.0+4"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7a5780a0d9c6864184b3a2eeeb833a0c871f00ab"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "0.1.6+4"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d713c1ce4deac133e3334ee12f4adff07f81778f"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2020.7.14+2"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "487da2f8f2f0c8ee0e83f39d13037d6bbf0a45ab"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.0.0+3"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─5928498a-cee5-4e32-a2c6-d46d3b9dcca4
# ╟─77f286eb-1757-4044-86dd-6a550affea75
# ╟─79a8d881-0e07-44b1-a4cd-783b97db124e
# ╟─b691771d-5ce5-4937-a249-cb01da6243ac
# ╟─30c6cfaa-abca-4099-8bee-c6cc582ef3b2
# ╟─2d9f3647-327d-4b3f-b002-a4e1e9dbd59c
# ╟─5437877a-5788-463c-97b3-3a203b878917
# ╟─90db9f5e-908a-45ca-acb5-7fded8d28183
# ╟─cadb7973-85ae-431e-9a61-448a2243b72e
# ╟─d1013f55-85ea-4e0b-9499-830511054374
# ╟─e8cb1b44-72df-4544-8c60-fa9413699391
# ╟─ca122aa2-cc1a-45c5-bf96-024eeca7b34e
# ╟─9d1bddd0-45bb-4789-a388-46719159f497
# ╟─fa1b0c79-0c40-4e88-ba65-e2565fc9b604
# ╟─0f489937-3266-4d83-880d-c4ec0ebd2d34
# ╟─392780ee-6c31-42fb-b311-a7e7df34b1dc
# ╟─751bff35-70ab-4d5a-909c-692157146226
# ╟─509d1cbe-f592-4768-a29c-03881e34441c
# ╟─506abf7e-04b1-4220-9e10-bc9e01383ec3
# ╟─3e2e4f04-b687-47e7-96d8-6774724a9207
# ╟─330410cd-40e8-4398-938d-ee5e4be6f650
# ╟─b7452e63-29ae-4c5e-87fa-5faecedf986c
# ╟─0de4e7c3-216a-432a-b998-493a963901cf
# ╟─8451fe29-cc80-426e-bc65-43ca5eafb5d3
# ╟─4844492f-27b9-4350-bbd6-66d97b87e909
# ╟─a3f94e01-c503-40f7-a4b2-d40d25401fe3
# ╟─52f1f384-ba7f-4a1d-b1fb-311a10ea21a9
# ╟─f185e1e7-fc9e-43d2-9f41-43be18c444c3
# ╟─c6ffb2c9-3b92-4032-ab96-1804b11290f5
# ╟─62afd674-74c5-429d-9d21-1673fdb7b385
# ╟─4c4bd58e-04d5-4da0-8943-6d128a6b749f
# ╟─db44a1a7-705d-4fc9-a16b-cd3032353b3a
# ╟─17ab5bf0-65b0-4079-99a0-ccd19721d5c6
# ╟─07f6f180-ac84-4952-a24f-234990c35191
# ╟─b942c102-5698-4086-8056-a053b0ca9d43
# ╟─f84606c4-3de8-4156-a4ea-3df383974de5
# ╟─f389ba80-79ea-45a8-a8cb-3fdbdcfb1971
# ╟─c3222e2a-cf12-42ba-bc9a-0f914dbfb529
# ╟─6363ee6a-df6d-4150-b4c3-5e959bdf8c6c
# ╟─60db1dcb-ef1e-4a4b-b4ad-e2d5415601d8
# ╟─ca86c2c9-7416-4d45-bcec-a644c151a43f
# ╟─efd8e89b-6f7d-44d2-b798-5e0c06fd3622
# ╟─fba21f02-4aa0-46a1-b900-386647c20450
# ╟─ed3590a1-b43a-4212-81a8-bed8b21f0858
# ╟─0615d26b-858d-4a03-8595-e85da72adca7
# ╟─36380af4-f8de-4b7e-97c8-65673754ff76
# ╟─5febf63a-7834-4248-8884-7be76bb08202
# ╟─d1b8b4eb-e856-4745-9b83-808e2c55b49f
# ╟─b0a55996-3e23-4696-9a5c-94284b982155
# ╟─02561b84-2447-4462-91e9-4872c44e1ebe
# ╟─a0f19e1e-2ee6-440c-8c39-9e1f0d07b2d5
# ╟─9badddfe-d618-421e-89c9-e54639b94423
# ╟─b0a38158-bfdb-4b6b-bf70-7d6405da0c2b
# ╟─506a2816-5d1f-4792-93ed-cbf526888897
# ╟─6bc1fefb-be6c-4e38-9299-13a632d2ee34
# ╟─de3153ea-e97d-4fda-b70f-bd23b49cd8b7
# ╟─4958f085-527e-4ea7-b004-82d1317e1ed7
# ╟─5aa87065-e94c-425a-bf24-d197aaacafd8
# ╟─ea5a1ed0-4494-4850-8943-18d8029cf356
# ╠═00f35bc7-5456-4a97-a474-166bde07dea8
# ╟─a53e17e1-bb2d-4d41-b8e0-1df381601af5
# ╟─fb03d35d-9a94-4824-bb95-5c48a9a6ad7f
# ╟─7fbbe932-5b35-4fa5-a649-46cba8b5df47
# ╟─5ab19faa-f029-48d4-b9d3-7cef666acd1e
# ╟─4fdfad2e-e860-4358-b888-b6642fa0b6c3
# ╟─3aa6fe06-ae82-48dc-9999-877533fa21ae
# ╟─b4305759-a665-429a-a888-fadf5280b6ba
# ╟─dd993ad7-758c-4cee-93ad-f9d62277e90b
# ╟─1e20d2a6-bea6-4c16-862c-eaec3bb99d9a
# ╟─ffd52892-0427-49b2-8b1e-891d64c6c688
# ╟─71e53648-d66d-400a-9f5c-1ee9a9b8b9c5
# ╟─b9c77a2e-c384-4313-bea8-2694ab4dfc56
# ╟─e997b0f3-98de-4fd6-8ea5-e4b3440f0f5a
# ╟─bbe13ec2-af6b-40fa-b4f7-74a68d59298d
# ╟─45aac5bb-e392-4eb0-8a5e-976fa6f17938
# ╟─904f1f79-4866-4837-b128-c8a30567fa76
# ╟─ce1a2204-adfd-4862-a5fc-8637dcc1ce9d
# ╟─81441e8f-253c-44a3-b160-c0ac350791ae
# ╟─ea174feb-97f8-497d-8a47-2fcdfeef29c8
# ╟─47751cc9-f4e0-400e-8b6e-5c37c954ce96
# ╟─a89ccf32-03b9-45c3-be77-7029ee8fdfee
# ╟─a0d30cd9-2fbf-45be-b20d-da3619755922
# ╟─b3f1155f-fc04-4a35-ab00-d9f7851c5b12
# ╟─2111a86d-051a-406f-9cf5-934f4f35e2c3
# ╟─05359c98-79e7-4c57-be02-eeda4e63d171
# ╟─cff2af83-26f1-4cc7-b1c6-d0d9617496ad
# ╟─59058a3f-6571-4090-90bf-b2cbac2dc709
# ╟─0cf4ab34-6479-43ba-aed4-c300ad7cf8a4
# ╟─7fb106db-664f-4139-9d29-96cd855a0605
# ╟─bfdd5cc2-8da9-4530-b44a-b2b45ffa119b
# ╟─776b5a94-b578-4fbd-8645-c5af41371b0a
# ╟─e39cb5f0-7110-4e98-af12-ed35a7b80f81
# ╟─364a5da1-fb7f-4878-98b3-318935195882
# ╟─86bf7675-87f8-4b4a-ba4d-99c3ba88c02e
# ╟─b71ec033-596d-4b1a-bfb1-aeb6e5e21914
# ╟─e80855f9-d666-46cc-a96d-49fada0341cb
# ╟─9c5f7821-d77f-40ae-bcc9-e558727c034f
# ╟─c41fae51-0d1b-4470-a2d9-507609b485f6
# ╟─295399f1-14b3-4cf4-83ab-349a92600337
# ╟─7b5c964d-c364-46b2-ab08-c817fb0bb08e
# ╠═3307e144-c20f-41a6-bba4-25b899079780
# ╟─0c20d396-1a35-4aef-8c1b-866c9496b4af
# ╠═5867632c-fff5-11eb-3a19-2f309efd424a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
