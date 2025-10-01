# Simulating a Wormhole

This page is a compilation of my notes and derivations for rendering, navigating around, and traversing a physically accurate Morrise-Thorn wormhole.

The wormhole implemented here was first described in 1987 in this paper:
[American Journal of Physics Volume 56 issue 5 1988 [doi 10.1119%2F1.15620] Morris, Michael S. -- Wormholes in spacetime and their use for interstellar travel- A tool for teaching general relativity.pdf (PDFy mirror) : Free Download, Borrow, and Streaming : Internet Archive](https://archive.org/details/pdfy--5NLGQyfpB61dmyG/mode/2up)

We will start by a derivation of the null geodesic equation and parallel transport equations in polar coordinates as well as a description of an implementation for a 2D observer. Once that is covered, we will go over a full 3D derivation and implementation.

## Tensor Calculus
Lets do a quick overview of some important notes regarding tensors.

$$
\begin{array}{|c|c|c|c|}
\hline
\textbf{Expression} & \textbf{Description} & \textbf{Tensor Type} & \textbf{Result Type} \\
\hline
v^\mu & \text{Contravariant vector (tangent)} & (1,0) & \text{Vector} \\
v_\mu & \text{Covariant vector (covector)} & (0,1) & \text{Covector} \\
g_{\mu\nu} & \text{Metric tensor} & (0,2) & \text{Covariant 2-tensor} \\
g^{\mu\nu} & \text{Inverse metric tensor} & (2,0) & \text{Contravariant 2-tensor} \\
M^\mu_{\ \nu} & \text{Linear map / mixed tensor} & (1,1) & \text{Transformation} \\
e_\mu = \partial_\mu & \text{Coordinate basis vector} & (1,0) & \text{Tangent basis} \\
dx^\mu & \text{Dual basis covector} & (0,1) & \text{Cotangent basis} \\
J^\mu_{\ \bar{\nu}} = \partial x^\mu / \partial \bar{x}^{\bar{\nu}} & \text{Jacobian (change of coordinates)} & (1,1) & \text{Basis transform} \\
J^{\bar{\mu}}_{\ \nu} = \partial \bar{x}^{\bar{\mu}} / \partial x^\nu & \text{Inverse Jacobian} & (1,1) & \text{Basis transform} \\
e^\mu_a & \text{Tetrad (global to local)} & (1,1) & \text{Change of basis} \\
e^a_\mu & \text{Inverse tetrad (local to global)} & (1,1) & \text{Change of basis} \\
\hline
u^\mu a^\nu & \text{Outer product of two vectors} & (1,0) \otimes (1,0) & (2,0)\ \text{Tensor} \\
u_\mu a^\nu & \text{Mixed outer product} & (0,1) \otimes (1,0) & (1,1)\ \text{Tensor} \\
u^\mu a_\nu & \text{Mixed outer product} & (1,0) \otimes (0,1) & (1,1)\ \text{Tensor} \\
u_\mu a^\mu & \text{Inner product (index contraction)} & \text{Scalar} & (0,0)\ \text{Scalar} \\
g_{\mu\nu} a^\nu & \text{Lower index of } a^\nu & (1,0) \otimes (0,2) & (0,1)\ \text{Covector} \\
g^{\mu\nu} a_\nu & \text{Raise index of } a_\nu & (0,1) \otimes (2,0) & (1,0)\ \text{Vector} \\
\hline
\end{array}
$$


## Polar Coordinates

In order to leverage the symetry of a wormhole we will be working in **Polar Coordinates**.
The 2D implementation of our wormhole will be using this coordinate system: $(r, \theta )$. The parametric equations are shown below.

$$\begin{align*}
x &= r \cos\theta \\
y &= r \sin\theta
\end{align*}$$

## Wormhole Geometry

The geometry of the wormhole is described using a very simple expression that is only in terms of the radial component. $$r(l) = \sqrt{ l^2 + r_0^2 } \quad \text{where } r_0 \text{ is the radius of the throat}$$

This is called the **shape function**. In this setup, the wormhole geometry is described using the coordinates $( l, \theta )$, where $l$ is the proper radial distance. The standard polar coordinates $(r, \theta)$ are modified so that the radial coordinate is now a function of $l$, i.e., $r = r(l)$
This ensures the coordinate transformation depends on the geometry of the wormhole. Unlike standard spherical coordinates where $r$ is independent, here $r$ is a function of the proper radial distance $l$, which encodes the spatial curvature of the wormhole throat and flare-out shape. We will call this new coordinate system **wormhole coordinates** $(r(l), \theta)$.


## Metric Tensor 

The first thing we must do is define is our Jacobian. The Jacobian is a Tensor Field that lets us transform vectors at some point P from cartesian $(x,y)$ to our wormhole coordinates $(r(l),\theta)$.

We first define our transformation:

$$
\begin{align*}
x &= r(l) \cos\theta \\
y &= r(l) \sin\theta
\end{align*}
$$

Then we derive the Jacobian. To compute $g_{uv}$ for our modified Polar coordinates we need the Jacobian that maps to Cartesian coordinates.

$$J: (l, \theta) \mapsto (x, y) = 
\begin{bmatrix}
{\frac{\partial x}{\partial l}} & {\frac{\partial x}{\partial \theta}} \\
{\frac{\partial y}{\partial l}} & {\frac{\partial y}{\partial \theta}} \\
\end{bmatrix}$$

$$\begin{align*}
{\frac{\partial x}{\partial l}} = r'(l)cos{\theta}     &\quad\quad {\frac{\partial x}{\partial \theta}} = -r(l)sin{\theta} \\
{\frac{\partial y}{\partial l}} = r'(l)sin{\theta}     &\quad\quad {\frac{\partial y}{\partial \theta}} = r(l)cos{\theta}
\end{align*}$$

And now we derive the Metric Tensor $g_{uv}$.
The **metric tensor** is a rank-2 tensor that defines how distances and angles are measured on a curved space or manifold. It generalizes the dot product to arbitrary coordinates and geometries, allowing us to compute quantities like lengths, angles, and volumes.
The metric $ds^2 = g_{\mu\nu}dx^{\mu}dx^{\upsilon}$ takes two vectors and returns a scalar: their inner product. It effectively tells us **how the basis vectors themselves are stretched, skewed, or rotated** at each point in space.

$$g_{uv} = J^TJ = 
\begin{bmatrix}
{\frac{\partial x}{\partial l}}^2 {\frac{\partial y}{\partial l}}^2 &
{\frac{\partial x}{\partial \theta}} {\frac{\partial x}{\partial l}} +  {\frac{\partial y}{\partial \theta}} {\frac{\partial y}{\partial l}}\\
{\frac{\partial x}{\partial \theta}} {\frac{\partial x}{\partial l}} +  {\frac{\partial y}{\partial \theta}} {\frac{\partial y}{\partial l}} &
{\frac{\partial y}{\partial l}}^2 {\frac{\partial x}{\partial l}}^2
\end{bmatrix} \\ = 
\begin{bmatrix}
r'(l)cos^2(\theta) + r'(l)sin^2(\theta) & -r(l)cos{\theta}sin{\theta} + r(l)cos{\theta}sin{\theta}\\
-r(l)cos{\theta}sin{\theta} + r(l)cos{\theta}sin{\theta} & r'(l)r(l)^{2}cos^2(\theta) + r'(l)r(l)^{2}sin^2(\theta)
\end{bmatrix} =
\begin{bmatrix}
r'(l) & 0 \\
0 & r'(l)^2r(l)^2
\end{bmatrix} $$

But wait! We defined $l$ as **proper radial distance**. It is defined such that:
$$dl = dr = 1$$

$$\dot{.\hspace{.05in}.} \quad g =\begin{bmatrix}
1 & 0 \\
0 & r(l)^2
\end{bmatrix} $$

I still dont fully understand the justification for dropping $r(l)$ from the components of $g_{uv}$ considering its in our jacobian, tbh.

## Christoffel Symbols

Now that we have a **metric tensor** we can describe how vectors stretch and rotate at any point in space. But this is not the full picture. If we want to traverse the curved space our wormhole produces, we need a way of tracking how our basis vectors change as we travel. This is where **Christoffel Symbols** (Connection Coefficients) come in. They are not tensors because they don't transform like tensors, but they are rich mathematical objects that describe how the coordinate basis vectors themselves change from point to point in a curved space. Unfortunately, they are quite tedious to calculate.

$$\Gamma^\mu_{\alpha\beta} = \frac{1}{2}\mathfrak{g}^{\mu\lambda}(\frac{\partial g_{\lambda\alpha}}{\partial x^\beta} + \frac{\partial g_{\lambda\beta}}{\partial x^\alpha} + \frac{\partial g_{\alpha\beta}}{\partial x^\lambda})$$

$$\text{where} \quad \mathfrak{g}=g^{-1}$$

For a 2D system like ours, there are a total of 6 christoffel symbols. Because our metric tensor is a diagonal matrix we can assume many of them will be zero.

$$\begin{align*}
&\Gamma^l_{ll} = \frac{1}{2}\cdot1(\frac{\partial}{\partial l}1 + \frac{\partial}{\partial l}1 - \frac{\partial}{\partial l}1) = 0 \\
&\Gamma^l_{\theta l} = \frac{1}{2}\cdot1(0 + 0- 0) = 0 = \Gamma^l_{l\theta} \\
&\Gamma^l_{\theta \theta} = \frac{1}{2}\cdot1( 0 + 0 - 2r(l)r'(l)) = -r(l)r'(l) \\
&\Gamma^\theta_{l l} = \frac{1}{2}\cdot \frac{1}{r(l)^2}(0 + 0- 0) = 0 \\
&\Gamma^\theta_{\theta l} = \frac{1}{2}\cdot \frac{1}{r(l)^2}(2r(l)r'(l) + 0- 0) = \frac{r(l)}{r'(l)} = &\Gamma^\theta_{l \theta} \\
&\Gamma^\theta_{\theta \theta} = \frac{1}{2}\cdot \frac{1}{r(l)^2}( 0 + 0 - 0 ) = 0 \\
\end{align*}$$

We have a total of three non-zero christoffel symbols:

$$\begin{align*}
&\Gamma^l_{\theta \theta} = -r(l)r'(l) \\
&\Gamma^\theta_{\theta l} = \Gamma^\theta_{l \theta} = \frac{r(l)}{r'(l)}
\end{align*}$$

## Null Geodesics
We now have everything we need to compute our geodesics! But first, lets describe what they are.
**Geodesics** are the "straightest possible" paths through curved spacetime — meaning a particle or ray moving along one experiences **no proper acceleration**. In other words, geodesics are the paths that extremize the interval $ds^2=0$, and they are solutions to the geodesic equation involving the metric and Christoffel symbols.

A **null geodesic** is a special type of geodesic where the spacetime interval is zero. These paths are followed by **massless particles** like light. Null geodesics trace how light propagates through spacetime, bending in response to curvature caused by mass and energy.

$$\frac{d^2 x^\alpha}{d\lambda^2 } + \Gamma^\alpha_{\beta\gamma} \frac{dx^\beta}{d\lambda} \frac{dx^\gamma}{d\lambda} = 0$$

$$\begin{align*}
&\frac{d^2l}{d\lambda^2} = r(l)(\frac{d\theta}{d\lambda})^2 = r(l)(\theta')^2 \\
&\frac{d^2 \theta}{d\lambda^2} = -\frac{2}{r(l)}\frac{dl}{d\lambda}\frac{d\theta}{d\lambda} = -\frac{2}{r(l)}r'(l)\theta'
\end{align*}$$

## Implementation
Our null geodesics are a system of 2nd order ordinary differential equations. Solving these analytically as very difficult and in many cases impossible. Thankfully we don't need a general solution to our geodesic equation. We can use a technique called **numerical integration**. The idea is simple, instead of trying to find a parametric equation where we can calculate the coordinate at param 500 along a path, we can just take 500 tiny steps along the path.
The simplest form of numerical integration is called **Eulers method**:

$$\begin{align*}
x_{n+1} = x_n + v_n \cdot \Delta \lambda \\
v_{n+1} = v_n + a_n \cdot \Delta \lambda
\end{align*}$$

It updates **position and velocity** using current values at time step $n$. The drawbacks is that its not energy conserving and it accumulates error quickly.

**Symplectic Euler** is a small but important variation of Euler’s method. It updates **velocity first**, then use the **new velocity** to update position.

$$\begin{align*}
v_{n+1} = v_n + a_n \cdot \Delta \lambda \\
x_{n+1} = x_n + v_{n+1} \cdot \Delta \lambda
\end{align*}$$

The last numerical integrator we will look at is **Runge–Kutta method**.
RK4 improves accuracy by sampling the derivative **at multiple points within each time step**, rather than just once at the beginning like Euler's method.

Given an ODE;
$$\frac{dx}{d\lambda} = f(\lambda,x), \quad y(\lambda_0)=y_0$$

Let $h = \Delta\lambda$ (time step):

$$\begin{align*}
&k_1 = f(\lambda_n, x_n) \\
&k_2 = f(\lambda_n + \frac{h}{2}, x_n + \frac{h}{2}k_1) \\
&k_3 = f(\lambda_n + \frac{h}{2}, x_n + \frac{h}{2}k_2) \\
&k_3 = f(\lambda_n + h, x_n + hk_3)
\end{align*}$$

$$x_{n+1} = x_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$

Each $k_n$​ is a derivative estimate. RK4 combines them to get a better average slope for the whole step.

We will be using either **Symplectic Euler** or **RK4** to solve our null geodesics (and other systems of ODEs) depending on the situation.

Once a ray has been marched remember to convert from our new polar coordinate system $(l,\theta)$ to the classical polar coordinate $(r,\theta)$ by applying our shape function to $l$.

## Tracking an Observer
We now have everything we need to trace rays through curved 2D space. A ray starting at $(l,\theta)$ that goes through the throat of the wormhole and end up at $(-l,\theta)$. This negative radial component implies the 3D space the 2D wormhole is embedded in.
But this is not going to be enough to traverse the wormhole. We are missing a crucial ingredient: **time**.
Even in a static simulation where the camera is being moved by fixed increments every frame, time is still essential to properly track observer position and orientation around and through the wormhole.

Our simulation can be split into two phases: Initialization and Runtime:
- During **Initialization** we set the initial position and orientation of the observer.
- During **Runtime** we track mouse movements to update observer orientation and key presses to offset the observer along one of the vectors of the local tangent frame. Once we offset the observer, his tangent frame must be updated for the new position in curved space he occupies.

## Tetrads
Moving forward we are going to be working with space time metrics, and positions and vectors with a time component. In our 2D polar implementation this means that the metric tensor gets upgraded to a 3x3 matrix.

$$g_{\mu\nu} =
\begin{bmatrix}
-1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & r(l)^2
\end{bmatrix} $$

The time component will always be the first column and negative. This is mainly to distinguish the temporal component from the spatial ones.

One of the core pillars of General Relativity is the notion that space is locally flat (Minkowskian) and globally curved. 

$$\eta_{\mu\nu} =
\begin{bmatrix}
-1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1
\end{bmatrix} $$

The **Minkowskian Metric** describes local flat space. It is what we will use to store the observer's orientation driven by the user (mouse move).

The observer's local tangent frame must be transformed to this globally curved space in order to be offset along the geodesic. This is where **tetrads** come in.
Tetrads will transform vectors from locally flat cartesian space to our globally curved wormhole space and back.

$$g_{\mu\nu} = \eta_{\alpha\beta}e^\alpha_\mu e^\beta_\nu$$

Or in matrix form:

$$g = e^T \cdot \eta \cdot e$$

Actually, building a tetrad is quite easy. First we observer that spatially it performs a similar function to our Jacobian from earlier. However, the Jacobian describe the differentials along each basis vector. It is not a tangent basis in it of itself. It must be orthonormalized into order to be a valid tangent frame. Naive orthonormalization with cross products may work in flat Euclidean space, but we have a curved metric to respect! We need to use **Relativistic Gramm-Schmidt** to orthonormalize our Jacobian.
I recomment reviewing regular Gramm-Schmidt from linear algebra.

$$\begin{align*}
\\
proj_{u^\mu}(v^\mu) &= \frac{g_{\mu\nu} v^\mu u^\nu}{g_{\alpha\beta} u^\alpha u^\beta} \equiv \frac{v^\mu u_\nu}{u^\alpha u_\alpha} \\
normalise(u^\mu) &= \frac{u^\mu}{g_{\alpha\beta} u^\alpha u^\beta} \equiv \frac{u^\mu}{u^\alpha u_\alpha}\\
\\
u_0 &= v_0\\
u_1 &= v_1 - proj_{u_0}(v_1) \\
\\
u_2 &= v_2 - proj_{u_0}(v_2) \\
u_2 &= u_2 - proj_{u_1}(u_2) \\
\\
u_3 &= v_3 - proj_{u_0}(v_3) \\
u_3 &= u_3 - proj_{u_1}(u_3) \\
u_3 &= u_3 - proj_{u_2}(u_3) \\
\\
u_0 &= normalise(u_0)\\
u_1 &= normalise(u_1)\\
u_2 &= normalise(u_2)\\
u_3 &= normalise(u_3)\\
\end{align*}$$

```c++	
	void BuildInitialOrientedTetrad()
	{
		//conver cartesian basis to polar basis. Here we use invverse Polar->Cartesian
		//Jacobian so we can preserver the negative r component if its present.
		Mat2 mInvJacobian = PolarToCartesianJacobianCurved( m_vPolarPosition );
		Mat2 mJacobian;
		mInvJacobian.Inverse( &mJacobian );

		//orthonormalize with respect to metric
		Mat3 mMetric = GetMetricTensor( m_vPolarPosition );
		Vec3 e0( 1.0f, 0.0f, 0.0f );
		Vec3 e1 = Vec3( 0.0f, mJacobian.GetColVec( 0 ) );
		Vec3 e2 = Vec3( 0.0f, mJacobian.GetColVec( 1 ) );
		m_mTetrad = OrthonormalizeWithMetric( e0, e1, e2, mMetric );
	}
```

## Parallel Transport
So now we have an oriented tetrad for the starting position of our observer that will transport his cartesian orientation to a globally curved spacetime. How do we move the tetrad with the observer?
Easy! We slide the columns of the tetrad with observer and re-orthonormalize with respect to the metric at the new location!

This operation is called **parallel transport** and it is captured by a system of first order ODEs. Thes must be numerically integrated just like our geodesic equations.

$$\begin{align*}
\frac{dv^\mu}{d\lambda} = -\Gamma^\mu_{\nu\rho} v^\nu \frac{dx^\rho}{d\lambda} \\
\quad dv = \text{geodesic velocity} \\
\quad v = \text{pre-transport vec} \\
\quad x = \text{position derivative} \\
\end{align*}$$

Recall the **Christoffel Symbols** we calculated earlier.

$$\begin{align*}
&\Gamma^l_{\theta \theta} = -r(l)r'(l) \\
&\Gamma^\theta_{\theta l} = \Gamma^\theta_{l \theta} = \frac{r(l)}{r'(l)}
\end{align*}$$

$$\begin{align*}
&\frac{dv^l}{d\lambda} = - \Gamma^l_{\theta\theta} v^\theta \theta' = r'(l)r(l) \cdot v^\theta \theta' \\
&\frac{dv^\theta}{d\lambda} = - \Gamma^\theta_{\theta l} v^\theta r'(l) - \Gamma^\theta_{l \theta} v^l \theta' = -\frac{r'(l)}{r(l)}(r'(l) v^\theta + \theta' v^l)
\end{align*}$$

```c++
	void TransportTetrad( const Vec2& vPolarStart, const Vec2& vPolarDir, const Vec2& vPolarEnd, const float fOffset )
	{
		// Extract basis vectors
		Vec3 e0 = m_mTetrad.GetColVec( 0 );
		Vec3 e1 = m_mTetrad.GetColVec( 1 );
		Vec3 e2 = m_mTetrad.GetColVec( 2 );

		// Parallel transport all 3 vectors
		ParallelTransportPolar( vPolarStart, vPolarDir, fOffset, e0 );
		ParallelTransportPolar( vPolarStart, vPolarDir, fOffset, e1 );
		ParallelTransportPolar( vPolarStart, vPolarDir, fOffset, e2 );

		// Re-orthonormalize with respect to the metric at the new location
		const Mat3 mMetric = GetMetricTensor( vPolarEnd );
		m_mTetrad = OrthonormalizeWithMetric( e0, e1, e2, mMetric );
	}
```

## Covariant Derivative
Before we move on, lets dig a little deaper about **Parallel Transport**.

The **covariant derivative** allows us to take derivatives of tensor fields in curved space. This ammounts to the partial derivatives of a vector filed with a **Connection Correction** provided by our christoffel symbols.

Let $\gamma(\tau)$ be a timelike curve parameterized by proper time $\tau$. A vector is being parallel transported along $\gamma$ if the covariant derivative is 0.

$$
\frac{Dv^\mu}{d\lambda} = 0
$$

$$
\text{where} \quad \frac{Dv^\mu}{d\lambda} = \frac{dv^\mu}{d\lambda} + \Gamma^\mu_{uv} \frac{dx^u}{d\lambda} v^v
$$

So to get the parallel transport equations we solve for the directional derivative $\frac{dv^\mu}{d\lambda}$.

But we can also get the **Null Geodesic Equations** from this.
The parallel transport equations arise from setting the covariant derivative of some vector $V^u$ on a curve $\gamma(\tau)$ to 0. For null geodesics we do this for velocity $u^\mu$ whic is already a derivative of position. So this covariant derivative is the *second directional derivative* + *connection correction* terms. Like with prallel transport, we solve fo the directional derivative.

$$
\frac{Du^\mu}{d\lambda} = \frac{d^2x^\mu}{d\lambda^2} + \Gamma^\mu_{uv} \frac{dx^u}{d\lambda} \frac{dx^v}{d\lambda} = 0
$$

## 3D Implementation
In principal, bringing this to 3D is a matter of rederiving everything in spherical coordinates. While true, there are quite a few gotchas to avoid due to some quirks about the spherical coordinate system and the extra degree of freedom.

The first think we decide is which coordinate system we are using.
In cartesian space (for rendering in OpenGL) we will be using a right-handed y-up coordinate system. We must keep this in mind during our derivation.

$$\begin{align*}
&x = r(l)sin{\theta}cos(\phi \\
&y = r(l)cos{\theta} \\
&z = r(l)sin{\theta}sin{\phi}
\end{align*}$$

Here, $\theta$ is the polar angle ranging from $[0,\pi]$, and $\phi$ is the azimuthal angle ranging from $[0, 2\pi]$ which rotates around the vertical y-axis.

$$J: (l, \theta, \phi) \mapsto (x, y, z) = 
\begin{bmatrix}
{\frac{\partial x}{\partial l}} & {\frac{\partial x}{\partial \theta}} & {\frac{\partial x}{\partial \phi}} \\
{\frac{\partial y}{\partial l}} & {\frac{\partial y}{\partial \theta}} & {\frac{\partial y}{\partial \phi}} \\
{\frac{\partial z}{\partial l}} & {\frac{\partial z}{\partial \theta}} & {\frac{\partial z}{\partial \phi}}
\end{bmatrix}$$

$$\begin{align*}
&{\frac{\partial x}{\partial l}} = r'(l)sin{\theta}cos{\phi} \quad\quad
&{\frac{\partial x}{\partial \theta}} = r(l)cos{\theta}cos{\phi} \quad\quad
&{\frac{\partial x}{\partial \phi}} = -r(l)sin{\theta}sin{\phi} \\
&{\frac{\partial y}{\partial l}} = r'(l)cos{\theta} \quad\quad
&{\frac{\partial y}{\partial \theta}} = -r(l)sin{\theta} \quad\quad
&{\frac{\partial y}{\partial \phi}} = 0.0 \\
&{\frac{\partial z}{\partial l}} = r'(l)sin{\theta}sin{\phi} \quad\quad
&{\frac{\partial z}{\partial \theta}} = r(l)cos{\theta}sin{\phi} \quad\quad
&{\frac{\partial z}{\partial \phi}} = r(l)sin{\theta}cos{\phi} \\
\end{align*}$$

This kinda blows up the complexity of our metric tensor. Instead of expressing it as a matrix, here is the formula for each component using eistein notation:

$$g_\mu\nu = \frac{\partial\vec{x}}{\partial q^\mu} \frac{\partial\vec{x}}{\partial q^\nu}$$

Before we start, keep in mind that we defined $l$ as **proper radial distance**. It is defined such tha $dl = dr = 1$.

$$
\begin{align*}
g_{ll} &= \frac{\partial\vec{x}}{\partial{l}} \cdot \frac{\partial\vec{x}}{\partial{l}} =
\begin{bmatrix} \sin{\theta}\cos{\phi} \\ \cos{\theta} \\ \sin{\theta}\sin{\phi} \end{bmatrix} \cdot
\begin{bmatrix} \sin{\theta}\cos{\phi} \\ \cos{\theta} \\ \sin{\theta}\sin{\phi} \end{bmatrix} \\
&= \sin^2(\theta)\cos^2(\phi) + \cos^2(\theta) + \sin^2(\theta)\sin^2(\phi) \\
&= \sin^2(\theta)(\cos^2(\phi) + \sin^2(\phi)) + \cos^2(\theta) \\
&= \sin^2(\theta) + \cos^2(\theta) = 1
\end{align*}
$$


$$\begin{align*}
g_{l\theta} &= \frac{\partial\vec{x}}{\partial{l}} \cdot \frac{\partial\vec{x}}{\partial{\theta}} =
\begin{bmatrix} sin{\theta}cos{\phi} \\ cos{\theta} \\ sin{\theta}sin{\phi}\end{bmatrix}
\cdot
\begin{bmatrix} r(l)cos{\theta}cos{\phi} \\ -r(l)sin{\theta} \\ r(l)cos{\theta}sin{\phi}\end{bmatrix} \\
&= r(l)sin{\theta}cos{\theta}cos^2(\phi) - r(l)sin{\theta}cos{\theta} + r(l)sin{\theta}cos{\theta}sin^2(\phi) \\
&= r(l)sin{\theta}cos{\theta}(cos^2(\phi) + sin^2(\phi)) - r(l)sin{\theta}cos{\theta} \\
&= r(l)sin( \theta )cos( \theta ) - r(l)sin{\theta}cos{\theta} \\
&= 0 = g_{\theta l}
\end{align*}$$


$$\begin{align*}
g_{l\phi} &= \frac{\partial\vec{x}}{\partial{l}} \cdot \frac{\partial\vec{x}}{\partial{\phi}} =
\begin{bmatrix} sin{\theta}cos{\phi} \\ cos{\theta} \\ sin{\theta}sin{\phi}\end{bmatrix}
\cdot
\begin{bmatrix} -r(l)sin{\theta}sin{\phi} \\ 0 \\ r(l)sin{\theta}cos{\phi}\end{bmatrix} \\
&=-r(l)sin^2(\theta)cos{\phi}sin{\phi} + r(l)sin^2{\theta}sin{\phi}cos{\phi} \\
&= 0 = g_{\phi l}
\end{align*}$$

$$\begin{align*}
g_{\theta\theta} &= \frac{\partial\vec{x}}{\partial{\theta}} \cdot \frac{\partial\vec{x}}{\partial{\theta}} = 
\begin{bmatrix} r(l)cos{\theta}cos{\phi} \\ -r(l)sin{\theta} \\ r(l)cos{\theta}sin{\phi}\end{bmatrix}
\cdot
\begin{bmatrix} r(l)cos{\theta}cos{\phi} \\ -r(l)sin{\theta} \\ r(l)cos{\theta}sin{\phi}\end{bmatrix} \\
&= r(l)^2cos^2(\theta)cos^2{\phi} + r(l)^2cos^2( \theta )sin^2(\phi) + r(l)^2sin^2(\theta) \\
&= r(l)^2cos^2(\theta)(cos^2(\phi) + sin^2(\phi)) + r(l)^2sin^2(\theta) \\
&= r(l)^2(cos^2(\theta)+sin^2(\theta)) \\
&= r(l)^2
\end{align*}$$


$$\begin{align*}
g_{\theta\phi} &= \frac{\partial\vec{x}}{\partial{\theta}} \cdot \frac{\partial\vec{x}}{\partial{\phi}} = 
\begin{bmatrix} r(l)cos{\theta}cos{\phi} \\ -r(l)sin{\theta} \\ r(l)cos{\theta}sin{\phi}\end{bmatrix}
\cdot
\begin{bmatrix} r(l)sin{\theta}sin{\phi} \\ 0 \\ r(l)sin{\theta}cos{\phi}\end{bmatrix} \\
&= -r(l)^2sin{\theta}cos{\theta}sin{\phi}cos{\phi} + r(l)^2sin{\theta}cos{\theta}sin{\phi}cos{\phi} \\
&= 0 = g_{\phi\theta}
\end{align*}$$

$$\begin{align*}
g_{\phi\phi} &= \frac{\partial\vec{x}}{\partial{\phi}} \cdot \frac{\partial\vec{x}}{\partial{\phi}} = 
\begin{bmatrix} -r(l)sin{\theta}sin{\phi} \\ 0 \\ r(l)sin{\theta}cos{\phi}\end{bmatrix}
\cdot
\begin{bmatrix} -r(l)sin{\theta}sin{\phi} \\ 0 \\ r(l)sin{\theta}cos{\phi}\end{bmatrix} \\
&= r(l)^2sin^2(\theta)sin^2(\phi) + r(l)sin^2(\theta)cos^2(\phi) \\
&= r(l)^2sin^2(\theta)(sin^2(\phi) + cos^2(\phi)) \\
&= r(l)^2sin^2(\theta)
\end{align*}$$

$$\dot{.\hspace{.05in}.} \quad g =\begin{bmatrix}
1 & 0 & 0 \\
0 & r(l)^2 & 0 \\
0 & 0 & r(l)^2sin{\theta}
\end{bmatrix} $$

$$
\mathfrak{g}=g^{-1} =
\begin{bmatrix}
1 & 0 & 0 \\
0 & \frac{1}{r(l)^2} & 0 \\
0 & 0 & \frac{1}{r(l)^2sin{\theta}}
\end{bmatrix}
$$

## Christoffel Symbols
That was brutal. But that is nothing compared to calculating all the non-zero christoffel symbols. In 3D there can be up to 18 of them! Since our metric tensor is a diagonal matrix, we can expect to only get 9 of them.

Recall from earlier:

$$\Gamma^\mu_{\alpha\beta} = \frac{1}{2}\mathfrak{g}^{\mu\lambda}(\frac{\partial g_{\lambda\alpha}}{\partial x^\beta} + \frac{\partial g_{\lambda\beta}}{\partial x^\alpha} + \frac{\partial g_{\alpha\beta}}{\partial x^\lambda})$$

$$
\begin{align*}
&\Gamma^l_{ll} = \frac{1}{2}(1)(0 + 0 - 0) = 0 \\
&\Gamma^l_{l \theta} = \frac{1}{2}(1)(0 + 0 - 0) = 0 = \Gamma^l_{\theta l} \\
&\Gamma^l_{l \phi} = \frac{1}{2}(1)(0 + 0 - 0) = 0 = \Gamma^l_{\phi l} \\
&\Gamma^\theta_{ll} = \frac{1}{2}(\frac{1}{r(l)^2})(0 + 0 - 0) = 0 \\
&\Gamma^\theta_{l\theta} = \frac{1}{2}(\frac{1}{r(l)^2})((\frac{\partial}{\partial l}(r(l)^2) + 0 - 0 ) + 0 + 0) = \frac{2r'(l)r(l)}{2r(l)^2} = \frac{r'(l)}{r(l)} \\
&\Gamma^\theta_{l\phi} = 0 = \Gamma^\theta_{\phi l} \\
&\Gamma^l_{\theta\theta} = \frac{1}{2}(1)((\partial_\theta0 + \partial_\theta0 - \partial_lr(l)^2 )+0+0) = \frac{1}{2}(-2r'(l)r(l)) = -r'(l)r(l) \\
&\Gamma^l_{\phi\phi} = \frac{1}{2}(1)(\partial_lr(l)^2sin^2(\theta)) = \frac{1}{2}(-2r'(l)r(l)sin^2(\theta)) = -r'(l)r(l)sin^2(\theta) \\
&\Gamma^\theta_{\theta\theta} = \frac{1}{2}(\frac{1}{r(l)^2}(0)) = 0 \\
&\Gamma^\theta_{\theta\phi} = \frac{1}{2}(\frac{1}{r(l)^2}(0)) = 0 = \Gamma^\theta_{\phi\theta} \\
&\Gamma^\phi_{l l} = \frac{1}{2}(\frac{1}{r(l)^2sin^2(\theta)}(0)) = 0 \\
&\Gamma^\phi_{l \theta} = \frac{1}{2}(\frac{1}{r(l)^2sin^2(\theta)}(0)) = 0 = \Gamma^\phi_{\theta l} \\
&\Gamma^\phi_{l \phi} = \frac{1}{2}(\frac{1}{r(l)^2sin^2(\theta)})(\partial_l(r(l)^2sin^2(\theta))) = \frac{2r'(l)r(l)sin^2(\theta)}{2r(l)^2sin^2(\theta)} = \frac{r'(l)}{r(l)} = \Gamma^\phi_{\phi l}\\
&\Gamma^\phi_{\theta \theta} = \frac{1}{2}(\frac{1}{r(l)^2sin^2(\theta)}(0)) = 0 \\
&\Gamma^\phi_{\theta \phi} = \frac{1}{2}(\frac{1}{r(l)^2\sin^2(\theta)})(\partial_\theta(r(l)^2\sin^2(\theta))) = \frac{2r(l)^2\sin{\theta}cos{\theta}}{2r(l)^2\sin^2(\theta)} = \cot \theta = \Gamma^\phi_{\phi\theta} \\
&\Gamma^\theta_{\phi\phi} = \frac{1}{2}(\frac{1}{r(l)^2})(-\partial_\theta(r(l)^2sin^2(\theta))) = \frac{-2r(l)^2sin{\theta}cos{\theta}}{2r(l)^2} = -sin{\theta}cos{\theta} \\
\end{align*}
$$

And finally we have our non-zero christoffel symbols:

$$\begin{align*}
&\Gamma^l_{\theta\theta} = -r'(l)r(l) \quad
&\Gamma^l_{\phi\phi} = -r'(l)r(l)sin^2(\theta) \\
&\Gamma^\theta_{\phi\phi} = -sin{\theta}cos{\theta} \quad
&\Gamma^\theta_{l\theta} = \Gamma^\theta_{\theta l} = \frac{r'(l)}{r(l)} \\
&\Gamma^\phi_{l\phi} = \Gamma^\phi_{\phi l} = \frac{r'(l)}{r(l)} \quad
&\Gamma^\phi_{\theta\phi} = \Gamma^\phi_{\phi\theta} = \frac{r'(l)}{r(l)}
\end{align*}$$

## Null Geodesic Equations

Recal the null geodesic formula from earlier:

$$\frac{d^2 x^\alpha}{d\lambda^2 } + \Gamma^\alpha_{\beta\gamma} \frac{dx^\beta}{d\lambda} \frac{dx^\gamma}{d\lambda} = 0$$

$$
\begin{align*}
&\frac{d^2l}{d\lambda^2} + \Gamma^l_{\theta\theta}(\frac{d\theta}{d\lambda}\frac{d\theta}{d\lambda^2}) + \Gamma^l_{\phi\phi}(\frac{d\phi}{d\lambda}\frac{d\phi}{d\lambda}) = 0 \\
\Rightarrow &\frac{d^2l}{d\lambda^2} = r(l)'r(l)\theta'' + r'(l)r(l)sin^2(\theta)(\phi')^2 = \\
&r'(l)r(l)(\theta'' + sin^2(\theta)(\phi')^2) \\\\
&\frac{d^2\theta}{d\lambda^2} + 2\frac{r'(l)}{r(l)}l'\theta' - sin{\theta}cos{\theta}(\phi')^2 = 0 \\
\Rightarrow&\frac{d^2\theta}{d\lambda^2} = -2\frac{r'(l)}{r(l)}l'\theta' + sin{\theta}cos{\theta}(\phi')^2 \\\\
&\frac{d^2\phi}{d\lambda^2} + \frac{r'(l)}{r(l)}l'\phi' - \frac{r'(l)}{r(l)}\phi'l' + cot{\theta}\theta'\phi' + cot{\theta}\theta'\phi' = 0 \\
\Rightarrow&\frac{d^2\phi}{d\lambda^2} = -2\frac{r'(l)}{r(l)}l'\phi' - 2cot{\theta}\theta'\phi'
\end{align*}
$$

## Parallel Transport 

Recall parallel transport formula from earlier:

$$\frac{dv^\mu}{d\lambda} = -\Gamma^\mu_{\nu\rho} v^\nu \frac{dx^\rho}{d\lambda}$$

$$
\begin{align*}
\frac{dv^l}{d\lambda} &= -\Gamma^l_{\theta\theta}v^\theta\theta' - \Gamma^l_{\phi\phi}v^\phi\phi' \\
&= r'(l)r(l)v^\theta\theta' + r'(l)r(l)sin{\theta}v^\phi\phi' \\
&= r'(l)r(l)(v^\theta\theta' + sin^2(\theta)v^\phi\phi')
\\\\
\frac{dv^\theta}{d\lambda} &= -\Gamma^\theta_{l\theta}v^l\theta' - \Gamma^\theta_{\theta l}v^\theta l' - \Gamma^\theta_{\phi\phi}v^\phi\phi' \\
&= -\frac{r'(l)}{r(l)}(v^l\theta' + v^\theta l') + sin{\theta}cos{\phi}v^\phi\phi'
\\\\
\frac{dv^\phi}{d\lambda} &= -\Gamma^\phi_{l\phi}v^l\phi' - \Gamma^\phi_{\phi l}v^\phi l' - \Gamma^\phi_{\theta\phi}v^\theta\phi' - \Gamma^\phi_{\phi\theta}v^\phi\theta' \\
&= -\frac{r'(l)}{r(l)}(v^l\phi' + v^\phi l') + cot{\theta}(v^\theta\phi + v^\phi\theta')
\end{align*}
$$

## Tetrads

Like in the 2D case, we need to build a minkowskian tangent frame for our observers orientation in flat space, and a tetrad to bring the frame into global curved space. This is achieved in very much the same way.

```c++
	void BuildInitialOrientedTetrad()
	{
		//convert cartesian basis to polar basis. Here we use inverse Spherical->Cartesian Jacobian so we can preserve the negative r component if its present.
		Mat3 mInvJacobian = SphericalToCartesianJacobianCurved( m_vSphericalPosition );
		Mat3 mJacobian;
		mInvJacobian.Inverse( &mJacobian );

		//promote jacobian with time component
		const Vec4 e0( 1.0f, 0.0f, 0.0f, 0.0f );
		const Vec4 e1( 0.0f, mJacobian.GetColVec( 0 ) );
		const Vec4 e2( 0.0f, mJacobian.GetColVec( 1 ) );
		const Vec4 e3( 0.0f, mJacobian.GetColVec( 2 ) );

		//orthonormalize with respect to metric
		const Mat4 mMetric = GetMetricTensor( m_vSphericalPosition );
		m_mTetrad = OrthonormalizeWithMetric( e0, e1, e2, e3, mMetric );
	}
```

And just like in the 2D case, we parallel transport the columns of the tetrad as the camera moved through curved space and re-orthonormalize with repsect to the metric using **Reletavistic Gram-Schmidt**.


## Rendering
At scene initialization, a camera is created with a position, minkowskian frame, and default tetrad. Every mouse move creates a euclidean rotation matrix which rotates the minkowski basis. When a move key is triggered a local space vector is passed to the camera's offset function. This vector gets transformed into global wormhole space with the tetrad. We integrate the geodesic equations and nudge the camera to its new position. We then parallel transport the tetrad to the new camera location.

A new frame is dispatched to draw. The new camera position and tetrad is passed to the GPU (strictly speaking, this is not a tetrad, but the tetrad * the orientated minkowski frame ).

The vertex shader generays view space rays by converting TexCoords to normalized device coordinates. It uses the inverse projection matrix to bring these into homogenious coordinates and get view space vectors by performing a perspective divide. In the fragment shader, we multiply by our tetrad to get the ray in wormhole space. We are now ready to raymarch.



## Variable Step Size
After a lot of experimentation I landed on this code for setting the step size during integration in the fragment shader:

```glsl
float ComputeStepSize( State s )
{
	// Max change rate — the fastest-changing velocity coordinate sets the limit
	float maxDeriv = max( abs( s.dl ), max( abs( s.dtheta ), abs( s.dphi ) ) );
	maxDeriv = max( maxDeriv, 1e-6 ); // prevent division by zero

	// Scale step size inversely to that
	float baseStep = 0.1 * fWormholeRadius;
	float step = baseStep / maxDeriv;

	// Add adaptive clamping near poles or throat
	float normDist = abs( s.l ) / fWormholeRadius;
	float t = smoothstep( 0.0, 10.0, normDist );
	float distScale = mix( flmin, 1.0, t );

	float poleDist = min(s.theta, PI - s.theta);  // Distance from either pole
	float poleT = smoothstep(0.01, 0.15, poleDist);  // Tightened range
	float poleScale = mix(fthetamin, 1.0, poleT);

	step *= distScale * poleScale;

	// Final clamp
	float minClamp = 0.001 * fWormholeRadius;
	float maxClamp = 0.5 * fWormholeRadius;
	return clamp( step, minClamp, maxClamp );
}
```

# Lagrangian

At this point we have everything we need to make the wormhole demo work. In fact, its what I ended up using. But I did A LOT of exploration to find different or better ways to do things.

From a physics perspective a **Lagrangian** is a function for how a system evolves such that when you integrate over time (calculate the action) it is stationary. There are many different functions that model any type of system. Our system is a massless particle traveling the shortest path between two points. Minimizing the action of the particle will give us the null geodesic.

To do this we will use the squared distance lagrangian: $g_{\mu\nu} q'^\mu q'^\nu$.

$$
\mathcal{L}(l, \theta, \phi, l', \theta', \phi') = g_{\mu\nu} q'^\mu q'^\nu = l'^2 + r(l)^2\theta'^2 + r(l)^2 sin^2\theta \phi'^2
$$

Now we can apply Euler-Lagrange to get geodesic equations.

$$
\frac{d}{d\lambda}(\frac{\partial\mathcal{L}}{\partial q'^\mu}) - \frac{\partial\mathcal{L}}{\partial q^\mu} = 0
$$

$$
\frac{\partial\mathcal{L}}{\partial r'} = 2r' \quad\quad \frac{\partial\mathcal{L}}{\partial r} = 2r\theta'^2
$$

$$
\dot{.\hspace{.05in}.} \quad 2r'' - 2r\theta' = 0 \Rightarrow r'' = r\theta'^2
$$

$$
\frac{\partial\mathcal{L}}{\partial \theta'} = 2\theta'r^2 \quad\quad \frac{\partial\mathcal{L}}{\partial \theta} = 0
$$

$$
\dot{.\hspace{.05in}.} \quad 4\theta'rr' + 2r^2\theta'' = 0 \Rightarrow \theta'' = -2\frac{r'}{r}\theta'
$$

Way easier than calculating from christoffel symbols explicitly!

# Hamiltonian

This is the point where my understanding starts to falter. The hamiltonain is similar to the lagrangian in that is models the evolution of a system, except the hamiltonian models the sum of potential and kenetic energy where the lagrandian models the difference of potential and kenetic. What this means in practice is beyond my understanding at the moment.

What I do know is that we currently have a system of second order ordinary differential equations to describe the path our massless particle takes along null geodesics. The hamiltonian equations of motion allow us to do the same in a system of first order ODEs. It achieves this by integrating in a *phase space* whatever that means. All I know is that instead of integrating velocity $v^\mu$ we integrate cannonical momenta (momentum) $v_\mu$ which is a covector.

First lets define the cannonical momenta for our Hamiltonian:

$$
\begin{align*}
&p_{l} = \frac{\partial\mathcal{L}}{\partial l'} = 2l' \\
&p_{\theta} = \frac{\partial\mathcal{L}}{\partial \theta'} = 2r(l)^2 \theta' \\
&p_{\phi} = \frac{\partial\mathcal{L}}{\partial \phi'} = 2r(l)^2sin^2\theta\phi' \theta' 
\end{align*}
$$

Now we express velocity in terms of momenta:

$$
\begin{align*}
&l' = \frac{p_l}{2} \\
&\theta' = \frac{p_\theta}{2r(l)^2} \\
&\phi' = \frac{p_\phi}{2r(l)^2sin^2\theta}
\end{align*}
$$

Define the Hamiltonian for the system:

$$
\begin{align*}
\mathcal{H} &= p_q q' - \mathcal{L} \\
\mathcal{H} &= \frac{1}{4}( p_l^2 + \frac{p_\theta^2}{r(l)^2} + \frac{p_\phi^2}{r(l)^2sin^2\theta} )
\end{align*}
$$

Now evolve the system using the Hamiltonian equations of motion:

$$
q'^i = \frac{\partial\mathcal{H}}{\partial p_i} \quad\quad p'^i = -\frac{\partial\mathcal{H}}{\partial q_i}
$$

$$
\begin{align*}
&p'_l = \frac{r'(l)}{2r^(l)^3}(p^2_\theta + \frac{p^2_\phi}{sin^2\theta}) \\
&p_\theta' = \frac{cos\theta}{2r(l)^2sin^3\theta}(p^2_\phi) \\
&p_\phi' = 0
\end{align*}
$$

While these first order derivatives allow for better numerical stability, once big catch is that they are not compatible with variable step sizes which is a pretty big performance hit. Furthurmore, variable step sizes during numerical integration combined with better numerical integrators give you the same stability that these equations offer.

A big downside of the Spherical Coordinate system is the presense of coordinate singularities where $\theta = 0$ and $\pi$.
One option to avoid them is to convert to Cartesian Coordinates. Doing this with the 2nd order geodesic equations is very difficult and yeilds monster expressions, but doing it with the Hamiltonian equations seemed more doable. The derivation below converts a 2D wormhole system from polar to cartesian. Polar coordinates have a coordinate singularity at the origin. Converting to cartesian would theoretically address this issue.

$$
\begin{align*}
\text{Shape Function} \quad\quad r(l) = \sqrt{l^2 + r_0^2} \\
\text{Metric} \quad\quad g_{\mu,\nu}(l, \theta) = \begin{bmatrix}
1 & 0 \\
0 & r(l)^2
\end{bmatrix} \\
\text{Inverse Metric} \quad\quad \mathfrak{g}^{\mu\nu} =
\begin{bmatrix}
1 & 0 \\
0 & \frac{1}{r(l)^2}
\end{bmatrix} \\
\text{Hamiltonian} \quad\quad \mathcal{H} = \frac{1}{2}({p_l^2 + \frac{p_\theta^2}{r(l)^2}}) \\
\text{Hamiltons Equations} \quad \frac{dx^\mu}{d\lambda} = \frac{\partial \mathcal{H}}{\partial p_\mu} \quad,\quad \frac{dp_\mu}{d\lambda} = -\frac{\mathcal{H}}{\partial x^\mu}
\end{align*}
$$

Now we conver $\mathcal{H}$ to Cartesian.

$$
\begin{align*}
&x = r(l)cos\theta \\
&y = r(l)sin\theta
\end{align*}
$$

Let's define

$$
\vec{x} = \begin{bmatrix} x \\ y \end{bmatrix}
$$

and re-express $\mathcal{H}(x, y, p_x, p_y)$.

$$
\hat{n} = \begin{bmatrix} cos{\theta} \\ sin{\theta} \end{bmatrix} \\
J = \begin{bmatrix}
&r'cos\theta & -rsin\theta \\
&r'sin\theta & rcos\theta
\end{bmatrix} \\
\begin{bmatrix} p_l \\ p_\theta\end{bmatrix} = J^{-1} \begin{bmatrix} p_x \\ p_y\end{bmatrix}
$$

which gives:

$$
\begin{align*}
&p_l = r'( p_x cos{\theta} + p_y sin\theta ) = r'(\hat{n} \cdot \vec{p}) \\
&p_\theta = -rsin\theta p_x + rcos\theta p_y = r(\hat{n}_\perp \cdot \vec{p})
\end{align*}
$$

$$
\mathcal{H} = \frac{1}{2}(p_l^2 + \frac{p_\theta^2}{r^2}) =
\frac{1}{2}(r'^2(\vec{p} \cdot \hat{n})^2 + (\hat{n_\perp \cdot \vec{p}})^2) \\
$$

Now evolve the rays wit $\frac{dx}{d\lambda} = \frac{\partial\mathcal{H}}{\partial p_x}$ and $\frac{dp_x}{d\lambda} = -\frac{\partial\mathcal{H}}{\partial x}$.

$$
\begin{align*}
\text{define:} \quad\quad &p_\parallel = \vec{p} \cdot \hat{n} = r'(p_x x + p_y y) \\
&p_\perp = \vec{p} \cdot \hat{n}_\perp = r'( -p_x y + p_y x )
\end{align*}
$$

$$
\dot{.\hspace{.05in}.} \quad \mathcal{H} = \frac{1}{2} ( (\frac{r'}{r}p_\parallel )^2 + p_{\perp}^2 )
$$

$$
\begin{align*}
&\frac{p_{\parallel}}{\partial p_x} = \frac{x}{r} \quad &\frac{\partial p_\parallel}{\partial p_y} = \frac{y}{r} \\
&\frac{p_{\perp}}{\partial p_x} = -\frac{y}{r} \quad &\frac{\partial p_\perp}{\partial p_y} = \frac{x}{r}
\end{align*}
$$

$$
\begin{align*}
&\frac{\partial\mathcal{H}}{\partial p_x} = (\frac{r'}{r})^2p_{\parallel} \cdot \frac{\partial p_\parallel}{\partial p_x} + p_\perp \cdot \frac{\partial p_\perp}{\partial p_x} = (\frac{r'}{r})^2 p_\parallel \frac{x}{r} + p_\perp \frac{-y}{r} = \frac{dx}{d\lambda} \\
&\frac{\partial\mathcal{H}}{\partial p_y} = (\frac{r'}{r})^2 p_\parallel \frac{y}{r} + p_\perp \frac{x}{r} = \frac{dy}{d\lambda}
\end{align*}
$$

$\frac{\partial\mathcal{H}}{\partial x}$ and $\frac{\partial\mathcal{H}}{\partial y}$ are harder because $r$ and $r'$ depend on $x$ and $y$ and so do $p_\parallel$ and $p_\perp$.

$$
r = \sqrt{x^2 + y^2} \Rightarrow \frac{\partial r}{\partial x} = \frac{x}{r} \quad\text{and}\quad \frac{\partial r}{\partial y} = \frac{y}{r} \\
A = (\frac{r'}{r})^2 = 1 - (\frac{r_0}{r})^2
$$

$$
\frac{\partial A}{\partial x} = \frac{dA}{dr} \cdot \frac{\partial r}{\partial x} = \frac{2r_0^2}{r^3} \cdot \frac{x}{r} = \frac{2r_0^2 x}{r^4} \\
\frac{\partial A}{\partial y} = \frac{2r_0^2 y}{r^4}
$$

$$
\begin{align*}
&\frac{\partial P_\parallel}{\partial x} = \frac{1}{r}p_x - \frac{x}{r^3}(p_x x + p_y y) = \frac{1}{r}p_x - \frac{x}{r^2}P_{\parallel} \\
&\frac{\partial P_\perp}{\partial y} = \frac{1}{r}p_y - \frac{x}{r^2}P_{\perp}
\end{align*}
$$

$$
\begin{align*}
\frac{\partial p_x}{\partial \lambda} = -\frac{\partial\mathcal{H}}{\partial x} &= -\frac{1}{2}(\frac{2r_0^2 x}{r^4}P_{\parallel}^2 + 2AP_{\parallel}(\frac{1}{r}p_x - \frac{x}{r^2}P_{\parallel}) + 2P_{\perp}(\frac{1}{r}p_y - \frac{x}{r^2}P_{\perp})) \\
&= -( \frac{r_0^2x}{r^4} P_{\parallel}^2 + AP_{\parallel}(\frac{1}{r}p_x - \frac{x}{r^2}P_{\parallel}) + P_{\perp}(\frac{1}{r}p_y - \frac{y}{r^2}P_{\perp}) )
\end{align*}
$$

$$
\begin{align*}
\frac{\partial p_y}{\partial \lambda} = -\frac{\partial\mathcal{H}}{\partial y}
&= -( \frac{r_0^2y}{r^4} P_{\parallel}^2 + AP_{\parallel}(\frac{1}{r}p_y - \frac{y}{r^2}P_{\parallel}) + P_{\perp}(\frac{1}{r}p_y - \frac{y}{r^2}P_{\perp}) )
\end{align*}
$$

Doing this highlights a major downside of Cartesian Coordinates. It has no way of tracking geodesics through the throat. Polar Coordinates express the higher dimensional nature of a wormhole by allowing the radial coordinate to become negative. There is no mechanism like that in Cartesian. Thus, geodesics degenrate past the throat. Similar stuff would likely happen if I went ahead and converted the 3D Hamiltonian to Cartesian. This combined with the fact that the Hamiltonian equations dont allow for variable step size during integrations makes this idea a non-starter.
