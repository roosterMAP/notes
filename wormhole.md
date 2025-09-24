# Simulating a Wormhole

This page is a compilation of my notes and derivations for rendering, navigating around, and traversing a physically accurate Morrise-Thorn wormhole.

The wormhole implemented here was first described in 1987 in this paper:
[American Journal of Physics Volume 56 issue 5 1988 [doi 10.1119%2F1.15620] Morris, Michael S. -- Wormholes in spacetime and their use for interstellar travel- A tool for teaching general relativity.pdf (PDFy mirror) : Free Download, Borrow, and Streaming : Internet Archive](https://archive.org/details/pdfy--5NLGQyfpB61dmyG/mode/2up)

We will start by a derivation of the null geodesic equation and parallel transport equations in polar coordinates as well as a description of an implementation for a 2D observer. Once that is covered, we will go over a full 3D derivation and implementation.

## Polar Coordinates

In order to leverage the symetry of a wormhole we will be working in **Polar Coordinates**.
The 2D implementation of our wormhole will be using this coordinate system: $(r, \theta )$. The parametric equations are shown below.
$$\begin{align*}
x &= r \cos\theta \\
y &= r \sin\theta
\end{align*}$$

## Wormhole Geometry

The geometry of the wormhole is described using a very simple expression that is only in terms of the radial component. $$r(l) = \sqrt{ l^2 + r_0^2 } \quad \text{where } r_0 \text{ is the radius of the throat}$$

In this setup, the wormhole geometry is described using the coordinates $( l, \theta )$, where $l$ is the proper radial distance. The standard polar coordinates $(r, \theta)$ are modified so that the radial coordinate is now a function of $l$, i.e., $r = r(l)$
This ensures the coordinate transformation depends on the geometry of the wormhole. Unlike standard spherical coordinates where $r$ is independent, here $r$ is a function of the proper radial distance $l$, which encodes the spatial curvature of the wormhole throat and flare-out shape.


## Metric Tensor 

The first thing we must define is our Jacobian. The Jacobian is a Tensor Field that lets us transform vectors at some point P from cartesian $(x,y)$ to our new coordinate system $(l,\theta)$.

We first define our transformation:
$$\begin{align*}
x &= r(l) \cos\theta \\
y &= r(l) \sin\theta
\end{align*}$$

Then we derive the Jacobian. To compute $g$ for our modified Polar coordinates we need the Jacobian that maps to Cartesian coordinates.
$$J: (l, \theta) \mapsto (x, y) = 
\begin{bmatrix}
{\frac{\partial x}{\partial l}} & {\frac{\partial x}{\partial \theta}} \\
{\frac{\partial y}{\partial l}} & {\frac{\partial y}{\partial \theta}} \\
\end{bmatrix}$$

$$\begin{align*}
{\frac{\partial x}{\partial l}} = r'(l)cos(\theta)     &\quad\quad {\frac{\partial x}{\partial \theta}} = -r(l)sin(\theta) \\
{\frac{\partial y}{\partial l}} = r'(l)sin(\theta)     &\quad\quad {\frac{\partial y}{\partial \theta}} = r(l)cos(\theta)
\end{align*}$$

And now we derive the Metric Tensor $g$.
The **metric tensor** is a rank-2 tensor that defines how distances and angles are measured on a curved space or manifold. It generalizes the dot product to arbitrary coordinates and geometries, allowing us to compute quantities like lengths, angles, and volumes.
The metric $ds^2 = g_{\mu\nu}dx^{\mu}dx^{\upsilon}$ takes two vectors and returns a scalar: their inner product. It effectively tells us **how the basis vectors themselves are stretched, skewed, or rotated** at each point in space.

$$g = J^TJ = 
\begin{bmatrix}
{\frac{\partial x}{\partial l}}^2 {\frac{\partial y}{\partial l}}^2 &
{\frac{\partial x}{\partial \theta}} {\frac{\partial x}{\partial l}} +  {\frac{\partial y}{\partial \theta}} {\frac{\partial y}{\partial l}}\\
{\frac{\partial x}{\partial \theta}} {\frac{\partial x}{\partial l}} +  {\frac{\partial y}{\partial \theta}} {\frac{\partial y}{\partial l}} &
{\frac{\partial y}{\partial l}}^2 {\frac{\partial x}{\partial l}}^2
\end{bmatrix} \\= 
\begin{bmatrix}
r'(l)cos^2(\theta) + r'(l)sin^2(\theta) & -r(l)cos(\theta)sin(\theta) + r(l)cos(\theta)sin(\theta)\\
-r(l)cos(\theta)sin(\theta) + r(l)cos(\theta)sin(\theta) & r'(l)r(l)^{2}cos^2(\theta) + r'(l)r(l)^{2}sin^2(\theta)
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

$$g_{\mu\nu} = \begin{bmatrix}
-1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & r(l)^2
\end{bmatrix} $$

The time component will always be the first column and negative. This is mainly to distinguish the temporal component from the spatial ones.

One of the core pillars of General Relativity is the notion that space is locally flat (Minkowskian) and globally curved. 

$$\eta_{\mu\nu} =\begin{bmatrix}
-1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1
\end{bmatrix} $$
The **Minkowskian Metric** describes local flat space. It is what we will use to store the observer's orientation driven by the user (mouse move).

The observer's local tangent frame must be transformed to this globally curved space in order to be offset along the geodesic. This is where **tetrads** come in.
Tetrads will transform vectors from locally flat cartesian space to our globally curved polar space and back.
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
