# Rendering A Black Hole


In the previous project, we rendered a simplified version of a morris-Thorne Wormhole using differential geometry. We started with a metric in modified spherical coordinates and derive Christoffell Sybols, then used those to build the first order system of Parallel Transport equations and a second order system of Null Goedesic equations. We also build the Hamiltonian for the system to get a first order system of null geodesic equations that evolve in a phase space. The observer was tracked along null geodesics and the tetrad frame was build by simply transforming the minkowskian basis of the observer with the jacobian of the modified spherical coordinate system we defined, then orthonormalizing with respect to the metric. While this works from a rendering perspective, many laws of physics were ignored. The observer was assumed to have no mass and no velocity. Instead of evolving the observes position and inertialm frame through the curved space-time of the wormhole, we were simply teleporting the player to different points along a null geodesic.

In this project we will faithfully simulate a spinning black hole and the physical timelike path an observer would take while traversing its gravitational field. We will start by exploring the Schwarzschild metric and the coordinate singularity that appears at the event horizon, motivating the need for horizon-penetrating coordinate systems. Then we will analyize the Kerr-Newmann metric and weigh the difficulties of solving geodesic equations with regular geometric methods. This will take us into a deep dive of the work of Brandon Carter who discovered how to leverage symmetris in the Kerr-Newmann system to elegantly solve both null and timelike geodesics. This seminal paper will be the backbone of our implementation. Once done, we will shift our focus to our innertial observe and do a deep dive of tetrad construction and parallel transport. We will analyze both numerical and analytical methods. Tieing these together will give us a faithful simulation of an observer in a Kerr gravitational field. To complete the picture, we will also briefly explore the physics of accretion flows, relativistic aberration, gravitational redshift, and radiative transfer. Finally, we will discuss analytical solutions of the null geodesic equations, the evaluation of elliptic integrals, and several other techniques used to render black holes in real time.


## What is a Black Hole?

Before jumping into the Schwarzschild metric, lets start with the older Newtonian intuition that led people toward the idea of black holes in the first place. A black hole is not merely an object with a large escape velocity; in general relativity, it is a region of spacetime whose causal structure prevents anything inside a boundary from escaping to infinity. But historically, the first hint of such an object came from a much simpler question: what happens if a star is so massive and compact that even light cannot escape?

In 1687, Isaac Newton published *Mathematical Principles of Natural Philosophy* which laid the groundwork for classical mechanics and calculus.

In it, he proposed a fundamental law of Gravity:

By started with Newtonian total mechanical enery and his law of gravitation:

$$
E = K + U
$$

$$
F = G\frac{Mm}{r^2} \quad\text{where G is a constant}
$$

By integrating the force law wrt $r$ we can get the gravitational potential for a test mass $m$ moving away from a larger mass $M$:
$$
K = \frac{1}{2} m v^2
\\
U = -\frac{GMm}{r}
$$

By setting $E = 0$, we define the minimum condition for $m$ to escape the gravitational pull of $M$.

Solving for $v$ gives us

$$
v = \sqrt{\frac{2GM}{r}}
$$

Setting the escape velocity to the speed of light we get:

$$
c = \sqrt{\frac{2GM}{r}}
$$

Btw, the speed of light was first estimated by Ole Romer in 1676 by measuring the timing of Io's eclipses given Earth's position relative to Jupiter. He measured approximately 226,663 km/s which is about %75 of our modern measurment.

Solving for r gives

$$
r = \frac{2GM}{c^2}
$$

This is exactly what amateur astronomer Reverend John Michell realized in 1783. For a sufficiently massive and compact object, the escape velocity would be greater than the speed of light. If this massive, compact object was a star emitting light, the light would not be able to escape, as it would always be pulled back towards the star. Thus he called in a *Dark Star*.

This is not yet the modern definition of a black hole. Michell's dark star is a Newtonian object: light fails to escape because gravity pulls it back, much like a cannonball launched below escape velocity. A relativistic black hole is more subtle. The event horizon is not a material surface, and light does not fail because it is too slow to climb out in the ordinary Newtonian sense. Instead, spacetime itself is curved so severely that all future-directed paths inside the horizon lead inward.

Still, the Newtonian argument is a useful historical doorway. It gives us the same characteristic radius,

$$
r = \frac{2GM}{c^2},
$$

which will reappear as the Schwarzschild radius. To understand why this radius is not merely an escape-velocity threshold, but a geometric boundary in spacetime, we need the Schwarzschild metric.

## Schwarzschild metric


Before we introduce the metric lets go over what they are and where they come from.

$$
R_{\mu\nu}
-\frac{1}{2}R\,g_{\mu\nu}
+\Lambda g_{\mu\nu}
=
\frac{8\pi G}{c^4}T_{\mu\nu}
$$

Above is the Einstein Field Equation. At a high level, the Einstein Field Equations are an equality between spacetime geometry and matter. The stress-energy tensor $T_{\mu\nu}$ describes how mass and energy are distributed in space, while the metric tensor $g_{\mu\nu}$ describes the resulting curvature of spacetime. It nice to think of matter telling spacetime how to curve, and spacetime telling matter how to move. Once a metric is known, all of the local geometry follows from it. The Christoffel symbols, geodesic equations, parallel transport equations, and curvature tensors can all be derived from $g_{\mu\nu}$. In principle, this is enough to simulate motion through Schwarzschild spacetime directly.


Using geometrized units where \(G = c = 1\), the line element for the Schwarzschild metric is

$$
ds^2 =
-\left(1 - \frac{2m}{r}\right)dt^2
+
\left(1 - \frac{2m}{r}\right)^{-1}dr^2
+
r^2 d\theta^2
+
r^2\sin^2\theta \, d\phi^2
$$

Out coordinate system $x_\mu = (t, r, \theta, \phi)$. This coordinate system represents a stationary observer infinetly far from the black hole. We can show this by taking the limit of the line element as $r \rightarrow \infty$ which gives us the Minkowski line element in spherical coordinates: $ds^2 \rightarrow -dt^2 + dr^2 + r^2d\theta^2 + r^2\sin^2\theta\,d\phi^2$.

The same metric in matrix form is

$$
g_{\mu\nu}
=
\begin{pmatrix}
-\left(1 - \frac{2m}{r}\right) & 0 & 0 & 0 \\
0 & \left(1 - \frac{2m}{r}\right)^{-1} & 0 & 0 \\
0 & 0 & r^2 & 0 \\
0 & 0 & 0 & r^2\sin^2\theta
\end{pmatrix}
$$

### Temporal Components

Lets have a close look at key components of this metric:

$$
g_{tt} = -\left(1 - \frac{2m}{r}\right)
$$

As $r \rightarrow \infty $,  $g_{tt} \rightarrow -1$, which is the minkowski metric.

A stationary observer ( $dr = d\theta = d\phi = 0$ ) experiences $ds^2 = -d\tau^2$. where $\tau$ is proper time.

$$
-d\tau^2 = -\left(1 - \frac{2m}{r}\right)dt^2
$$

$$
d\tau = \sqrt{1 - \frac{2m}{r}}\,dt
$$

As the observer falls deeped into the gravitational field, his clock ticks more slowly. By graphing $d\tau$ we can see that as r decreases from $+\infty$, $d\tau$ decreases and hits 0 at $r=2m$. This is the event horizon.

A clock at $r=4m$ has $d\tau = \sqrt\frac{1}{2}dt$. This clock runs about 70% slower than a clock thats in assymptotally flat space. 

### Radial Components


First, lets graph the experience of an infalling observer by setting, $ds^2 = d\theta = d\phi = 0$. Technically, this is an infalling photon. An infalling observer would have $ds^2<0$.

$$
0 =
-\left(1 - \frac{2m}{r}\right)dt^2
+
\left(1 - \frac{2m}{r}\right)^{-1}dr^2
$$

We want to see how t changes wrt r:

$$
\frac{dt}{dr}
=
\pm
\left(1 - \frac{2m}{r}\right)^{-1}
$$

- \+ branch corresponds to outgoing ray
- \- branch corresponds to ingoing ray
- vertical asymptote at &r=2m& is the event horizon where coordinate time t becomes singular.

Graphing this shows us that as our photon approaches the horizon at 2m, it takes longer and longer to reach it. In fact, it never reaches it as coordinate time t explodes to $\infty$. This is, of course, from the perspective of a distant observer.

By integrating $\frac{dt}{dr}$ we can graph how the t and r coordinate behave as the photon falls in.

$$
-\int
\left(1 - \frac{2m}{r}\right)^{-1} dr
=
-(r + 2m \ln\left|\frac{r}{2m} - 1\right|) + C
$$

Graphing multiple curves each with a different $C$ gives us an array of photon paths falling into the black hole. We can see that as r decreases, t increases, and the curve approaches a slope of nearly -1 far from the black hole. But t diverges as it approaches the horizon at $r=2m$. But on the other side of the horizon it seems to progress as normal.

It almost makes you want to pull t down from infinity to make this curve continuous through the horizon. The good news is we can do just that. It turns out that this singularity at $r=2m$ is actually a coordinate singularity and it can be removed using simple coordinate transformation.

Lets take our radial integral and use it to define a new term $r_\ast$. As we obvserved in our graph earlier, $r_\ast$ contains exactly the logarithmic divergence needed to cancel the divergence in t. This is exactly the regularizing behavior we want. As I stated above in the most handwavy way possible, we can "pull down" t from infinity by subtracting the divergence in $r_\ast$ (for infalling rays). We do the reverse for outgoing rays.

Thus we can define two new temporal coordinates: advanced and retarded null coordinates are:

$$
v = t + r_\ast
$$

$$
u = t - r_\ast
$$

For an infalling particle or ingoing light ray, use advanced *Eddington-Finkelstein* time:
$$
v = t + r_\ast
$$

This coordinate transformation was initially disovered by Arthur Eddington in 1924 back befor the modern interpretation of the horizon as a one-way causal boundary. Finkelstein’s 1958 contribution was the physical reinterpretation that showed the Schwarzchild surface $r=2m$ as a one-way point of no-return.

By solving for t and substituting into the metric we get the Schwarzschild metric in ingoing Eddington-Finkelstein $x^\mu = (v, r, \theta, \phi)$.


$$
ds^2 =
-\left(1 - \frac{2m}{r}\right)dv^2
+
2\,dv\,dr
+
r^2 d\theta^2
+
r^2\sin^2\theta\,d\phi^2
$$

And in matrix form:

$$
g_{\mu\nu}
=
\begin{pmatrix}
-\left(1 - \frac{2m}{r}\right) & 1 & 0 & 0 \\
1 & 0 & 0 & 0 \\
0 & 0 & r^2 & 0 \\
0 & 0 & 0 & r^2\sin^2\theta
\end{pmatrix}
$$

From here you can see that there is nononger divergent behavior at the horizon. We still have divergence then $r=0$, but we are not going into proving that this is a real curvature singularity and not a coordinate singularity.

In theory, we could calculate christoffel symbols and solve parallel transport and geodesic equations.

One last thing, you may have noticed that we have non-diagonal terms at $g_{vr}$ and $g_{rv}$. This is perfectly fine, it only means that our basis vector are non-orthogonal. Its perfectly fine to work in a coordinate system like this.

## Kerr-Newmann metric

The Kerr metric was discovered by Roy Kerr in 1963 as the exact solution describing a rotating black hole, and Ezra Newman and collaborators generalized it in 1965 to include electric charge, producing the Kerr-Newman metric.

$$
x^\mu = ( t, r, \theta, \phi )
$$

$$
\begin{aligned}
ds^2 ={}&
-\left(1 - \frac{2mr - e^2}{\rho^2}\right)dt^2
-\frac{2a(2mr - e^2)\sin^2\theta}{\rho^2}\,dt\,d\phi
+\frac{\rho^2}{\Delta}\,dr^2
+\rho^2\,d\theta^2 \\
&+
\left(
r^2 + a^2
+
\frac{a^2(2mr - e^2)\sin^2\theta}{\rho^2}
\right)
\sin^2\theta\,d\phi^2 .
\end{aligned}
$$

$$
\rho^2 = r^2 + a^2\cos^2\theta
\\
\Delta = r^2 - 2mr + a^2 + e^2
$$

$$
g_{\mu\nu}
=
\begin{pmatrix}
-\left(1 - \frac{2mr - e^2}{\rho^2}\right)
&
0
&
0
&
-\frac{a(2mr - e^2)\sin^2\theta}{\rho^2}
\\

0
&
\frac{\rho^2}{\Delta}
&
0
&
0
\\

0
&
0
&
\rho^2
&
0
\\

-\frac{a(2mr - e^2)\sin^2\theta}{\rho^2}
&
0
&
0
&
\left(
r^2 + a^2
+
\frac{a^2(2mr - e^2)\sin^2\theta}{\rho^2}
\right)\sin^2\theta
\end{pmatrix}
$$