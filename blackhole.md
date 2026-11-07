# Rendering A Black Hole

In the previous project, we rendered a simplified version of a Morris-Thorne Wormhole using differential geometry. We started with a metric in modified spherical coordinates, derived Christoffel Symbols, then used those to build a first order system of parallel transport equations and a second order system of geodesic equations. We also constructed the Hamiltonian for the system to get a first order set of geodesic equations that evolve in a phase space. The observer was tracked along geodesics and his tetrad frame was built by transforming local basis of the observer with the jacobian of the modified spherical coordinate system we defined, then orthonormalizing with respect to the metric. While this works from a rendering perspective, many laws of physics were ignored. Observer mass and velocity was never meaningfully considered. Instead of evolving the observer's position and inertial frame through the curved space-time of the wormhole, we were simply teleporting the player to different points along a geodesic.

In this project we will faithfully simulate a spinning black hole and the physical time-like path an observer would take while traversing its gravitational field. We will start by exploring the Schwarzschild metric and the coordinate singularity at the event horizon, motivating the need for a horizon-penetrating coordinate system. Then we will analyze the Kerr-Newman metric and weigh the difficulties of solving geodesic equations with regular geometric methods. This will take us into a deep dive of the work of Brandon Carter who discovered how to leverage symmetries in the Kerr-Newman system to elegantly solve both null and time-like geodesics. This seminal paper will be the backbone of our implementation. Once done, we will shift our focus to the inertial observer and do a deep dive of tetrad construction and parallel transport. We will analyze both numerical and analytical methods. Tying these together will give us a faithful simulation of an observer in a Kerr-Newman gravitational field. To complete the picture, we will also briefly explore the physics of accretion flows, relativistic aberration, gravitational redshift, and radiative transfer. Finally, we will discuss analytical solutions of the null geodesic equations, the evaluation of elliptic integrals, and several other techniques used to render black holes in real time.

# Geometry

Before jumping into the Schwarzschild metric, lets start with the older Newtonian intuition that led people toward the idea of black holes in the first place. A black hole is not merely an object with a large escape velocity; in general relativity, it is a region of spacetime whose causal structure prevents anything inside a boundary from escaping to infinity. But historically, the first hint of such an object came from a much simpler question: what happens if a star is so massive and compact that even light cannot escape?

In 1687, Isaac Newton published *Mathematical Principles of Natural Philosophy* which laid the groundwork for classical mechanics and calculus.

In it, he proposed a fundamental law of gravity. Let's start with the Newtonian total mechanical energy $E$ and his law of gravitation $F$:

$$
E = K + U \quad\text{where}\quad K = \frac{1}{2} m v^2
$$

$$
F = G\frac{Mm}{r^2} \quad\text{where G is a constant}
$$

By integrating the force law $F$ wrt $r$ (from $r$ out to $\infty$, taking $U(\infty)=0$) we can get the gravitational potential $U$ for a test mass $m$ moving away from a larger mass $M$:

$$
U = -\frac{GMm}{r}
$$

By setting $E = 0$, we define the minimum condition for $m$ to escape the gravitational pull of $M$:

$$
\frac{1}{2}mv^2 - \frac{GMm}{r} = 0
$$

Solving for $v$ gives us (note that $m$ cancels, so the escape velocity does not depend on the test mass)

$$
v = \sqrt{\frac{2GM}{r}}
$$

Setting the escape velocity to the speed of light $c$ we get:

$$
c = \sqrt{\frac{2GM}{r}}
$$

Solving for $r$ gives

$$
r = \frac{2GM}{c^2}
$$

This is exactly what amateur astronomer Reverend John Michell realized in 1783. For a sufficiently massive and compact object, the escape velocity would be greater than the speed of light. If this massive, compact object was a star emitting light, the light would not be able to escape, as it would always be pulled back towards the star. Thus he called it a *Dark Star*.

This is not yet the modern definition of a black hole. Michell's dark star is a Newtonian object: light fails to escape because gravity pulls it back, much like a cannonball launched below escape velocity. A relativistic black hole is more subtle. The event horizon is not a material surface, and light does not fail because it is too slow to climb out in the ordinary Newtonian sense. Instead, spacetime itself is curved so severely that all future-directed paths inside the horizon lead inward.

Still, the Newtonian argument is a useful historical doorway. It gives us the same characteristic radius,

$$
r = \frac{2GM}{c^2}
$$

which will reappear as the Schwarzschild radius. To understand why this radius is not merely an escape-velocity threshold, but a geometric boundary in spacetime, we need the Schwarzschild metric.

## Schwarzschild metric

Before we introduce the metric lets go over what they are and where they come from.

$$
R_{\mu\nu} -
\frac{1}{2}Rg_{\mu\nu} +
\Lambda g_{\mu\nu} =
\frac{8\pi G}{c^4}T_{\mu\nu}
$$

Above is the Einstein Field Equation. At a high level, the Einstein Field Equations are an equality between spacetime geometry and matter. The stress-energy tensor $T_{\mu\nu}$ describes how mass and energy are distributed in space, while the metric tensor $g_{\mu\nu}$ describes the resulting curvature of spacetime. It's nice to think of matter telling spacetime how to curve, and spacetime telling matter how to move. Once a metric is known, all of the local geometry follows from it. The Christoffel symbols, geodesic equations, parallel transport equations, and curvature tensors can all be derived from $g_{\mu\nu}$. In principle, this is enough to simulate motion through Schwarzschild spacetime directly.

Using geometrized units where $G = c = 1$ (and writing the black hole mass as $m$ from here on, so $2GM/c^2 \rightarrow 2m$), the line element for the Schwarzschild metric is

$$
ds^2 =
-\left(1 - \frac{2m}{r}\right)dt^2 +
\left(1 - \frac{2m}{r}\right)^{-1}dr^2 +
r^2 d\theta^2 +
r^2\sin^2\theta  d\phi^2
$$

This coordinate system $x^\mu = (t, r, \theta, \phi)$ is adapted to a stationary observer infinitely far from the black hole: $t$ is the time measured by that observer's clock. We can show this by taking the limit of the line element as $r \rightarrow \infty$. Substituting gives us the Minkowski line element in spherical coordinates: $ds^2 \rightarrow -dt^2 + dr^2 + r^2d\theta^2 + r^2\sin^2\theta d\phi^2$.

The same metric in matrix form is

$$
g_{\mu\nu} =
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
d\tau = \sqrt{1 - \frac{2m}{r}}dt
$$

As the observer falls deeper into the gravitational field, his clock ticks more slowly (as measured by external observers at $r=\infty$). By [graphing](https://www.desmos.com/calculator/uiwjothwn3) $d\tau$ we can see that as r decreases from $+\infty$, $d\tau$ decreases and hits 0 at $r=2m$. This is the event horizon.

A clock at $r=4m$ has $d\tau = \sqrt{\frac{1}{2}}dt \approx 0.707\,dt$. This clock runs at about 71% of the rate (about 29% slower than) a clock that's in asymptotically flat space. 

### Radial Components

First, lets consider radial motion by setting $ds^2 = d\theta = d\phi = 0$. Technically, this is a radially moving photon (ingoing or outgoing). A massive infalling observer would have $ds^2<0$.

$$
0 =
-\left(1 - \frac{2m}{r}\right)dt^2 +
\left(1 - \frac{2m}{r}\right)^{-1}dr^2
$$

We want to see how t changes wrt r:

$$
\frac{dt}{dr} =
\pm
\left(1 - \frac{2m}{r}\right)^{-1}
$$

- \+ branch corresponds to outgoing ray
- \- branch corresponds to ingoing ray
- vertical asymptote at $r=2m$ is the event horizon where coordinate time t becomes singular.

[Graphing](https://www.desmos.com/calculator/lnzoz9ivdm) this shows us that as our photon approaches the horizon at 2m, it takes longer and longer to reach it. In fact, it never reaches it as coordinate time t explodes to $\infty$. This is, of course, from the perspective of a distant observer.

By integrating $\frac{dt}{dr}$ we can graph how the t and r coordinate behave as the photon falls in.

$$
-\int \left(1 - \frac{2m}{r}\right)^{-1} dr =
-(r + 2m \ln\left|\frac{r}{2m} - 1\right|) + C
$$

[Graphing](https://www.desmos.com/calculator/0vj8rdktb9) multiple curves each with a different $C$ gives us an array of photon paths falling into the black hole. We can see that as r decreases, t increases, and the curve approaches a slope of nearly $-1$ far from the black hole. It diverges as it approaches the horizon at $r=2m$. On the other side of the horizon the curve comes back down from infinity and is perfectly well behaved again, ending at a finite $t$ at $r=0$.

It almost makes you want to pull the divergent region down from infinity to make this curve continuous through the horizon. The good news is we can do just that. It turns out that this singularity at $r=2m$ is actually a coordinate singularity and it can be removed using simple coordinate transformation.

Lets take our radial integral and use it to define a new term, the *tortoise coordinate* $r_\ast$:

$$
r_\ast = \int \left(1 - \frac{2m}{r}\right)^{-1} dr = r + 2m \ln\left|\frac{r}{2m} - 1\right|
$$

As we observed in our graph earlier, $r_\ast$ contains exactly the logarithmic divergence needed to cancel the divergence in t. As I stated above in the most hand-wavy way possible, we can "pull down" $t$ from infinity by subtracting the divergence in $r_\ast$ (for infalling rays). We do the reverse for outgoing rays.

Thus we can define two new temporal coordinates: *advanced* and *retarded* null coordinates:

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

This coordinate transformation was initially discovered by Arthur Eddington in 1924 back before the modern interpretation of the horizon as a one-way causal boundary. Finkelstein’s 1958 contribution was the physical reinterpretation that showed the Schwarzschild surface $r=2m$ as a one-way point of no-return.

Along an ingoing ray $t = -r_\ast + C$, so $v = C$ is constant: each ingoing light ray is simply a line of constant $v$.

By solving for $t$ and substituting into the metric we get the Schwarzschild metric in ingoing Eddington-Finkelstein $x^\mu = (v, r, \theta, \phi)$. Explicitly, $dt = dv - \left(1 - \frac{2m}{r}\right)^{-1}dr$, so

$$
-\left(1 - \frac{2m}{r}\right)dt^2 =
-\left(1 - \frac{2m}{r}\right)dv^2 + 2dvdr - \left(1 - \frac{2m}{r}\right)^{-1}dr^2
$$

and the last term exactly cancels the $g_{rr}$ term of the original metric.

$$
ds^2 =
-\left(1 - \frac{2m}{r}\right)dv^2 +
2dvdr +
r^2 d\theta^2 +
r^2\sin^2\theta d\phi^2
$$

And in matrix form:

$$
g_{\mu\nu} =
\begin{pmatrix}
-\left(1 - \frac{2m}{r}\right) & 1 & 0 & 0 \\
1 & 0 & 0 & 0 \\
0 & 0 & r^2 & 0 \\
0 & 0 & 0 & r^2\sin^2\theta
\end{pmatrix}
$$

From here you can see that there is no longer divergent behavior at the horizon. We still have divergence at $r=0$, but we are not going into proving that this is a real curvature singularity and not a coordinate singularity.

From here, we can easily calculate Christoffel symbols and solve parallel transport and geodesic equations like we did with the Wormhole project.

One last thing, you may have noticed that we have non-diagonal terms at $g_{vr}$ and $g_{rv}$. This is perfectly fine, it only means that our basis vectors are non-orthogonal. Also note $g_{rr}=0$: the coordinate basis vector $\partial_r$ is itself null (it points along the ingoing light rays).

## Kerr-Newman metric

The Kerr metric was discovered by Roy Kerr in 1963 as the exact solution describing a rotating black hole, and Ezra Newman and collaborators generalized it in 1965 to include electric charge, producing the Kerr-Newman metric.

This new metric takes into consideration two additional properties of the black hole:
- spin $a$
- charge $e$

We are including $e$ for completeness, but in nature we generally assume that no black holes have charge simply because if they did they would quickly attract matter of the opposite charge and return to neutral immediately. We also don't consider the other Standard Model gauge charges such as color charge or weak isospin even though there are fun theories about them in some Big Bang models.

$$
x^\mu = ( t, r, \theta, \phi )
$$

$$
\begin{aligned}
ds^2 ={}&
-\left(1 - \frac{2mr - e^2}{\rho^2}\right)dt^2
-\frac{2a(2mr - e^2)\sin^2\theta}{\rho^2}dtd\phi
+\frac{\rho^2}{\Delta}dr^2
+\rho^2d\theta^2 \\
&+
\left(
r^2 + a^2 +
\frac{a^2(2mr - e^2)\sin^2\theta}{\rho^2}
\right)
\sin^2\theta d\phi^2 .
\end{aligned}
$$

The metric uses two new terms: $\rho^2 = r^2 + a^2\cos^2\theta$ and $\Delta = r^2 - 2mr + a^2 + e^2$.
- $\rho^2$ is the equivalent of $r^2$ in Schwarzschild. Because Kerr rotates about its spin axis, its geometry is flattened/oblate. $\rho^2$ appropriately scale based off the spin $a$ and latitudinal angle of approach $\theta$.
- Horizon(s) occur when $\Delta = 0$. You can solve this using the quadratic equation: $r_\pm = m \pm \sqrt{m^2 - a^2 - e^2}$. [Graphing](https://www.desmos.com/calculator/gurvrqqy9s) this will give you two circles. When $a=0$ the inner horizon collapses into a point and the system reduces down to the Schwarzschild case. When $a=m$ we get an extremal black hole where the inner and outer horizon perfectly overlap (assuming $e=0$; in general the extremal condition is $a^2+e^2=m^2$). If $a^2+e^2>m^2$ there is no horizon at all (a naked singularity).

Before we dive in, it is useful to check how Kerr-Newman reduces to simpler metrics.
- If $e=0$, the metric becomes the Kerr metric: a rotating, uncharged black hole.
- If $a=0$, the metric becomes the Reissner-Nordström metric: a charged, non-rotating black hole.
- If both $a=0$ and $e=0$, then $\rho^2 = r^2$ and $\Delta = r^2 - 2mr = r(r-2m)$, so the metric reduces to Schwarzschild.

Here is the same metric in matrix form:

$$
\rho^2 = r^2 + a^2\cos^2\theta
\\
\Delta = r^2 - 2mr + a^2 + e^2
$$

$$
g_{\mu\nu} =
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
\left( r^2 + a^2 + \frac{a^2(2mr - e^2)\sin^2\theta}{\rho^2} \right)\sin^2\theta
\end{pmatrix}
$$

### The Ergosphere

Lets analyze components of the metric like we did for Schwarzschild.

Remember $ds^2 = -d\tau^2$:

Lets consider a static observer where $dr=d\theta=d\phi=0$.

$$
d\tau = \sqrt{1-\frac{2mr-e^2}{\rho^2}}dt
$$

[Graphing](https://www.desmos.com/calculator/0y724ov3c9) looks very similar to Schwarzschild except $g_{tt}$ depends both on $r$ and $\theta$. So, gravitational time dilation is no longer purely radial. It depends on the latitude relative to the black hole's spin axis. Another thing that's strange is that when $a\ne0$, $d\tau$ hits $0$ outside the outer horizon (at $r=2m$ on the equator when $e=0$) rather than at the outer horizon.
Setting $g_{tt}=0$ and solving for $r$ gives us:

$$
r_{\text{ergo}} = m + \sqrt{m^2 - a^2\cos^2{\theta}-e^2}
$$

We can add this to our [graph](https://www.desmos.com/calculator/ofqdzk3z9a) of the inner and outer horizon. We get a third, larger, oblate surface. It is not a circle: it touches the outer horizon at the poles ($\cos^2\theta=1$) and bulges out to $r=2m$ at the equator (for $e=0$).
Lets dig deeper into the nature of this new *ergosphere*.

Consider an observer trying to stay fixed at $r$ and $\theta$: $dr = d\theta = 0$.

$$
\Omega = \text{angular velocity} = \frac{d\phi}{dt}
\\
d\phi = \Omega dt
$$

Substituting into the metric gives us:

$$
ds^2 = (g_{tt} + 2g_{t\phi}\Omega + g_{\phi\phi}\Omega^2)dt^2
$$

For a massive observer, the path must be timelike:

$$
g_{tt} + 2g_{t\phi}\Omega + g_{\phi\phi}\Omega^2 < 0
$$

Once again we can solve for $\Omega$ using the quadratic equation, giving us two limiting angular velocities.

$$
\Omega_\pm = \frac{-g_{t\phi} \pm \sqrt{g_{t\phi}^2 - g_{tt}g_{\phi\phi}}}{g_{\phi\phi}}
$$

Earlier, we [graphed](https://www.desmos.com/calculator/0y724ov3c9) $d\tau$ where $dr=d\theta=d\phi=0$. We could see the curve was discontinuous at $r=2m$. This is because our solution became complex for $0 < r < 2m$ (on the equator, with $e=0$). Now we know why. Inside the ergosphere $g_{tt} > 0$, meaning the vector $\partial_t$ (a static observer's worldline) is spacelike. Setting $\Omega=0$ leaves $ds^2 = g_{tt}dt^2 > 0$, so zero angular velocity is not allowed for a physical observer. Since $g_{tt}>0$ there, both roots $\Omega_\pm$ are positive (same sign as the spin). Thus we define the ergosphere as a region where all observers must corotate with the black hole. The effect of a spinning black hole giving static observers angular velocity via its rotating spacetime is called *frame dragging* because the inertial frame of the observer is being dragged along the $\phi$ direction.

### Principal Null Directions

In the Schwarzschild case we were able to set $ds^2 = d\theta = d\phi = 0$ to track a radially infalling photon. This gave us a pretty simple expression we could integrate. The direction of these radial light rays is called a *Principal Null Direction* (PND). In Schwarzschild the PND photons are simply the photons with zero angular momentum.

In Kerr, frame dragging makes this subtler. Setting $d\phi=0$ no longer gives a geodesic, and the PND photons are *not* zero angular momentum photons. Instead they are the photons that stay at constant $\theta$ while being dragged in $\phi$; they carry angular momentum $\Phi = aE\sin^2\theta$ and Carter constant $K=0$ (both defined later). We will verify this once we have Carter's equations (see the note at the end of the Jacobi Separation section). For now, we will just take it for granted below (+ outgoing, − ingoing):

$$
\frac{dt}{dr} = \pm\frac{r^2 + a^2}{\Delta}
$$

$$
\frac{d\phi}{dr} = \pm\frac{a}{\Delta}
$$

Integrating $\frac{dt}{dr}$ gives us (using $\Delta = (r-r_+)(r-r_-)$, writing $r^2+a^2 = \Delta + 2mr - e^2$, and splitting the remainder into partial fractions with $r_\pm^2 + a^2 = 2mr_\pm - e^2$):

$$
\int\frac{r^2 + a^2}{\Delta}dr = r +
\frac{r_+^2+a^2}{r_+-r_-}\ln|r-r_+| -
\frac{r_-^2+a^2}{r_+-r_-}\ln|r-r_-| + C
$$

When you [graph](https://www.desmos.com/calculator/xnqtojq2im) this you can clearly see the divergence of geodesics in $t$ as $r$ approaches the inner and outer horizons.

Like with Schwarzschild, this will allow us to define a coordinate transformation that will remove the divergence at the horizons. This integral plays the role of the Kerr-Newman tortoise coordinate $r_\ast$, and similarly we define $r_\sharp = \int \frac{a}{\Delta}dr$ for the angle.

### Regularizing Coordinate Transformation

We start by defining our advanced/retarded coordinates. Advanced coordinates regularize the horizon from ingoing geodesics, and retarded coordinates regularize for outgoing ones. 

But why would we care about outgoing geodesics if nothing can escape a black hole? The answer is that the structure of a rotating black hole is far more nuanced and complex than one might think and there is actually a lot of utility for such a coordinate patch.

$$
\text{retarded (outgoing)}\quad\quad
\begin{aligned}
du &= dt - (r^2 + a^2)\Delta^{-1}dr
\\
d\phi_u &= d\phi - a\Delta^{-1}dr
\end{aligned}
$$

$$
\text{advanced (ingoing)}\quad\quad
\begin{aligned}
dv &= dt + (r^2 + a^2)\Delta^{-1}dr
\\
d\phi_v &= d\phi + a\Delta^{-1}dr
\end{aligned}
$$

Integrated, these are $v = t + r_\ast$ and $\phi_v = \phi + r_\sharp$, with $r_\ast$ and $r_\sharp$ the integrals from the previous section. Along an ingoing principal null ray both $v$ and $\phi_v$ are constant. Note that the angle has to be shifted too, not just time; otherwise the dragging term $d\phi/dr = a/\Delta$ would still diverge at the horizon. From here on, whenever we work in EF coordinates we drop the subscript and simply write $\phi$ for $\phi_v$.

Substituting advanced time into the metric gives us:

$$
\begin{aligned}
ds^2 ={}&
-\left(1-\frac{2mr-e^2}{\rho^2}\right)dv^2
+2dvdr
+\rho^2d\theta^2
-\frac{2a(2mr-e^2)\sin^2\theta}{\rho^2}dvd\phi
\\
&-2a\sin^2\theta drd\phi +
\frac{
\left(r^2+a^2\right)^2 -
a^2\Delta\sin^2\theta
}{
\rho^2
}
\sin^2\theta d\phi^2
\end{aligned}
$$

In matrix form, the coordinate transformation to advance time is:

$$
J^\mu_\nu =
\frac{\partial x^{\mu}_{\mathrm{EF}}}{\partial x^{\nu}_{\mathrm{BL}}} =
\begin{pmatrix}
1 & \dfrac{r^2+a^2}{\Delta} & 0 & 0
\\
0 & 1 & 0 & 0
\\
0 & 0 & 1 & 0
\\
0 & \dfrac{a}{\Delta} & 0 & 1
\end{pmatrix}
$$

Vectors transform with $J$ ($V_{\mathrm{EF}} = J V_{\mathrm{BL}}$). The metric is a covariant tensor, so it transforms with the inverse Jacobian: $g_{\mathrm{EF}} = (J^{-1})^T g_{\mathrm{BL}} J^{-1}$. Conveniently, $J^{-1}$ is just $J$ with the two off-diagonal terms negated. Multiplying against the Boyer-Lindquist Metric gives us:

$$
g_{\mu\nu} =
\begin{pmatrix}
-\left(1-\frac{2mr-e^2}{\rho^2}\right)
&
1
&
0
&
-\frac{a(2mr-e^2)\sin^2\theta}{\rho^2}
\\
1
&
0
&
0
&
-a\sin^2\theta
\\
0
&
0
&
\rho^2
&
0
\\
-\frac{a(2mr-e^2)\sin^2\theta}{\rho^2}
&
-a\sin^2\theta
&
0
&
\frac{ \left(r^2+a^2\right)^2 - a^2\Delta\sin^2\theta}{\rho^2}\sin^2\theta
\end{pmatrix}
$$

Following Brandon Carter's paper can be a bit confusing because older Kerr literature does not always match the modern $u$/$v$ naming convention. In modern notation, advanced or ingoing time is usually written as

$$
v = t + r_\ast
$$

while retarded or outgoing time is written as

$$
u = t - r_\ast
$$

Carter instead uses the symbol $u$ and refers to it as retarded time, but the sign of his coordinate transformation matches what we would now call the ingoing or advanced coordinate. To avoid confusion, these notes will also use the word *ingoing* and the variable $v$ throughout.

Also, Carter does not refer to this as an Eddington-Finkelstein coordinate transformation. The modern terminology was not yet standard in the Kerr literature, so the transformation appears as part of the coordinate construction rather than under the name *Eddington-Finkelstein*.

### Evaluating the ingoing Eddington-Finkelstein Metric

A quick look at our new metric in matrix form makes it clear that it is regular at the horizon. Unlike the Boyer-Lindquist metric, $\Delta$ never appears in the denominator of any term and $\rho^2$ is well behaved.

Lets evaluate what it would take to naively derive our second order geodesic equations from the Christoffel symbols like we did for the Morris-Thorne wormhole:

In four dimensions, our Christoffel symbol $\Gamma_{\alpha\beta}^\mu$ has 4x4x4=64 raw components. Since the lower indices are symmetric, only 40 are independent. For this metric, many are nonzero, and expanding the lower-index symmetry gives dozens of terms in the geodesic equations.

$$
\frac{d^2 x^\alpha}{d\lambda^2 } + \Gamma^\alpha_{\beta\mu} \frac{dx^\beta}{d\lambda} \frac{dx^\mu}{d\lambda} = 0
$$

This becomes four coupled second-order equations, each containing a pile of terms with up to 10 velocity-pair terms per equation.

The good news is that the EF Christoffel symbols should also be regular at the horizon, because the EF metric and inverse metric are regular there. The bad news is that they are a steaming pile of algebraic shit.

We need a better way. By leveraging the symmetries in Kerr-Newman, we can dramatically simplify the algebraic complexity of our equations.

# Jacobi Separation

### Kerr Symmetries

The equations of motion of a test particle of mass $\mu$ and charge $\epsilon$ may be derived from the Lagrangian:

$$
L = \frac{1}{2} g_{ij}\dot{x}^i\dot{x}^j + \epsilon A_i\dot{x}^i
$$

where the dot over a symbol denotes ordinary differentiation wrt an affine parameter $\lambda$. We are free to choose the scale of $\lambda$, so we choose it such that it is related to proper time by $d\tau = \mu\, d\lambda$.

Remember our line element for a time-like observer:

$$
ds^2 = g_{ij} dx^i dx^j \quad\text{and}\quad ds^2 = -d\tau^2
$$

Dividing both by $d\lambda^2$ and using $d\tau/d\lambda = \mu$ we obtain:

$$
g_{ij}\dot{x}^i\dot{x}^j = \frac{ds^2}{d\lambda^2} = -\left(\frac{d\tau}{d\lambda}\right)^2 = -\mu^2
$$

The nice thing about this normalization is that photons are included for free: they simply have $\mu = 0$.

In order to support the electromagnetic field of the black hole with charge $e$, we define $A$ as the covariant vector potential:

$$
A = e \rho^{-2}r(dv-a \sin^2\theta\, d\phi)
$$

or in components: $A_v = \dfrac{er}{\rho^2}$, $A_\phi = -\dfrac{aer\sin^2\theta}{\rho^2}$, $A_r = A_\theta = 0$. (If you transform the Boyer-Lindquist potential into EF coordinates you also pick up an $A_r = -er/\Delta$ term. Since it only depends on $r$ it is a pure gradient, which is a gauge transformation, so it can be dropped.)

Honestly, I know little about electromagnetism and have only included charge for completeness. It is enough to know that the vector potential $A$ only shifts the canonical momentum for test particles with charge. Setting $\epsilon = 0$ (or $e = 0$) removes it everywhere.

In order to transform to a Hamiltonian formulation we introduce the canonical momenta:

$$
p_i = \frac{\partial L}{\partial \dot{x}^i} = g_{ij}\dot{x}^j + \epsilon A_i
$$

Inverting this gives the velocity in terms of the momenta, which we will use a lot:

$$
p_i - \epsilon A_i = g_{ij}\dot{x}^j
\quad\Longrightarrow\quad
\dot{x}^i = g^{ij}(p_j - \epsilon A_j)
$$

The Hamiltonian is the Legendre transform $H = p_i\dot{x}^i - L$:

$$
\begin{aligned}
H &= (g_{ij}\dot{x}^j + \epsilon A_i)\dot{x}^i - \frac{1}{2} g_{ij}\dot{x}^i\dot{x}^j - \epsilon A_i\dot{x}^i
= \frac{1}{2} g_{ij}\dot{x}^i\dot{x}^j
\end{aligned}
$$

Written in terms of the momenta (which is what a Hamiltonian must be a function of), this is:

$$
H = \frac{1}{2}g^{ij}(p_i - \epsilon A_i)(p_j - \epsilon A_j)
$$

Substituting $p_i - \epsilon A_i = g_{ik}\dot{x}^k$ back in, and using $g^{ij}g_{ik} = \delta^j_k$, we can check its value along a trajectory:

$$
H = \frac{1}{2}g^{ij} g_{ik}\dot{x}^k g_{jl}\dot{x}^l = \frac{1}{2}\delta^j_k\dot{x}^k g_{jl}\dot{x}^l = \frac{1}{2}g_{kl} \dot{x}^k \dot{x}^l = -\frac{1}{2} \mu^2
$$

The motion is generated by Hamilton's equations:

$$
\dot{x}^i = \frac{\partial H}{\partial p_i} = g^{ij}(p_j - \epsilon A_j)
\qquad
\dot{p}_i = -\frac{\partial H}{\partial x^i}
$$

Since $H$ has no explicit dependence on $\lambda$, $H$ itself is conserved. That is our first constant of motion: the mass $\mu$.

Next, if we assume that the black hole is eternal, does not change size, and spins at a constant rate, then the metric and the potential are *stationary* (independent of $v$) and *axisymmetric* (independent of $\phi$). Coordinates that do not appear in $H$ are called *cyclic*, and Hamilton's equations immediately give us

$$
\dot{p}_v = -\frac{\partial H}{\partial v} = 0
\qquad
\dot{p}_\phi = -\frac{\partial H}{\partial \phi} = 0
$$

Geometrically these correspond to the Killing vectors $\partial_v$ (time translation, giving conservation of energy) and $\partial_\phi$ (rotation about the spin axis, giving conservation of angular momentum about that axis).

Thus we define our next two constants of motion:

$$
E = -p_v
\\
\Phi = p_\phi
$$

The minus sign makes $E$ positive for future-directed particles (in flat space $p_t = -E$). Also, because $v$ and $t$ (and the two $\phi$'s) only differ by functions of $r$, the vector $\partial_v$ at fixed $(r,\theta,\phi)$ is the same vector as $\partial_t$ in Boyer-Lindquist. So $E$ and $\Phi$ take the same values in both coordinate systems. For a massive particle, $E/\mu$ is the energy per unit mass measured at infinity ($E/\mu > 1$ means unbound).

### Jacobi Separation

Instead of building a set of coupled ordinary differential equations to express the motion of a test particle, Carter uses the Hamilton-Jacobi equation to express the motion as a single first order partial differential equation. But this only helps if that equation turns out to be separable.

The idea is to look for a single function $S(\lambda, x)$, the action, whose gradient gives the momentum everywhere:

$$
p_i = \frac{\partial S}{\partial x^i}
$$

Plugging this into the Hamiltonian gives the Hamilton-Jacobi equation:

$$
\frac{\partial S}{\partial \lambda} + H\!\left(x, \frac{\partial S}{\partial x}\right) = 0
\quad\Longrightarrow\quad
\frac{\partial S}{\partial \lambda} = -\frac{1}{2}g^{ij} ( \frac{\partial S}{\partial x^i} - \epsilon A_i ) ( \frac{\partial S}{\partial x^j} - \epsilon A_j )
$$

We already know a lot about $S$:

- $H = -\frac{1}{2}\mu^2$ is constant, so $\partial S/\partial\lambda = \frac{1}{2}\mu^2$.
- $p_v = -E$ and $p_\phi = \Phi$ are constant, so $S$ must be linear in $v$ and $\phi$.
- Whatever is left can only depend on $r$ and $\theta$. Carter's guess is that it splits into a sum of a function of $r$ and a function of $\theta$.

So if there is a separable solution it must take the form:

$$
S = \frac{1}{2}\mu^2\lambda - Ev + \Phi\phi + S_\theta(\theta) + S_r(r)
$$

with $p_\theta = dS_\theta/d\theta$ and $p_r = dS_r/dr$.

Inserting this into the Hamiltonian means replacing each momentum $p_i$ with $\partial S/\partial x^i$. If the resulting Hamilton-Jacobi equation separates into an $r$-dependent piece and a $\theta$-dependent piece, then the motion admits an additional conserved quantity. This separation constant is the Carter constant.

Lets actually do it. We need the inverse of the ingoing EF metric. You can get it either by inverting the matrix above, or by transforming the (much easier to invert) Boyer-Lindquist inverse metric as a contravariant tensor, $g^{-1}_{\mathrm{EF}} = J\, g^{-1}_{\mathrm{BL}} J^T$. In the order $(v, r, \theta, \phi)$:

$$
g^{ij} =
\frac{1}{\rho^2}
\begin{pmatrix}
a^2\sin^2\theta & r^2+a^2 & 0 & a
\\
r^2+a^2 & \Delta & 0 & a
\\
0 & 0 & 1 & 0
\\
a & a & 0 & \dfrac{1}{\sin^2\theta}
\end{pmatrix}
$$

Notice that every term shares a common factor of $\rho^{-2}$, and there is no $\Delta$ in any denominator. Both facts are what make everything below work.

To keep things short, define the *kinetic* momenta $\pi_i = p_i - \epsilon A_i$:

$$
\pi_v = -E - \frac{\epsilon e r}{\rho^2}
\qquad
\pi_\phi = \Phi + \frac{\epsilon a e r\sin^2\theta}{\rho^2}
\qquad
\pi_r = \frac{dS_r}{dr}
\qquad
\pi_\theta = \frac{dS_\theta}{d\theta}
$$

Multiplying $2H = g^{ij}\pi_i\pi_j$ through by $\rho^2$ and grouping terms:

$$
2\rho^2 H =
\underbrace{\Delta\,\pi_r^2 + 2\left[(r^2+a^2)\pi_v + a\pi_\phi\right]\pi_r}_{\text{radial}}
+
\underbrace{\pi_\theta^2 + a^2\sin^2\theta\,\pi_v^2 + 2a\pi_v\pi_\phi + \frac{\pi_\phi^2}{\sin^2\theta}}_{\text{polar}}
$$

The radial bracket simplifies using $r^2 + a^2 - a^2\sin^2\theta = \rho^2$:

$$
(r^2+a^2)\pi_v + a\pi_\phi
= -E(r^2+a^2) + a\Phi - \frac{\epsilon e r}{\rho^2}\left(r^2+a^2-a^2\sin^2\theta\right)
= -\left(E(r^2+a^2) - a\Phi + \epsilon e r\right)
\equiv -P(r)
$$

The last three polar terms are a perfect square, and the charge terms cancel out of it entirely:

$$
a^2\sin^2\theta\,\pi_v^2 + 2a\pi_v\pi_\phi + \frac{\pi_\phi^2}{\sin^2\theta}
= \left(a\sin\theta\,\pi_v + \frac{\pi_\phi}{\sin\theta}\right)^2
= \left(aE\sin\theta - \frac{\Phi}{\sin\theta}\right)^2
$$

So

$$
2\rho^2 H = \Delta\left(\frac{dS_r}{dr}\right)^2 - 2P\frac{dS_r}{dr} + \left(\frac{dS_\theta}{d\theta}\right)^2 + \left(aE\sin\theta - \frac{\Phi}{\sin\theta}\right)^2
$$

Now set $2H = -\mu^2$, so the left side becomes $-\mu^2\rho^2 = -\mu^2 r^2 - \mu^2 a^2\cos^2\theta$. This is the step where the magic happens: the factor of $\rho^2$ that was coupling $r$ and $\theta$ splits cleanly into an $r$ piece and a $\theta$ piece. Moving every $\theta$ term to the left and every $r$ term to the right:

$$
(\frac{dS_\theta}{d\theta})^2 + a^2 \mu^2 \cos^2\theta + (aE\sin\theta - \Phi \sin^{-1}\theta)^2 =
-\Delta (\frac{dS_r}{dr})^2 + 2( E( r^2 + a^2 ) - a \Phi + \epsilon e r )(\frac{dS_r}{dr}) - \mu^2r^2
$$

The left side depends only on $\theta$ and the right side depends only on $r$. Since $r$ and $\theta$ can be varied independently, the only way the two sides can be equal everywhere is if both are equal to the same constant. Thus, both sides are equal to a new constant of motion $K$.

$$
K = p_\theta^2 + a^2 \mu^2 \cos^2\theta + (aE\sin\theta - \Phi \sin^{-1}\theta)^2
\\
K = -\Delta p_r^2 + 2( E( r^2 + a^2 ) - a \Phi + \epsilon e r )p_r - \mu^2r^2
$$

This new constant of motion is harder to intuit. The easiest way in is to set $a = 0$: then $K = p_\theta^2 + \Phi^2/\sin^2\theta$, which is exactly the total angular momentum squared $L^2$. In Schwarzschild that is conserved because of spherical symmetry. Kerr is only axisymmetric, so $L^2$ is no longer conserved, but $K$ survives as its generalization. Unlike $E$ and $\Phi$, $K$ is quadratic in the momenta and does not come from any coordinate symmetry; it comes from a *hidden* symmetry (a Killing tensor, which we will meet again in the Marck section).

Also notice what separability buys us: $p_\theta$ depends only on $\theta$ and $p_r$ depends only on $r$. The $r$ and $\theta$ motions are no longer coupled to each other except through the overall factor $\rho^2$ that converts between $\lambda$ and the parameter they are each "naturally" evolving in. This will become explicit when we introduce Mino time below.

With $K$ we can now solve for $\frac{dS_\theta}{d\theta}$ and $\frac{dS_r}{dr}$ at any point during the system's evolution. Here we define two new functions $\Theta(\theta)$ and $R(r)$.

For $\theta$, the first equation gives directly $p_\theta^2 = K - a^2\mu^2\cos^2\theta - (aE\sin\theta - \Phi/\sin\theta)^2 \equiv \Theta$. To get it into Carter's form, expand the square using $\sin^2\theta = 1 - \cos^2\theta$ and $1/\sin^2\theta = 1 + \cos^2\theta/\sin^2\theta$:

$$
\left(aE\sin\theta - \frac{\Phi}{\sin\theta}\right)^2
= a^2E^2\sin^2\theta - 2aE\Phi + \frac{\Phi^2}{\sin^2\theta}
= (\Phi - aE)^2 - a^2E^2\cos^2\theta + \frac{\Phi^2\cos^2\theta}{\sin^2\theta}
$$

so, defining $Q = K - (\Phi - aE)^2$, every remaining term is proportional to $\cos^2\theta$.

For $r$, the second equation is a quadratic in $p_r$:

$$
\Delta p_r^2 - 2P p_r + (\mu^2 r^2 + K) = 0
\quad\Longrightarrow\quad
p_r = \frac{P \pm \sqrt{P^2 - \Delta(\mu^2r^2 + K)}}{\Delta}
$$

We absorb the $\pm$ into the sign of $\sqrt{R}$. Putting it together:

$$
\frac{dS_\theta}{d\theta} = \sqrt{\Theta}
$$

$$
\frac{dS_r}{dr} = \frac{P + \sqrt{R}}{\Delta}
$$

where

$$
\begin{aligned}
\Theta(\theta) &= Q - \cos^2\theta \left[ a^2(\mu^2 - E^2) + \frac{\Phi^2}{\sin^2\theta} \right]
\\
R(r) &= P^2 - \Delta \left(\mu^2 r^2 + K\right)
\\
P(r) &= E(r^2+a^2)-a\Phi+\epsilon e r
\\
Q &= K-(\Phi-aE)^2
\end{aligned}
$$

For $a=0$, $Q = L^2 - \Phi^2 = L_x^2 + L_y^2$: the part of the angular momentum that tilts the orbit out of the equatorial plane. An equatorial orbit has $Q = 0$. In the literature you will see both $K$ and $Q$ called "the Carter constant", so check which one a paper means.

Thus, the final solution for the Jacobi action is:

$$
S = \frac{1}{2} \mu^2 \lambda -
E v + \Phi\phi +
\int^\theta \sqrt{\Theta}d\theta +
\int^r \frac{P}{\Delta}dr +
\int^r \frac{\sqrt{R}}{\Delta}dr
$$

An important thing to keep in mind is that the sign of the two square roots are independent of each other. In fact, managing these signs becomes critical in correctly modeling the motion of a particle in the system.

The integrated forms of the geodesic and orbit equations can now be obtained by using the fact that partial derivatives of the Jacobi action wrt the constants of motion are themselves constant. Thus, by differentiating wrt $K$, $\mu$, $E$, and $\Phi$, we obtain:

$$
\begin{aligned}
& \int^\theta \frac{d\theta}{\sqrt{\Theta}} = \int^r \frac{dr}{\sqrt{R}}
\\
\lambda &= \int^\theta \frac{a^2\cos^2{\theta}}{\sqrt{\Theta}}d\theta +
\int^r \frac{r^2}{\sqrt{R}}dr
\\
v &= \int^\theta \frac{-a(aE\sin^2{\theta}-\Phi)}{\sqrt{\Theta}}d\theta +
\int^r \frac{r^2+a^2}{\Delta}(1+\frac{P}{\sqrt{R}})dr
\\
\phi &= \int^\theta \frac{-(aE-\Phi \sin^{-2}{\theta})}{\sqrt{\Theta}}d\theta + \int^r \frac{a}{\Delta}(1+\frac{P}{\sqrt{R}})dr
\end{aligned}
$$

We can re-express this in terms of the first-order differential system:

$$
\begin{aligned}
\rho^2 \dot{\theta} &= \sqrt{\Theta}
\\
\rho^2 \dot{r} &= \sqrt{R}
\\
\rho^2 \dot{v} &= -a(aE\sin^2\theta - \Phi) + (r^2+a^2)\Delta^{-1}(P+\sqrt{R})
\\
\rho^2 \dot{\phi} &= -(aE -\Phi\sin^{-2}{\theta}) + a\Delta^{-1}(P+\sqrt{R})
\end{aligned}
$$

Going from the integrated form to this differential form: along the trajectory, the first equation says $d\theta/\sqrt{\Theta} = dr/\sqrt{R}$, and the $\lambda$ equation then says $d\lambda = (a^2\cos^2\theta + r^2)\,dr/\sqrt{R} = \rho^2\,dr/\sqrt{R}$, i.e. $\rho^2\dot{r} = \sqrt{R}$. Similarly $\rho^2\dot{\theta} = \sqrt{\Theta}$, and the $v$ and $\phi$ equations follow by multiplying their integrands by $d\theta/d\lambda$ or $dr/d\lambda$.

We can also cross-check against Hamilton's equations $\dot{x}^i = g^{ij}\pi_j$ using the EF inverse metric from before. For $\theta$: $\rho^2\dot{\theta} = \pi_\theta = \sqrt{\Theta}$. For $r$: $\rho^2\dot{r} = \Delta\pi_r + (r^2+a^2)\pi_v + a\pi_\phi = \Delta\frac{P+\sqrt{R}}{\Delta} - P = \sqrt{R}$, using the same simplification of the radial bracket as before. For $v$:

$$
\rho^2\dot{v} = (r^2+a^2)\pi_r + a^2\sin^2\theta\,\pi_v + a\pi_\phi
= (r^2+a^2)\frac{P+\sqrt{R}}{\Delta} - a^2E\sin^2\theta + a\Phi
$$

(the charge terms cancel again, by the same $r^2 + a^2 - a^2\sin^2\theta = \rho^2$ trick), which matches the $\dot{v}$ equation above. The $\dot{\phi}$ equation follows the same way from the last row of $g^{ij}$.

A quick sanity check on the signs: a massive particle at rest far from the hole has $E \approx \mu$, $\Phi \approx 0$, $K \approx 0$, so $P \approx \mu r^2$, $\sqrt{R} \approx 0$, and $\rho^2\dot{v} \approx r^2 \cdot \mu r^2 / r^2$, i.e. $\dot{v} \approx \mu$. Time runs forward at the expected rate ($dv/d\tau = 1$).

**Principal null directions, revisited.** Now we can verify the PND claim from earlier. Take a photon ($\mu = 0$) with $\Phi = aE\sin^2\theta_0$ and $K = 0$ (and $\epsilon = 0$). Then $P = E(r^2 + a^2 - a^2\sin^2\theta_0) = E\rho^2$ at $\theta_0$, $R = P^2$, so $\sqrt{R} = \pm P$, and $\Theta(\theta_0) = 0$ with $\Theta'(\theta_0) = 0$, so the photon stays at $\theta = \theta_0$. In Boyer-Lindquist coordinates the $t$ and $\phi$ equations are $\rho^2\dot{t} = -a(aE\sin^2\theta - \Phi) + (r^2+a^2)P/\Delta$ and $\rho^2\dot{\phi} = -(aE - \Phi\sin^{-2}\theta) + aP/\Delta$. The angular terms vanish, leaving $\rho^2\dot{t} = (r^2+a^2)P/\Delta$, $\rho^2\dot{\phi} = aP/\Delta$, $\rho^2\dot{r} = \pm P$. Dividing gives exactly $dt/dr = \pm(r^2+a^2)/\Delta$ and $d\phi/dr = \pm a/\Delta$.

### Turning Points

An important property of the separated potentials $R(r)$ and $\Theta(\theta)$ is that their zeros define turning points of the motion. Simply put, a turning point is a point where one component of the geodesic velocity shrinks to zero and then reverses sign.

For example, imagine a massive test particle in a scattering orbit through the black hole's gravitational field. As the particle falls toward the black hole, $r$ decreases and $dr/d\lambda$ is negative. But if the particle does not plunge, it eventually reaches a minimum radius, slingshots around the black hole, and escapes back toward infinity. At that minimum radius,

$$
R(r)=0
$$

and therefore

$$
\frac{dr}{d\lambda}=0.
$$

After this point, the radial branch changes sign and $r$ begins increasing again. This is a radial turning point.

The same idea applies to polar motion. A particle in a non-equatorial orbit oscillates above and below the equatorial plane as it winds around the black hole. The peaks and valleys of this oscillation occur where

$$
\Theta(\theta)=0
$$

and therefore

$$
\frac{d\theta}{d\lambda}=0.
$$

These are latitudinal turning points.

In Carter's separated equations, we account for turning points by tracking the signs of $\sqrt{R}$ and $\sqrt{\Theta}$. When the particle reaches a radial or latitudinal turning point, the corresponding sign must flip:

$$
\sqrt{R} \rightarrow -\sqrt{R}
$$

or

$$
\sqrt{\Theta} \rightarrow -\sqrt{\Theta}.
$$

Ordinary turning points correspond to simple roots of $R$ or $\Theta$. These are not singularities; the motion reverses in finite affine parameter. Double roots are different. They describe limiting orbits that are approached asymptotically, such as unstable spherical photon orbits, and the affine parameter can diverge. This is covered in the next section.

These turning points are also the basis for analytically solving the system. The roots of $R(r)$ and $\Theta(\theta)$ determine the allowed regions of motion and eventually lead to the elliptic integrals used in the analytic geodesic solutions.


### Geodesic Completeness

Now we consider if we can extend the geodesics to unbounded values of $\lambda$. We can easily see that any geodesic can be extended indefinitely unless it reaches the singularity at $\rho^2=0$ or unless one of the integrals diverges.

We can show that $d\lambda$ diverges if $\Theta$ and $R$ are zeros by taking our $\dot{\theta}$ and $\dot{r}$ equations and solving for $d\lambda$:

$$
\begin{aligned}
d\lambda &= \rho^2(\frac{d\theta}{\sqrt{\Theta}})
\\
d\lambda &= \rho^2(\frac{dr}{\sqrt{R}})
\end{aligned}
$$

At first glance, this makes it look like zeros of $R$ or $\Theta$ are singular. Usually they are not. A simple zero of $R$ or $\Theta$ represents a turning point where the corresponding component of motion momentarily vanishes and then reverses sign. The affine-parameter integral remains finite.

A divergence only occurs when the zero is repeated; a *double zero*. For example, if

$$
R(r) \sim (r-r_0)^2,
$$

then

$$
\int \frac{dr}{\sqrt{R(r)}}
\sim
\int \frac{dr}{|r-r_0|}
$$

diverges. This describes a geodesic that asymptotically approaches a limiting orbit, such as the spherical photon orbit mentioned in the section above.

Thus, the dangerous surfaces are not ordinary zeros of $R$ and $\Theta$, but true singularities such as $\rho^2=0$, or coordinate singularities such as $\Delta=0$ when using a chart that is not regular at the horizon (like using retarded time for an ingoing path).

Lets shift our focus to the coordinate singularity at $\Delta=0$ in the $\dot{v}$ and $\dot{\phi}$ equations. You might be thinking "why did we bother transforming to ingoing Eddington-Finkelstein coordinates if it doesn't remove divergence at the horizons? Wasn't that the point?" The truth is that making the metric regular at the horizons does not guarantee the geodesic equations are also regular *as written*. The divergence here is only apparent: the combination $(P + \sqrt{R})/\Delta$ is really just $p_r$, and we will see that it stays finite for ingoing motion.

Thankfully the fix is quite simple. The trick is to multiply the numerator and denominator by the "conjugate" $P - \sqrt{R}$:

$$
\frac{P + \sqrt{R}}{\Delta}
= \frac{(P + \sqrt{R})(P - \sqrt{R})}{\Delta\,(P - \sqrt{R})}
= \frac{P^2 - R}{\Delta\,(P - \sqrt{R})}
$$

From the definition of $R$, $P^2 - R = \Delta(\mu^2r^2 + K)$, so the $\Delta$ cancels *exactly*:

$$
\frac{P + \sqrt{R}}{\Delta} = \frac{\mu^2 r^2 + K}{P - \sqrt{R}}
$$

This is an exact identity, not an approximation. Now look at what happens at a horizon. There $\Delta = 0$ so $R = P^2$ and $\sqrt{R} = \pm P$. There are two cases:

- $\sqrt{R} \rightarrow -P$ (the radial velocity has the opposite sign to $P$): the denominator $P - \sqrt{R} \rightarrow 2P$, which is finite. The ingoing equations are regular.
- $\sqrt{R} \rightarrow +P$: the denominator goes to zero and the ingoing chart genuinely fails. This is a particle crossing the horizon in the "wrong" direction for this chart.

For a future-directed particle crossing the outer (future) horizon we always have $P > 0$ and $\dot{r} < 0$, so we are in the first case. Substituting the identity gives us the horizon-regular ingoing equations:

$$
\begin{aligned}
\rho^2 \dot{v} &= -a(aE\sin^2\theta - \Phi) + \frac{(r^2+a^2)(\mu^2r^2 + K)}{P - \sqrt{R}}
\\
\rho^2 \dot{\phi} &= -(aE -\Phi\sin^{-2}{\theta}) + \frac{a(\mu^2r^2 + K)}{P - \sqrt{R}}
\end{aligned}
$$

You will often see this written using a series expansion in $\Delta$ instead, $\sqrt{R} = -P\left(1 - \frac{\Delta(\mu^2r^2+K)}{2P^2} + O(\Delta^2)\right)$, which gives $P - \sqrt{R} \approx 2P$ and replaces the denominator above with $2P$. That form is only correct to leading order near the horizon, so don't use it for the whole trajectory. The exact form costs nothing extra.

Implementation note: the two forms are algebraically identical but have different numerical behavior. When $\sqrt{R}$ and $P$ have opposite signs, $(P + \sqrt{R})/\Delta$ subtracts two nearly equal numbers near the horizon, so use $(\mu^2r^2+K)/(P - \sqrt{R})$. When they have the same sign, $P - \sqrt{R}$ is the one that cancels, so use $(P + \sqrt{R})/\Delta$ (this case only blows up if you are actually crossing a horizon outward, where the chart is genuinely bad).

Now we can see that an ingoing geodesic can smoothly traverse the future outer horizon.

What happens next is where Kerr gets strange. Between the horizons ($r_- < r < r_+$) everything must move inward, so the particle reaches the inner horizon $r_-$. Whether the ingoing chart is regular there depends on the sign of $P(r_-)$, by the same two cases above: depending on the particle's $E$ and $\Phi$, it may cross through one part of the inner horizon or another, and only one of those is covered smoothly by the ingoing chart. Past $r_-$ things get even stranger. Unlike Schwarzschild, most infalling observers do *not* hit the singularity. The ring singularity at $\rho^2 = 0$ ($r = 0$, $\theta = \pi/2$) can only be hit by equatorial motion; almost every other geodesic reaches a radial turning point inside $r_-$ and heads back *outward* through $r_-$, which is again the bad case for the ingoing chart. Some even pass through the disk inside the ring ($r = 0$, $\theta \ne \pi/2$) into a negative-$r$ region of the extended Kerr-Newman geometry.

For our renderer, we stop the simulation at (or just before) the inner horizon $r_-$. This is not only a convenience. The inner horizon of the exact Kerr-Newman solution is known to be unstable: any infalling radiation is infinitely blueshifted there ("mass inflation"), so the smooth region beyond $r_-$ is not believed to exist inside a real black hole.

For completeness, here is how Carter handles the rest of the maximal extension. He introduces additional coordinate patches to describe horizon crossings with the opposite null orientation. The starting point is the remark that the ingoing transformation can be mirrored by reversing the sign of $t$ and $\phi$, giving a new time and angle coordinate $w$ and $\tilde{\phi}$:

$$
\begin{aligned}
dt &= -dw + (r^2+a^2)\Delta^{-1}dr
\\
d\phi &= -d\tilde{\phi} + a\Delta^{-1}dr
\end{aligned}
$$

That is, $w = -t + r_\ast$ and $\tilde{\phi} = -\phi + r_\sharp$, a time-reversed retarded time. Since the Boyer-Lindquist metric is unchanged by $(t, \phi) \rightarrow (-t, -\phi)$, substituting into it gives the Kerr-Newman metric in exactly the same form as the ingoing EF metric, just with $(w, \tilde{\phi})$ in place of $(v, \phi_v)$. The transformation between these two forms can be given directly as:

$$
\begin{aligned}
dv + dw &= 2(r^2 + a^2)\Delta^{-1}dr
\\
d\phi_v + d\tilde{\phi} &= 2a\Delta^{-1}dr
\end{aligned}
$$

This means transforming from an ingoing coordinate system to a (time-reversed) outgoing one, or vice-versa. It is convenient that this transformation is its own inverse.

Carrying the same steps through ($\dot{w} = -\dot{t} + (r^2+a^2)\dot{r}/\Delta$, and similarly for $\tilde{\phi}$), the $w$ chart equations are:

$$
\begin{aligned}
\rho^2 \dot{w} &= a(aE\sin^2\theta - \Phi) + (r^2+a^2)\frac{\sqrt{R} - P}{\Delta}
= a(aE\sin^2\theta - \Phi) - \frac{(r^2+a^2)(\mu^2r^2 + K)}{P + \sqrt{R}}
\\
\rho^2 \dot{\tilde{\phi}} &= (aE -\Phi\sin^{-2}{\theta}) + a\frac{\sqrt{R} - P}{\Delta}
= (aE -\Phi\sin^{-2}{\theta}) - \frac{a(\mu^2r^2 + K)}{P + \sqrt{R}}
\end{aligned}
$$

These are regular exactly in the case where the ingoing chart fails ($\sqrt{R} \rightarrow +P$). So at every horizon crossing, one of the two charts is regular, and a simulation can hop between them using the transformation above. (The $r$ and $\theta$ equations are the same in both charts.)

Finally, we have obtained our maximally extended geodesics. So long as our observer/photon doesn't hit the singularity at $\rho^2=0$, its affine parameter $\lambda$ can run all the way from $-\infty$ to $+\infty$.

### Mino Time

We can make our lives easier by defining a new affine-like parameter $\gamma$, commonly called Mino time, by

$$
d\lambda = \rho^2d\gamma,
$$

or equivalently

$$
d\gamma = \frac{d\lambda}{\rho^2}.
$$

Then

$$
\frac{dx^\mu}{d\gamma} =
\rho^2\frac{dx^\mu}{d\lambda}.
$$

This removes the explicit $\rho^2$ factors from the separated equations. Notice that $r$ now evolves using only $r$, and $\theta$ using only $\theta$: in Mino time the radial and polar motions are completely decoupled. (Here we use the horizon-regular ingoing form of the $v$ and $\phi$ equations.) To recover proper time, integrate $d\tau = \mu\, d\lambda = \mu\rho^2 d\gamma$ alongside.

$$
\begin{aligned}
\frac{d\theta}{d\gamma} &= \sqrt{\Theta}
\\
\frac{dr}{d\gamma} &= \sqrt{R}
\\
\frac{dv}{d\gamma} &= -a(aE\sin^2\theta - \Phi) + \frac{(r^2+a^2)(\mu^2r^2 + K)}{P - \sqrt{R}}
\\
\frac{d\phi}{d\gamma} &= -(aE -\Phi\sin^{-2}{\theta}) + \frac{a(\mu^2r^2 + K)}{P - \sqrt{R}}
\end{aligned}
$$

### Implementation notes

Numerically, crossing turning points by tracking the sign of the potential functions can be a little unreliable if we simply wait for $R$ or $\Theta$ to cross zero. With finite step sizes, the integrator may step slightly past $R=0$ or $\Theta=0$, causing the square root to become invalid before we ever had the chance to flip the sign. A more robust approach is to integrate the radial and polar motion as a second-order system. This removes the need to manually flip signs at turning points. The velocity $dr/d\gamma$ or $d\theta/d\gamma$ naturally passes through zero and reverses direction under the acceleration term.

To get the second-order equations, differentiate $(dr/d\gamma)^2 = R(r)$ wrt $\gamma$. The chain rule gives $2\frac{dr}{d\gamma}\frac{d^2r}{d\gamma^2} = \frac{dR}{dr}\frac{dr}{d\gamma}$, and dividing out $\frac{dr}{d\gamma}$ leaves $\frac{d^2r}{d\gamma^2} = \frac{1}{2}\frac{dR}{dr}$. The same holds for $\theta$. Note the sign of the square root has disappeared: it is now carried by the velocity. Expanding the derivatives with $\frac{dP}{dr} = 2Er + \epsilon e$ and $\frac{d\Delta}{dr} = 2r - 2m$:

$$
\begin{aligned}
\frac{d^2r}{d\gamma^2} &= \frac{1}{2}\frac{dR}{dr} = \frac{dP}{dr}P - \frac{1}{2} \frac{d\Delta}{dr}(\mu^2r^2+K) - \Delta\mu^2r
\\
\frac{d^2\theta}{d\gamma^2} &= \frac{1}{2}\frac{d\Theta}{d\theta} = -\left(aE\sin\theta-\frac{\Phi}{\sin\theta}\right)
\left(aE\cos\theta+\frac{\Phi\cos\theta}{\sin^2\theta}\right)
+
a^2\mu^2\sin\theta\cos\theta
\end{aligned}
$$


# Tracking the Observer

## Types of Observers

Something that will be very important this time around is the distinction between time-like, null, and space-like paths through spacetime.

Throughout these notes, we use the metric signature

$$
(-,+,+,+)
$$

With this convention, the squared norm of a vector $v^\mu$ is

$$
g_{\mu\nu}v^\mu v^\nu
$$

Massive observers follow time-like curves. If the curve is parameterized by proper time $\tau$, then the observer's 4-velocity is

$$
v^\mu = \frac{dx^\mu}{d\tau}
$$

and its normalization is

$$
g_{\mu\nu}v^\mu v^\nu = -1
$$

More generally, if we use an affine parameter $\lambda$ related to proper time by $\tau=\mu\lambda$, then

$$
g_{\mu\nu}\dot{x}^\mu\dot{x}^\nu = -\mu^2
$$

Massless particles, such as photons, follow null curves. Since photons experience no proper time, it is better to describe them using a null 4-velocity. Their defining condition is

$$
g_{\mu\nu}v^\mu v^\nu = 0
$$

Finally, space-like curves have positive norm:

$$
g_{\mu\nu}v^\mu v^\nu > 0
$$

A physical observer cannot follow a space-like worldline; doing so would require faster-than-light motion and would break causality. If a massive observer's 4-velocity ever becomes space-like in the simulation, then something has gone wrong.

It is important to remain conscious of the type of observer or test particle being modeled. Massive observers require time-like worldlines and local orthonormal tetrads. Photons follow null geodesics and are usually evolved using null tangent vectors or momenta. Keeping these cases separate is crucial for deriving equations, constructing frames, and debugging the simulation.


## Tetrads

This time around we are going to dive much deeper into tetrads and attempt to intuit them a bit better.

A good way to think of tetrads is to think of them as the inertial frame of the observer. Lets say we have an observer moving through space. This observer has both a time-like velocity and an orientation. This is actually all you need to describe an observer and its all encapsulated in this tetrad tensor.

Keep in mind, we are working with tetrads for massive observers. Tetrads for massless particles like photons need their own special null tetrad which we will not cover in these notes.

The first column of the tensor encodes the 4-velocity, and the other three store the bases vectors. Thats it!

More formally, a tetrad is defined as a set of four linearly independent, orthonormal vector fields $e_\mu^\alpha$ that map a local inertial frame (flat Minkowski space) to the tangent space of the curved manifold.

$$
g_{\mu\nu}e_{a}{}^\mu e_{b}{}^\nu =
\eta_{ab}
$$

If our observer emits a ray in his local Minkowski space, we can transform that ray into global curved space with the tetrad.

### ZAMO Tetrad

The first tetrad we can build is the ZAMO tetrad, which stands for *zero angular momentum observer*.

It uses one more shorthand, which must not be confused with $\rho^2$:

$$
A = (r^2+a^2)^2 - \Delta a^2\sin^2\theta
$$

(this is the same combination that appears in $g_{\phi\phi}$ of the EF metric).

$$
e_{(a)}{}^\mu =
\left[
\begin{array}{cccc}
\vert & \vert & \vert & \vert \\
e_{(0)}^\mu & e_{(1)}^\mu & e_{(2)}^\mu & e_{(3)}^\mu \\
\vert & \vert & \vert & \vert
\end{array}
\right] =
\begin{pmatrix}
\sqrt{\dfrac{A}{\rho^2\Delta}} & 0 & 0 & 0
\\
0 & \sqrt{\dfrac{\Delta}{\rho^2}} & 0 & 0
\\
0 & 0 & -\dfrac{1}{\sqrt{\rho^2}} & 0
\\
\dfrac{(2mr-e^2)a}{\sqrt{A\rho^2\Delta}} & 0 & 0 &
\dfrac{1}{\sin\theta}\sqrt{\dfrac{\rho^2}{A}}
\end{pmatrix}
$$

In the previous project we covered orthonormalizing wrt a metric. This ZAMO tetrad is not what you get by starting from the coordinate basis with $e_0 = \partial_t = (1,0,0,0)$ and Gram-Schmidt orthonormalizing wrt the metric. That would give the *static* observer, which only exists outside the ergosphere (inside it $\partial_t$ is spacelike and can't be normalized to $-1$). Instead, the ZAMO's $e_0$ is the unit normal to the surfaces of constant $t$, which works out to $e_0 \propto \partial_t + \omega\,\partial_\phi$ with $\omega = -g_{t\phi}/g_{\phi\phi} = a(2mr - e^2)/A$ (check: the ratio of the $\phi$ and $t$ entries of the first column is exactly $\omega$). If you look at $e_0$ (the first column), you will see that it has been finely tuned to perfectly account for the frame dragging of the black hole: the ZAMO rotates at exactly the frame-dragging rate $\omega$, which is why it has zero angular momentum ($\Phi = 0$).

Another note, this ZAMO tetrad is in Boyer-Lindquist coordinates and will have to be transformed to ingoing Eddington-Finkelstein using the Jacobian provided in a previous section. This also makes this tetrad construction singular at the horizon: it contains $\sqrt{\Delta}$, so it only exists outside $r_+$ (the ZAMO's lapse goes to zero there). But we will tackle that problem in a later section.

### Constructing a tetrad from an orientation and velocity

Suppose we have an observer at some spacetime point $x^\mu$, and we know the metric $g_{\mu\nu}(x)$ at that point. We want to construct a tetrad whose time axis is the observer's 4-velocity and whose spatial axes encode the observer's orientation.

The cleanest way to think about this is in two stages.

1. Construct a static tetrad $E$.
2. Apply a local Lorentz boost to encode velocity.

Lets start on the first step:

We construct a local rest-frame tetrad $E_{(a)}{}^\mu$ satisfying $g_{\mu\nu}E_{(a)}{}^\mu E_{(b)}{}^\nu = \eta_{ab}$.
Here $E_{(0)}{}^\mu $ is the time axis of the local rest frame, while $E_{(1)}{}^\mu, E_{(2)}{}^\mu, E_{(3)}{}^\mu$
are the oriented spatial axes. In practice, these spatial axes can come from the observer's orientation matrix, transformed into the coordinate basis using the appropriate Jacobian, and then orthonormalized with respect to the metric.

Once we have this rest-frame tetrad, we can encode the observer's velocity as a local Minkowski-space 3-velocity
$\vec{v} = (v^{(1)},v^{(2)},v^{(3)}),$
with
$|\vec{v}|^2 < 1.$

The corresponding Lorentz factor is
$\gamma =\frac{1}{\sqrt{1-|\vec{v}|^2}}$
and the observer's local 4-velocity is

$$
\gamma
\begin{pmatrix}
1\\
v^{(1)}\\
v^{(2)}\\
v^{(3)}
\end{pmatrix}
$$

The second step is to construct a local Lorentz boost to encode this velocity into the static tetrad $E_a$. Define $q_i = \gamma v^{(i)}$ and $q^2 = q_1^2+q_2^2+q_3^2$ (when $q=0$ the boost is just the identity, so special-case it to avoid dividing by zero).

$$
\Lambda^{a}_{b} =
\begin{pmatrix}
\gamma & q_1 & q_2 & q_3
\\
q_1 & 1+\dfrac{\gamma-1}{q^2}q_1^2 &
\dfrac{\gamma-1}{q^2}q_1q_2 &
\dfrac{\gamma-1}{q^2}q_1q_3
\\
q_2 &
\dfrac{\gamma-1}{q^2}q_2q_1 &
1+\dfrac{\gamma-1}{q^2}q_2^2 &
\dfrac{\gamma-1}{q^2}q_2q_3
\\
q_3 &
\dfrac{\gamma-1}{q^2}q_3q_1 &
\dfrac{\gamma-1}{q^2}q_3q_2 &
1+\dfrac{\gamma-1}{q^2}q_3^2
\end{pmatrix}
$$

Finally, we commit this Lorentz boost to the static tetrad to obtain our final boosted tetrad:

$$
e = E \Lambda
$$

The important idea is that orientation is handled by the rest-frame tetrad $E_{(a)}{}^\mu$, while velocity is handled by a Lorentz boost in the local Minkowski frame. This keeps the construction conceptually clean: first build a local inertial frame, then boost it.

### The Symmetric Carter Tetrad

In a separate publication in 1968, Brandon Carter published another paper describing the construction of a Boyer-Lindquist tetrad from the symmetries we used earlier. This tetrad is the first step towards analytical parallel transport. This tetrad is called *symmetric*, not because it's a symmetric matrix, but because it treats the ingoing and outgoing directions symmetrically: $e_{(0)} \pm e_{(1)}$ point exactly along the two principal null directions we met earlier, $\left(\frac{r^2+a^2}{\Delta}, \pm 1, 0, \frac{a}{\Delta}\right)$ (up to scale). This alignment with the principal null directions is what makes Carter's separable structure show up so cleanly in its components. [3]

$$
e_{(a)}{}^\mu =
\left[
\begin{array}{cccc}
\vert & \vert & \vert & \vert \\
e_{(0)}^\mu & e_{(1)}^\mu & e_{(2)}^\mu & e_{(3)}^\mu \\
\vert & \vert & \vert & \vert
\end{array}
\right] =
\begin{pmatrix}
\dfrac{r^2+a^2}{\rho\sqrt{\Delta}} & 0 & 0 & -\dfrac{a\sin\theta}{\rho}
\\
0 & \dfrac{\sqrt{\Delta}}{\rho} & 0 & 0
\\
0 & 0 & \dfrac{1}{\rho} & 0
\\
\dfrac{a}{\rho\sqrt{\Delta}} & 0 & 0 & -\dfrac{1}{\rho\sin\theta}
\end{pmatrix}
$$

## Analytic Parallel Transport

If we follow J.-A. Marck's 1983 paper, we have what looks like a great option for analytically transporting our observer tetrad.
By starting with Carters Symmetric tetrad, we can construct a local Lorentz transform that will correctly transport our tetrad to the desired point on the timelike geodesic.

While the paper works with a tetrad in Boyer-Lindquist coordinates, the transformation occurs in the tangent space of the observer's orthonormal frame, so it will work just as well with the same tetrad in ingoing Eddington-Finkelstein coordinates. More specifically, Marck constructs a velocity-adapted orthonormal frame and then applies a residual spatial rotation.

Marck's construction is not just a brute-force integration of the parallel transport equation. It exploits a special hidden symmetry of Kerr spacetime: the existence of the Killing-Yano tensor. This tensor is closely related to the Carter constant and is sometimes described as the geometric reason Kerr geodesic motion is separable.

In this section we shall construct an orthonormal tetrad that is parallel transported along an arbitrary timelike geodesic. Let $\Lambda = ( \lambda_0, \lambda_1, \lambda_2, \lambda_3 )$ be the desired tetrad. Each $\lambda_A^{(a)}$ below is written in components with respect to the Carter tetrad above (the $(a)$ index), not the coordinate basis.

Two conventions to keep in mind. First, Marck normalizes $\mu = 1$, so the dot here means $d/d\tau$ and $E$, $\Phi$, $K$ are per unit mass. For a general $\mu$, divide $E$, $\Phi$ by $\mu$ and $K$ by $\mu^2$ first. The normalization matters: the orthonormality of $\tilde{\lambda}_1$ and $\tilde{\lambda}_3$ below relies on it. Second, Marck's paper treats the uncharged case; with charge the $E(r^2+a^2) - a\Phi$ terms would become $P$. Also note the sign of $e_{(3)}$ in the Carter tetrad above: it is chosen so that the components below come out with Marck's signs. With the opposite sign of $e_{(3)}$, every $(3)$ component below flips sign. You can check it directly: projecting the geodesic velocity onto $e_{(3)}$ gives $(aE\sin\theta - \Phi/\sin\theta)/\rho$, matching $\lambda_0^{(3)}$.

Lets start with the timelike vector which is tangent to the geodesic, so it is automatically parallel transported.

$$
\begin{aligned}
\lambda_{0}^{(0)}
&=
\frac{1}{\sqrt{\Delta \rho^2}}
( E(r^{2}+a^{2}) - a\Phi)
\\
\lambda_{0}^{(1)}
&=
\sqrt{\frac{\rho^2}{\Delta}}\dot{r}
\\
\lambda_{0}^{(2)}
&=
\sqrt{\rho^2}\dot{\theta}
\\
\lambda_{0}^{(3)}
&=
\frac{1}{\sqrt{\rho^2}}
\left( aE\sin\theta - \frac{\Phi}{\sin\theta} \right)
\end{aligned}
$$

The clever part is that the Killing-Yano tensor can be contracted with the geodesic velocity to produce another vector that is also parallel transported. In other words, two of the four tetrad legs are essentially fixed by the geodesic and Kerr's hidden symmetry. I dont know much about the Killing-Yano tensor other than it is responsible for Carter’s hidden constant $K$ and it is being exploited to obtain $\lambda_{2}$.

$$
\begin{aligned}
\lambda_{2}^{(0)}
&=
\sqrt{\frac{\rho^2}{K\Delta}} a \cos\theta \dot{r}
\\
\lambda_{2}^{(1)}
&=
\frac{a\cos\theta}{\sqrt{K\rho^2\Delta}}( E(r^{2}+a^{2}) - a\Phi)
\\
\lambda_{2}^{(2)}
&=
-\frac{r}{\sqrt{K\rho^2}}\left( aE\sin\theta - \frac{\Phi}{\sin\theta} \right)
\\
\lambda_{2}^{(3)}
&=
r \sqrt{\frac{\rho^2}{K}} \dot{\theta}
\end{aligned}
$$

This leaves only the two-dimensional plane orthogonal to both $\lambda_{0}$ and $\lambda_{2}$. Within that plane, the only remaining freedom is an ordinary spatial rotation. Marck chooses two convenient intermediate basis vectors, $\tilde{\lambda}_{1}$ and $\tilde{\lambda}_{3}$, spanning this plane. These intermediate vectors are orthonormal, but they are not individually parallel transported. Instead, as the observer moves along the geodesic, they rotate within their own two-plane.

$$
\begin{aligned}
{\tilde{\lambda}}_{1}^{(0)}
&=
\alpha \sqrt{\frac{\rho^{2}}{K\Delta}}\, r\dot{r}
\\
{\tilde{\lambda}}_{1}^{(1)}
&=
\alpha \frac{r}{\sqrt{K\rho^{2}\Delta}}
\left( E(r^{2}+a^{2}) - a\Phi \right)
\\
{\tilde{\lambda}}_{1}^{(2)}
&=
\beta \frac{a\cos\theta}{\sqrt{K\rho^{2}}}
\left( aE\sin\theta - \frac{\Phi}{\sin\theta} \right)
\\
{\tilde{\lambda}}_{1}^{(3)}
&=
-\beta \sqrt{\frac{\rho^{2}}{K}}\, a\cos\theta\,\dot{\theta}
\end{aligned}
$$

and

$$
\begin{aligned}
{\tilde{\lambda}}_{3}^{(0)}
&=
\alpha \frac{1}{\sqrt{\Delta\rho^{2}}}
\left( E(r^{2}+a^{2}) - a\Phi \right)
\\
{\tilde{\lambda}}_{3}^{(1)}
&=
\alpha \sqrt{\frac{\rho^{2}}{\Delta}}\,\dot{r}
\\
{\tilde{\lambda}}_{3}^{(2)}
&=
\beta \sqrt{\rho^{2}}\,\dot{\theta}
\\
{\tilde{\lambda}}_{3}^{(3)}
&=
\beta \frac{1}{\sqrt{\rho^{2}}}
\left( aE\sin\theta - \frac{\Phi}{\sin\theta} \right)
\end{aligned}
$$

where

$$
\begin{aligned}
\alpha &=
\frac{\sqrt{K-a^{2}\cos^{2}\theta}}{\sqrt{r^{2}+K}}
\\
\beta &= \alpha^{-1}
\end{aligned}
$$

The failure of $\tilde{\lambda}_{1}$ and $\tilde{\lambda}_{3}$ to be parallel transported is therefore captured by a single angular velocity $\dot{\Psi}$. By accumulating this angle and applying the rotation we cancel that residual rotation. The result is a full orthonormal tetrad $\Lambda$ that is parallel transported along the timelike geodesic.

$$
\dot{\Psi} =
\frac{\sqrt{K}}{\rho^2}
\left(
\frac{E(r^{2}+a^{2})-a\Phi}{r^{2}+K} +
a\frac{\Phi-aE\sin^{2}\theta}{K-a^{2}\cos^{2}\theta}
\right)
$$

$$
\begin{pmatrix}
\lambda_{1} \\
\lambda_{3}
\end{pmatrix} =
\begin{pmatrix}
\cos\Psi & -\sin\Psi \\
\sin\Psi & \cos\Psi
\end{pmatrix}
\begin{pmatrix}
\tilde{\lambda}_{1} \\
\tilde{\lambda}_{3}
\end{pmatrix}
$$

As stated earlier, this Lorentz transformation $\Lambda_A^{(a)}$ (which in this case is a tetrad in and of itself) occurs in the tangent space of the observer tetrad $e_{(a)}^\mu$. Thus, it must be committed to the base tetrad:

$$
E_A^\mu = e_{(a)}^\mu \Lambda_A^{(a)}
$$

Unfortunately, there is a practical problem. Marck's construction is tied to Carter's symmetric tetrad, which is naturally written in Boyer--Lindquist coordinates. Although we can transform this tetrad into ingoing Eddington--Finkelstein coordinates using a Jacobian, this does not necessarily make the tetrad itself horizon-regular. The coordinate system is regular at the outer horizon, but the Carter frame being used as the local reference basis is still a Boyer--Lindquist-adapted frame. Near the horizon it becomes singular or infinitely boosted relative to a regular infalling frame.

## Numeric Parallel Transport

The first step is to compute the christoffel symbols. Instead of calculating them symbolically, we are going to compute them numerically using *finite difference method*. More specifically, we will use the *central finite difference* method. We sample the metric a little bit behind and a little bit in front of the observer.

Here is the Christoffel symbol formula:

$$
\Gamma^\mu_{\alpha\beta} = \frac{1}{2}\mathfrak{g}^{\mu\lambda}(\frac{\partial g_{\lambda\alpha}}{\partial x^\beta} + \frac{\partial g_{\lambda\beta}}{\partial x^\alpha} - \frac{\partial g_{\alpha\beta}}{\partial x^\lambda})
$$

The inverse metric $\mathfrak{g}$ is evaluated at the current point. The metric derivatives $\frac{\partial g}{\partial x}$ are estimated with a *central finite difference*. They compute the partial derivative of the metric $g$:

$$
\frac{\partial g_{\alpha\beta}}{\partial x^\sigma}
\approx
\frac{
g_{\alpha\beta}(x+\epsilon \hat{x}_\sigma) -
g_{\alpha\beta}(x-\epsilon \hat{x}_\sigma)
}{2\epsilon}
$$

First you compute all metric derivatives using *central finite difference* then compute and cache all 64 Christoffel symbol components for parallel transport. The parallel transport equations are below:

$$
\begin{aligned}
\frac{d v^\mu}{d\lambda} &= -\Gamma^\mu_{\nu\rho} v^\nu u^\rho,
\\
\quad
u^\rho &= \frac{dx^\rho}{d\lambda}
\end{aligned}
$$

where $v^\mu$ is the vector being transported and $u^\rho$ is the tangent vector to the observer's worldline.

So, before an integration step, we recompute Christoffel symbols, step the observer to the next location on the worldline, recompute the metric at the new position, transport each arm of the tetrad, and orthonormalize wrt the new metric.

# References

[1]
B. Carter,
Global Structure of the Kerr Family of Gravitational Fields,
Physical Review 174, 1559–1571 (1968).

[2]
J.-A. Marck,
Solution to the equations of parallel transport in Kerr geometry; tidal tensor,
Proceedings of the Royal Society A 385, 431–438 (1983).

[3]
B. Carter,
Hamilton-Jacobi and Schrödinger separable solutions of Einstein's equations,
Communications in Mathematical Physics 10, 280–310 (1968).
