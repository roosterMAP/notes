# Rendering A Black Hole


In the previous project, we rendered a simplified version of a morris-Thorne Wormhole using differential geometry. We started with a metric in modified spherical coordinates and derive Christoffell Sybols, then used those to build the first order system of Parallel Transport equations and a second order system of Null Goedesic equations. We also build the Hamiltonian for the system to get a first order system of null geodesic equations that evolve in a phase space. The observer was tracked along null geodesics and the tetrad frame was build by simply transforming the minkowskian basis of the observer with the jacobian of the modified spherical coordinate system we defined, then orthonormalizing with respect to the metric. While this works from a rendering perspective, many laws of physics were ignored. The observer was assumed to have no mass and no velocity. Instead of evolving the observes position and inertialm frame through the curved space-time of the wormhole, we were simply teleporting the player to different points along a null geodesic.

In this project we will faithfully simulate a spinning black hole and the physical timelike path an observer would take while traversing its gravitational field. We will start by exploring the Schwarzschild metric and the coordinate singularity that appears at the event horizon, motivating the need for horizon-penetrating coordinate systems. Then we will analyze the Kerr-Newman metric and weigh the difficulties of solving geodesic equations with regular geometric methods. This will take us into a deep dive of the work of Brandon Carter who discovered how to leverage symmetris in the Kerr-Newman system to elegantly solve both null and timelike geodesics. This seminal paper will be the backbone of our implementation. Once done, we will shift our focus to our innertial observe and do a deep dive of tetrad construction and parallel transport. We will analyze both numerical and analytical methods. Tieing these together will give us a faithful simulation of an observer in a Kerr gravitational field. To complete the picture, we will also briefly explore the physics of accretion flows, relativistic aberration, gravitational redshift, and radiative transfer. Finally, we will discuss analytical solutions of the null geodesic equations, the evaluation of elliptic integrals, and several other techniques used to render black holes in real time.


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


Using geometrized units where $G = c = 1$, the line element for the Schwarzschild metric is

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

## Kerr-Newman metric

The Kerr metric was discovered by Roy Kerr in 1963 as the exact solution describing a rotating black hole, and Ezra Newman and collaborators generalized it in 1965 to include electric charge, producing the Kerr-Newman metric.

This new metric takes into consideration two additional properties of the black hole:
- spin $a$
- charge $e$

We are including $e$ for completness, but in nature we generally assume that no blackholes have charge simply because if they did they would quicly attract matter of the opposite charge and return to neutral immidiatly. We also dont consider the other Standard Model gauge charges such as color charge or weak isospin even though there are fun theories about them in some Big Bang models.

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

The metric uses two new terms: $\rho^2$ and $\Delta$.
- $\rho^2$ is the equivelant of $r^2$ in Schwarzschild. Because Kerr rotates about its splin axis, its geometry is flattened/oblate. $\rho^2$ appropriatly scale based off the spin $a$ and latitudinal angle of approach $\theta$.
- Horizon(s) occure when $\Delta = 0$. You can solve this using the quadratic equation: $r_\pm = m \pm \sqrt{m^2 - a^2 - e^2}$. [Graphing](https://www.desmos.com/calculator/qbt6lvyepg) this will give you two circles. When $a=0$ the inner horizon collapses into a point and the system reduces down to the Schwarzschild case. When $a=m$ we can an extremal black hole where the inner and outer horizon perfectly overlap (assuming $e=0$).

Before we dive in, it is useful to check how Kerr-Newman reduces to simpler metrics.
- If $e=0$, the metric becomes the Kerr metric: a rotating, uncharged black hole.
- If $a=0$, the metric becomes the Reissner-Nordström metric: a charged, non-rotating black hole.
- If both $a=0$ and $e=0$, then $rho^2 = r^2$ and\Delta = r^2 - 2mr = r(r-2m), so the metric reduces to Schwarzschild.

Here is the same metric in matrix form:

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

### The Ergosphere

Lets analyze components of the metric like we did for Schwarzschild.

Remember $ds^2 = -d\tau^2$:

Lets consider a static observer where $dr=d\theta=d\phi=0$.

$$
d\tau = \sqrt{1-\frac{2mr-e^2}{\rho^2}}dt
$$

[Graphing](https://www.desmos.com/calculator/iw5vagheof) this looks very similar to Schwarzchild except $g_{tt}$ depends both on $r$ and $\theta$. So, gravitational time dilation is nolonger purely radial. It depends on the latitude relative to the back holes spin axis.
Setting $g_{tt}=0$ and solving for $r$ gives us:

$$
r_{\text{ergo}} = m + \sqrt{m^2 - a^2\cos^2{\theta}-e^2}
$$

We can add this to our [graph](https://www.desmos.com/calculator/ofqdzk3z9a) of the inner and outer horizon. We get a third larger circle.
Lets dig deeper into the nature of this new ergosphere.

Consider an observe trying to stay fixed at $r$ and $\theta$: $dr = d\theta = 0$.

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

Earlier, we [graphed](https://www.desmos.com/calculator/iw5vagheof) the $d\tau$ where $dr=d\theta=d\phi=0$. We could see the curve was discontinuous at $r=2m$. This is because our solution became complex for $0 < r < 2m$. Now we know why. Having zero angular velocity past this point is not allowed for a physical observer because $g_{tt}$ becomes spacelike past the ergosphere. Thus we define the ergosphere as a region where all observers must corotate with the black hole because past that point they must have some angular velocity. The effect of a spinning black hole giving static observers angular velocity via its rotating spacetime is called *frame dragging* because the inertial frame of the observer is being dragged along the $\phi$ direction.


### Principal Null Directions

In the Schwarzchild case we were able to set $ds^2 = d\theta = d\phi = 0$ to track a radially infalling photon. This gave us a pretty simple expression we could itegrate. This expression is called a *Principal Null Direction*: the path of a photon with no angular momentum falling from infinity. This is not so easy to derive in Kerr because of the frame dragging. We will just have to take this for granted bellow:

$$
\frac{dt}{dr} = \frac{r^2 + a^2}{\Delta}
$$
$$
\frac{d\phi}{dr} = \frac{a}{\Delta}
$$

Integrating $\frac{dt}{dr}$ gives us:

$$
\int\frac{r^2 + a^2}{\Delta}dr = r +
\frac{r_+^2+a^2}{r_+-r_-}\ln|r-r_+| -
\frac{r_-^2+a^2}{r_+-r_-}\ln|r-r_-| + C
$$

When you [graph](https://www.desmos.com/calculator/eruofhwhbc) this you can clearly see the divergence at the inner and outer horizons.

Like with Schwarzchild, this will allow us to define a coordinate transformation that will remove the divergence at the horizons.

### Regularizing Coordinate Transformation

Like with Schwarzchild, we define our advanced/retarded coordinates:

$$
du = dt + (r^2 + a^2)\Delta^{-1}
\\
d\phi_u = d\phi + a\Delta^{-1}
$$

$$
dv = dt - (r^2 + a^2)\Delta^{-1}
\\
d\phi_v = d\phi - a^2\Delta^{-1}
$$

Substituting into the metric gives us:

$$
\begin{aligned}
ds^2 ={}&
-\left(1-\frac{2mr-e^2}{\rho^2}\right)dv^2
+2\,dv\,dr
+\rho^2\,d\theta^2
-\frac{2a(2mr-e^2)\sin^2\theta}{\rho^2}\,dv\,d\phi
\\
&-2a\sin^2\theta\,dr\,d\phi
+
\frac{
\left(r^2+a^2\right)^2
-
a^2\Delta\sin^2\theta
}{
\rho^2
}
\sin^2\theta\,d\phi^2 .
\end{aligned}
$$


In matrix form, the coordinate transformation is:

$$
J^\mu_\nu
=
\frac{\partial x^{\mu}_{\mathrm{EF}}}{\partial x^{\nu}_{\mathrm{BL}}}
=
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

Multiplying against the Boyer-Lindquist Metric gives us:

$$
g_{\mu\nu}
=
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
\frac{
\left(r^2+a^2\right)^2
-
a^2\Delta\sin^2\theta
}{
\rho^2
}
\sin^2\theta
\end{pmatrix}
$$

Moving forward, we will use the ingoing Eddington-Finkelstein time coordinate $v$. This can be confusing when following Brandon Carter's paper, because older Kerr literature does not always match the modern $u$/$v$ naming convention. In modern notation, advanced or ingoing time is usually written as

$$
v = t + r_\ast,
$$

while retarded or outgoing time is written as

$$
u = t - r_\ast.
$$

Carter instead uses the symbol $u$ and refers to it as retarded time, but the sign of his coordinate transformation matches what we would now call the ingoing or advanced coordinate. To avoid confusion, these notes will use the word *ingoing* and the variable $v$ throughout.

Also, Carter does not refer to this as an Eddington-Finkelstein coordinate transformation. The modern terminology was not yet standard in the Kerr literature, so the transformation appears as part of the coordinate construction rather than under the name *Eddington-Finkelstein*.


### Evaluating the ingoing Eddington-Finkelstein Metric

A quick look at our new metric in matrix form makes it clear that it is regular at the horizon. Unlike the Boyer-Lindquist metric, $\Delta$ never appears in the denominator of any term and $\rho^2$ is well behaved.

Lets evaluate what it would take to naivley derive our second order geodesic equations from the Christoffel-symbols like we did for the Morris-Thorne wormhole:

In four dimensions, our Christoffel-symbol $\Gamma_{\alpha\beta}^\mu$ has 4x4x4=64 raw components. Since the lower indices are symmetric, only 40 are independent. For this metric, many are nonzero, and expanding the lower-index symmetry gives dozens of terms in the geodesic equations.

$$
\frac{d^2 x^\alpha}{d\lambda^2 } + \Gamma^\alpha_{\beta\mu} \frac{dx^\beta}{d\lambda} \frac{dx^\mu}{d\lambda} = 0
$$

becomes four coupled second-order equations, each containing a pile of terms with up to 10 velocity-pair terms per equation.

The good news is that the EF Christoffels should also be regular at the horizon, because the EF metric and inverse metric are regular there. The bad news is that they are a steaming pile of algebraic shit.


### Kerr Symmetries

The equations of motion of a test particle of mass $\mu$ and charge $\epsilon$ may be derived by the Lagrangian:

$$
L = \frac{1}{2} g_{ij}\dot{x}^i\dot{x}^j + \epsilon A_i\dot{x}^i
$$

where the dot over a symbol denotes ordinary differentiation wrt an affine parameter $\lambda$. In order to obtain the equations of motion for our test particle, we can relate to propter time by using $\tau = \mu \lambda$. 

Remember our line element for a timelike observer:

$$
ds^2 = g_{ij} dx^i dx^j \quad\text{and}\quad ds^2 = -d\tau^2
$$

By substituting $\mu = \tau/\lambda$ we obtain:

$$
g_{ij}\dot{x}^i\dot{x}^j = -\mu^2
$$


In order to support the electromagnetic field of the black hole with charge $e$, we define $A$ as the covariant vector potential:

$$
A = e \rho^{-2}r(dv-a sin^2\theta d\phi)
$$

Honestly, I know little about electromagnetism and have only included charge for completness. It is enough to know that the vector potential $A$ only shifts the canonical momentum.

In order to transform to a Hamiltonian formulation we introduce the momenta to obtain:

$$
p_i = g_{ij}\dot{x}^j + \epsilon A_i
$$

and thus obtaining the Hamiltonian:

$$
H = \frac{1}{2}g^{ij}(p_i - \epsilon A_i)(p_j - \epsilon A_j))
$$

From our definition of the canonical momentum we get:

$$
p_i - \epsilon A_i = g_{ij}\dot{x}^j
$$

Substituting into out definition of the hamiltonian we get:

$$
H = \frac{1}{2}g^{ij} g_{ik}\dot{x}^k g_{jl}\dot{x}^l = \frac{1}{2}g_{kl} \dot{x}^k \dot{x}^l = -\frac{1}{2} \mu^2
$$


When defining our initial conditions we calculate our canonical momenta. But some of these momenta are conserved quantities that arrise from various symmetries in the system since they correspond to constants in the covector field. If we assume that the black hole is eternal and does not change size and splins at a constant rate, we can assume a conservation of energy symmetry and a conservation of angular momentum symmetry.

Thus we define out first two constants of motion:

$$
E = -p_v
\\
\Phi = p_\phi
$$


### Jacobi Separation

Instead of building a set of coupled ordinary differential equations to express the motion of a test particle, Carter uses the Hamilton-Jacobi equation to express the motion as a single first order partial differential equation. But this can only happen if the system is fully seperable.

The general form of the Hamilton-Jacobi equation is

$$
\frac{\partial S}{\partial \lambda} = \frac{1}{2}g^{ij} ( \frac{\partial S}{\partial x^i} - \epsilon A_i ) ( \frac{\partial S}{\partial x^j} - \epsilon A_j )
$$

If there is a seperable solution it must take the form:

$$
S = -\frac{1}{2}\mu^2\lambda - Ev + \Phi\phi + S_\theta(\theta) + S_r(r)
$$

Where S is the action. The canonical momenta are recovered from $S$ by

$$
p_i = \frac{\partial S}{\partial x^i}
$$


Inserting this into the Hamiltonian means replacing each momentum $p_i$ with $\partial S/\partial x^i$. If the resulting Hamilton-Jacobi equation separates into an $r$-dependent piece and a $\theta$-dependent piece, then the motion admits an additional conserved quantity. This separation constant is the Carter constant.

$$
(\frac{dS_\theta}{d\theta})^2 + a^2 \mu^2 \cos^2\theta + (aE\sin\theta - \Phi \sin^{-1}\theta)^2 =
-\Delta (\frac{dS_r}{dr})^2 + 2( E( r^2 + a^2 ) - a \Phi + \epsilon e r )(\frac{dS_r}{dr}) - \mu^2r^2
$$

Thus, both sides are equal to a new constant of motion $K$.

$$
K = p_\theta^2 + a^2 \mu^2 \cos^2\theta + (aE\sin\theta - \Phi \sin^{-1}\theta)^2
\\
K = -\Delta p_r^2 + 2( E( r^2 + a^2 ) - a \Phi + \epsilon e r )p_r - \mu^2r^2
$$

With $K$ can also solve for $\frac{dS_\theta}{d\theta}$ and $\frac{dS_r}{dr}$ at any point during the systems evolution. Here we define two new functions $\Theta(\theta)$ and $R(r)$.

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
\\[6pt]
R(r) &= P^2 - \Delta \left(\mu^2 r^2 + K\right)
\\[6pt]
P(r) &= E(r^2+a^2)-a\Phi+\epsilon e r
\\[6pt]
Q &= K-(\Phi-aE)^2
\end{aligned}
$$


Thus, the final solution for the Jacobi action is:

$$
S = -\frac{1}{2} \mu^2 \lambda -
E u + \Phi\phi +
\int^\theta \sqrt{\Theta}d\theta +
\int^r \frac{P}{\Delta}dr +
\int^r \frac{\sqrt{R}}{\Delta}dr
$$

An important thing to keep in mind is that the sign of the two square roots are independend of eachother.

The integrated forms of the geodesic and orbit equations can now be obtained by using the fact that partial derivatives of the Jacobi actin wrt the constants of motion are themselves constant. Thus, but differentiating wrt $K$, $\mu$, $E$, and $\Phi$, we obtain:

$$
\begin{aligned}
& \int^\theta \frac{d\theta}{\sqrt{\Theta}} = \int^r \frac{dr}{\sqrt{R}}
\\[6pt]
\lambda &= \int^\theta \frac{a^2\cos^2{\theta}}{\sqrt{\Theta}}d\theta +
\int^r \frac{r^2}{\sqrt{R}}dr
\\[6pt]
u &= \int^\theta \frac{-a(aE\sin^2{\theta}-\Phi)}{\sqrt{\Theta}}d\theta +
\int^r \frac{r^2+a^2}{\Delta}(1-\frac{P}{\sqrt{R}})dr
\\[6pt]
\phi &= \int^\theta \frac{-(aE-\Phi \sin^{-2}{\theta})}{\sqrt{\Theta}}d\theta + \int^r \frac{a}{\Delta}(1-\frac{P}{\sqrt{R}})dr
\end{aligned}
$$

We can re-express this in terms of the first-order differential system by either directly integrating the form above or directly from the deffinition of the constants:

$$
\begin{aligned}
\rho^2 \dot{\theta} &= \sqrt{\Theta}
\\[6pt]
\rho^2 \dot{r} &= \sqrt{R}
\\[6pt]
\rho^2 \dot{u} &= -a(aE\sin^2\theta - \Phi) + (r^2+a^2)\Delta^{-1}(\sqrt{R}-P)
\\[6pt]
\rho^2 \dot{\phi} &= -(aE -\Phi\sin^{-2}{\theta}) + a\Delta^{-1}(\sqrt{R}-P)
\end{aligned}
$$


### Geodesic Completness

The separated equations are

Now we consider if we can extend the geodesics to unbounded values of $\lambda$. We can easily see that any geodesic can be extended indefinitely unless is reaches the singularity at $\rho^2=0$ or unless the one the integrals diverge.

We can show that $d\lambda$ diverges if $\Theta$ and $R$ are zeros by taking our $\dot{\theta}$ and $\dot{r}$ equations and solving for $d\lambda$:

$$
\begin{aligned}
d\lambda &= \rho^2(\frac{d\theta}{\sqrt{\Theta}})
\\[6pt]
d\lambda &= \rho^2(\frac{dr}{\sqrt{R}})
\end{aligned}
$$


At first glance, this makes it look like zeros of $R$ or $\Theta$ are singular. Usually they are not. A simple zero of $R$ or $\Theta$ represents a turning point where the corresponding component of motion momentarily vanishes and then reverses sign. The affine-parameter integral remains finite.

A divergence only occurs when the zero is repeated. For example, if

$$
R(r) \sim (r-r_0)^2,
$$

then

$$
\int \frac{dr}{\sqrt{R(r)}}
\sim
\int \frac{dr}{|r-r_0|}
$$

diverges. This describes a geodesic that asymptotically approaches a limiting orbit, such as a spherical photon orbit.

Thus, the dangerous surfaces are not ordinary zeros of $R$ and $\Theta$, but true singularities such as $\rho^2=0$, or coordinate singularities such as $\Delta=0$ when using a chart that is not regular at the horizon.

Lets shift out focus to the coordinate singularity at $\Delta=0$ in the $\dot{u}$ and $\dot{\phi}$ equations. You might be thinking "why did we bother transforming to ingoing Eddington-Finkelstein cordinates if it doesnt remove divergence at the horizons? Wasnt that the point?" The truth is that making the metric regular at the horizons does not guarantee the geodesic equations are also regular.

Thankfully the fix is quite simple. First we can re-express  $\dot{u}$ and $\dot{\phi}$ as:

$$
\begin{aligned}
\rho^2 \dot{u} &= -a(aE\sin^2\theta - \Phi) + (r^2+a^2)\Delta^{-1}(1+\frac{P}{\sqrt{R}})
\\[6pt]
\rho^2 \dot{\phi} &= -(aE -\Phi\sin^{-2}{\theta}) + a\Delta^{-1}(1+\frac{P}{\sqrt{R}})
\end{aligned}
$$

Then, provided P is nonzero, we obtain the expansion:

$$
\frac{P}{\sqrt{R}} = \pm[1\pm\frac{\mu^2r^2 + K}{2P}(\frac{\Delta}{P}) + O(\frac{\Delta}{P})^2]
$$

The sign depends on which choise fo $\sqrt{R}$ is under consideration.

Substituting with the correct sign will cancel out the $\Delta$ terms and give us:

$$
\begin{aligned}
\rho^2 \dot{u} &= -a(aE\sin^2\theta - \Phi) + \frac{(r^2+a^2)(\mu^2r^2 + K)}{2P^2}
\\[6pt]
\rho^2 \dot{\phi} &= -(aE -\Phi\sin^{-2}{\theta}) + \frac{a(\mu^2r^2 + K)}{2P^2}
\end{aligned}
$$


Now we can see that an ingoing geodesic can smoothly traverse the future horizon.
This will be enough for the case of our renderer. However, this does not cover the entire maximally extended global structure. If a geodesic avoids the ring singularity at $\rho^2=0$, it may pass through $r=0$ into a negative-$r$ region of the extended Kerr-Newman geometry. Carter introduces additional coordinate patches to describe horizon crossings with the opposite null orientation and to cover more of the maximal extension. For the purposes of an infalling observer that terminates at the singularity, this extra patching is not necessary.

The starting point for this transformation is the remark that the invertible form can be extended in a symmetric manner in an inverted direction in terms of a new time and angle coordinate $w$ and $\tilde{\phi}$:

$$
\begin{aligned}
d\hat{t} &= -dw + (r^2+a^2)\Delta^{-1}dr
\\
d\hat{\phi} &= -d\tilde{\phi} + a\Delta^{-1}dr
\end{aligned}
$$

Substituting into the Boyer-Lindquist metric gives us a new form of the Kerr-Newman metric. This transformation between these two forms can be given directly as:

$$
\begin{aligned}
du + dw &= 2(r^2 + a^2)\Delta^{-1}dr
\\
d\phi + d\tilde{\phi} &= 2a\Delta^{-1}dr
\end{aligned}
$$

We can substitute into our $\dot{v}$ and $\dot{\phi}$ equations to get a new set of $\dot{w}$ and $\dot{\tilde{\phi}}$ equations. From there we can once again cancel out the divergent $\Delta$ terms via the series expansion, taking care to respect the sign of $\sqrt(R)$.

Finally, we have obtained out maximally extended geodesics. Solong as our observer/photon doesn't hit the singularity at $\rho^2=0$, he can go from $+\infty$ to $-\infty$.


### Mino Time


We can make our lives easier by defining a new affine-like parameter $\gamma$, commonly called Mino time, by

$$
d\lambda = \rho^2\,d\gamma,
$$

or equivalently

$$
d\gamma = \frac{d\lambda}{\rho^2}.
$$

Then

$$
\frac{dx^\mu}{d\gamma}
=
\rho^2\frac{dx^\mu}{d\lambda}.
$$

This removes the explicit $\rho^2$ factors from the separated equations.

$$
\begin{aligned}
\frac{d\theta}{d\gamma} &= \sqrt{\Theta}
\\[6pt]
\frac{dr}{d\gamma} &= \sqrt{R}
\\[6pt]
\frac{du}{d\gamma} &= -a(aE\sin^2\theta - \Phi) + (r^2+a^2)\Delta^{-1}(\sqrt{R}-P)
\\[6pt]
\frac{d\phi}{d\gamma} &= -(aE -\Phi\sin^{-2}{\theta}) + a\Delta^{-1}(\sqrt{R}-P)
\end{aligned}
$$
