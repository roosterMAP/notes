# Rigid Body Physics

## _Book One_

Notes from [GamePhysicsWeekend](<https://gamephysicsweekend.github.io/>) series.


# Eulers Method

$$x(t) = x_0 + \int_{t_0}^t \left( v_0 + \int_{\tau_0}^\tau a(s)\, ds \right) d\tau$$
$$\text{where } t \text{ and } \tau \text{ are time variables, } v_0 \text{ is initial velocity, and } a \text{ is acceleration.}$$

$$\dot{.\hspace{.05in}.} \quad x(t) = \int_{t_0}^t v(\tau) \, d\tau \quad \text{and} \quad  v(t) = \int_{t_0}^t a(t) \, dt$$

We discritize using Eulers Method:
$$x_{n+1} = x_n + v(t_n) \cdot \Delta t$$
$$v_{n+1} = v_n + a(t_n) \cdot \Delta t$$

Eulers method only works for small and regular time intervals with no repid variations in acceleration. As an alternative to Euler’s method, the Runge-Kutta 4th order method (RK4) achieves much higher accuracy by combining multiple slope estimates within each timestep. But it is more computationally expensive and is not directly compatible with linear solvers for constraints.

Acceleration due to gravity is $$9.8m/s^2$$.
```
v += |0, 0, -9.8| * fDeltaTime
p += v * fDeltaTime
```


# Generalize Gravty as a Force
$$F = m \cdot a \quad \text{Newton's Second Law}$$
$$\dot{.\hspace{.05in}.} \quad F_g = m \cdot -9.8{m/s^2}$$

$$F = dp/dt \quad \text{Force in terms of Momentum (p)}$$
$$p = m \cdot v \quad \text{Momentum is Velocity scaled by Mass}$$

We handle Force indirectly by applying it as an __Impulse__, which is simply Force integrated over time.
So, and __Impulse__(J) is the total Force accumulated over a period of time (change in momentum).
$$J = F \cdot dt$$

To calculate the impulse due to gravity we scale the Force by the time delta.
To apply the impulse, it must get commited to object velocity.
$$J = \Delta p = mv_\text{after} - mv_\text{before} = m \cdot \Delta v$$
$$\dot{.\hspace{.05in}.} \quad \Delta v = J / m$$

1. Force acts over a short time generating an __Impulse__.
2. The __Impulse__ causes a change in Momentum.
3. Momentum causes a change in Velocity.
4. Velocity causes a change in Position.

```
J = |0, 0, -9.8| * m * fDeltaTime
v += J / m
```


# Collision Between Spheres

Two spherescolide if the distance between their centers is less than the sum of their radiuses.
Becuase we test at discrete time intervals, collisions are detected late, so spheres are always interpenetrating. We need to move them so that they are exactly kissing. But the move must be proportional to their masses.

$$x_\text{cm} = \frac{ \sum_i x_i \cdot m_i}{\sum_i m_i} \quad \text{center of mass between i bodies}$$
$$x_\text{cm} = \frac{x_A \cdot m_A + x_B \cdot m_B}{m_A + m_B}$$

We want to solve for $$x'_\text{cm}$$ such that $$d \equiv x'_B - x'_A = 0$$.

$$x_\text{cm} = \frac{x'_A \cdot m_A + x'_B \cdot m_B }{ m_A + m_B }$$
$$x_A \cdot m_A + x_B \cdot m_B = x'_A \cdot m_A + x'_B \cdot m_B$$
$$\dot{.\hspace{.05in}.} \quad x'_A = x_A + \frac{d \cdot m_B}{m_A + m_B}$$
$$\dot{.\hspace{.05in}.} \quad x'_B = x_B + \frac{d \cdot m_A}{m_A + m_B}$$


# Conservation of Momentum

In order for two colliding object to bounce off eachother we need to generate an impulse for each body at the point of contact, in the direction of the normal of the contact, and with the magnitude considering mass, elasticity, and velocity of the two bodies.

Conservation of Momentum means that the total momentum in a system needs to be the same before and after a collision.
$$x_A \cdot m_A + x_B \cdot m_B = x'_A \cdot m_A + x'_B \cdot m_B$$

Conservation of Energy is similar and those must also be equal.
$$T = \frac{m \cdot v^2}{2} \quad \text{where T is Kenetic Energy}$$
$$\frac{m_A \cdot v_A^2}{2} + \frac{m_B \cdot v_B^2}{2} = \frac{m_A \cdot v_A^{'2}}{2} + \frac{m_B \cdot v_B^{'2}}{2}$$

We want to solve for velocities post collision. Since swe have two linear equation we can solve for the primes after rearanging.

$$\begin{cases}
m_A(v_A - v'_A) = -m_B(v_B - v'_B) \\
m_A(v^2_A - v^{'2}_A) = -m_B(v^2_B - v^{'2}_B)
\end{cases}$$

Remember $$J = m(v' - v)$$.Plugging in our solutions for $$v'_A$$ and $$v'_B$$ and solving for $$J_A$$ and $$J_B$$ gives us:
$$J_A = \frac{2m_Am_B(v_B-v_A)}{m_A-m_B}$$
$$J_B = \frac{2m_Am_B(v_B-v_A)}{m_B-m_A} = -J_A$$

This only solves this system in 1D. To extend to 3D we simply scale the normal fo the contact point to $$J_A$$ and $$J_B$$ respectively.

Btw, we are using inverse mass, so:
$$J_A = \frac{2(v_B-v_A)}{m^{-1}_A-m^{-1}_B}$$

# Elasticity

Elasticity captures how much kenetic energy is lost during a collision. It dampens the resulting impulses.

$$\epsilon = -\frac{v'_B - v'_A}{v_B - v_A} \quad \text{}(2)$$

Remember $$J = m(v'-v)$$
$$\dot{.\hspace{.05in}.} \quad v'_A = v_A + J / m_A \quad \text{and} \quad v'_B = v_B - J/m_B$$

subtracting the two gives:
$$v'_B - v'_A = v_B - v_A - J/m_B - J / m_A$$
$$-\epsilon(v_B-v_A) = v_B - v_A - J/m_B - J / m_A \quad \text{substitute } (v'_B-v'_A) \text{ using (2)}$$
$$\implies J = ( 1 - \epsilon ) \cdot \frac{v_B-v_A}{m^{-1}_A+m^{-1}_B}$$


# Angular Velocity

So far we are not representing rotations. Thankfully there are rotation analogs to velocity, acceleration, and momentum.

| Linear Motion | Rotational Motion |
|----------------|----------------|
| Offset _x_                   | Angular Offset $$\theta$$ |
| Velocity _v_                 | Angular Velosity $$\omega$$ |
| Acceleration _a_              | Angular Acceleration $$\alpha$$ |
| Mass _m_                      | Moment of Inertia _I_ (Inertia Tensor) |
| Force _F_                     | Torque $$\tau$$ |
| Momentum $$\rho = m \cdot v$$ | Angular Momentum $$L = I\omega$$ |
| Impuse $$J = \Delta p$$       | Angular Impulse $$\tau \Delta t = \Delta L$$ |
| Newtons 2nd Law $$F=ma$$      | Rotational $$\tau = I \alpha$$

$$\tau = I \cdot \alpha$$
$$L = I \cdot \omega$$
$$d \theta = \omega \cdot dt$$
$$d \omega = \alpha \cdot dt = I^{-1} \cdot J \quad \text{where J is the angular impuse} J \Delta t$$

In the same way that mass is an objects resistance to being moved, Moment of Inertia (Expressed via the Inertia Tensor) is an objects resistance to being rotated. It captures the mass distribution of the body. Rather than a scalar its a 3x3 matrix.

$$\mathbf{I} = \frac{2}{5} m R^2
\begin{bmatrix}
1 & 0 & 0 \\
0 & 1 & 0 \\
0 & 0 & 1 \\
\end{bmatrix} \quad \text{Inertia Tensfor for a sphere with radius R}$$

Like with velocity, we want to track and upadate the angular velocity of a body.

```
w += I_inv * J
```

Currently, we model collisions as if they pass through the center of mass, inducing NO angula momentum. This is generally not going to happen. If an object is hit off center it should generate both an impulse proportional to mass AND an andulare impulse proportional to moment of inertia.
$$L = I \cdot \omega = \vec{r} \times \vec{p} \quad \text{where } \vec{r} = \text{impulsePoint - centerOfMass}$$
$$dL = I \cdot d\omega = r \times J_\text{linear}$$
$$\implies J_\text{angular}=r \times J$$

Committing angular velocity to orientation is more complicated than commiting velocity to position.
If the body is asymetric the object will precess and cause internal torque on itself (like tossing a tennis raquette).

$$\tau = \omega \times I \cdot \omega \quad \text{and} \quad \tau = I \cdot \alpha$$
$$\dot{.\hspace{.05in}.} \quad \alpha = I^{-1} \cdot (\omega \times I \cdot \omega)$$

# Conservation of Angular Momentum

Lets update our equation for Impulse _J_ to consider conservation of angular momentum.
_See Gregs book for derivation_

$$J = ( 1 + \epsilon ) \cdot \frac{v_B - v_A}{m^{-1}_A + m^{-1}_B + (I^{-1}_A( r_A \times n ) \times r_A + I^{-1}_B( r_B \times n ) \times r_B ) \cdot n }$$


# Frication

Static friction prevents a body from sliding.
Kenetic friction resists motion once started.
In our model we will treat both as one ins a single __Friction__ coefficient.

$$F_\text{friction} = \mu N$$
$$\quad \text{where } \mu \text{is friction coefficient}$$
$$\quad \text{where } N \text{is normal force}$$

Friction force is always going to be in the oposite direction of motion.
$$\quad v_t \quad \text{is tangent velocity}$$
The normal force $$N$$ is just gonna be our contact normal $$n$$.

To apply friction we're gonna apply an impulse as a tangent collision that removes some energy from the system each frame.

$$J = \frac{\mu - v_t}{m^{-1}_A + m^{-1}_B + (I^{-1}_A( r_A \times v_t ) \times r_A + I^{-1}_B( r_B \times v_t ) \times r_B ) \cdot v_t }$$


# Continuouse Collision Detection

Rather than testing collisions in discrete time steps, we can sweep the spheres in the direction of motion at the length alloted by the time ste. This is effectively a collision between two capsuls.

We can simplify by making the motion relative to one of the bodies. This produces a collision between a static sphere and a capsule.

By subtracting the radius of the capsule and adding it to the sphere we can simplify the intersection furthur into a collision between a line and a sphere. Now we can solve for the point and time of impact.

$$d = B - A \quad \text{where A is start of line and B is center of sphere}$$
$$r(t) = A + d(t) \quad \text{line extending from A to distance allowed by time}$$
$$s(t) = B^2 - r^2_B \quad \text{sphere equation}$$

$$0 = r^2_B - s \cdot s$$
$$\quad = r^2_B - B^2 - r^2_B + 2 \cdot B \cdot r_B$$
$$\quad = r^2_B - B^2 - A^2 - A \cdot d \cdot t - d^2 \cdot t^2 \quad \text{this is a quadratic equation}$$

$$t = \frac{ b \pm \sqrt{4ac} }{2a}$$
$$a = -d^2$$
$$b = -A \cdot d$$
$$c = r^2_B - B^2 - A^2$$

This gives us the time of impact (toi). Now we can calculate _toi_ between frames, but we must properly resolve collision.


# Time of Impact

Say that we have three spheres that may collide. The sphere sweep tells use 3 collisions occure. But _toi_ between A and B is way earlier than with C. A and B collide and ricochet and C never collides. But how do we caputre this?

Solution:
1. sort collisions by _toi_
2. resolve contacts in temporal _toi_ order
3. update positions
4. update positions on all bodies from time remaining between last _toi_ and time for entire frame


# Broad Phase

1. Expand bounding boxes of all bodies by a dist they would travel in a single frame.
2. Project bounds onto the |1, 1, 1| diagonal world axis.
3. Sort bodies by position on sorting axis.
4. Create collision pairs from overlapping bodies.
5. Iterate through pairs and resolve collisions.

More info here: https://www.youtube.com/watch?v=MKeWXBgEGxQ