# Parrameterizing the 4D Hyperboloid

A hyperboloid is a generalization of a hyperbola in two dimensions. It’s basically a hyperbola spun around its principal axis forming two 2D surfaces embedded in 3D space. The first thing we must do is map Euclidean space onto one of the sheets of the hyperboloid. This is done by parametrizing the hyperboloid. By doing this, we are saying that for every point $(u,v)$ on the 2D Euclidean plane, there is a point $(x,y,z)$ on the surface of the hyperboloid. A simple google search will get you the parametrization you need to do this. At this point, you may think you’re off the hook, and everything is fine and dandy. But wait! All we can do is map a 2D point on the plane to a 3D point on the surface of the hyperboloid. That’s not very exciting. We want to map fully 3-dimensional objects to this hyperboloid, damn it!

Fortunately, this can be achieved by raising the dimension of the hyperboloid. Sadly, this means we’re going to have to do the parametrization by hand. Google can't save you from math this time!
To parameterize, we will be heavily relying on the following hyperbolic trig identity:

$${cosh^2(x) - sinh^2(x) = 1}$$

Since we want to map 3-dimensional Euclidean space, we must use a 4-dimensional hyperboloid.

$${a^2 + b^2 + c^2 - d^2 = -1}$$
$${t^2 - c^2 = 1 \qquad t^2 = d^2 - a^2 - b^2}$$
$${t^2 - c^2 = 1 \qquad s^2 - a^2 = t^2 \qquad s^2 = d^2 - b^2}$$

Solve for c:

$${t^2 - c^2 = 1}$$
$${c = sinh(v) \qquad t = cosh(v)}$$

Solve for a:

$${s^2 - a^2 = t^2}$$
$${(s/t)^2 - (a/t)^2 = 1}$$
$${s/t = cosh(u) \qquad a/t = sinh(u)}$$
$${a = tsinh(u) \qquad s = tcosh(u)}$$

Solve for b and d:

$${s^2 = d^2 - b^2}$$
$${(d/s)^2 - (b/s)^2 = 1}$$
$${d/s = cosh(w) \qquad b/s = sinh(w)}$$
$${d = scosh(w) \qquad b = ssinh(w)}$$

Solve fo a, b, c, and d and get rid of s and t:

$${c = sinh(v)}$$
$${a = tsinh(u) = cosh(v)sinh(u)}$$
$${b = ssinh(w) = tcosh(u)sinh(w) = cosh(v)cosh(u)sinh(w)}$$
$${d = scosh(w) = tcosh(u)cosh(w) = cosh(v)cosh(u)cosh(w)}$$

And that’s it! Now we can easily map a set of 3D points of to the surface of hyperboloid. Don’t worry too much about the fourth term, for the purpose of mapping the vertices of a 3d model to the surface of the hyperboloid you can drop it entirely.
