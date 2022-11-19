# 2D deacaying turbulence

## initial condition
stream function is given.

$$ \psi \propto \sin(2\pi k x)\sin(2\pi k y).$$

Initial velocity is set from $\psi$.

$$ v_x = \partial_y \psi, v_y = -\partial_x \psi, $$

## analysis

Animation of vorticity will be made.

$$ \omega =  \partial_x v_y - \partial_y v_x $$

The spectrum of kinetic energy and enstropy is calculated. First 2D Fourier transformation is given.

$$\hat{E}_{{\rm 2D},c}(k_x,k_y) = \iint dx dy \rho v^2 \cos(2\pi (k_x x+k_y y))$$

$$\hat{E}_{{\rm 2D},s}(k_x,k_y) = \iint dx dy \rho v^2 \sin(2\pi (k_x x+k_y y))$$

$$\hat{V}_{{\rm 2D},c}(k_x,k_y) = \iint dx dy \omega^2 \cos(2\pi (k_x x+k_y y))$$

$$\hat{V}_{{\rm 2D},s}(k_x,k_y) = \iint dx dy \omega^2 \sin(2\pi (k_x x+k_y y))$$

Then 1D Fourier transformation is given.

$$\hat{E}_{\rm 1D} dk = \sqrt{ \hat{E}_{{\rm 2D},c}^2+\hat{E}_{{\rm 2D},s}^2} dk_x dk_y $$
