# 2D hydrodynamic turbulence

[Go to top](./README.md)  

## initial condition
stream function is given.

$$ \psi \propto \sin(2\pi k x)\sin(2\pi k y).$$

Initial velocity is set from $\psi$.

$$ v_x = \partial_y \psi, v_y = -\partial_x \psi, $$


## Results

https://user-images.githubusercontent.com/20675833/202837591-93161c38-5829-420b-810d-da3190799c21.mp4

Animation of vorticity is shown.

$$ \omega =  \partial_x v_y - \partial_y v_x $$

The spectrum of kinetic energy and enstropy is calculated. First 2D Fourier transformation is given.

$$\hat{E}_{{\rm 2D},c}(k_x,k_y) = \iint dx dy \rho v^2 \cos(2\pi (k_x x+k_y y))$$

$$\hat{E}_{{\rm 2D},s}(k_x,k_y) = \iint dx dy \rho v^2 \sin(2\pi (k_x x+k_y y))$$

$$\hat{V}_{{\rm 2D},c}(k_x,k_y) = \iint dx dy \omega^2 \cos(2\pi (k_x x+k_y y))$$

$$\hat{V}_{{\rm 2D},s}(k_x,k_y) = \iint dx dy \omega^2 \sin(2\pi (k_x x+k_y y))$$

Then 1D Fourier transformation is given.

$$\hat{E}_{\rm 1D} dk = \sqrt{ \hat{E}_{{\rm 2D},c}^2+\hat{E}_{{\rm 2D},s}^2} dk_x dk_y $$


![k-E_k](https://user-images.githubusercontent.com/20675833/202843901-1f5c51f6-7bf0-43d6-8f4e-30ea8ae5d9ca.png)
![k-V_k](https://user-images.githubusercontent.com/20675833/202843904-e5db78cc-5ffb-4a68-af75-1f97bd21b2dd.png)




