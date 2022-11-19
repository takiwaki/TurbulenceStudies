# 2D hydrodynamic turbulence

[Go to summary](./README.md)  

## initial condition
stream function is given.

$$ \psi \propto \sin(2\pi k x)\sin(2\pi k y).$$

Initial velocity is set from $\psi$.

$$ v_x = \partial_y \psi, v_y = -\partial_x \psi, $$


## Results

https://user-images.githubusercontent.com/20675833/202837591-93161c38-5829-420b-810d-da3190799c21.mp4

Animation of vorticity is shown.

$$ \omega =  \partial_x v_y - \partial_y v_x $$

The spectrum of 2D variables are calculated as follows.

$$\hat{X}_{{\rm 2D},c}(k_x,k_y) = \iint dx dy X \cos(2\pi (k_x x+k_y y))$$

$$\hat{X}_{{\rm 2D},s}(k_x,k_y) = \iint dx dy X \sin(2\pi (k_x x+k_y y))$$

From it, then 1D Fourier transformation is given.

$$\hat{X}_{\rm 1D} dk = \sqrt{ \hat{X}_{{\rm 2D},c}^2+\hat{X}_{{\rm 2D},s}^2} dk_x dk_y $$



![k-E_k](https://user-images.githubusercontent.com/20675833/202852849-9acb138b-085e-490a-b166-75b8d98fb3d4.png)

![k-V_k](https://user-images.githubusercontent.com/20675833/202852852-4ae2517c-b293-46d9-bdaf-b6d87c0ac2e8.png)

