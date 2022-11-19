# 2D hydrodynamic deacaying turbulence

[Go to top](../README.md)  

## How to run
To run the code, you just type `make`.
    
    make
    
Then `Simulation.x`is compiled and automatically executed.
The simulation data is saved in `bindata/`.

The data is binary file, to make figures analyis is done by `Analysis.x`.

## Initial condition
stream function is given.

$$ \psi \propto \sin(2\pi k x)\sin(2\pi k y).$$

Initial velocity is set from $\psi$.

$$ v_x = \partial_y \psi, v_y = -\partial_x \psi, $$

## Analysis

Animation of vorticity will be made.

$$ \omega =  \partial_x v_y - \partial_y v_x $$

The spectrum of variable, $X$ ,is calculated as follows. First 2D Fourier transformation is given.

$$\hat{X}_{{\rm 2D},c}(k_x,k_y) = \iint dx dy X \cos(2\pi (k_x x+k_y y))$$

$$\hat{X}_{{\rm 2D},s}(k_x,k_y) = \iint dx dy X\sin(2\pi (k_x x+k_y y))$$

Then 1D Fourier transformation is given.

$$\hat{X}_{\rm 1D} dk = \sqrt{ \hat{X}_{{\rm 2D},c}^2+\hat{X}_{{\rm 2D},s}^2} dk_x dk_y $$
