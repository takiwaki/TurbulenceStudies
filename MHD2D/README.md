# 2D magneto-hydrodynamic deacaying turbulence

[Go to top](../README.md)  

## How to run
To run the code, you just type `make`.
    
    make
    
Then `Simulation.x`is compiled and automatically executed.
The simulation data is saved in `bindata/`.

The data is binary file, to make figures analyis is done by `Analysis.x`.

## initial condition
stream function is given.

$$ \psi \propto \sin(2\pi k x)\sin(2\pi k y).$$

Initial velocity and magnetic fields are set from $\psi$.

$$ v_x = v_0 \partial_y \psi, v_y = - v_0 \partial_x \psi, $$

$$ B_x = B_0 \partial_y \psi, B_y = - B_0 \partial_x \psi, $$

## analysis

Animation of vorticity and current density are made.

$$ \omega =  \partial_x v_y - \partial_y v_x $$

$$ j =  \partial_x B_y - \partial_y B_x $$

The spectrum of variable, $X$ ,is calculated as follows. First 2D Fourier transformation is given.

$$\hat{X}_{{\rm 2D},c}(k_x,k_y) = \iint dx dy X \cos(2\pi (k_x x+k_y y))$$

$$\hat{X}_{{\rm 2D},s}(k_x,k_y) = \iint dx dy X\sin(2\pi (k_x x+k_y y))$$

Then 1D Fourier transformation is given.

$$\hat{X}_{\rm 1D} dk = \sqrt{ \hat{X}_{{\rm 2D},c}^2+\hat{X}_{{\rm 2D},s}^2} dk_x dk_y $$



