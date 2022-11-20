# 2D hydrodynamic deacaying turbulence

[Go to top](../README.md)  

## How to run

### compile 
To run the code, you need to compile 'Simulation.f90'.
    
    make simulation-code
    
Then `Simulation.x`is made in this directory.

### run
Let's run the code.
    
    make run
    
The simulation data is saved in `bindata/`.

### analysis
To analyze the data, let us make `Analysis.x`.
    
    make analysis-code
    
Now you have many time-snapshots of data. To count it, use a script.
    
    make count-number
   
See the file, `cat control.dat`. You can know the number of files.
Then preparation is done. Run the analyis.
    
    make run-analysis
    
The output is saved in `output/`.
### 2D plots and animation.
If you need 2D snapshots. 
    
    make 2Dsnaps
   
Using `output/vor*.dat`, image files are made and save as `figures/vor*.png`.
To make movie from the files. Type as follows.

    make movie
   
The movie files in saved in `movie/anivor`.

### spectrum
To obtain the spectrum
   
      make spectrum
      
### Do all of them
To do all in one command, you just type `make` or `make all`.
   
      make all
      
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
