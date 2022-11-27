# 2D hydrodynamic deacaying turbulence

[Go to top](../README.md)  

## Setups and Results
[Go to Notion page.](https://www.notion.so/Turbulence-Studies-e4836ad642684f8f992d54a1f7e22635#97db0fffb85541a891157c14669bd36e)

## How to run

### compile 
To run the code, you need to compile 'Simulation.f90' in GPU server.
    
    make Simulation.x
    
Then `Simulation.x`is made in this directory.

### run
Let's run the code.
    
    sbatch sj_g00.sh
    
The simulation data is saved in `bindata/`.

### analysis
To analyze the data, let us make `Analysis.x` in GPU server..
    
    make Analysis.x
    
Now you have many time-snapshots of data. To count it, use a script.
    
    CountBindata.sh
   
See the file, `cat control.dat`. You can know the number of files.
Then preparation is done. Run the analyis.
    
    sbatch sj_g00_ana.sh
    
The output is saved in `output/`.
### 2D plots and animation.
If you need 2D snapshots. cp `bindata` and output` in GPU server to 'HYD2D' in analyis server.
    
    cp -r ???/bindata .
    cp -r ???/output .
    make 2Dsnaps
   
Using `output/vor*.dat`, image files are made and save as `figures/vor*.png`.
To make movie from the files. Type as follows.

    make movie
   
The movie files in saved in `movie/anivor`.

### spectrum
If you need 2D snapshots. cp `bindata` and output` in GPU server to `HYD2D` in analyis server.
To obtain the spectrum
   
    cp -r ???/bindata .
    cp -r ???/output .
    make spectrum
      
### Do all of them
To do all in one command, you just type `make` or `make all`.
   
      make all
      
