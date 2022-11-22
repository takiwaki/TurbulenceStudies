# 2D hydrodynamic deacaying turbulence

[Go to top](../README.md)  

## Setups and Results
[Go to Notion page.](https://www.notion.so/Turbulence-Studies-e4836ad642684f8f992d54a1f7e22635#0c0a5a621c4f42cca05c08ec140ee73f)

## How to run

### compile 
To run the code, you need to compile 'Simulation.f90'.
    
    make Simulation.x
    
Then `Simulation.x`is made in this directory.

### run
Let's run the code.
    
    ./Simulation.x
    
The simulation data is saved in `bindata/`.

### analysis
To analyze the data, let us make `Analysis.x`.
    
    make Analysis.x
    
Now you have many time-snapshots of data. To count it, use a script.
    
    CountBindata.sh
   
See the file, `cat control.dat`. You can know the number of files.
Then preparation is done. Run the analyis.
    
    ./Analyis.x
    
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
      
