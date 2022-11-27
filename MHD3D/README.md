# 3D magneto-hydrodynamic deacaying turbulence
[Go to top](../README.md)  

## Setups and Results
[Go to Notion page.](https://www.notion.so/Turbulence-Studies-e4836ad642684f8f992d54a1f7e22635#9f4fd3d31ffc4343aaa5b0c82e430090)

## How to run

### compile 
To run the code, you need to compile 'Simulation.f90'.
    
    make Simulation.x
    
Then `Simulation.x`is made in this directory.

### run
Let's run the code.
    
    ./Simulation.x
    
The simulation data is saved in `bindata/`.

### visualization
To visulaize the data, let us make `Visulatization.x`.
    
    make Visulatization.x
    
Now you have many time-snapshots of data. To count it, use a script.
    
    CountBindata.sh
   
See the file, `cat control.dat`. You can know the number of files.
Then preparation is done. Run the analyis.
    
    ./Visulatization.x
    
The output is saved in `output/`.
### 2D plots and animation.
If you need 2D snapshots. 
    
    make 2Dsnaps
   
### 3D plots
The output is saved in `vstdata/`. Using VisIt, you can make 3Dfigures.

### analysis
To analyze the data, let us make `Analysis.x`.
    
    make Analysis.x
    
Now you have many time-snapshots of data. To count it, use a script.
    
    CountBindata.sh
   
See the file, `cat control.dat`. You can know the number of files.
Then preparation is done. Run the analyis.
    
    ./Analyis.x
    
The output is saved in `output/`.

### spectrum
To obtain the spectrum
   
      make spectrum
      
## Do all of them
To do all in one command, you just type `make` or `make all`.
   
      make all
