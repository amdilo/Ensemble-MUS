# Ensemble-MUS

Ensemble generation for Made-Up Sensor from a flat truth image and a simple error structure: relative error at pixel level, e.g. random error of Earth counts; relative error at scanline level, e.g. combination of random errors of PRTs' temperature and ICT counts. 

The script gens_v0.py generates FCDR ensemble of a Made-Up Sensor with random errors at pixel and scanline level for which the relative standard uncertainty is given as a constant. 

The other scripts aim the same for a more complex error structure and are organised in modules with classes and/or functions offering a similar functionality. 
