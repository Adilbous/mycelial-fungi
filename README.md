# Mycelial Fungi 

This project aims to model the growth of mycelium and the active and inacative propagation of glucosis into the mycelial network. 
Active propagation of nutrients in the mycelium, also called translocation, translates the ability of fungi to move nutrients from an area to another to sustain the growth of the mycilum in nutrients-poor area (or polluted areas, for example). 
Hence it is a remarkable property of fungi, making them able to grow in areas where other organism (plants) can't.
The equations and numerical resolution approach are similar to those exposed by BOSWELL ET AL. (2002, 2003). 


## Content 

- function.R 
  - Biological parameters of the model 
  - Utils function necessary to resolve the model

- solver.R
  - Numerical resolution of the equations

- plot.R 
  - Model's "front-end"
  - Allows the user to solve the equations and to plot a graphical representation of the main system variables
  
## To use this project 

- Open solver.R, edit 'setwd()' and indicate the path to 'function.R' script
- Open plot.R, edit 'setwd()' and indicate the path to 'solver.R' script

- Execute plot.R : U = solver() resolves the equations and saves the solution matrix U

- Use plot_mmp(U), plot_p(U), plot_si(U) et plot_se(U) to display graphical representation of system variables
