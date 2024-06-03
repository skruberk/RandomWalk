This code simulates a random walk in 2D with a circular boundary on a grid to model the way an actin filament tip interacts with polymerases on a membrane surface. Determines whether the walker has hit a target by checking distance vs circle perimeter and then sending the walker back to a random point after a hit. Circular boundary wraps around to the other side and in earlier versions was reflective. The targets also move and their walk is a diffusing point tied by an elastic spring to the center, so the degree of their movement can be controlled with a spring constant. Any code labeled "runs" (RandomWalkRuns.m) allows the code to run multiple times and generates a histogram while other code (RandomWalk.m) generates a plot of the random walk and shows the targets moving over time. For plotting the random walk, each run is displayed as a different color. The data can be fit with an exponential (FitExponential.m) and then the output data can be cleaned and plotted in R (CleanMatlabData.R)
More descriptions can be found in the code comments.

Variables and Parameters to Set  \
radius     : filament movement radius variable  \
num_discs   : Number of discs or binding points  \
target   : Membrane disc variable size  \
tstep  : number of time steps in random walk \ 
num_runs : how many times to run code  \
spring_target : controls how much target is tied to the center  \  
tau : time increment between steps units of seconds  \
D : diffusion coefficient for the walk units of um^2/sec  \
delta : distance increment (in units of um)  \
