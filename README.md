This code simulates a random walk in 2D with a circular boundary on a grid to model the way an actin filament tip interacts with polymerases on a membrane surface. Determines whether the walker has hit a target by checking distance vs circle perimeter and then sending the walker back to a random point after a hit. Circular boundary wraps around to the other side and in earlier versions was reflective. The targets also move and  their walk is a diffusing point tied by an elastic spring to the center, so the degree of their movement can be controlled with a spring constant. Any code labeled "runs" allows the code to run multiple times and generates a histogram while other code generates a plot of the random walk and shows the targets moving over time. For plotting the random walk, each run is displayed as a different color.Elastic spring modeled in files named "elasticspring" for a plot of the randomwalk and "elasticspringruns" to generate the histogram. 
More descriptions can be found in the code comments.

Variables and Parameters to Set
radius     : filament movement radius variable \n
num_discs   : Number of discs or binding points \n
target   : Membrane disc variable size  \n
tstep  : number of time steps in random walk \n
num_runs : how many times to run code
spring_target : controls how much target is tied to the center  
tau : time increment between steps units of seconds 
D : diffusion coefficient for the walk units of um^2/sec
delta : distance increment (in units of um)
