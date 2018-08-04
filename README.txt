Barrett Hudecek
CSC 422
Parallel Project


The task for this project was to implement both a sequential and parallelized version of the n-body problem, which simulates the interactions of particles in a closed space. The purpose was to better understand the way parallelization affects (and doesn’t affect) performance.


I wrote my programs in Java. The sequential and parallelized programs are (respectively) Nbodies.java and NbodiesParallel.java. Both programs simulate the problem in 2 dimensions within a 500 x 500 square. Bodies which approach the edges of the simulation area are bounced back. Both programs have the same required and optional arguments.


Required Arguments:
- numWorkers (int): The number of threads to employ for solving the problem. For the sake of consistency, both versions take this argument, but the procedural version ignores it completely.
- numBodies (int): The number of bodies to simulate. These will be generated pseudo-randomly unless detailed information is provided about them in the optional input file.
- mass (double): The mass of each body for the sake of calculating its gravitational effects on other bodies. For the sake of simplicity, all bodies have the same mass.
- time (int): The number of discrete steps to take in the simulation. If infinite mode is enabled via optional arguments, this will be ignored.

Optional Arguments:
- input filename (String): This is used to open a text file containing any of the below options to change the simulation or its output. Since the input file format is rather specific, several example input files are provided.

Input File Options:
- radius (double): The radius of each body. For the sake of simplicity, all bodies have the same radius.
- delta time (double): A multiplier for the length of each discrete increment. Values larger than 1 will cause the simulation to take larger steps and thus move faster, but also tend to result in strange bugs (such as bodies overlapping and escaping from the frame).
- gravity multiplier (double): The gravitational constant is multiplied by this number before the simulation begins. Values higher than 1 will thus cause bodies to be more attracted to each other.
- maximum initial velocity (int): Bodies which are randomly generated will not start with a velocity higher than this along either single axis (x or y). Higher values tend to result in faster, more chaotic simulations.
- seed (long): The seed which is used to generate initial position and velocities of bodies, as well as colors for them to be animated with. For the sake of consistent testing, the default is a generic value (12345678).
- graphics mode (boolean): False by default. If set to true, the simulation will be visible as an animation.
- infinity mode (boolean): False by default. Causes the limit on increments to be ignored (ie the program will continue until manually shut down). Mainly added because the animation can be fun to watch.
- nudge mode (boolean): True by default. Causes bodies which are overlapping to be artificially pushed away from each other. This prevents some overlapping, but in cramped simulations tends to result in jittery behavior by bodies.
- output filename (String): By default, detailed information about each increment of the simulation is printed to “output.txt”, but a preferred file can be provided here. If no output file is desired, the value “none” can be used.

Known Bugs:
- When using the command “make all” to compile both programs, subsequent runs of Nbodies will sometimes return a “NoSuchMethodError”. I believe this has something to do with both programs having identical copies of several sub-classes, but I am not very experienced with Java Makefiles. Compiling just one program or the other at a time resolves the issue however.
- Extreme values, particularly of delta time and the gravity multiplier, will sometimes result in all bodies becoming trapped on top of each other in one corner of the simulation. I have not been able to figure why this happens.


Most of my verifications of correctness were done via visual inspections of the animated mode. While bugs undoubtedly remain, I am reasonably sure that gravity and collisions between bodies will work correctly in most situations. The following runs were helpful for confirming this. Though I am writing them with Nbodies, it should work the same when NbodiesParallel is substituted.

java Nbodies 4 4 100000 100 4stillBodies.txt

This run has gravity and the mass of bodies turned up quite high. 4 equidistant bodies with no starting velocity are eventually pulled together by gravitational forces. Afterwards, the gravity is so high that they bounce and curve quite chaotically.

java Nbodies 5 5 100 100 infBig.txt

	A few random large bodies start off bouncing quickly, but eventually lose speed and begin to clump together.

java Nbodies 50 50 100 100 infLittle.txt

	Similar to the previous case, but with a larger number of small bodies.

java Nbodies 1 15 50 1 infinite.txt

  This one is just fun to watch. :)
