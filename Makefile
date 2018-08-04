JFLAGS = -g
JC = javac

all: Nbodies NbodiesParallel

Nbodies: Nbodies.class

NbodiesParallel: NbodiesParallel.class

Nbodies.class: Nbodies.java
	$(JC) $(JFLAGS) Nbodies.java

NbodiesParallel.class: NbodiesParallel.java
	$(JC) $(JFLAGS) NbodiesParallel.java

clean:
	$(RM) *.class
