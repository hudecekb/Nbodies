/*----------------------------------------------------------
 * NbodiesParallel
 * by Barrett Hudecek
 * 
 * A parallelized simulation of the n-bodies problem.
 * 
 * Arguments:
 *  - numProcs
 *  - numBodies
 *  - mass
 *  - time
 *  - input filename (optional)
 * 
 ---------------------------------------------------------*/
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Toolkit;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;
import java.util.concurrent.Semaphore;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

public class NbodiesParallel {

  // numWorkers, numBodies, time, max initial velocity for generated bodies
  static int workern;
  static int n;
  static int t;
  private int init = 25;
  static int[] collisions;
  // seed for random generation
  private long seed = 12345678;
  // radius, mass, length of increment, gravity multiplier
  static double r = 10.0;
  static double dt = 1;
  static double m;
  static double gm = 1.0;
  static double G = 6.67e-11;
  static Pnt[] p; // points, velocities, forces
  static Pnt[] v;
  static Pnt[][] f;
  static Pnt[] np, nv; // holds the updated values
  static Semaphore[][] sems;
  static Boolean[][] c; // whether a collision occured last increment
  // whether to display gui, whether to ignore time limit, whether to nudge overlapping bodies
  static Boolean  writingOut = false, graphics = false, infinite = false, nudge = true;
  static double[][] dist;
  static String input, filename = "output.txt";

  static int height = 500;
  static int width = 500;
  static JFrame frame;
  
  /*-----------------------
   *   main
   ----------------------*/
  public static void main(String[] args) {
    
    // create the simulation with initial values
    NbodiesParallel simulation = new NbodiesParallel(args);

    // create gui to show simulation
    SwingUtilities.isEventDispatchThread();
    frame = new JFrame();
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    frame.setBounds(50, 50, width, height);
    if (simulation.getGraphicsMode()) {
      BodiesPanel panel = new BodiesPanel(simulation, width, height);
      frame.add(panel);
      frame.pack();
      frame.setVisible(true);
    }
    
    // read the time
    long startTime = System.currentTimeMillis();
    
    // start threads
    int workernAdj = workern;
    if (writingOut || graphics) {
      workernAdj--;
    }
    int bodiesPer = n / workernAdj;
    int leftovers = n % workernAdj;
    Worker workers[] = new Worker[workern];
    for (int i = 0; i < workern; i++) {
      int myBodies = bodiesPer;
      int startBody = (i*bodiesPer + Math.min(i,leftovers));
      if (i < leftovers) { myBodies++; }
      workers[i] = new Worker(workern, n, i, t, p, v, f, np, nv, dist, c, collisions,
          graphics, writingOut, infinite, nudge, filename, sems, frame, r, G, m, dt,
          myBodies, startBody);
      workers[i].start();
    }
    
    // join threads
    for (int i = 0; i < workern; i++) {
      try {
        workers[i].join();
      } catch (InterruptedException e) {
        System.err.println("Thread " + i + " join failed.");
        e.printStackTrace();
        return;
      }
    }
    
    // print final time
    long endTime = System.currentTimeMillis();
    long startSeconds = startTime / 1000;
    long endSeconds = endTime / 1000;
    long seconds = endSeconds - startSeconds;
    long startMillis = startTime % 1000;
    long endMillis = endTime % 1000;
    long millis = endMillis - startMillis;
    if (endMillis < startMillis) {
      millis += 1000;
      seconds--;
    }
    System.out.println("computation time: "+seconds+" seconds "+millis+" milliseconds");
    int totalC = 0;
    for (int i = 0; i < workern; i++) {
      totalC += collisions[i];
    }
    System.out.println("total collisions: "+totalC+"\n");
    
  }
    
  /*-----------------------
   *   constructor
   ----------------------*/
  public NbodiesParallel(String[] args) {
    
    if (args.length < 4) {
      System.out.println("Required arguments: numProcs numBodies mass time");
    }
    
    // read the command line arguments
    workern = Integer.parseInt(args[0]);
    n = Integer.parseInt(args[1]);
    m = Double.parseDouble(args[2]);
    t = Integer.parseInt(args[3]);
    
    // read optional arguments from input file
    int count = 0;
    if (args.length > 4) {
      try {
        input = args[4];
        BufferedReader buffRead = new BufferedReader(new FileReader(new File(input)));
        String line = buffRead.readLine();
        if (line == null) {
          System.err.println("Failed to read from file. Quitting.");
          buffRead.close();
          return;
        }
        
        // skip first comment block in input file
        while (line.charAt(0) == '-') {
          line = buffRead.readLine();
        }
      
        // radius
        if (!line.equals("default")) {
          r = Double.parseDouble(line);
        }
      
        // deltaT
        line = buffRead.readLine();
        if (!line.equals("default")) {
          dt = Double.parseDouble(line);
        }
        
        line = buffRead.readLine();
        if (!line.equals("default")) {
          gm = Double.parseDouble(line);
        }
        G *= gm;

        line = buffRead.readLine();
        if (!line.equals("default")) {
          init = Integer.parseInt(line);
        }
        
        line = buffRead.readLine();
        if (!line.equals("default")) {
          seed = Long.parseLong(line);
        }
        
        line = buffRead.readLine();
        if (!line.equals("default")) {
          graphics = Boolean.parseBoolean(line);
        }
        
        line = buffRead.readLine();
        if (!line.equals("default")) {
          infinite = Boolean.parseBoolean(line);
        }
        
        line = buffRead.readLine();
        if (!line.equals("default")) {
          nudge = Boolean.parseBoolean(line);
        }
        
        line = buffRead.readLine();
        if (!line.equals("default")) {
          filename = line;
        }
        writingOut = !filename.equals("none");
        
        // skip second comment block in input file
        while (line != null && line.charAt(0) != '(') {
          line = buffRead.readLine();
        }
        
        // initialize arrays
        p = new Pnt[n];
        v = new Pnt[n];
      
        while (line != null && line.charAt(0) == '(') {
          line = line.substring(1, line.length()-1);
          String[] nums = line.split("\\)\\(");
          String[] pos = nums[0].split(",");
          String[] vel = nums[1].split(",");
          p[count] = new Pnt(Double.parseDouble(pos[0]), Double.parseDouble(pos[1]));
          v[count] = new Pnt(Double.parseDouble(vel[0]), Double.parseDouble(vel[1]));
          count++;
          line = buffRead.readLine();
        }
        
        buffRead.close();
      }
      catch (FileNotFoundException e) {
        e.printStackTrace();
        return;
      }
      catch (IOException e) {
        e.printStackTrace();
        return;
      }
    } else { // no optional input file provided
      // initialize arrays
      p = new Pnt[n]; // starting location (initialized as 0, randomly generated below)
      v = new Pnt[n]; // velocities of bodies (initialized to 0)
    }

    f = new Pnt[n][workern]; // forces acting on each body (initialized to 0)
    np = new Pnt[n]; // arrays for threads to work in without stepping on each other
    nv = new Pnt[n];
    c = new Boolean[n][n]; // whether a collision has occured last iteration
    int semDepth = (int) Math.ceil(Math.sqrt(workern));
    sems = new Semaphore[workern][semDepth]; // semaphore array for barrier
    for (int i = 0; i < workern; i++) {
      for (int j = 0; j < semDepth; j++) {
        sems[i][j] = new Semaphore(0, true);
      }
    }
    collisions = new int[workern]; //array for each thread to record its collisions separately
    for (int i = 0; i < n; i++) {
      if (i < workern) {
        collisions[i] = 0;
      }
      for (int j = 0; j < workern; j++) {
        f[i][j] = new Pnt();
      }
      for (int j = 0; j < n; j++) {
        c[i][j] = false;
      }
    }
    dist = new double[n][n]; // distances between bodies (calculated later by calcForce())
    
    // randomly generate any unprovided initial positions and velocities
    // if two bodies are overlapping, re-generate
    // if the frame gets crowded with bodies, this will get very inefficient
    if (count < n) {
      Random rand = new Random(seed);
      for (int i = count; i < n; i++) {
        Boolean done = false;
        while (!done) {
          p[i] = new Pnt(rand.nextInt(width-(4*(int)Math.ceil(r)))+r,
                  rand.nextInt(height-(4*(int)Math.ceil(r)))+r);
          v[i] = new Pnt(rand.nextInt(init*2)-init,
              rand.nextInt(init*2)-init);
          np[i] = p[i];
          nv[i] = v[i];
          done = true;
          for (int j = 0; j < i; j++) {
            double dist = Math.sqrt( Math.pow((p[i].x-p[j].x), 2) 
                                   + Math.pow((p[i].y-p[j].y), 2) );
            if (dist < 2*r) {
              done = false;
              break;
            }
          } 
        }
      }
    }
    
    for (int i = 0; i < n; i++) {
      np[i] = new Pnt(p[i].x, p[i].y);
      nv[i] = new Pnt(v[i].x, v[i].y);
    }
    
  } // constructor
  
  // many get methods
  public Pnt[] getPoints() {
    return p;
  }
  
  public double getRadius() {
    return r;
  }
  
  public Boolean getGraphicsMode() {
    return graphics;
  }
  
  public Boolean getInfiniteMode() {
    return infinite;
  }
  
  public Boolean getNudgeMode() {
    return nudge;
  }
  
} // NbodiesParallel

/*----------------------------------------------------
 *  class Worker
 *  
 *  a Thread responsible for some number of bodies
 ----------------------------------------------------*/
class Worker extends Thread {
  
  // how many bodies this thread is responsible for,
  //   and which one is its first
  int myBodies, startBody;
  // numWorkers, numBodies, iteration, workerID
  int workern, n, iter, id, t;
  int[] collisions;
  double r, G, m, dt;
  Semaphore[][] sems;
  Pnt[] p, v; // points, velocities
  Pnt[][] f; // forces
  Pnt[] np, nv; // holds the updated values
  double[][] dist;
  Boolean[][] c;
  String filename;
  Boolean writingOut, graphics, infinite, nudge;
  OutputWriter writer;
  JFrame frame;
  
  // constructor
  // assign a whole boatload of local references
  // i feel like this is horribly messy, but i really struggled finding a better way to get
  //    all of this information to the threads. i could've passed the NbodiesParallel object,
  //    but then there would just be this many get methods :/
  public Worker(int workern, int n, int id, int t,  Pnt[] p, Pnt[] v, Pnt[][] f, Pnt[] np,
      Pnt[] nv, double[][] dist, Boolean[][] c, int[] collisions, Boolean graphics,
      Boolean writingOut, Boolean infinite, Boolean nudge, String filename, Semaphore[][] sems,
      JFrame frame, double r, double G, double m, double dt, int myBodies, int startBody) {
    this.workern = workern;
    this.n= n;
    iter = 0;
    this.id = id;
    this.t= t;
    this.p = p;
    this.v = v;
    this.f= f;
    this.np = np;
    this.nv = nv;
    this.dist = dist;
    this.c = c;
    this.collisions = collisions;
    this.filename = filename;
    this.graphics = graphics;
    this.infinite = infinite;
    this.writingOut = writingOut;
    this.nudge = nudge;
    this.sems = sems;
    this.frame = frame;
    this.r = r;
    this.G = G;
    this.m = m;
    this.dt = dt;
    this.myBodies = myBodies;
    this.startBody = startBody;
  }
  
  // where the thread actually starts
  @Override
  public void run() {
    
    Boolean time = iter < t;
    while (time) {
      
      // last worker handles the gui and writing to output file
      if (id == workern-1 && (writingOut || graphics)) {
        if (iter == 0) {
          myBodies = 0;
          startBody = -1;
          if (writingOut) {
            writer = new OutputWriter(filename);
            writer.write(0);
          }
        } else {
          if (writingOut || graphics) {
            writeAndDraw();
          }
        }
        barrier(); // just to let the other threads move on
      }
     
      // other workers handle the actual bodies
      if (id != workern-1 || (!writingOut && !graphics)) {
        // clear out this thread's part of f
        for (int i = startBody; i < startBody+myBodies; i++) {
          f[i][0] = new Pnt();
        }
        calcForce();
        barrier(); // all forces need to be in f before we can sum them
        // sum all forces on our bodies recorded by other threads
        // (also clear out those parts of f while we're at it)
        for (int i = startBody; i < startBody+myBodies; i++) {
          for (int j = 1; j < workern; j++) {
            f[i][0].x += f[i][j].x;
            f[i][0].y += f[i][j].y;
            f[i][j] = new Pnt();;
          }
        }
        moveBodies();
        collide();
      }

      // iterate, but don't actually check the value if infinite mode is on
      iter++;
      if (!infinite) time = iter < t;
      
      // wait for everyone to finish
      barrier();

      // copy the point and velocity values back from the updated arrays
      for (int i = startBody; i < startBody+myBodies; i++) {
        p[i] = new Pnt(np[i].x, np[i].y);
        v[i] = new Pnt(nv[i].x, nv[i].y);
      }
      // wait again
      barrier();
      
    }
    
    // close output file
    if (id == workern-1 && writingOut) writer.close();
  }
  
  /*------------------------------------------------------
   *   writeAndDraw
   *   
   *   manages the gui and file output
   -----------------------------------------------------*/
  public void writeAndDraw() {

    if (graphics) {
      frame.repaint();
    }
    if (writingOut) {
      writer.write(iter);
    }
    
  }

  /*--------------------------------------------------
   *   calcForce
   *   
   *   calculates forces between this thread's bodies,
   *   and all of the following bodies
   *   
   *   forces calculated for other threads are left in
   *   f[bodyNum][id] so that we don't run into any
   *   other threads in f
   --------------------------------------------------*/
  public void calcForce() {
    
    double mag; // magnitude
    Pnt dir;  // direction
    for (int i = startBody; i < (startBody+myBodies); i++) {
      for (int j = i+1; j < n; j++) {
        dist[i][j] = Math.sqrt( Math.pow((p[i].x-p[j].x), 2) 
                              + Math.pow((p[i].y-p[j].y), 2) );
        dist[j][i] = dist[i][j];
        mag = ( G * Math.pow(m, 2) / Math.pow(dist[i][j], 2) );
        dir = new Pnt(p[j].x-p[i].x, p[j].y-p[i].y);
        f[i][0].x += mag*dir.x/dist[i][j];
        f[j][id].x -= mag*dir.x/dist[i][j];
        f[i][0].y += mag*dir.y/dist[i][j];
        f[j][id].y -= mag*dir.y/dist[i][j];
      }
    }
    
  }

  /*------------------------------------------------------
   *   moveBodies
   *   
   *   calculates new velocity and position for this
   *   thread's bodies
   ------------------------------------------------------*/
  public void moveBodies() {
    
    Pnt deltav = new Pnt(), deltap = new Pnt();
    for (int i = startBody; i < (startBody+myBodies); i++) {
      deltav.set( (f[i][0].x/m*dt), (f[i][0].y/m*dt) );
      deltap.set( (v[i].x + deltav.x/2) * dt,
                  (v[i].y + deltav.y/2) * dt);
      nv[i].x += deltav.x;
      nv[i].y += deltav.y;
      np[i].x += deltap.x;
      np[i].y += deltap.y;

      // detect collisions with walls of frame
      if (np[i].x <= r) {
        np[i].x = r;
        nv[i].x *= -1;
      }
      if (np[i].x >= (frame.getWidth()-r)) {
        np[i].x = frame.getWidth()-r;
        nv[i].x *= -1;
      }
      if (np[i].y <= r) {
        np[i].y = r;
        nv[i].y *= -1;
      }
      if (np[i].y >= (frame.getHeight()-r)) {
        np[i].y = frame.getHeight()-r;
        nv[i].y *= -1;
      }
    }
  
  }
  
  /*------------------------------------------------------
   *   collide
   *   
   *   checks for collisions and calculates new velocities
   ------------------------------------------------------*/
  public void collide() {
    
    for (int i = startBody; i < (startBody+myBodies); i++) {
      for (int j = i+1; j < n; j++) {
        if (dist[i][j] <= 2*r) {
          if (!c[i][j]) {
            c[i][j] = true;
            collisions[id]++;
          }
          if (nudge) {
            // nudge overlapping bodies away from each other
            double overlap = (2*r)/dist[i][j];
            if (overlap > 1.2) {
              double xnudge = (np[j].x-np[i].x)*((overlap-1)/2);
              double ynudge = (np[j].y-np[i].y)*((overlap-1)/2);
              np[i].x += -1*xnudge;
              np[i].y += -1*ynudge;
              np[j].x += xnudge;
              np[j].y += ynudge;
            }
          }
          // recalculate velocities for each body
          // this is a nightmare to read, but it seems to work
          double dx = Math.pow(p[i].x - np[j].x, 2);
          double dy = Math.pow(p[i].y - np[j].y, 2);
          double dx2 = Math.pow(dx, 2);
          double dy2 = Math.pow(dy, 2);
          // equation A
          nv[i].x = ( (nv[j].x*dx2 + nv[j].y*dx*dy +
                      nv[i].x*dy2 - nv[i].y*dx*dy)
                      / (dx2 + dy2) );
          // equation B
          nv[i].y = ( (nv[j].x*dx*dy + nv[j].y*dy2 -
                      nv[i].x*dy*dx + nv[i].y*dx2)
                      / (dx2 + dy2) );
          // equation C
          nv[j].x = ( (nv[i].x*dx2 + nv[i].y*dx*dy +
                      nv[j].x*dy2 - nv[j].y*dx*dy)
                      / (dx2 + dy2) );
          // equation D
          nv[j].y = ( (nv[i].x*dx*dy + nv[i].y*dy2 -
                      nv[j].x*dy*dx + nv[j].y*dx2)
                      / (dx2 + dy2) );
        } else {
          c[i][j] = false;
        }
      }
    }   
  } // collide
  
  /*--------------------------------------------------------------
   *   barrier
   *   
   *   dissemination barrier to keep threads on the same iteration
   --------------------------------------------------------------*/
  private void barrier() {

    int round = 0, target_distance = 1;
    while (target_distance < workern) {
      // V on our semaphore
      sems[id][round].release();
      // P on target
      int target = (target_distance + id) % workern;
      try {
        sems[target][round].acquire();
      } catch (InterruptedException e) {
        e.printStackTrace();
      }
      target_distance *= 2;
      round++;
    }

  } // barrier()
  
  /*----------------------------------------------------
   *  class OutputWriter
   *  
   *  writes to an output file
   ----------------------------------------------------*/
  private class OutputWriter {
    
    BufferedWriter writer;
    
    public OutputWriter(String filename) {
      // don't open the file if we aren't going to use it
      if (filename.equals("none")) return;
      
      // open the file
      try {
        writer = new BufferedWriter(new FileWriter(filename));
      } catch (IOException e) {
        e.printStackTrace();
        return;
      }
    }
    
    // writes the current details of the simulation to the output file
    public void write(int increment) {
      try {
        writer.write("------ Increment " + (increment) + " ------\n");
        for (int i = 0; i < n; i++) {
          writer.write("Body " + (i+1) + ":");
          writer.write("\n\tPosition - (" + String.format("%.2f",p[i].x) + ","
              + String.format("%.2f",p[i].y) + ")\n");
          writer.write("\tVelocity - (" + String.format("%.2f",v[i].x) + "," 
              + String.format("%.2f",v[i].y) + ")\n");
          writer.write("\n");
        }
      } catch (IOException e) {
        e.printStackTrace();
        System.err.println("Increment " + increment + " failed to write to output file.");
      }
    }
    
    // closes the file stream
    public void close() {
      try {
        writer.close();
      } catch (IOException e) {
        e.printStackTrace();
      }
    }
    
  } // OutputWriter
  
} // class Worker

/*----------------------------------------------------
 *  class BodiesPanel
 *  
 *  JPanel to animate the bodies on
 ----------------------------------------------------*/
class BodiesPanel extends JPanel {

  static final long serialVersionUID = 1L;
  private int height, width;
  private double r;
  private Pnt[] points;
  private Color[] colors;
  
  public BodiesPanel(NbodiesParallel sim, int width, int height) {
    points = sim.getPoints();
    r = sim.getRadius();
    this.width = width;
    this.height = height;
    
    // generate random colors for the bodies
    Random rand = new Random();
    colors = new Color[points.length];
    for (int i = 0; i < colors.length; i++) {
      float r = rand.nextFloat();
      float g = rand.nextFloat();
      float b = rand.nextFloat();
      colors[i] = new Color(r, g, b);
    }
  }
  
  public Dimension getPreferredSize() {
    return new Dimension(width, height);
  }
  
  public void paintComponent(Graphics g) {
    super.paintComponent(g);
 
    // draw all the bodies
    for (int i = 0; i < points.length; i++) {
      int x = (int) (points[i].x-r);
      int y = (int) (points[i].y-r);
      g.setColor(colors[i]);
      g.fillOval(x, y, (int) (2*r), (int) (2*r));
      g.setColor(Color.BLACK);
      g.drawOval(x, y, (int) (2*r), (int) (2*r));
    }
    Toolkit.getDefaultToolkit().sync();
  }

}

/*----------------------------------------------------
 *  class Pnt
 *  
 *  stores a Pnt as two public doubles
 ----------------------------------------------------*/
class Pnt {
  
  public double x;
  public double y;
  
  public Pnt() {
    x = 0;
    y = 0;
  }
  
  public Pnt(double x, double y) {
    this.x = x;
    this.y = y;
  }
  
  public void set(double x, double y) {
    this.x = x;
    this.y = y;
  }
} // Pnt
