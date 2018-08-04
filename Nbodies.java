/*----------------------------------------------------------
 * Nbodies
 * by Barrett Hudecek
 * 
 * A procedural simulation of the n-bodies problem.
 * 
 * Arguments:
 *  - numProcs (no effect because procedural)
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

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

public class Nbodies {

  //numBodies, time, max initial velocity for generated bodies
  private int n, t, init = 25;
  private static int collisions = 0;
  //seed for random generation
  private long seed = 12345678;
  //radius, mass, length of increment, gravity multiplier
  private double r = 10.0,  dt = 1, m, gm = 1.0;
  private Pnt[] p, v, f; // points, velocities, forces
  private Boolean[][] c; // whether a collision occured last increment
  // whether to display gui, whether to ignore time limit, whether to nudge overlapping bodies
  private Boolean graphics = false, infinite = false, nudge = true;
  private double[][] dist;
  private String input, filename = "output.txt";
  private double G = 6.67e-11;

  private static int height = 500;
  private static int width = 500;
  private static JFrame frame;
  
  /*-----------------------
   *   main
   ----------------------*/
  public static void main(String[] args) {
    
    // create the simulation with initial values
    Nbodies simulation = new Nbodies(args);
    
    // read the time
    long startTime = System.currentTimeMillis();
    
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
    
    // run it
    simulation.simulate();
    
    // print final time
    long endTime = System.currentTimeMillis();
    long seconds = (endTime - startTime) / 1000;
    long millis = (endTime % 1000) - (startTime % 1000);
    if ((endTime % 1000) < (startTime % 1000)) {
      millis += 1000;
      seconds--;
    }
    System.out.println("computation time: "+seconds+" seconds "+millis+"milliseconds");
    System.out.println("total collisions: "+collisions+"\n");
    
  }
    
  /*-----------------------
   *   constructor
   ----------------------*/
  public Nbodies(String[] args) {
    
    if (args.length < 4) {
      System.out.println("Required arguments: numProcs numBodies mass time");
    }
    
    // read the command line arguments
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

    f = new Pnt[n]; // forces acting on each body (initialized to 0)
    c = new Boolean[n][n];
    for (int i = 0; i < n; i++) {
      f[i] = new Pnt();
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
    
  } // constructor
  
  /*------------------------------------------------------
   *   simulate
   *   
   *   manages the calculating methods and the file output
   -----------------------------------------------------*/
  public void simulate() {

    // print initial status
    Boolean writingOut = !filename.equals("none");
    OutputWriter writer = new OutputWriter(filename);
    if (writingOut) {
      writer.write(0);
    }
    
    // loop through each step
    int i = 0;
    Boolean time = i < t;
    while (time) {
      for (int j = 0; j < n; j++) {
        f[j] = new Pnt();
      }
      calcForce();
      moveBodies();
      collide();
      if (graphics) {
        frame.repaint();
      }
      if (writingOut) writer.write(i+1);
      i++;
      if (!infinite) time = i < t;
    }
    
    // close output file
    if (writingOut) writer.close();
    
  }
  
  /*--------------------------------------------------
   *   calcForce
   *   
   *   calculates forces between every pair of bodies
   --------------------------------------------------*/
  public void calcForce() {
  
    double mag; // magnitude
    Pnt dir;  // direction
    for (int i = 0; i < n; i++) {
      for (int j = i+1; j < n; j++) {
        dist[i][j] = Math.sqrt( Math.pow((p[i].x-p[j].x), 2) 
                              + Math.pow((p[i].y-p[j].y), 2) );
        dist[j][i] = dist[i][j];
        mag = ( G * Math.pow(m, 2) / Math.pow(dist[i][j], 2) );
        dir = new Pnt(p[j].x-p[i].x, p[j].y-p[i].y);
        f[i].x += mag*dir.x/dist[i][j];
        f[j].x -= mag*dir.x/dist[i][j];
        f[i].y += mag*dir.y/dist[i][j];
        f[j].y -= mag*dir.y/dist[i][j];
      }
    }

  }
  
  /*------------------------------------------------------
   *   moveBodies
   *   
   *   calculates new velocity and position for each body
   ------------------------------------------------------*/
  public void moveBodies() {
  
    Pnt deltav = new Pnt(), deltap = new Pnt();
    for (int i = 0; i < n; i++) {
      deltav.set( (f[i].x/m*dt), (f[i].y/m*dt) );
      deltap.set( (v[i].x + deltav.x/2) * dt,
                  (v[i].y + deltav.y/2) * dt);
      v[i].x += deltav.x;
      v[i].y += deltav.y;
      p[i].x += deltap.x;
      p[i].y += deltap.y;
      
      // detect collisions with walls of frame
      if (p[i].x <= r) {
        p[i].x = r;
        v[i].x *= -1;
      }
      if (p[i].x >= (width-r)) {
        p[i].x = width-r;
        v[i].x *= -1;
      }
      if (p[i].y <= r) {
        p[i].y = r;
        v[i].y *= -1;
      }
      if (p[i].y >= (height-r)) {
        p[i].y = height-r;
        v[i].y *= -1;            
      }
    }
  
  }
  
  /*------------------------------------------------------
   *   collide
   *   
   *   checks for collisions and calculates new velocities
   ------------------------------------------------------*/
  public void collide() {
    
    for (int i = 0; i < n; i++) {
      for (int j = i+1; j < n; j++) {
        if (dist[i][j] <= 2*r) {
          if (!c[i][j]) {
            c[i][j] = true;
            collisions++;
          }
          if (nudge) {
            // nudge overlapping bodies away from each other
            double overlap = (2*r)/dist[i][j];
            if (overlap > 1.2) {
              double xnudge = (p[j].x-p[i].x)*((overlap-1)/2);
              double ynudge = (p[j].y-p[i].y)*((overlap-1)/2);
              p[i].x += -1*xnudge;
              p[i].y += -1*ynudge;
              p[j].x += xnudge;
              p[j].y += ynudge;
            }
          }
          // recalculate velocities for each points
          // this is a nightmare to read, but it seems to work
          double dx = Math.pow(p[i].x - p[j].x, 2);
          double dy = Math.pow(p[i].y - p[j].y, 2);
          double dx2 = Math.pow(dx, 2);
          double dy2 = Math.pow(dy, 2);
          // equation A
          v[i].x = ( (v[j].x*dx2 + v[j].y*dx*dy +
                      v[i].x*dy2 - v[i].y*dx*dy)
                      / (dx2 + dy2) );
          // equation B
          v[i].y = ( (v[j].x*dx*dy + v[j].y*dy2 -
                      v[i].x*dy*dx + v[i].y*dx2)
                      / (dx2 + dy2) );
          // equation C
          v[j].x = ( (v[i].x*dx2 + v[i].y*dx*dy +
                      v[j].x*dy2 - v[j].y*dx*dy)
                      / (dx2 + dy2) );
          // equation D
          v[j].y = ( (v[i].x*dx*dy + v[i].y*dy2 -
                      v[j].x*dy*dx + v[j].y*dx2)
                      / (dx2 + dy2) );
        } else {
          c[i][j] = false;
        }
      }
    }   
  } // collide
  
  // get methods
  public Pnt[] getPoints() {
    return p;
  }
  
  public double getRadius() {
    return r;
  }
  
  public Boolean getGraphicsMode() {
    return graphics;
  }

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
   
} // Nbodies

/*----------------------------------------------------
 *  class BodiesPanel
 *  
 *  JPanel to animate the bodies on
 ----------------------------------------------------*/
class BodiesPanel extends JPanel {

  private static final long serialVersionUID = 1L;
  private int height, width;
  private double r;
  private Pnt[] points;
  private Color[] colors;
  
  public BodiesPanel(Nbodies sim, int width, int height) {
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
