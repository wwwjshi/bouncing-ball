// room depth, z
int depth = 400;

// ball properties
ArrayList<PImage> textures = new ArrayList<PImage>();
ArrayList<Float> densities = new ArrayList<Float>();

// balls in room
ArrayList<Ball> balls = new ArrayList<Ball>();


void setup() {
  size(600, 400, P3D);
  
  // read all texture imgs and load as PImage
  File f = new File(dataPath("")); // data folder
  String imgnames[] = f.list();
  for (String name : imgnames){
    PImage img = loadImage(name);
    textures.add(img);
    densities.add(random(5,10)); // generate random density associated to each texture
  }
}


void draw() {
  background(0);
  // display walls
  walls(); 
  
  // display balls
  noStroke();
  for(Ball b : balls){
    b.display();
    b.move();
    for(Ball otherBall: balls){
      if(otherBall != b){
        b.ballCollision(otherBall);
      }
    }
  }
 
}


// shoots a new ball when mouse clicked
void mouseClicked(){
  int randIndex = int(random(0, textures.size()));
  PImage img = textures.get(randIndex); 
  float density = densities.get(randIndex);
  
  // at a click, pass in a ball at location of mouse with random texture and associated density   
  Ball b = new Ball(mouseX, mouseY, 0, img, density);
  balls.add(b);
}
 
  
void walls(){
  fill(255);
  stroke(255);
  line(0, 0, 0, 0, 0, -depth); 
  line(0, 0, -depth, 0, height, -depth); 
  line(0, height, -depth, 0, height, 0); 
  line( width, 0, 0,  width, 0, -depth); 
  line( width, 0, -depth,  width, height, -depth); 
  line( width, height, 0,  width, height, -depth); 
  line(0, 0, -depth,  width, 0, -depth); 
  line(0, height, -depth,  width, height, -depth); 
}
  
  
  
/*========================================================================
==            Ball class and Force calculation related variables                                 
==========================================================================*/

// properties/variables/coeficents for force calculation
float rho = 1.2; // density of air in temp=20 degree  [kg/m^3]
float mu = 1.8 * pow(10,-5); // viscosity of air in temp=20 degree  [kg/m*s^2]
// PVector gravity = new PVector(0, 9.8, 0); // acceleration due to gravity [m^2/s]
PVector gravity = new PVector(0, 0.98, 0); // less gravity for better visualisation
float dragCoefSphere = 0.47;
float frictionCoef = 0.1;


// ball class
class Ball{
  // size, mass
  PShape shape;
  float radius;
  float mass;
  
  // location, velocity, acceleration, rotation
  PVector location;
  PVector prevLocation;
  PVector velocity;
  PVector acceleration;
  PVector rotation;
  
  // Forces acting on a ball  N = kg*m/s^2
  PVector gravityForce;
  PVector dragForce;
  PVector buoyancyForce;
  PVector collideForce; // force exert by other ball when collided
  PVector frictionForce;
  PVector netForce; // overall
  
  // variable to check 'bad collision' with other ball
  /* bad collision means that the location correction after collision 
     results in the ball been collided with another ball or wall
     And so location correction is countinously applied, making the ball 'trembling'
  */
  int collided; 

  
  // Constructor
  Ball(float x, float y, float z, PImage textureimg, float density){
    // random radius between 10 to ~50
    this.radius = random(10,50);
    this.mass = 4/3 * PI * pow(radius, 3) * density; // mass of ball  = volume * density
    
    // sphere shape and apply texture
    this.shape = createShape(SPHERE, radius);
    this.shape.setTexture(textureimg);

    // location, intial x,y = where mouse is clicked, z = screen wall
    location = new PVector(x, y, z-radius); 
    prevLocation = location.copy();

    // veloctiy, v = u + a*t, initially, ball is thrown at some random velocity & direction
    velocity = new PVector(random(-100, 100), random(-100, 100), random(-100, 100));
    
    // acceleration, determined by the forces acting on the ball
    acceleration = new PVector(0, 0, 0);
    
    // rotation/spinning 
    rotation = new PVector(0, 0, 0);

    // Force
    netForce = new PVector(0, 0, 0);
    frictionForce = new PVector(0, 0, 0);
    collideForce = new PVector(0, 0, 0);

    collided = 0;
  }
  
  
  // Display ball object 
  void display(){
    noStroke();
    pushMatrix();
    translate(this.location.x, this.location.y, this.location.z); 
    
    rotateX(rotation.x);
    rotateY(rotation.y);
    rotateZ(rotation.z);
    
    shape(shape);
    popMatrix();
  }
  
  
  // move ball object
  void move(){
    // if ball is badly collided with other ball, make it stationary
    if(collided > 1 && netForce.mag() < gravity.mag()){
      velocity.set(0,0,0);
    // otherwise move the ball as it is supposed to
    } else {
      // update force, acceleration, velocity, location
      updateNetForce(); // netForce at this frame and so updates acceleration
      velocity.add(acceleration); // acculmuate velocity v=u+at
      prevLocation = location.copy();
      location.add(velocity); // acculmuate location r = sum(v*t)
      
      // check for collision with walls and rebound
      bounceWall();
      
      // check for stationary balls and 
      PVector locChange = PVector.sub(prevLocation, location);
      if(locChange.mag() == 0){
        velocity.set(0,0,0);
      }
      
      // make ball rotation wrt to amount of change in displacement
      locChange.x = locChange.x/width * 2*PI;
      locChange.y = locChange.y/height * 2*PI;
      locChange.z = locChange.z/depth * 2*PI;
      rotation.add(locChange);
    }
    netForce.set(0,0,0);
    acceleration.set(0,0,0); // don't accumulate acceleration as it changes each time/frame
  }
  
   
  // function to calculate netForce and updates the acceleration 
  void updateNetForce(){
    // add force due to collision if any
    netForce.add(collideForce);
    collideForce.set(0,0,0); // collision is an impulse, only exert once, set to zero once added
  
    // gravity fg = m*g 
    gravityForce = PVector.mult(gravity, mass);
    netForce.add(gravityForce);
    
    // drag 
    dragForce = dragForce(rho, mu, velocity);
    netForce.add(dragForce);
    
    // buoyancy
    buoyancyForce = PVector.mult(gravity, -4/3 * PI * pow(radius,3) * rho);
    netForce.add(buoyancyForce);
    
    // friction, if rolling on ground (xz-plane of y=0)
    frictionForce.set(0,0,0);
    if(location.y >= height-radius) {
      frictionForce = velocity.copy().normalize().mult(-1); // acts in opposite position of velocity
      frictionForce.mult(netForce.y*frictionCoef); // f = normal * mu, since rolling on flat ground, normal = netForce.y
      frictionForce.y = 0; // rolling on ground, no y direction friction
      netForce.add(frictionForce);
    }
    
    // get acceleration due to net force
    acceleration = PVector.div(netForce, mass);
  }
  

  // check collide with wall, and bounce back
  void bounceWall(){
    // check left/right wall
    if(location.x > width-this.radius){
      velocity.x *= -1;
      location.x = width-this.radius;
    }
    if(location.x < this.radius){
      velocity.x *= -1;
      location.x = this.radius;
    }
    // check bot/top wall
    if(location.y > height-this.radius){
      velocity.y *= -1;
      location.y = height-this.radius;
    }
    if(location.y < this.radius){
      velocity.y *= -1;
      location.y = this.radius;
    }
    // check front/back wall
    if(location.z < -depth+this.radius){
      velocity.z *= -1;
      location.z = -depth+this.radius;
    }
    if(location.z > this.radius){
      velocity.z *= -1;
      location.z = this.radius;
    }
  }

  
  // function to check collision with another ball
  void ballCollision(Ball other){
    // get distance vector between current and other ball
    PVector distVect = PVector.sub(this.location, other.location); //  sub(curr, other) => result vector point from other to curr
    
    // get magnitude of distance
    float distMag = distVect.mag();
    
    // allowed distance between two ball 
    float allowedDist = this.radius + other.radius;
    
    // check collision and reaction
    if (distMag < allowedDist) {
      
      // other ball
      PVector otherDir = distVect.copy().normalize().mult(-1);
      other.velocity = otherDir.copy().setMag(other.velocity.mag()*0.8);
      other.location.add(other.velocity);
      other.bounceWall();
      // current exert force to other ball
      if(other.velocity.mag() == 0){
        other.collideForce = otherDir.copy().setMag(this.mass);
      }

      // current ball 
      // change movement direction
      this.velocity = distVect.copy().normalize().setMag(velocity.mag()*0.8);
      this.location.add(velocity);
      bounceWall();
      
      
      distVect = PVector.sub(this.location, other.location);
      distMag = distVect.mag();
      allowedDist = this.radius + other.radius;
      
      if (distMag < allowedDist) {
        this.location.sub(velocity);
        collided += 1;
      } 
      else {
        if(collided > 1){
          collided -= 1;
        }
      }

    }
  }
  
  
  // function to calculate drag Force
  PVector dragForce(float rho, float mu, PVector velocity){
    float dragMag;
    // get Reynolds number
    float re = rho * 2*radius * velocity.mag() / mu;
    // calculate drag base on Re
    if(re < 1){
      // Fd = 6*pi*mu*r*v
      dragMag = 6 * PI * mu * radius * velocity.mag();
    } else {
      // Fd = 1/2*rho*Cd*A*v^2
      dragMag = 0.5 * rho * dragCoefSphere * PI*pow(radius,2) * pow(velocity.mag(),2);
    }
    PVector dragForce = velocity.copy().normalize().mult(-1); // opposite direction of velocity
    dragForce.setMag(dragMag);
    return dragForce;
  }  
}
