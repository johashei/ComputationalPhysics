#include "physics_simulator.hpp"

PhysicsSimulator::PhysicsSimulator(){
  m_N_objects = 0; // initialize number of objects to simulate
  m_N_fixed = 0;  // initialize number of fixed objects
  m_gen_rel = 0;  // default: no relativistic correction
  m_CoM = true;   // default: simulate in center of mass system
}

void PhysicsSimulator::add_fixed_object(PhysicsObject object){
  m_N_fixed++;
  m_CoM = false; // Stay in the input coordinate system
  m_fixed.push_back(object);
}

void PhysicsSimulator::add_object(PhysicsObject object){
  m_N_objects++;
  m_initial.push_back(object);
}

void PhysicsSimulator::add_objects(unsigned int N, const char* posvelFilename, const char* massFilename){
  /* Read mass and initial conditions from files, generate the corresponding
  PhysicsObject and add them to m_initial. */
  double x,y,z,vx,vy,vz,m;
  FILE *fp_init = fopen(posvelFilename, "r");
  FILE *fp_mass = fopen(massFilename, "r");

  for(unsigned i=0;i<N;i++){
    fscanf(fp_init, "%lf %lf %lf %lf %lf %lf",&x,&y,&z,&vx,&vy,&vz);
    fscanf(fp_mass, "%lf", &m);
    vec3 pos(x,y,z);
    vec3 vel(vx,vy,vz);
    PhysicsObject planet(pos,vel*365.25,m/2e30);
    add_object(planet);
  }
  fclose(fp_init);
  fclose(fp_mass);
}

void PhysicsSimulator::set_parameters(double simulation_time, unsigned int steps, bool gen_rel){
  m_simulation_time = simulation_time;
  m_N_steps = steps;
  m_h = m_simulation_time/(m_N_steps-1); // actually m_N_{data points}
  m_objects = new PhysicsObject[m_N_objects*m_N_steps];
  if(gen_rel){
    if(m_N_objects!=1 || m_N_fixed!=1){
      cout << "Sorry, general relativity is only implemented for one moving and one fixed object"<<endl;
    }
    else{
      vec3 r = m_initial[0].Get_pos() - m_fixed[0].Get_pos();
      vec3 l = m_initial[0].Get_pos().cross(m_initial[0].Get_vel());
      double c = 63198; // AU/yr
      m_gen_rel = 3*l.lengthSquared() / (r.lengthSquared()*c*c);
      cout << "r x v = " << r <<" x "<< m_initial[0].Get_vel() << " = " << l << endl;
      cout << "l^2 = " << l.lengthSquared() << endl;
      cout << "r^2 = " << r.lengthSquared() << endl;
      cout << "general relativity correction :" << m_gen_rel << endl;
    }
  }

}


void PhysicsSimulator::get_shape(unsigned int shape[2]){
  shape[0] = m_N_steps;
  shape[1] = m_N_objects;
}

void PhysicsSimulator::run(string method){
  set_initial_conditions();
  clock_t start,end;
  if(method=="velocityVerlet"){
    start = clock();
    velocityVerlet();
    end = clock();
  }
  else if(method=="EulerCromer"){
    start = clock();
    EulerCromer();
    end = clock();
  }
  else if(method=="Euler"){
    start = clock();
    Euler();
    end = clock();
  }
  else{cout<< method<<" is not implemented.";exit(1);}
  double algotime = (double) (end-start)/CLOCKS_PER_SEC;
  cout << "Simulation completed in "<< algotime << " s."<< endl;
}

void PhysicsSimulator::write_to_file(string filename,unsigned step){
cout << "write to file" << endl;
  unsigned kN, idx;
  ofstream outfile;
  outfile.open(filename);
  outfile << "["<<m_N_steps/step << ", " << m_N_objects << "]\n";
  for(unsigned k=0;k<m_N_steps/step;k++){
    kN = k*m_N_objects*step;
    outfile << k*m_h << ";  ";
    for(unsigned i=0;i<m_N_objects;i++){
      idx = kN + i;
      outfile << m_objects[idx].Get_pos() << " " << m_objects[idx].Get_vel() << "   ";
    }
    outfile <<"\n";
  }
  outfile.close();
}

vec3 PhysicsSimulator::perihelion(unsigned int k, unsigned int i){
  unsigned j = (k-1)*m_N_objects + i;
  double rr[3];
  rr[0] = m_objects[j+m_N_objects].Get_pos().length();
  rr[1] = m_objects[j].Get_pos().length();
  rr[2] = m_objects[j-m_N_objects].Get_pos().length();

  while(!(rr[0] > rr[1] && rr[1] < rr[2])){
    j = j-m_N_objects;
    rr[0] = rr[1];
    rr[1] = rr[2];
    rr[2] = m_objects[j-m_N_objects].Get_pos().length();
  }
  cout << rr[1]-rr[0]<<","<<rr[1]-rr[2]<<endl;
  return(m_objects[j].Get_pos());
}

void PhysicsSimulator::set_initial_conditions(){
  vec3 rCoM, vCoM; // Stay 0 if m_CoM == false.
  if(m_CoM){ // Calculate Center of Mass
    double M = 0;
    double m;
    for(unsigned i=0;i<m_N_objects;i++){
      m = m_initial[i].Get_mass();
      rCoM = rCoM + m_initial[i].Get_pos()*m;
      vCoM = vCoM + m_initial[i].Get_vel()*m;
      M = M + m;
    }
    rCoM = rCoM/M;
    vCoM = vCoM/M;
    cout << "CoM = "<<rCoM<< " "<<vCoM << endl;
  }
  // Move to CoM system
  for(unsigned i=0;i<m_N_objects;i++){
    m_objects[i].Set_pos(m_initial[i].Get_pos() - rCoM);
    m_objects[i].Set_vel(m_initial[i].Get_vel() - vCoM);
    m_objects[i].Set_mass(m_initial[i].Get_mass());
  }
}

void PhysicsSimulator::velocityVerlet(){
  cout << "Simulating with velocity Verlet." << endl;
//cout<<"Values to set outside loop:"<< endl;
  double h5 = m_h*0.5;
  double hh5 = m_h*m_h*0.5;
//cout << "Values to be updated each iteration:" << m_N_objects << endl;
  vec3 acc[m_N_objects], accNew[m_N_objects]; // Store all accelerations
  for(unsigned i=0;i<m_N_objects;i++){// initialize acceleration
    acc[i] = acceleration(0,i);
  }
  unsigned j,kN,jn;

  // Simulation loop
//cout << "Starting loop" << endl;
  for(unsigned k=0;k<m_N_steps-1;k++){
    kN = k*m_N_objects;
    for(unsigned i=0;i<m_N_objects;i++){ // Calculate all new positions
      j = kN + i; // index of the object (k,i)
      jn = j + m_N_objects; // index of object (k+1,i)
      m_objects[jn].Set_pos( m_objects[j].Get_pos() + m_h*m_objects[j].Get_vel() + hh5*acc[i]);
      m_objects[jn].Set_mass(m_objects[j].Get_mass());
    }
    for(unsigned i=0;i<m_N_objects;i++){  // Calculate new accelerations
      accNew[i] = acceleration(k+1,i);
    }
    for(unsigned i=0;i<m_N_objects;i++){ // Calculate new velocities
      j = kN + i; // index of the object (k,i)
      jn = j + m_N_objects; // index of object (k+1,i)
      m_objects[jn].Set_vel( m_objects[j].Get_vel() + h5*(acc[i] + accNew[i]));
      acc[i] = accNew[i];
    }
  }
//cout << "done" << endl;
}

void PhysicsSimulator::Euler(){
  cout << "Simulating with Euler." << endl;
  unsigned j,kN,jn;
  for(unsigned k=0;k<m_N_steps-1;k++){
    kN = k*m_N_objects;
    for(unsigned i=0;i<m_N_objects;i++){
      j = kN + i; // index of object (k,i)
      jn = j + m_N_objects; // index of object (k+1,i)
      m_objects[jn].Set_vel(m_objects[j].Get_vel() + acceleration(k,i)*m_h);
      m_objects[jn].Set_pos(m_objects[j].Get_pos() + m_objects[j].Get_vel()*m_h);
      m_objects[jn].Set_mass(m_objects[j].Get_mass());
    }
  }
}

void PhysicsSimulator::EulerCromer(){
  cout << "Simulating with EulerCromer." << endl;
  unsigned j,kN,jn;
  for(unsigned k=0;k<m_N_steps-1;k++){
    kN = k*m_N_objects;
    for(unsigned i=0;i<m_N_objects;i++){
      j = kN + i; // index of object (k,i)
      jn = j + m_N_objects; // index of object (k+1,i)
      m_objects[jn].Set_vel(m_objects[j].Get_vel() + acceleration(k,i)*m_h);
      m_objects[jn].Set_pos(m_objects[j].Get_pos() + m_objects[jn].Get_vel()*m_h);
      m_objects[jn].Set_mass(m_objects[j].Get_mass());
    }
  }
}

vec3 PhysicsSimulator::acceleration(unsigned k,unsigned i){
  unsigned kN = k*m_N_objects;

  double G = 4*M_PI*M_PI; // AU^2 yr^-2 Mo^-1
  vec3 ri = m_objects[kN+i].Get_pos();

  vec3 a, rj, dr;
  double mj, d;
//cout << "moving objects loop" << endl;
  for(unsigned j=0;j<m_N_objects;j++){
    if(j!=i){ // Calculateing if false gives division by zero
      mj = m_objects[kN+j].Get_mass(); // These are inlined functions
      rj = m_objects[kN+j].Get_pos();
      dr = ri-rj;
      d = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
      a = a + mj*dr/pow(d,3);
    }
  }
//cout << a << endl;
//cout << "fixed objects loop" << endl;
  vec3 a_from_fixed;
  if(!m_CoM){ // if there are fixed objects
//cout << m_N_fixed << endl;
    for(unsigned j=0;j<m_N_fixed;j++){
//cout << j << endl;
      mj = m_fixed[j].Get_mass();
      rj = m_fixed[j].Get_pos();
      dr = ri-rj;
      d = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
      a_from_fixed = a_from_fixed + mj*dr/pow(d,3);
//cout << a_from_fixed << endl;
    }
  }
//cout << a_from_fixed << endl;

  a = -G*(a + a_from_fixed)*(1 + m_gen_rel);

  return(a);

}
