#include "physics_simulator.hpp"

PhysicsSimulator::PhysicsSimulator(){
  m_N_objects = 0;
}

void PhysicsSimulator::add_object(PhysicsObject object){
  m_N_objects++;
  m_initial.push_back(object);
}

void PhysicsSimulator::add_objects(const char* posvelFilename, const char* massFilename){
  /* Read mass and initial conditions from files, generate the corresponding
  PhysicsObject and add them to m_initial. */
  double x,y,z,vx,vy,vz,m;
  FILE *fp_init = fopen(posvelFilename, "r");
  FILE *fp_mass = fopen(massFilename, "r");

  for(unsigned i=0;i<9;i++){
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
  m_gen_rel = gen_rel;
  m_h = m_simulation_time/(m_N_steps-1); // actually m_N_{data points}
  m_objects = new PhysicsObject[m_N_objects*m_N_steps];
}

PhysicsObject PhysicsSimulator::get_object(unsigned k, unsigned i){
  return(m_objects[k*m_N_objects + i]);
}

void PhysicsSimulator::get_shape(unsigned int shape[2]){
  shape[0] = m_N_steps;
  shape[1] = m_N_objects;
}

void PhysicsSimulator::run(){
  set_initial_conditions();
  velocityVerlet();
}

void PhysicsSimulator::write_to_file(string filename){
  unsigned kN, idx;
  ofstream outfile;
  outfile.open(filename);
  outfile << "["<<m_N_steps << ", " << m_N_objects << "]\n";
  for(unsigned k=0;k<m_N_steps;k++){
    kN = k*m_N_objects;
    outfile << k*m_h << ";  ";
    for(unsigned i=0;i<m_N_objects;i++){
      idx = kN + i;
      outfile << m_objects[idx].Get_pos() << " " << m_objects[idx].Get_vel() << "   ";
    }
    outfile <<"\n";
  }
  outfile.close();
}

void PhysicsSimulator::set_initial_conditions(){
  // Calculate Center of Mass
  vec3 rCoM, vCoM;
  double M,m;
  for(unsigned i=0;i<m_N_objects;i++){
    m = m_initial[i].Get_mass();
    rCoM = rCoM + m_initial[i].Get_pos()*m;
    vCoM = vCoM + m_initial[i].Get_vel()*m;
    M = M + m;
  }
  rCoM = rCoM/M;
  vCoM = vCoM/M;
  cout << "CoM = "<<rCoM<< " "<<vCoM << endl;
  // Move to CoM system
  for(unsigned i=0;i<m_N_objects;i++){
    m_objects[i].Set_pos(m_initial[i].Get_pos() - rCoM);
    m_objects[i].Set_vel(m_initial[i].Get_vel() - vCoM);
    m_objects[i].Set_mass(m_initial[i].Get_mass());
  }
}

void PhysicsSimulator::velocityVerlet(){

  cout<<"Values to set outside loop:"<< endl;
  double h5 = m_h*0.5;
  double hh5 = m_h*m_h*0.5;
  cout << "Values to be updated each iteration:" << m_N_objects << endl;
  vec3 acc[m_N_objects], accNew[m_N_objects]; // Store all accelerations
  for(unsigned i=0;i<m_N_objects;i++){// initialize acceleration
    acc[i] = acceleration(0,i);
  }
  unsigned idx,kN,j;

  // Simulation loop
  cout << "Starting loop" << endl;
  for(unsigned k=0;k<m_N_steps;k++){
    // cout << k << " => " << flush;
    kN = k*m_N_objects;
    for(unsigned i=0;i<m_N_objects;i++){ // Calculate all new positions
      idx = kN + i; // index of the object (k,i)
      j = idx + m_N_objects; // index of object (k+1,i)
      m_objects[j].Set_pos( m_objects[idx].Get_pos() + m_h*m_objects[idx].Get_vel() + hh5*acc[i]);
      m_objects[j].Set_mass(m_objects[idx].Get_mass());
    }
    for(unsigned i=0;i<m_N_objects;i++){  // Calculate new accelerations
      accNew[i] = acceleration(k+1,i);
    }
    for(unsigned i=0;i<m_N_objects;i++){ // Calculate new velocities
      idx = kN + i; // index of the object (k,i)
      j = idx + m_N_objects; // index of object (k+1,i)
      m_objects[j].Set_vel( m_objects[idx].Get_vel() + h5*(acc[i] + accNew[i]));
      acc[i] = accNew[i];
    }
    // cout << endl;
  }
  cout << "done" << endl;
}

vec3 PhysicsSimulator::acceleration(unsigned k,unsigned i){
//  cout << "acceleration" << endl;
  unsigned kN = k*m_N_objects;

  double G = 4*M_PI*M_PI; // AU^2 yr^-2 Mo^-1
  vec3 ri = m_objects[kN+i].Get_pos();
//   cout << i << ri<<endl;

  vec3 a, rj, dr;
  double mj, d;
  bool noti; // avoid if test inside loop

  //cout << "start loop" << endl;
  for(unsigned j=0;j<m_N_objects;j++){
    //cout << "iteration" << j << endl;
    //noti = j!=i;
    if(j!=i){ // Calculateing if false gives division by zero
      mj = m_objects[kN+j].Get_mass(); // These are inlined functions
      rj = m_objects[kN+j].Get_pos();
      dr = ri-rj;
      d = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
      a = a + mj*dr/(d*d*d);
  //    cout <<j<<  rj<<mj<<" ";
    }
  }
//  cout << endl;
  a = -G*a;
//  cout << a << endl;

  return(a);

}
