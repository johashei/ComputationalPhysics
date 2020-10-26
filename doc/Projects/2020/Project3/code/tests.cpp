#include "tests.hpp"

Tests::Tests(PhysicsSimulator& system){
cout << "constructor"<<endl;
  m_system = &system;
  m_system->get_shape(m_shape);
}

int Tests::test_energy_conservation(double tolerance){
  double E0 = kinetic_E(0) + potential_E_dynamic(0) + potential_E_fixed(0);
  double EN = kinetic_E(m_shape[0]-1) + potential_E_dynamic(m_shape[0]-1) + potential_E_fixed(m_shape[0]-1);
  double DE = EN - E0;
  cout << "DE/E = " << DE<<"/"<<E0 <<" = "<<DE/E0<<  " Mo AU^2 yr^-2." << endl;
  if(abs(DE)<=tolerance){
    cout << "Energy is conserved within specified tolerance :" << tolerance << " Mo AU^2 yr^-2." << endl;
    return(1);
  }
  else{
    cout << "Energy is not conserved." << endl;
    return(0);
  }
}

int Tests::test_angular_momentum_conservation(double tolerance){
  vec3 L0 = angular_momentum(0);
  vec3 LN = angular_momentum(m_shape[0]-1);
  vec3 DL = L0 - LN;
  cout << "DL/L = "<< DL<<"/"<<L0<<" = "<<DL/L0 << " Mo Au^2 yr^-1." << endl;
  if(DL.length()<=tolerance){
    cout << "Angular momentum is conserved within specified tolerance :" << tolerance << " Mo Au^2 yr^-1." << endl;
    return(1);
  }
  else{
    cout << "Angular momentum is not conserved."<< endl;
    return(0);
  }
}

double Tests::kinetic_E(unsigned k){
//cout<<"K"<<endl;
  /* Return the kinetic energy of the system K = sum(Ki) at time t_k*/
  double K = 0;

  for(unsigned i=0;i<m_shape[1];i++){ // kinetic energy of each body : 0.5*m*v^2
//cout <<i<<" "<<flush;
    K = K + m_system->get_object(k,i).Get_mass()*m_system->get_object(k,i).Get_vel().lengthSquared();
  }
//cout<<endl;
  K = 0.5*K;
  return(K);
}

double Tests::potential_E_dynamic(unsigned k){
//cout<<"V"<<endl;
  /* Return potential energy of the system V = -G*sum( mj*mi/|ri-rj| ) at time t_k  */
  if(m_shape[1]<2){return 0;} // If only one moving object
  double V = 0;
  double d;
  vec3 dr;
  double Mm;
  for(unsigned i=0;i<m_shape[1]-1;i++){
    for(unsigned j=i+1;j<m_shape[1];j++){
      dr = m_system->get_object(k,i).Get_pos() - m_system->get_object(k,j).Get_pos(); // vector between objects i and j
      d = sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]); // distance between objects i and j
      Mm = m_system->get_object(k,i).Get_mass() * m_system->get_object(k,j).Get_mass();
      V = V + Mm/d;
    }
  }
  V = -4*M_PI*M_PI * V;
  return(V);
}

double Tests::potential_E_fixed(unsigned k){
  if(m_system->m_N_fixed==0){return 0;} // If no fixed objects
  double V = 0;
  double Mm;
  for(unsigned i=0;i<m_shape[1];i++){
    for(unsigned j=0;j<m_system->m_N_fixed;j++){
      Mm = m_system->get_object(k,i).Get_mass() * m_system->m_fixed[j].Get_mass();
      V = V + Mm/m_system->get_object(k,i).Get_pos().length();
    }
  }
  V = -4*M_PI*M_PI * V;
  return(V);
}

vec3 Tests::angular_momentum(unsigned k){
//cout<<"L"<<endl;
  vec3 L(0,0,0);
  vec3 r,p;
  for(unsigned i=0;i<m_shape[1];i++){
    r = m_system->get_object(k,i).Get_pos();
    p = m_system->get_object(k,i).Get_vel() * m_system->get_object(k,i).Get_mass();
    L = L + r.cross(p);
  }
  return(L);
}
