# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cmath>

using namespace std;

int main ( int argc, char *argv[] );
void compute ( int np, int nd, double pos[], double vel[], 
  double mass, double f[], double &pot, double &kin );
double cpu_time ( );
double dist ( int nd, double r1[], double r2[], double dr[] );
void initialize ( int np, int nd, double pos[], double vel[], double acc[] );
void r8mat_uniform_ab ( int m, int n, double a, double b, int &seed, double r[] );
void timestamp ( );
void update ( int np, int nd, double pos[], double vel[], double f[], 
  double acc[], double mass, double dt );

//****************************************************************************80

int main ( int argc, char *argv[] ){
  freopen ("serialBig1.txt","w",stdout);
  cout<<"This work is contributed by Himangshu, Adarsh, Avaneesh, Abhinav"<<endl;
  double *acc;
  double ctime;
  double dt;
  double e0;
  double *force;
  double kinetic;
  double mass = 1.0;
  int nd;
  int np;
  double *pos;
  double potential;
  int step;
  int step_num;
  int step_print;
  int step_print_index;
  int step_print_num;
  double *vel;

  timestamp ( ); // Just print the start time
  cout << "\n";
  cout << " Serial version\n";
  cout << "  A molecular dynamics program.\n";
//
//  Get the spatial dimension.
//
  if ( 1 < argc )
  {
    nd = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter ND, the spatial dimension (2 or 3).\n";
    cin >> nd;
  }
//
//  Get the number of particles.
//
  if ( 2 < argc )
  {
    np = atoi ( argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter NP, the number of particles (500, for instance).\n";
    cin >> np;
  }
//
//  Get the number of time steps.
//
  if ( 3 < argc )
  {
    step_num = atoi ( argv[3] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter STEP_NUM, the number of time steps (500 or 1000, for instance).\n";
    cin >> step_num;
  }
//
//  Get the time step.
//
  if ( 4 < argc )
  {
    dt = atof ( argv[4] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter DT, the time step size (0.1, for instance).\n";
    cin >> dt;
  }
//
//  Report.
//
  cout << "\n";
  cout << "  ND, the spatial dimension, is " << nd << "\n";
  cout << "  NP, the number of particles in the simulation is " << np << "\n";
  cout << "  STEP_NUM, the number of time steps, is " << step_num << "\n";
  cout << "  DT, the size of each time step, is " << dt << "\n";
//
//  Allocate memory.
//
  acc = new double[nd*np];
  force = new double[nd*np];
  pos = new double[nd*np];
  vel = new double[nd*np];
//
//  This is the main time stepping loop:
//    Compute forces and energies,
//    Update positions, velocities, accelerations.
//
  cout << "\n";
  cout << "  At each step, we report the potential and kinetic energies.\n";
  cout << "  The sum of these energies should be a constant.\n";
  cout << "  As an accuracy check, we also print the relative error\n";
  cout << "  in the total energy.\n";
  cout << "\n";
  cout << "      Step      Potential       Kinetic        (P+K-E0)/E0\n";
  cout << "                Energy P        Energy K       Relative Energy Error\n";
  cout << "\n";

  step_print = 0;
  step_print_index = 0;
  step_print_num = 10;
  
  ctime = cpu_time ( );

  for ( step = 0; step <= step_num; step++ )
  {
    if ( step == 0 )
    {
      initialize ( np, nd, pos, vel, acc );
    }
    else
    {
      update ( np, nd, pos, vel, force, acc, mass, dt );
    }

    compute ( np, nd, pos, vel, mass, force, potential, kinetic );

    if ( step == 0 )
    {
      e0 = potential + kinetic;
    }
    if ( step == step_print )
    {
      cout << "  " << setw(8) << step
           << "  " << setw(14) << potential
           << "  " << setw(14) << kinetic
           << "  " << setw(14) << ( potential + kinetic - e0 ) / e0 << "\n";
      step_print_index = step_print_index + 1;
      step_print = ( step_print_index * step_num ) / step_print_num;
    }

  }
//
//  Report timing.
//
  ctime = cpu_time ( ) - ctime;
  cout << "\n";
  cout << "  Elapsed cpu time " << ctime << " seconds.\n";
//
//  Free memory.
//
  delete [] acc;
  delete [] force;
  delete [] pos;
  delete [] vel;
//
//  Terminate.
//
  cout << "\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void compute ( int np, int nd, double pos[], double vel[], double mass, 
  double f[], double &pot, double &kin ){
  double d;
  double d2;
  int i;
  int j;
  int k;
  double PI2 = 3.141592653589793 / 2.0;
  double rij[3];

  pot = 0.0;
  kin = 0.0;

  for ( k = 0; k < np; k++ )
  {
//
//  Compute the potential energy and forces.
//
    for ( i = 0; i < nd; i++ )
    {
      f[i+k*nd] = 0.0;
    }

    for ( j = 0; j < np; j++ )
    {
      if ( k != j )
      {
        d = dist ( nd, pos+k*nd, pos+j*nd, rij );
//
//  Attribute half of the potential energy to particle J.
//
        if ( d < PI2 )
        {
          d2 = d;
        }
        else
        {
          d2 = PI2;
        }

        pot = pot + 0.5 * pow ( sin ( d2 ), 2 );

        for ( i = 0; i < nd; i++ )
        {
          f[i+k*nd] = f[i+k*nd] - rij[i] * sin ( 2.0 * d2 ) / d;
        }
      }
    }
//
//  Compute the kinetic energy.
//
    for ( i = 0; i < nd; i++ )
    {
      kin = kin + vel[i+k*nd] * vel[i+k*nd];
    }
  }

  kin = kin * 0.5 * mass;
  
  return;
}
//****************************************************************************80

double cpu_time ( ){
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
//****************************************************************************80

double dist ( int nd, double r1[], double r2[], double dr[] ){
  // Dist between (x1,y1,z1) and (x1,y1,z1) where nd = 3 ==> sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2)
  double d;
  int i;

  d = 0.0;
  for ( i = 0; i < nd; i++ )
  {
    dr[i] = r1[i] - r2[i];
    d = d + dr[i] * dr[i];
  }
  d = sqrt ( d );

  return d;
}
//****************************************************************************80

void initialize ( int np, int nd, double pos[], double vel[], double acc[] ){
  int i;
  int j;
  int seed;
//
//  Set the positions.
//
  seed = 123456789;
  r8mat_uniform_ab ( nd, np, 0.0, 10.0, seed, pos );
//
//  Set the velocities.
//
  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      vel[i+j*nd] = 0.0;
    }
  }
//
//  Set the accelerations.
//
  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      acc[i+j*nd] = 0.0;
    }
  }
  return;
}
//****************************************************************************80

void r8mat_uniform_ab ( int m, int n, double a, double b, int &seed, double r[] ){
  int i;
  const int i4_huge = 2147483647;
  int j;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + i4_huge;
      }

      r[i+j*m] = a + ( b - a ) * ( double ) ( seed ) * 4.656612875E-10;
    }
  }

  return;
}
//****************************************************************************80

void timestamp ( ){
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

void update ( int np, int nd, double pos[], double vel[], double f[], 
  double acc[], double mass, double dt ){
  int i;
  int j;
  double rmass;

  rmass = 1.0 / mass;

  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      pos[i+j*nd] = pos[i+j*nd] + vel[i+j*nd] * dt + 0.5 * acc[i+j*nd] * dt * dt;
      vel[i+j*nd] = vel[i+j*nd] + 0.5 * dt * ( f[i+j*nd] * rmass + acc[i+j*nd] );
      acc[i+j*nd] = f[i+j*nd] * rmass;
    }
  }

  return;
}