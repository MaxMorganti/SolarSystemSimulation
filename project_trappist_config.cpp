
//code fully generalized to multiple bodies 

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>


using namespace std;

inline double sqr(double num)
{return num*num;}

double rhs (double t, double data[][2][2], int body_id, int coord_id, int i, double mass_array[]);

void euler (double t, double h, double data[][2][2], int body_id, int coord_id, double mass_array[],
		double (*f)(double t, double data[][2][2], int body_id, int coord_id, int i, double mass_array[]));
	
void update_2D (double time, ofstream& output, double data[][2][2], int body_id);
		


const int body_num = 7;     //number of bodies in system - can be freely changed but must also change arrays of initial conditions

int
main (void)
{
   
   double h = 2;  //mesh spacing for Euler alg. - use bigger size for many bodies or longer simulation times to conserve space
                  //but this will give larger error
   
   //astronomical quantities defined - numbers given in meters, kg, seconds
   double earth_rad = 6300000.0;
   double AU = 1.496e11;
   double earth_mass = 5.972e24;
   double moon_mass = 7.347e22;
   double solar_mass = 1.989e30;
  
   //declare main data array: main_array[body_num][dim][N]   :  key: dimension 0=x,1=y and N+1 is deriv order
   //could potentially easily expand to 3 dimensions by changing the second array dimension size to be 3
   //init. conditions format key: {{{body1 x, body1 x'}, {body1 y, body1 y'} }, {{body2 x, body2 x'}, {body2 y, body2 y'}}, ...etc.}
   //initial configuration given is for Trappist-1 system assuming all circular orbits and starting in alignment
   double main_array[body_num][2][2] = {{{0, 0},{0, 0}},
                                        {{0, 80300},{.011*AU, 0}},
					{{0, -68000},{-.015*AU, 0}},
					{{.0214*AU, 0},{0, -57600}},
					{{-.028*AU, 0},{0, 50300}},
					{{0, 43800},{.037*AU, 0}},
					{{0, -39700},{-.045*AU, 0}}};				
						
   //most probable values of planetary masses 										
   double mass[body_num] = {.08*solar_mass, .85*earth_mass, 1.38*earth_mass, .41*earth_mass, .62*earth_mass, .68*earth_mass, 1.34*earth_mass}; 
   //alternative mass configuration below
   //double mass[body_num] = {.08*solar_mass, .2*earth_mass, 1.9*earth_mass, .65*earth_mass, .62*earth_mass, .68*earth_mass, 1.34*earth_mass};  
   //alternative mass config. 2
   //double mass[body_num] = {.08*solar_mass, .2*earth_mass, 4*earth_mass, 2*earth_mass, .62*earth_mass, .68*earth_mass, 1.34*earth_mass};
   
   //initiate clocks for returning runtime
   clock_t t1, t2;
   t1 = clock();
   
   double tmin = 0;
   double tmax = 15000000;    //10e6 seconds is over 100 days - enough for 10 periods of every body in Trappist-1 sytem
   
   
   //create arrays of ofstreams and their filename streams and then string arrays for final filenames
   ofstream out[body_num];
   ostringstream  file_name_stream[body_num];
   string file_name[body_num];
   
   
   //loop over individual file for each body
   for (int i = 0; i < body_num; i++)
   {
     //open each outfile - give name with ostringstream
     file_name_stream[i] << "Trappist_" << i << ".dat";
     file_name[i] = file_name_stream[i].str();
     out[i].open(file_name[i].c_str());
     
     //print header columns and initial conditions
     out[i] << "#        t          x(t)          x'(t)         y(t)         y'(t)" <<endl;
     update_2D (tmin, out[i], main_array, i);
   }
   
   //looping over whole simulation time
   for (double t = tmin; t <= tmax; t+=h) {
     for (int i = 0; i < body_num; i++) {
       //euler alg. call for each coord for each body
       euler (t, h, main_array, i, 0, mass, rhs);
       euler (t, h, main_array, i, 1, mass, rhs);
       
       //update each data file, one for each body	
       update_2D (t + h, out[i], main_array, i);
     }
     //progress counter
     if (fmod(t,(tmax/10)) == 0)
      {cout << (t/tmax)*100 << "% completion..."<< endl;}	 
   }
   
   //out messages for user
   t2 = clock();
   float diff((float)t2 - (float)t1);
   cout << "Runtime: " << diff/CLOCKS_PER_SEC << " seconds." << endl;
   
   
   for (int i = 0; i < body_num; i++) {
    cout << "Body " << i << " data stored in " << file_name[i] << endl;
    out[i].close();
   }
   
   return(0);

}


//*******************************************************************************************
//
//                                     *FUNCTIONS*
//
//*******************************************************************************************

//right hand side of diff eq. of order i+1, when i = 1 this returns right side of graviational accel. equation
double rhs (double t, double data[][2][2], int body_id, int coord_id, int i, double mass_array[])
{
  switch(i) {
    case 0:
      return data[body_id][coord_id][1];
    
    case 1:
      double output = 0;
      //summing over gravitational contributions from every other body 
      for (int j = 0; j < body_num; j++) {
        if (j != body_id) {
	
	  double dx = (data[body_id][0][0] - data[j][0][0]);  
          double dy = (data[body_id][1][0] - data[j][1][0]);
	  double rad = sqrt(sqr(dx) + sqr(dy));
	  double norm[2] = {(dx/rad), (dy/rad)};
	
	  output += (-6.67e-11*mass_array[j]*norm[coord_id])/sqr(rad);
	}
       }
       return output;
   }
   return(1);
}


//Euler algorithm taken from equation 6.45 in the class notes
void euler (double t, double h, double data[][2][2], int body_id, int coord_id, double mass_array[],
		double (*f)(double t, double data[][2][2], int body_id, int coord_id, int i, double mass_array[]))
{ 
  for (int deriv_order = 0; deriv_order < 2; deriv_order++)
      {data[body_id][coord_id][deriv_order] += h*f(t, data, body_id, coord_id, deriv_order, mass_array);}
}



//updates each column in the output file
void update_2D (double time, ofstream& output, double data[][2][2], int body_id)
{
  output << setw (10) << time << "  " 
         << setw (12) << data[body_id][0][0] << "  " 
         << setw (12) << data[body_id][0][1] << "  "
         << setw (12) << data[body_id][1][0] << "  " 
         << setw (12) << data[body_id][1][1] << "  " << endl;	 
}
   
   
