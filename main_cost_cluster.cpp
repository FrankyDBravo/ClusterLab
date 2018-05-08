#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cstddef>      // std::size_t
#include <cmath>        // std::pow
#include <math.h>      // M_PI
#include <set>
#include <chrono>
#include <ctime>
#include <string>
#include <iomanip> // setprecision
#include <sstream> // stringstream
#include <cstdlib> // for exit()
#include <algorithm>  // find
#include <unordered_map>

using namespace std;


// -------------- CLASSES ---------------------
struct waiting_time{
  vector <int> sources;
  vector <double> current_time;
  vector <double> cost;
};

class lattice
{
  public:
    // mutation_position for each dimension
    set<long int> mutation_position_1d;
    vector<long int> mutation_position_vec1;
    set<vector <double> > mutation_position_d;
    // mutation_number per time unit
    vector<long int> number_mutated_pert;
    vector<double> mean_r;   // will calculate the mean of the power law
    vector<double> t; // current time to keep track of
    vector<double> time_save; // time that will be saved
    vector<double> true_time_jump; // time that will be saved
    vector<long int> number_mutated_save;// #mut sites that will be saved
    unordered_map<long int, waiting_time> waiting_sites;  //current site waiting to be approved
    unordered_map<long int, double> old_site_cost;  //old site used  for time

    vector<vector<double>> true_evo_link; // keep track of source and destination at each time

    double u; // mutation rate of deme
    double m; // migration Parameters
    int d;
    // constructor
    lattice(int);
    // evolution by one time step
    void lattice_evo(double, double, double);
    // track of the mutated demes
    void print_mutated_number();
    int current_mutated_number();
    int current_true_size();
    void save_mutated_number(string);
    void save_true_time(string);
    void save_evolution(string);

    double return_meanr();

  private:
    vector <double> start_times={1,2,3,4,5,6,7,8,9,10,25,50,100,250}; //times(mutation number) needed at the beginning

  };

  // ---------------- USEFULL FUNCTIONS -----------------–
bool sort_vect_first(vector<double> a, vector<double> b)
{
  return a[0]<b[0];
}
bool sort_vect_last(vector<double> a, vector<double> b)
{
  return a[4]<b[4];
}
//---------------- Problem related FUNCTIONS -----------------------
double cost_f(double r, double prefactor, double exponent){
    //double cost=0;
    double cost=prefactor*pow(r,exponent);
    return cost;
}
int update_cost(vector <double>& costs, double dt)
{
  for(double& c : costs){
    c += -dt;     // reduce all cost by dt
  }
  vector<double>::iterator it2;
  it2=min_element(costs.begin(), costs.end());
  if(*it2<=0){
    int indice=distance(costs.begin(),it2);

    return indice+1;
  }
  else{
    return 0;
  }

}
//--------function for simulation ----------
void simulate(int d, double mu, double max_size,double prefactor,double exponent, string file_name,bool print_bool)
{

  lattice test(d);// default migr_rate=2 & muta_rate=2

  auto start = chrono::system_clock::now();
  //for(int i=0; i<1100; i++)
  //while (test.current_mutated_number()<max_size)
  while (test.current_true_size()<max_size-1)
  {



    test.lattice_evo(mu,prefactor,exponent);
  }



  auto end = chrono::system_clock::now();
  chrono::duration<double> elapsed_seconds = end-start;
  test.print_mutated_number();
  cout << "\n elaspsed time: "<< elapsed_seconds.count() <<"s \n";

  if (print_bool){
      cout <<"\n mean r:, "<< test.return_meanr();
      test.save_true_time(file_name);
    //  test.save_evolution(file_name);
  }
  test.save_mutated_number(file_name);

}


// ------------------ MAIN -------------------
int main(int argc, char *argv[])
{

  srand (1);
  // -      -      -      -      -    -
  //int d;

  bool to_print=true;

  // -- --  -- --   --  --

  if (argc <= 1)
  	{
  		// On some operating systems, argv[0] can end up as an empty string instead of the program's name.
  		// We'll conditionalize our response on whether argv[0] is empty or not.
  		if (argv[0])
  			std::cout << "Usage: " << argv[0] << " <number>" << '\n';
  		else
  			std::cout << "Usage: <program name> <number>" << '\n';

  		exit(1);
  	}

    int d=stoi(argv[1]);
    cout << "Got integer: " << d << '\n';

    int number_file=stoi(argv[2]);
    double mu=stod(argv[3]);
    cout << "Got double: " << mu << '\n';

    double max_size=stod(argv[4]);
    cout << "Got double: " << max_size << '\n';

    double prefactor=stod(argv[5]);

    double exponent=stod(argv[6]);

    string id=argv[7];


  // find folder name
  string file_name="./data/outbreak";

  file_name=file_name+"/d"+argv[1]+"u"+argv[3]+"c"+argv[5]+"r"+argv[6]+"t"+argv[4]+"/number";
  cout << file_name << endl;



    //add file number
    simulate(d,mu,max_size,prefactor,exponent,file_name+id,to_print);


  
  return 0;
}

// ------------------ END ---------------------
// -------------------------------------------


// constructor
lattice::lattice(int dim){
  d=dim;
  t.push_back(1);
  time_save.push_back(1);
  number_mutated_save.push_back(1);
  old_site_cost[0]=0;

  mutation_position_1d.insert(0);
  mutation_position_vec1.push_back(0);
  number_mutated_pert.push_back(mutation_position_1d.size());
}
//----------- main function for evolution--------------

void lattice::lattice_evo(double mu, double prefactor, double exponent)
{

  // stuff required by random
  random_device rd;
  mt19937 gen(rd());  // choosing mersenne_twister_engine with predefined parameters
  uniform_real_distribution<double> unif(0, 1);  // defining the distribution
  uniform_int_distribution<int>  dir_1d(0,1);

  // problem variable
  double n=-(mu+d);
  double x_0=1;
  long int new_pos; // will store new calculated position


  // --- for dimension 1 ---
  if (d==1){
    int current_size=(mutation_position_vec1).size();
    //keep track of time$
    double dt=1./(current_size);
    t.push_back(t.back()+dt);
    // choose random site
    int rand_site=floor(unif(gen)*(current_size));
    long int mutated_dem_pos=mutation_position_vec1.at(rand_site);
    // choose direction
    int x=dir_1d(gen)*2-1;
    double u=unif(gen);
    double r=pow(u,1./(n+1));
    double cost=cost_f(r,prefactor,exponent);
    // if x0=1 can be reduce topow((-pow(x_0,(n+1))*u+pow(x_0,(n+1))), (1./(n+1))); pow(u,1/(n+1));
    if (x==-1){
      r=r-1;
    }
    mean_r.push_back(cost);
    new_pos=floor(r*x+mutated_dem_pos);
    if (not mutation_position_1d.count(new_pos)){

      waiting_sites[new_pos].sources.push_back(mutated_dem_pos);
      waiting_sites[new_pos].current_time.push_back(t.back());
      waiting_sites[new_pos].cost.push_back(cost);

    }
    unordered_map<long int, waiting_time>  :: iterator it;
    unordered_map<long int, waiting_time>  :: iterator old_it;

    for(it = waiting_sites.begin(); it != waiting_sites.end();){
      long int indice=update_cost((*it).second.cost,dt);
      if (indice){ //return the indice+1(avoid 0 problem) of the one to update
        new_pos=(*it).first; //redefine new_pos
        mutation_position_1d.insert(new_pos);
        mutation_position_vec1.push_back(new_pos);
        current_size=current_size+1;

        old_site_cost[new_pos]=old_site_cost[(*it).second.sources[indice-1]]+t.back()-(*it).second.current_time[indice-1];
        double true_time=(*it).second.current_time[indice-1]-old_site_cost[(*it).second.sources[indice-1]];
                  


        true_time_jump.push_back(true_time);
        if (current_size<100)
        {
        vector<double> link;
        link.push_back(true_time); //true time
        link.push_back((*it).second.sources[indice-1]); //source
        link.push_back((*it).second.current_time[indice-1]);//source time
        link.push_back(new_pos);  // target
        link.push_back(t.back()); // target time

        true_evo_link.push_back(link);//with first time
        }
        if ((current_size % 250 == 0) || (find(start_times.begin(), start_times.end(), current_size) != start_times.end()))
        {
          //cout << "current size="<<current_size <<endl;
          number_mutated_save.push_back(current_size);
          time_save.push_back(t.back());
        }
        old_it=it;
        ++it;
        waiting_sites.erase(old_it);

      }
        else{
        ++it;
        }


    }
    // keep trackof the new number of mutated deme
    number_mutated_pert.push_back(current_size);



  }
  else{
    cout << "Wrong dimensionality" ;
  }

}

//-------------------------- other functions --------------------------
//r mean
double lattice::return_meanr()
{
  double average = accumulate( mean_r.begin(), mean_r.end(), 0.0)/ mean_r.size();
  return average;
}

// print number of mutated site at all t
void lattice::print_mutated_number()
{
  cout << endl;
  cout << "number of step: "<< number_mutated_pert.size() <<endl;
  cout << "final size: "<< mutation_position_vec1.size() << endl;
  //cout << "final size: "<< true_time_jump.size() << endl;
}
void lattice::save_true_time(string name)
{
  ofstream outputfile2;
  outputfile2.open(name+"_true_time.txt" );
  ofstream outputfile;
  outputfile.open(name+"_true.txt");
  vector<double> true_t_save;
  vector<double> true_mutated_save={1};
  int size;
  sort(true_time_jump.begin(), true_time_jump.end()); //inplace
  vector<double> :: iterator tt;
  for(tt=true_time_jump.begin();tt!=true_time_jump.end();++tt)
  {
    size=true_mutated_save.back()+1;
    true_mutated_save.push_back(size);
      if ((size % 200 == 0) || (find(start_times.begin(), start_times.end(), size) != start_times.end()))
      {
        outputfile2<<*tt<<"\n";
        outputfile<<size<<"\n";

      }
  }
  outputfile.close();
  outputfile2.close();
}
int lattice::current_mutated_number()
{
  return number_mutated_pert.back();
}

int lattice::current_true_size()
{
  return true_time_jump.size();
}

// save number of mutated site at all t
void lattice::save_mutated_number(string name)
{

  //-- time file --
  ofstream outputfile2;
  outputfile2.open(name+"_time.txt" ); // opens file named "filename" for output

  for(int i=0; i<time_save.size();i++)
  {
    outputfile2<<time_save[i]<<"\n";
  }
  outputfile2.close();

  //-- # mut file ---
  ofstream outputfile;
  outputfile.open(name+".txt"); // opens file named "filename" for output
  for(int i=0;i<number_mutated_save.size();i++)
  {
        outputfile<<number_mutated_save[i]<<"\n";
    }
  outputfile.close();

}

void lattice::save_evolution(string name)
{

  //-- true time---
  ofstream outputfile;
  sort (true_evo_link.begin(), true_evo_link.end(), sort_vect_first);

  outputfile.open(name+"_source_true.txt" ); // opens file nam
  for(int i=0;i<true_evo_link.size();i++)
  {
    outputfile<<true_evo_link[i][0]<<" "<<true_evo_link[i][1]<<" "<<true_evo_link[i][3]<<"\n";
  }
  //------ actual time ----
  outputfile.close();
  ofstream outputfile2;
  sort (true_evo_link.begin(), true_evo_link.end(), sort_vect_last);
  outputfile2.open(name+"_source_actual.txt" ); // opens file nam
  for(int i=0;i<true_evo_link.size();i++)
  {
    outputfile2<<true_evo_link[i][2]<<" "<<true_evo_link[i][1]<<" "<<true_evo_link[i][3]<<" "<<true_evo_link[i][4]<<"\n";
  }
  outputfile2.close();
}
