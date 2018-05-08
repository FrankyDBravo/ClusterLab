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
#include <algorithm>  // find
#include <unordered_map>

using namespace std;

// -------------- CLASSES ---------------------

class lattice
{
public:
  // mutation_position for each dimension
    set<long int> mutation_position_1d;

    //store number of new mutation at each time
    vector<long int> number_mutated_pert;
    long int pop_size;

    //data needed for tracking
    vector<long int> waiting_sources;
    vector<double> waiting_sources_cost_run;
    vector<long int> parent;              //parent in graph
    vector<long int> parent_link_num;       //# step in graph
    vector<vector<double>> true_evo_link;  // keep track of source and destination at each time


    int d;  //dimension
    int k;  //number of edges allowed

    // constructor
    lattice(int,int);
    // evolution by one time step
    void lattice_evo(double,double,double);
    // track of the mutated demes
    int current_mutated_number();
    void save_mutated_number(string);
    void save_evolution(string);
    void print_mutated_number();

  };

  // ---------------- USEFULL FUNCTIONS -----------------–
bool sort_vect_first(vector<double> a, vector<double> b)
{
  return a[0]<b[0];
}
void print_vector(vector<long int> path)
{
  for (vector<long int>::const_iterator i = path.begin(); i != path.end(); ++i)
  {
      std::cout << *i << ' ';
    }
    cout <<endl;
}
void print_vector_double(vector<double> path)
{
  for (vector<double>::const_iterator i = path.begin(); i != path.end(); ++i)
  {
      std::cout << *i << ' ';
    }
    cout <<endl;
}
bool sort_vect_last(vector<double> a, vector<double> b)
{
  return a[4]<b[4];
}

//---------------- Problem related FUNCTIONS -----------------------

double cost_f(double r,double prefactor, double exponent){
    //double cost=0;
    double cost=prefactor*pow(r,exponent);
    return cost;
}

void update_cost(vector <double>& costs, double dt)
{
  for(double& c : costs){
    c += -dt;     // reduce all cost by dt
  }
}
//--------function for simulation ----------
void simulate(int d, double mu, int num_edges, double max_size,double prefactor,double exponent, string file_name,bool print_bool)
{

  lattice test(d,num_edges);
  cout<<max_size<<endl;
  auto start = chrono::system_clock::now();
  while (test.current_mutated_number()<max_size)
  {
    test.lattice_evo(mu,prefactor,exponent);
  }


  auto end = chrono::system_clock::now();
  chrono::duration<double> elapsed_seconds = end-start;
  test.print_mutated_number();
  cout << "\n elaspsed time: "<< elapsed_seconds.count() <<"s \n";

  if (print_bool){
    //  test.save_evolution(file_name);
  }
  test.save_mutated_number(file_name);
}

// -------------------------------------------
// ------------------ MAIN -------------------
// -------------------------------------------
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
    double num_edges=stod(argv[7]);


      // find folder name
      string file_name="./data/network";

    //add file name
    file_name=file_name+"/d"+argv[1]+"u"+argv[3]+"c"+argv[5]+"r"+argv[6]+"t"+argv[4]+"k"+argv[7]+"/number";
    cout << file_name << endl;


     simulate(d,mu,num_edges,max_size,prefactor,exponent,file_name+argv[8],to_print);

     

      return 0;
    }

    // ------------------ End -------------------
    // -------------------------------------------

    // constructor
    lattice::lattice(int dim, int num_edges){
      d=dim;
      k=num_edges;
        
      pop_size=1;
      parent.push_back(0);
      parent_link_num.push_back(0);
      waiting_sources.push_back(0);
      waiting_sources_cost_run.push_back(1);
      number_mutated_pert.push_back(0);
    }


//----------- main function for evolution--------------
    void lattice::lattice_evo(double mu,double prefactor, double exponent)
    {
      // stuff required by random
      random_device rd;
      mt19937 gen(rd());  // choosing mersenne_twister_engine with predefined parameters
      uniform_real_distribution<double> unif(0, 1);  // defining the distribution
      uniform_int_distribution<int>  dir_1d(0,1);
      // problem variable
      double n=-(mu+d);
      double x_0=1;

      // --- for dimension 1 ---
      if (d==1){

        if (waiting_sources.empty())
        {
          cout << "plus de source.. fini"<< endl;
        }

        // choose site and delete it
        vector<double>::iterator it2;
        //find new elem to update
        it2=min_element(waiting_sources_cost_run.begin(), waiting_sources_cost_run.end());
        //keep track of #links
        int indice=distance(waiting_sources_cost_run.begin(),it2);

        long int link_num=parent_link_num[indice];

        //print_vector(waiting_sources);
        long int mutated_dem_pos=waiting_sources[indice];
        double cost=waiting_sources_cost_run[indice];
        //cout<<"source _mut: "<<mutated_dem_pos<<" indice "<<indice<<endl;
        long int mother=parent[indice];
        waiting_sources.erase(waiting_sources.begin()+indice);
        waiting_sources_cost_run.erase (it2);
        parent.erase(parent.begin()+indice);
        parent_link_num.erase(parent_link_num.begin()+indice);


        // prevent for re-chosing a already done node
        if (not mutation_position_1d.count(mutated_dem_pos))
        {

          while (number_mutated_pert.size()<=link_num)
          {
            number_mutated_pert.push_back(0);
          }
          number_mutated_pert[link_num]+=1; // stores at the right number of links
          pop_size=pop_size+1;
          mutation_position_1d.insert(mutated_dem_pos);

          // ---- storing for graph ---
          if (pop_size<100)
          {
            vector<double> link;
            link.push_back(1); //true time
            link.push_back(mother); //source
            link.push_back(link_num);//source time
            link.push_back(mutated_dem_pos);  // target
            link.push_back(link_num+1); // target time
            true_evo_link.push_back(link);//with first time
          }

          for(int i = 0; i < k; ++i)
          {
            // choose direction
            int x=dir_1d(gen)*2-1;
            double u=unif(gen);
            double r=pow(u,1./(n+1));
            double cost_new=cost_f(r,prefactor,exponent);
            // if x0=1 can be reduce topow((-pow(x_0,(n+1))*u+pow(x_0,(n+1))), (1./(n+1))); pow(u,1/(n+1));
            if (x==-1){
              r=r-1;
            }

            long int new_pos=floor(r*x+mutated_dem_pos);
            if (not mutation_position_1d.count(new_pos))
            {
              parent_link_num.push_back(link_num+1);
              parent.push_back(mutated_dem_pos);
              waiting_sources.push_back(new_pos);
              waiting_sources_cost_run.push_back(cost_new+cost+1);
            }
          }
        }

      }
      else{
        cout << "Wrong dimensionality" ;
      }

    }
 //-------------------------- other functions --------------------------


    // save number of mutated site at all t
    void lattice::save_mutated_number(string name)
    {
      //-- # mut file ---
      ofstream outputfile;
      outputfile.open(name+".txt"); // opens file named "filename" for output
      for(int i=0;i<number_mutated_pert.size();i++)
      {
            outputfile<<number_mutated_pert[i]<<"\n";
        }
      outputfile.close();

    }
    void lattice::print_mutated_number()
    {
      cout << endl;
      cout << "number of step: "<< number_mutated_pert.size() <<endl;
      cout << "final size: "<< pop_size <<endl;
    }
    void lattice::save_evolution(string name)
    {
      //------ actual time ----
      ofstream outputfile2;
      sort (true_evo_link.begin(), true_evo_link.end(), sort_vect_last);
      outputfile2.open(name+"_source_actual.txt" ); // opens file nam
      for(int i=0;i<true_evo_link.size();i++)
      {
        outputfile2<<true_evo_link[i][2]<<" "<<true_evo_link[i][1]<<" "<<true_evo_link[i][3]<<" "<<true_evo_link[i][4]<<"\n";
      }
      outputfile2.close();
    }


    int lattice::current_mutated_number()
    {
      return pop_size;
    }
