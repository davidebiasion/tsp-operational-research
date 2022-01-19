#ifndef VRP_H_  

#define VRP_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>  

#include <cplex.h>  
#include <pthread.h>  


#define VERBOSE				   60		// printing level  (=0 csv file, =10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)


//hard-wired parameters
#define XSMALL		  		  1e-5 		// 1e-4*	// tolerance used to decide ingerality of 0-1 var.s
#define EPSILON		  		  1e-9		// 1e-9		// very small numerical tolerance 
#define EPS 				  1e-5
#define TICKS_PER_SECOND 	  1000.0  	// cplex's ticks on Intel Core i7 quadcore @2.3GHZ
                   

//data structures  
typedef struct {   
	
	//input data
	int nnodes; 	
	double *demand;   
	double *xcoord;
	double *ycoord;
	int depot;
	double capacity; 
	int nveh;

	// parameters 
	int model_type; 
	int heuristic;
	int old_benders;
	int randomseed;
	int num_threads;
	double timelimit;						// overall time limit, in sec.s
	char input_file[1000];		  			// input file
	//char node_file[1000];		  			// cplex node file
	int available_memory;
	int max_nodes; 							// max n. of branching nodes in the final run (-1 unlimited)
	double cutoff; 							// cutoff (upper bound) for master
	int integer_costs;
	double tl;								// internal time limit

	//global data
	double	tstart;								
	double zbest;							// best sol. available  
	double tbest;							// time for the best sol. available  
	double *best_sol;						// best sol. available    
	double *last_sol;						 
	double	best_lb;						// best lower bound available  
	double *load_min;						// minimum load when leaving a node
	double *load_max;						// maximum load when leaving a node
	int nthread;							// number of threads running on this program
	double** patch_vec;						// thread safe storage of patching heuristics
	double* patch_cost;						// costs of the patching heuristcs (default -1)

	// model;     
	int xstart;
	int ustart;
	int qstart;
	int bigqstart;  						
	int sstart;
	int bigsstart;
	int ystart;
	int fstart;
	int zstart;
	int ncols;

	// plot
	int nplot;
} instance;        


//inline
inline int imax(int i1, int i2) { return ( i1 > i2 ) ? i1 : i2; } 
inline double dmin(double d1, double d2) { return ( d1 < d2 ) ? d1 : d2; } 
inline double dmax(double d1, double d2) { return ( d1 > d2 ) ? d1 : d2; } 


// tsp opt
int TSPopt(instance *inst);
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
void build_model_0(instance *inst, CPXENVptr env, CPXLPptr lp);
void build_model_1(instance *inst, CPXENVptr env, CPXLPptr lp);
void build_model_11(instance *inst, CPXENVptr env, CPXLPptr lp);
void build_model_2(instance *inst, CPXENVptr env, CPXLPptr lp);
void build_model_21(instance *inst, CPXENVptr env, CPXLPptr lp);
void build_model_22(instance *inst, CPXENVptr env, CPXLPptr lp);
void build_sol(const double* xstar, instance* inst, int* succ, int* comp, int* ncomp);
void build_sol_0(const double* xstar, instance* inst, int* succ, int* comp, int* ncomp);
void build_sol_1(const double* xstar, instance* inst, int* succ, int* comp, int* ncomp);
int addSEC(instance* inst, CPXENVptr env, CPXLPptr lp, int* succ, int* comp, int* ncomp);
void patching(instance* inst, double* xstar, int thread_idx);
int myseparation(instance* inst, double* xstar, CPXCENVptr env, void* cbdata, int wherefrom);
int callback_candidate(instance* inst, CPXCALLBACKCONTEXTptr context);
int callback_heuristic(instance* inst, CPXCALLBACKCONTEXTptr context);
int salvagnin_separation(instance* inst, double* xstar, CPXCALLBACKCONTEXTptr context, int thread_idx);


// tsp approx
int TSPapprox(instance *inst);
void heuristic(instance *inst, CPXENVptr env, CPXLPptr lp);
int hard_fixing_0(instance *inst, CPXENVptr env, CPXLPptr lp);
int hard_fixing_1(instance *inst, CPXENVptr env, CPXLPptr lp);
int hard_fixing_2(instance *inst, CPXENVptr env, CPXLPptr lp);
int local_branching(instance *inst, CPXENVptr env, CPXLPptr lp);
int nearest_neighbour_manager(instance* inst);
int grasp_manager(instance* inst);
int nearest_neighbour(instance* inst, int start, double* xstar);
int grasp(instance* inst, int start, double* xstar);
int closest_node(instance* inst, int node, int* visited);
int rand_closest_node(instance* inst, int node, int* visited);
int insertion(instance* inst, int node1, int node2);
void get_hab(instance* inst, int* visited, double* xstar, int* min_hab);
int rand_insertion_manager(instance* inst);
int rand_insertion(instance* inst, int node1, int node2, double* xstar);
void rand_get_hab(instance* inst, int* visited, double* xstar, int* min_hab);
int two_opt_move(instance * inst, double *xstar);
int vns(instance *inst);
int tabu(instance* inst);
int simulated_annealing(instance* inst);


// tsp utilities
void debug(const char *err);
void print_error(const char *err);
int xpos(int i, int j, instance *inst);  
int xpos_1(int i, int j, instance *inst);
int upos(int i, instance *inst);
int ypos(int i, int j, instance *inst);                              
double dist(int i, int j, instance *inst);
double cost(instance* inst, double* xstar);
int findNext(int x, int* succ, int* nodes, int size);
int contains(int x, int *vec, int size);
int is_in_tabu_list(int e1, int e2, int* tabu_list, int tabu_size);
double min(double x1, double x2);
void write_sol(instance *inst);
void plot_sol();
//void plot_sols(instance *inst);
//void write_sol_name(instance *inst, double *xstar, int n);
void read_input(instance *inst);
void parse_command_line(int argc, char** argv, instance *inst); 
void free_instance(instance *inst);


// chrono
double second();  


#endif   /* VRP_H_ */ 
