#include "tsp.h" 
#include <unistd.h>


/**************************************************************************************************************************/
void debug(const char *err) { printf("\nDEBUG: %s \n", err); fflush(NULL); }
/**************************************************************************************************************************/


/**************************************************************************************************************************/
void print_error(const char *err) { printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1); }
/**************************************************************************************************************************/


/***************************************************************************************************************************/
int xpos(int i, int j, instance *inst)                                          
/***************************************************************************************************************************/
{ 
	if ( i == j ) print_error(" i == j in xpos" );
	if ( i > j ) return xpos(j,i,inst);
	int pos = i * inst->nnodes + j - (( i + 1 ) * ( i + 2 )) / 2;
	return pos;
}


/**************************************************************************************************************************/
int xpos_1(int i, int j, instance *inst) { return inst->xstart + i * inst->nnodes + j; }  
/**************************************************************************************************************************/


/**************************************************************************************************************************/
int upos(int i, instance *inst) { return inst->ustart + i; }
/**************************************************************************************************************************/


/**************************************************************************************************************************/
int ypos(int i, int j, instance *inst) { return inst->ystart + i * inst->nnodes + j; }       
/**************************************************************************************************************************/


/**************************************************************************************************************************/
double dist(int i, int j, instance *inst)
/**************************************************************************************************************************/
{
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j]; 
	if ( !inst->integer_costs ) return sqrt(dx*dx+dy*dy);
	int dis = sqrt(dx*dx+dy*dy) + 0.499999999; // nearest integer 
	return dis+0.0;
}


/**************************************************************************************************************************/
double cost(instance* inst, double* xstar)
/**************************************************************************************************************************/
{
	double c = 0;

	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = i + 1; j < inst->nnodes; j++)
		{
			if (xstar[xpos(i, j, inst)] > 0.5)
				c = c + dist(i, j, inst);
		}
	}

	return c;
}


/**************************************************************************************************************************/
int contains(int x, int* vec, int size)
/**************************************************************************************************************************/
{
	for (int i = 0; i < size; i++)
	{
		if (x == vec[i])
			return x;
	}

	return -1;
}


/**************************************************************************************************************************/
int findNext(int x, int* succ, int* nodes, int size)
/**************************************************************************************************************************/
{
	int next;
	int temp = succ[x];
	while (1)
	{
		next = contains(temp, nodes, size);
		if (next >= 0)
			return  next;

		temp = succ[temp];
	}

	print_error("findNext");
	return -1;
}


/**************************************************************************************************************************/
int is_in_tabu_list(int e1, int e2, int* tabu_list, int size)
/**************************************************************************************************************************/
{
	for (int i = 0; i < size; i++)
	{
		if (tabu_list[i] == e1 || tabu_list[i] == e2)
			return 1;
	}

	return -1;
}


/**************************************************************************************************************************/
double min(double x1, double x2)
/**************************************************************************************************************************/
{
	return x1 < x2 ? x1 : x2;
}


/**************************************************************************************************************************/
void write_sol(instance *inst) 
/**************************************************************************************************************************/
{
	FILE *fp;
	fp = fopen("plot.dat", "w");

	if ( inst->model_type == 0 || 
		 inst->model_type == 3 ||
		 inst->model_type == 30 || 
		 inst->model_type == 31 ||
		 inst->model_type == 310) // loop and callback
	{
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			for ( int j = i+1; j < inst->nnodes; j++ )
			{
				if (inst->best_sol[xpos(i,j,inst)] > 0.5) {
					//fprintf(fp,"%f %f\n",inst->xcoord[i], inst->ycoord[i]);
					//fprintf(fp,"%f %f\n",inst->xcoord[j], inst->ycoord[j]);
					fprintf(fp,"%f %f %d\n", inst->xcoord[i], inst->ycoord[i], i);
					fprintf(fp,"%f %f %d\n", inst->xcoord[j], inst->ycoord[j], j);
					fprintf(fp, "\n");
				}
			}
		}
	}

	else if ( inst->model_type == 1 ||
	 		  inst->model_type == 11||
	 		  inst->model_type == 2 ||
	 		  inst->model_type == 21||
	 		  inst->model_type == 22 ) // MTZ and F1
	{
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			for ( int j = 0; j < inst->nnodes; j++ )
			{
				if (inst->best_sol[xpos_1(i,j,inst)] > 0.5) {
					fprintf(fp,"%f %f\n",inst->xcoord[i], inst->ycoord[i]);
					fprintf(fp,"%f %f\n",inst->xcoord[j], inst->ycoord[j]);
					fprintf(fp, "\n");
				}
			}
		}		
	}

	else
		print_error("printing model type unknown!");

   	fclose(fp);
}


/**************************************************************************************************************************/
//void write_sol_name(instance *inst, double *xstar, int n) 
/**************************************************************************************************************************/
/*
{
	FILE *fp;
	char file_name[40];
	sprintf(file_name, "plot/dat/plot_%d.dat", n);
	fp = fopen(file_name, "w");

	if ( inst->model_type == 0 || 
		 inst->model_type == 3 ||
		 inst->model_type == 30 || 
		 inst->model_type == 31 ||
		 inst->model_type == 310) // loop and callback
	{
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			for ( int j = i+1; j < inst->nnodes; j++ )
			{
				if (xstar[xpos(i,j,inst)] > 0.5) {
					//fprintf(fp,"%f %f\n",inst->xcoord[i], inst->ycoord[i]);
					//fprintf(fp,"%f %f\n",inst->xcoord[j], inst->ycoord[j]);
					fprintf(fp,"%f %f %d\n", inst->xcoord[i], inst->ycoord[i], i);
					fprintf(fp,"%f %f %d\n", inst->xcoord[j], inst->ycoord[j], j);
					fprintf(fp, "\n");
				}
			}
		}
	}

	else if ( inst->model_type == 1 ||
	 		  inst->model_type == 11||
	 		  inst->model_type == 2 ||
	 		  inst->model_type == 21||
	 		  inst->model_type == 22 ) // MTZ and F1
	{
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			for ( int j = 0; j < inst->nnodes; j++ )
			{
				if (xstar[xpos_1(i,j,inst)] > 0.5) {
					fprintf(fp,"%f %f\n",inst->xcoord[i], inst->ycoord[i]);
					fprintf(fp,"%f %f\n",inst->xcoord[j], inst->ycoord[j]);
					fprintf(fp, "\n");
				}
			}
		}		
	}

	else
		print_error("printing model type unknown!");

   	fclose(fp);
}
*/


/**************************************************************************************************************************/
void plot_sol() 
/**************************************************************************************************************************/
{
    FILE *fp;
    fp = popen("gnuplot -p","w");
    fprintf(fp, "load \"script.txt\"\n");
    pclose(fp);
}


/**************************************************************************************************************************/
//void plot_sols(instance *inst) 
/**************************************************************************************************************************/
/*
{
    FILE *fp;
   
   	for (int i = 0; i < inst->nplot; i++) {
   		fp = popen("gnuplot -p","w");
	   	fprintf(fp, "set terminal png medium\n");

	   	char cmd1[40];
	   	sprintf(cmd1, "set output \'plot/png/plot_%d.png\'\n", i);
	   	fprintf(fp, "%s", cmd1);
	   	
	   	char title[30];
	   	sprintf(title, "set title \"TSP GRAPH %d\"\n", i);
    	fprintf(fp, "%s", title);
	   	
	   	fprintf(fp, "set xlabel \"x\"\n");
	   	fprintf(fp, "set ylabel \"y\"\n");
	   	fprintf(fp, "set style line 1 \
				    linecolor rgb '#0060ad' \
				    linetype 1 linewidth 2 \
				    pointtype 7 pointsize 1.5\n");

	   	char cmd2[60];
	   	sprintf(cmd2, "plot \"plot/dat/plot_%d.dat\" with linespoints linestyle 1\n", i);
		fprintf(fp, "%s", cmd2);
		pclose(fp);
	}
}
*/


/**************************************************************************************************************************/
void read_input(instance *inst) // simplified CVRP parser, not all SECTIONs detected  
/**************************************************************************************************************************/
{
                            
	FILE *fin = fopen(inst->input_file, "r");
	if ( fin == NULL ) print_error("input file not found!");
	
	inst->nnodes = -1;
	inst->depot = -1;  
	inst->nveh = -1;

	char line[180];
	char *par_name;   
	char *token1;
	char *token2;
	
	int active_section = 0; // =1 NODE_COORD_SECTION, =2 DEMAND_SECTION, =3 DEPOT_SECTION 
	
	int do_print = ( VERBOSE >= 1000 );

	while ( fgets(line, sizeof(line), fin) != NULL ) 
	{
		if ( VERBOSE >= 2000 ) { printf("%s",line); fflush(NULL); }
		if ( strlen(line) <= 1 ) continue; // skip empty lines
	    par_name = strtok(line, " :");
		if ( VERBOSE >= 3000 ) { printf("parameter \"%s\" ",par_name); fflush(NULL); }

		if ( strncmp(par_name, "NAME", 4) == 0 ) 
		{
			active_section = 0;
			continue;
		}

		if ( strncmp(par_name, "COMMENT", 7) == 0 ) 
		{
			active_section = 0;   
			token1 = strtok(NULL, "");  
			if ( VERBOSE >= 10 ) printf(" ... solving instance %s with model %d\n\n", token1, inst->model_type);
			continue;
		}   
		
		if ( strncmp(par_name, "TYPE", 4) == 0 ) 
		{
			token1 = strtok(NULL, " :");  
			if ( strncmp(token1, "TSP",3) != 0 ) print_error("format error:  only TYPE == TSP implemented so far!!!!!!"); 
			active_section = 0;
			continue;
		}
		

		if ( strncmp(par_name, "DIMENSION", 9) == 0 ) 
		{
			if ( inst->nnodes >= 0 ) print_error("repeated DIMENSION section in input file");
			token1 = strtok(NULL, " :");
			inst->nnodes = atoi(token1);
			if ( do_print ) printf(" ... nnodes %d\n", inst->nnodes); 
			inst->demand = (double *) calloc(inst->nnodes, sizeof(double)); 	 
			inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double)); 	 
			inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double));    
			active_section = 0;  
			continue;
		}

		if ( strncmp(par_name, "CAPACITY", 8) == 0 ) 
		{
			token1 = strtok(NULL, " :");
			inst->capacity = atof(token1);
			if ( do_print ) printf(" ... vehicle capacity %lf\n", inst->capacity); 
			active_section = 0;
			continue;
		}


		if ( strncmp(par_name, "VEHICLES", 8) == 0 ) 
		{
			token1 = strtok(NULL, " :");
			inst->nveh = atoi(token1);
			if ( do_print ) printf(" ... n. vehicles %d\n", inst->nveh);  
			active_section = 0;
			continue;
		}


		if ( strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0 ) 
		{
			token1 = strtok(NULL, " :");
			if ( strncmp(token1, "EUC_2D", 6) != 0 ) print_error("format error:  only EDGE_WEIGHT_TYPE == EUC_2D implemented so far!!!!!!"); 
			active_section = 0;
			continue;
		}            
		
		if ( strncmp(par_name, "NODE_COORD_SECTION", 18) == 0 ) 
		{
			if ( inst->nnodes <= 0 ) print_error(" ... DIMENSION section should appear before NODE_COORD_SECTION section");
			active_section = 1;   
			continue;
		}
		
		if ( strncmp(par_name, "DEMAND_SECTION", 14) == 0 ) 
		{
			if ( inst->nnodes <= 0 ) print_error(" ... DIMENSION section should appear before DEMAND_SECTION section");
			active_section = 2;
			continue;
		}  

		if ( strncmp(par_name, "DEPOT_SECTION", 13) == 0 )  
		{
			if ( inst->depot >= 0 ) print_error(" ... DEPOT_SECTION repeated??");
			active_section = 3;   
			continue;
		}

		
		if ( strncmp(par_name, "EOF", 3) == 0 ) 
		{
			active_section = 0;
			break;
		}
		
			
		if ( active_section == 1 ) // within NODE_COORD_SECTION
		{
			int i = atoi(par_name) - 1; 
			if ( i < 0 || i >= inst->nnodes ) print_error(" ... unknown node in NODE_COORD_SECTION section");     
			token1 = strtok(NULL, " :,");
			token2 = strtok(NULL, " :,");
			inst->xcoord[i] = atof(token1);
			inst->ycoord[i] = atof(token2);
			if ( do_print ) printf(" ... node %4d at coordinates ( %15.7lf , %15.7lf )\n", i+1, inst->xcoord[i], inst->ycoord[i]); 
			continue;
		}    
		  
		if ( active_section == 2 ) // within DEMAND_SECTION
		{
			int i = atoi(par_name) - 1; 
			if ( i < 0 || i >= inst->nnodes ) print_error(" ... unknown node in NODE_COORD_SECTION section");     
			token1 = strtok(NULL, " :,");
			inst->demand[i] = atof(token1);
			if ( do_print ) printf(" ... node %4d has demand %10.5lf\n", i+1, inst->demand[i]); 
			continue;
		}  

		if ( active_section == 3 ) // within DEPOT_SECTION
		{
			int i = atoi(par_name) - 1; 
			if ( i < 0 || i >= inst->nnodes ) continue;
			if ( inst->depot >= 0 ) print_error(" ... multiple depots not supported in DEPOT_SECTION");     
			inst->depot = i;
			if ( do_print ) printf(" ... depot node %d\n", inst->depot+1); 
			continue;
		}  
		
		printf(" final active section %d\n", active_section);
		print_error(" ... wrong format for the current simplified parser!!!!!!!!!");     
		    
	}                

	fclose(fin);    
	
}

/**************************************************************************************************************************/
void parse_command_line(int argc, char** argv, instance *inst) 
/**************************************************************************************************************************/
{ 
	
	if ( VERBOSE >= 100 ) printf(" running %s with %d parameters \n", argv[0], argc-1); 
		
	// default   
	inst->model_type = 0;
	inst->heuristic = -1;
	inst->old_benders = 0;
	strcpy(inst->input_file, "NULL");
	inst->randomseed = 0; 
	inst->num_threads = 0;
	inst->timelimit = CPX_INFBOUND+0.0;
	inst->tl = 600;
	inst->cutoff = CPX_INFBOUND; 
	inst->integer_costs = 0;

	inst->available_memory = 12000;   			// available memory, in MB, for Cplex execution (e.g., 12000)
	inst->max_nodes = -1; 						// max n. of branching nodes in the final run (-1 unlimited)        

    int help = 0; if ( argc < 1 ) help = 1;	
	for ( int i = 1; i < argc; i++ ) 
	{ 
		if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-input") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-f") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 				// input file
		if ( strcmp(argv[i],"-time_limit") == 0 ) { inst->timelimit = atof(argv[++i]); continue; }		// total time limit
		if ( strcmp(argv[i],"-tl") == 0 ) { inst->tl = atof(argv[++i]); continue; }						// internal time limit		
		if ( strcmp(argv[i],"-model_type") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 	// model type
		if ( strcmp(argv[i],"-model") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 			// model type
		if ( strcmp(argv[i],"-m") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 				// model type
		if ( strcmp(argv[i],"-heur") == 0 ) { inst->heuristic = atoi(argv[++i]); continue; } 			// heuristic
		if ( strcmp(argv[i],"-old_benders") == 0 ) { inst->old_benders = atoi(argv[++i]); continue; } 	// old benders
		if ( strcmp(argv[i],"-seed") == 0 ) { inst->randomseed = abs(atoi(argv[++i])); continue; } 		// random seed
		if ( strcmp(argv[i],"-s") == 0 ) { inst->randomseed = abs(atoi(argv[++i])); continue; } 		// random seed
		if ( strcmp(argv[i],"-threads") == 0 ) { inst->num_threads = atoi(argv[++i]); continue; } 		// n. threads
		if ( strcmp(argv[i],"-memory") == 0 ) { inst->available_memory = atoi(argv[++i]); continue; }	// available memory (in MB)
		//if ( strcmp(argv[i],"-node_file") == 0 ) { strcpy(inst->node_file,argv[++i]); continue; }		// cplex's node file
		if ( strcmp(argv[i],"-max_nodes") == 0 ) { inst->max_nodes = atoi(argv[++i]); continue; } 		// max n. of nodes
		if ( strcmp(argv[i],"-cutoff") == 0 ) { inst->cutoff = atof(argv[++i]); continue; }				// master cutoff
		if ( strcmp(argv[i],"-int") == 0 ) { inst->integer_costs = 1; continue; } 						// inteher costs
		if ( strcmp(argv[i],"-help") == 0 ) { help = 1; continue; } 									// help
		if ( strcmp(argv[i],"--help") == 0 ) { help = 1; continue; } 									// help
		if ( strcmp(argv[i],"-h") == 0 ) { help = 1; continue; } 										// help		
		help = 1;
    }      

	if ( help || (VERBOSE >= 10) )		// print current parameters
	{
		printf("\n\nAVAILABLE PARAMETERS -------------------------------------------------------------------------\n");
		printf("-file %s\n", inst->input_file); 
		printf("-time_limit %lf\n", inst->timelimit); 
		printf("-tl %lf\n", inst->tl); 
		printf("-model_type %d\n", inst->model_type);
		printf("-heuristic %d\n", inst->heuristic); 
		printf("-old_benders %d\n", inst->old_benders); 
		printf("-seed %d\n", inst->randomseed); 
		printf("-threads %d\n", inst->num_threads);  
		printf("-max_nodes %d\n", inst->max_nodes); 
		printf("-memory %d\n", inst->available_memory); 
		printf("-int %d\n", inst->integer_costs); 
		//printf("-node_file %s\n", inst->node_file);
		printf("-cutoff %lf\n", inst->cutoff); 
		printf("\nenter -h, -help or --help for help\n");

		printf("\nAVAILABLE MODELS -----------------------------------------------------------------------------\n");
		printf("0:  loop\n");
		printf("1:  MTZ\n");
		printf("11: MTZ with lazy constraints\n");
		printf("2:  F1\n");
		printf("21: F1 with n-2\n");
		printf("22: F1 with n-2 and lazy constraints\n");
		printf("3:  callback\n");
		printf("30: callback + patching\n");
		printf("31: generic callback\n");		
		printf("310:generic callback + patching\n");
		printf("\nAVAILABLE HEURISTICS -------------------------------------------------------------------------\n");
		printf("0:  hard fixing 0 (fixed threshold)\n");
		printf("1:  hard fixing 1 (decreasing threshold)\n");
		printf("2:  hard fixing 2 (random threshold)\n");
		printf("3:  local branching\n");
		printf("4:  nearest neighbourhood\n");
		printf("5:  GRASP\n");
		printf("6:  insertion\n");
		printf("7:  random insertion\n");
		printf("8:  GRASP + 2-opt move\n");
		printf("9:  VNS\n");
		printf("10: tabu search\n");
		printf("11: simulated annealing\n");		
		printf("----------------------------------------------------------------------------------------------\n\n");
	}        
	
	if ( help ) exit(1);

}    

/**************************************************************************************************************************/
void free_instance(instance *inst)
/**************************************************************************************************************************/
{   
	free(inst->demand);	
	free(inst->xcoord);
	free(inst->ycoord);
	//free(inst->load_min);
	//free(inst->load_max);
	free(inst->best_sol);
	if (inst->model_type == 30 || inst->model_type == 310) {
		for (int i = 0; i < inst->nthread; i++)
			free(inst->patch_vec[i]);
		free(inst->patch_vec);
		free(inst->patch_cost);
	}
}
