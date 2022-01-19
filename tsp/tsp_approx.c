#include "tsp.h"
#include <time.h>


/**************************************************************************************************************************/
int TSPapprox(instance *inst)
/**************************************************************************************************************************/
{
	// open cplex model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP"); 
	if ( VERBOSE >= 60 ) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // cplex output on screen
	
	// cplex parameters
	CPXsetintparam(env, CPXPARAM_RandomSeed, inst->randomseed); // Cplex random seed
	CPXsetintparam(env, CPXPARAM_MIP_Strategy_RINSHeur, 10);

	// build model
	if(inst->model_type != 3 || inst->model_type != 30 || inst->model_type != 31 || inst->model_type != 310)
		inst->model_type = 310; //using generic callback
	build_model(inst, env, lp);

	inst->ncols = CPXgetnumcols(env, lp);
	inst->best_sol = (double*)calloc(inst->ncols, sizeof(double)); 	// all entries to zero  
	inst->best_lb = CPX_INFBOUND;
	inst->zbest = CPX_INFBOUND;

	srand(inst->randomseed);

	heuristic(inst, env, lp);
	
	/*
	// use solution
	if (inst->heuristic < 4)
		if ( CPXgetx(env, lp, inst->best_sol, 0, inst->ncols-1) ) print_error("CPXgetx");

	// find succ, comp and ncomp
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	int* comp = (int*)malloc(inst->nnodes * sizeof(int));
	int ncomp;
	build_sol(inst->best_sol, inst, succ, comp, &ncomp);

	if (ncomp > 1)
		if ( inst->model_type != 0 ) print_error("calling addSEC with model_type != 0");

	// write .dat file for gnuplot
	write_sol(inst);

	// print solution on screen
	if ( VERBOSE >= 70 )
	{
		int x = 0;
		for (int i = 0; i < inst->nnodes; i++)
		{
			printf("\nEdge for node %4d: ", x + 1);
			printf("[%4d, %4d] ", x + 1, succ[x] + 1);
			x = succ[x];
		}
		printf("\nThe number of connected components is: %d\n", ncomp);
	}
	
	free(succ);
	free(comp);
	*/

	// free pools and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 	
	
	return 0;
}


/***************************************************************************************************************************/
void heuristic(instance* inst, CPXENVptr env, CPXLPptr lp)
/***************************************************************************************************************************/
{
	switch (inst->heuristic) 
	{
		case 0: 
			hard_fixing_0(inst, env, lp);
			break;

		case 1:
			hard_fixing_1(inst, env, lp);
			break;

		case 2:
			hard_fixing_2(inst, env, lp);
			break;

		case 3:
			local_branching(inst, env, lp);
			break;

		case 4:
			nearest_neighbour_manager(inst);
			break;

		case 5:
			grasp_manager(inst);
			break;
		
		case 6:
			insertion(inst, 0, 1);
			break;
		
		case 7:
			rand_insertion_manager(inst);
			break;

		case 8:
			grasp_manager(inst);
			break;

		case 9:
			vns(inst);
			break;

		case 10:
			tabu(inst);
			break;

		case 11:
			simulated_annealing(inst);
			break;

		default:
			print_error("heuristic type unknown");
	}
}


/**************************************************************************************************************************/
int hard_fixing_0(instance *inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{
	double tl = inst->tl; 	// internal timelimit
	CPXsetdblparam(env, CPX_PARAM_TILIM, tl);
	double remaining_time = inst->timelimit;
	double ratio; 			// ratio best_lb/zbest
	double p; 				// probability to fix a variable
	int* index = (int*)malloc(inst->ncols * sizeof(int));
	double* bound0 = (double*)malloc(inst->ncols * sizeof(double));
	double* bound1 = (double*)malloc(inst->ncols * sizeof(double));
	char* lu = (char*)malloc(inst->ncols * sizeof(double));

	// add a heuristic starting solution
	double* value = (double*)malloc(inst->ncols * sizeof(double));
	inst->timelimit = tl;
	vns(inst);

	int izero = 0;
	int nnz = 0;
	for (int i = 0; i < inst->ncols; i++) {
		if (inst->best_sol[i] > 0.5) {
			index[nnz]   = i;
			value[nnz++] = 1;
		}
	}

	int effortlevel = 5;
	CPXaddmipstarts(env, lp, 1, nnz, &izero, index, value, &effortlevel, NULL);
	remaining_time -= tl;

	for (int i = 0; i < inst->ncols; i++)
	{
		index[i] = i;
		bound0[i] = 0.0;
		bound1[i] = 1.0;
	}

	//fix new variables
	for ( int i = 0; i < inst->ncols; i++ )  
	{
		if (inst->best_sol[i] > 0.5 ) 
		{
			p = (double) rand() / (RAND_MAX + 0.0);
			if ( p <= 0.9 )
				lu[i] = 'L';
			else
				lu[i] = 'U';
		}
		else
			lu[i] = 'U';
	}
	if ( CPXchgbds(env, lp, inst->ncols, index, lu, bound1) ) print_error("CPXchgbds");


	int t = 1; // iteration
	double t1 = second();
	double elapsed_time = 0;

	while ( remaining_time > 0 )
	{
		if (VERBOSE >= 60)
		{
			printf(" ... HARD FIXING 0 iteration %d ... \n", t+1);
			printf(" *** remaining time: %f\n", remaining_time);
		}

		// find the optimal solution
		if ( CPXmipopt(env, lp) ) print_error("CPXmipopt");
		elapsed_time = second() - t1;
		remaining_time -= elapsed_time;
		t1 = second();

		// get optimal solution found by cplex
		if ( CPXgetx(env, lp, inst->best_sol, 0, CPXgetnumcols(env, lp) - 1) ) 
		{
			if ( VERBOSE >= 60 )
				printf(" !!! Iteration %d didn't find an integer solution. Let's double the internal timelimit !!!\n", t+1);
			
			if ( remaining_time > 0)
			{
				if ( remaining_time <= tl )
				{
					printf(" !!! no solution found !!!\n");
					print_error("hard_fixing_0: try to increase tl\n");
				}
				else if ( remaining_time >= 2*tl )
					tl = 2 * tl;
				else
					tl = remaining_time;
			}
			CPXsetdblparam(env, CPX_PARAM_TILIM, tl);
			
			// unfix variables
			for ( int i = 0; i < inst->ncols; i++ ) 
				lu[i] = 'L';
			if ( CPXchgbds(env, lp, inst->ncols, index, lu, bound0) ) print_error("CPXchgbds");
			
			continue;
		}

		if ( CPXgetbestobjval(env, lp, &inst->best_lb) ) print_error("CPXgetbestobjval");
		if ( CPXgetobjval(env, lp, &inst->zbest) ) print_error("CPXgetobjval");

		if ( VERBOSE >= 60 )
			printf("--> iteration %d hard_fixing_0 with best_lb: %f, zbest: %f\n", t+1, inst->best_lb, inst->zbest);

		ratio = inst->best_lb / inst->zbest;
		if ( t == 0 && ratio > 1.0-EPS && ratio < 1.0+EPS )
		{
			if ( VERBOSE >= 60 )
				printf(" !!! Heuristic found optimal solution !!!\n");
			break;
		}
		
		// unfix variables
		for ( int i = 0; i < inst->ncols; i++ )
			lu[i] = 'L';
		if ( CPXchgbds(env, lp, inst->ncols, index, lu, bound0) ) print_error("CPXchgbds");
		
		//fix new variables
		for ( int i = 0; i < inst->ncols; i++ )  
		{
			if (inst->best_sol[i] > 0.5 ) 
			{
				p = (double) rand() / (RAND_MAX + 0.0);
				if ( p <= 0.9 )
					lu[i] = 'L';
				else
					lu[i] = 'U';
			}
			else
				lu[i] = 'U';
		}
		if ( CPXchgbds(env, lp, inst->ncols, index, lu, bound1) ) print_error("CPXchgbds");
		
		if ( VERBOSE >= 60 ) CPXwriteprob(env, lp, "model.lp", NULL);

		// last iteration
		if ( remaining_time < tl )
		{
			if ( remaining_time > 180 )
			{
				tl = remaining_time;
				CPXsetdblparam(env, CPX_PARAM_TILIM, tl);
			}
			else
				remaining_time = 0;
		}

		t++;
	}

	if ( VERBOSE >= 60 )
		printf("--> exiting hard_fixing_0 with zbest: %f\n", inst->zbest);

	free(index);
	free(value);
	free(bound0);
	free(bound1);
	free(lu);
	
	return 0;
}


/**************************************************************************************************************************/
int hard_fixing_1(instance *inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{
	double tl = inst->tl; 		// internal timelimit
	CPXsetdblparam(env, CPX_PARAM_TILIM, tl);
	double remaining_time = inst->timelimit;
	double ratio; 				// ratio best_lb/zbest
	double p; 					// probability to fix a variable
	int count = 0;				// non-improvement counter
	double last = CPX_INFBOUND;	// last zbest found
	double threshold = 0.9; 	// fixing threshold
	int* index = (int*)malloc(inst->ncols * sizeof(int));
	double* bound0 = (double*)malloc(inst->ncols * sizeof(double));
	double* bound1 = (double*)malloc(inst->ncols * sizeof(double));
	char* lu = (char*)malloc(inst->ncols * sizeof(double));

	// add a heuristic starting solution
	double* value = (double*)malloc(inst->ncols * sizeof(double));
	inst->timelimit = tl;
	vns(inst);

	int izero = 0;
	int nnz = 0;
	for (int i = 0; i < inst->ncols; i++) {
		if (inst->best_sol[i] > 0.5) {
			index[nnz]   = i;
			value[nnz++] = 1;
		}
	}

	int effortlevel = 5;
	CPXaddmipstarts(env, lp, 1, nnz, &izero, index, value, &effortlevel, NULL);
	remaining_time -= tl;
	
	for (int i = 0; i < inst->ncols; i++)
	{
		index[i] = i;
		bound0[i] = 0.0;
		bound1[i] = 1.0;
	}

	//fix new variables
	for ( int i = 0; i < inst->ncols; i++ )  
	{
		if (inst->best_sol[i] > 0.5 ) 
		{
			p = (double) rand() / (RAND_MAX + 0.0);
			if ( p <= 0.9 )
				lu[i] = 'L';
			else
				lu[i] = 'U';
		}
		else
			lu[i] = 'U';
	}
	if ( CPXchgbds(env, lp, inst->ncols, index, lu, bound1) ) print_error("CPXchgbds");


	int t = 1; // iteration
	double t1 = second();
	double elapsed_time = 0;

	while ( remaining_time > 0 )
	{
		if (VERBOSE >= 60)
		{
			printf(" ... HARD FIXING 1 iteration %d ... \n", t+1);
			printf(" *** remaining time: %f\n", remaining_time);
		}

		// find optimal solution
		if ( CPXmipopt(env, lp) ) print_error("CPXmipopt");
		elapsed_time = second() - t1;
		remaining_time -= elapsed_time;
		t1 = second();

		// get optimal solution found by cplex
		if ( CPXgetx(env, lp, inst->best_sol, 0, CPXgetnumcols(env, lp) - 1) ) 
		{
			if ( VERBOSE >= 60 )
				printf(" !!! Iteration %d didn't find an integer solution. Let's double the internal timelimit !!!\n", t+1);
			
			if ( remaining_time > 0)
			{
				if ( remaining_time <= tl )
				{
					printf(" !!! no solution found !!!\n");
					print_error(" !!! Try to increase tl !!!\n");
				}
				else if ( remaining_time >= 2*tl )
					tl = 2 * tl;
				else
					tl = remaining_time;
			}
			CPXsetdblparam(env, CPX_PARAM_TILIM, tl);
			
			// unfix variables
			for ( int i = 0; i < inst->ncols; i++ ) 
				lu[i] = 'L';
			if ( CPXchgbds(env, lp, inst->ncols, index, lu, bound0) ) print_error("CPXchgbds");
			
			continue;
		}

		if ( CPXgetbestobjval(env, lp, &inst->best_lb) ) print_error("CPXgetbestobjval");
		if ( CPXgetobjval(env, lp, &inst->zbest) ) print_error("CPXgetobjval");

		if ( VERBOSE >= 60 )
			printf("--> iteration %d hard_fixing_1 with best_lb: %f, zbest: %f\n", t+1, inst->best_lb, inst->zbest);

		ratio = inst->best_lb / inst->zbest;
		if ( (t == 0 || threshold < 0) && ratio > 1.0-EPS && ratio < 1.0+EPS )
		{
			if (VERBOSE >= 60)
				printf(" !!! Heuristic found optimal solution !!!\n");
			break;
		}

		if ( inst->zbest < last )
			count = 0;
		else
			count++;

		last = inst->zbest;

		if ( count == 5 )
		{
			threshold = threshold - 0.15;
			count = 0;
			if ( VERBOSE >= 60 ) printf(" *** threshold: %f\n", threshold);
		}
		
		// unfix variables
		for ( int i = 0; i < inst->ncols; i++ )
			lu[i] = 'L';
		if ( CPXchgbds(env, lp, inst->ncols, index, lu, bound0) ) print_error("CPXchgbds");
		
		// fix new variables
		for ( int i = 0; i < inst->ncols; i++ ) 
		{
			if ( inst->best_sol[i] > 0.5 ) 
			{
				p = (double) rand() / (RAND_MAX + 0.0);
				if ( p <= threshold ) 
					lu[i] = 'L';
				else
					lu[i] = 'U';
			}
			else
				lu[i] = 'U';
		}
		if ( CPXchgbds(env, lp, inst->ncols, index, lu, bound1) ) print_error("CPXchgbds");
		
		if ( VERBOSE >= 60 ) CPXwriteprob(env, lp, "model.lp", NULL);

		// last iteration
		if ( remaining_time < tl )
		{
			if ( remaining_time > 180 )
			{
				tl = remaining_time;
				CPXsetdblparam(env, CPX_PARAM_TILIM, tl);
			}
			else
				remaining_time = 0;
		}

		t++;
	}

	if ( VERBOSE >= 60 )
		printf("-->	exiting hard_fixing_1 with zbest: %f\n", inst->zbest);

	free(index);
	free(value);
	free(bound0);
	free(bound1);
	free(lu);
	
	return 0;
}


/**************************************************************************************************************************/
int hard_fixing_2(instance *inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{
	double tl = inst->tl; 	// internal timelimit
	CPXsetdblparam(env, CPX_PARAM_TILIM, tl);
	double remaining_time = inst->timelimit;
	double ratio; 			// ratio best_lb/zbest
	double p; 				// probability to fix a variable
	double threshold; 		// threshold for variable fixing
	int* index = (int*)malloc(inst->ncols * sizeof(int));
	double* bound0 = (double*)malloc(inst->ncols * sizeof(double));
	double* bound1 = (double*)malloc(inst->ncols * sizeof(double));
	char* lu = (char*)malloc(inst->ncols * sizeof(double));
	
	// add a heuristic starting solution
	double* value = (double*)malloc(inst->ncols * sizeof(double));
	inst->timelimit = tl;
	vns(inst);

	int izero = 0;
	int nnz = 0;
	for (int i = 0; i < inst->ncols; i++) {
		if (inst->best_sol[i] > 0.5) {
			index[nnz]   = i;
			value[nnz++] = 1;
		}
	}

	int effortlevel = 5;
	CPXaddmipstarts(env, lp, 1, nnz, &izero, index, value, &effortlevel, NULL);
	remaining_time -= tl;

	for (int i = 0; i < inst->ncols; i++)
	{
		index[i] = i;
		bound0[i] = 0.0;
		bound1[i] = 1.0;
	}

	//fix new variables
	for ( int i = 0; i < inst->ncols; i++ )  
	{
		if (inst->best_sol[i] > 0.5 ) 
		{
			p = (double) rand() / (RAND_MAX + 0.0);
			if ( p <= 0.9 )
				lu[i] = 'L';
			else
				lu[i] = 'U';
		}
		else
			lu[i] = 'U';
	}
	if ( CPXchgbds(env, lp, inst->ncols, index, lu, bound1) ) print_error("CPXchgbds");

	int t = 1; // iteration
	double t1 = second();
	double elapsed_time = 0;

	while ( remaining_time > 0 )
	{
		if (VERBOSE >= 60)
		{
			printf(" ... HARD FIXING 2 iteration %d ... \n", t+1);
			printf(" *** remaining time: %f\n", remaining_time);
		}

		// find optimal solution
		if ( CPXmipopt(env, lp) ) print_error("CPXmipopt");
		elapsed_time = second() - t1;
		remaining_time -= elapsed_time;
		t1 = second();

		// get optimal solution found by cplex
		if ( CPXgetx(env, lp, inst->best_sol, 0, CPXgetnumcols(env, lp) - 1) ) 
		{
			if ( VERBOSE >= 60 )
				printf(" !!! Iteration %d didn't find an integer solution. Let's double the internal timelimit !!!\n", t+1);
			
			if ( remaining_time > 0)
			{
				if ( remaining_time <= tl )
				{
					printf(" !!! no solution found !!!\n");
					print_error("hard_fixing_2: try to increase tl\n");
				}
				else if ( remaining_time >= 2*tl )
					tl = 2 * tl;
				else
					tl = remaining_time;
			}
			CPXsetdblparam(env, CPX_PARAM_TILIM, tl);

			// unfix variables
			for ( int i = 0; i < inst->ncols; i++ ) 
				lu[i] = 'L';
			if ( CPXchgbds(env, lp, inst->ncols, index, lu, bound0) ) print_error("CPXchgbds");
			
			continue;
		}

		if ( CPXgetbestobjval(env, lp, &inst->best_lb) ) print_error("CPXgetbestobjval");
		if ( CPXgetobjval(env, lp, &inst->zbest) ) print_error("CPXgetobjval");

		if ( VERBOSE >= 60 )
			printf("--> iteration %d hard_fixing_2 with best_lb: %f, zbest: %f\n", t+1, inst->best_lb, inst->zbest);

		ratio = inst->best_lb / inst->zbest;
		if ( (t == 0 || threshold < 0) && ratio > 1.0-EPS && ratio < 1.0+EPS )
		{
			if (VERBOSE >= 60)
				printf(" !!! Heuristic found optimal solution !!!\n");
			break;
		}
		
		// unfix variables
		for ( int i = 0; i < inst->ncols; i++ )
			lu[i] = 'L';
		if ( CPXchgbds(env, lp, inst->ncols, index, lu, bound0) ) print_error("CPXchgbds");
		
		// fix new variables
		for ( int i = 0; i < inst->ncols; i++ ) 
		{
			if ( inst->best_sol[i] > 0.5 ) 
			{
				p = (double) rand() / (RAND_MAX + 0.0);
				threshold = (double) rand() / (RAND_MAX + 0.0);
				if ( p <= threshold ) 
					lu[i] = 'L';
				else
					lu[i] = 'U';
			}
			else
				lu[i] = 'U';
		}
		if ( CPXchgbds(env, lp, inst->ncols, index, lu, bound1) ) print_error("CPXchgbds");
		
		if ( VERBOSE >= 60 ) CPXwriteprob(env, lp, "model.lp", NULL);

		// last iteration
		if ( remaining_time < tl )
		{
			if ( remaining_time > 180 )
			{
				tl = remaining_time;
				CPXsetdblparam(env, CPX_PARAM_TILIM, tl);
			}
			else
				remaining_time = 0;
		}

		t++;
	}

	if ( VERBOSE >= 60 )
		printf("--> exiting hard_fixing_2 with zbest: %f\n", inst->zbest);

	free(index);
	free(value);
	free(bound0);
	free(bound1);
	free(lu);
	
	return 0;
}


/**************************************************************************************************************************/
int local_branching(instance* inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{
	double tl = inst->tl;		//timelimit for each iteration
	CPXsetdblparam(env, CPX_PARAM_TILIM, tl);
	double remaining_time = inst->timelimit;
	double ratio; 				// ratio best_lb/zbest
	int count = 0; 				// non-improvement counter
	double last = CPX_INFBOUND;	// last zbest found
	int k = 2; 					// number of free edges
	int nnz = 0;
	int izero = 0;
	int* index = (int*)malloc(inst->ncols * sizeof(int));
	double* value = (double*)malloc(inst->ncols * sizeof(double));
	char** rname = (char**)calloc(1, sizeof(char*));
	rname[0] = (char*)calloc(100, sizeof(char));
	char sense = 'G';
	double rhs;
	int t = 0;					//number of iterations

	inst->timelimit = tl;
	vns(inst);

	for (int i = 0; i < inst->ncols; i++) {
		if (inst->best_sol[i] > 0.5) {
			index[nnz]   = i;
			value[nnz++] = 1;
		}
	}

	int effortlevel = 5;
	CPXaddmipstarts(env, lp, 1, nnz, &izero, index, value, &effortlevel, NULL);

	// add local branching
	nnz = 0;
	for (int i = 0; i < inst->ncols; i++)
	{
		if (inst->best_sol[i] > 0.5)
		{
			index[nnz] = i;
			value[nnz++] = 1.0;
		}
	}
	rhs = inst->nnodes - k;
	sprintf(rname[0], "localbranching(%d)", t);
	if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, rname)) print_error("wrong CPXaddrows() for local_branching()");

	if (VERBOSE >= 60)
	{
		printf(" ... LOCAL BRANCHING AFTER VNS ... \n");
		printf(" *** starting zbest: %10lf\n", inst->zbest);
	}

	remaining_time -= tl;

	t++;
	double t1 = second();
	double elapsed_time = 0;

	while (remaining_time > 0)
	{
		if (VERBOSE >= 60)
		{
			printf(" ... LOCAL BRANCHING iteration %d ... \n", t + 1);
			printf(" *** remaining time: %f\n", remaining_time);
			printf(" *** k: %d\n", k);
		}

		//find the optimal solution
		if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");
		elapsed_time = second() - t1;
		remaining_time -= elapsed_time;
		t1 = second();

		//take the optimal solution found by cplex
		if (CPXgetx(env, lp, inst->best_sol, 0, CPXgetnumcols(env, lp) - 1))
		{
			if (VERBOSE >= 60)
				printf(" !!! Iteration %d didn't find an integer solution. Let's double the internal timelimit !!!\n", t + 1);


			if (remaining_time <= tl)
			{
				printf(" !!! no solution found !!!\n");
				print_error(" !!! Try to increase tl !!!\n");
			}
			else if (remaining_time >= 2 * tl)
				tl = 2 * tl;
			else
				tl = remaining_time;
			CPXsetdblparam(env, CPX_PARAM_TILIM, tl);

			// remove local branching
			if (t > 0)
				if (CPXdelrows(env, lp, CPXgetnumrows(env, lp) - 1, CPXgetnumrows(env, lp) - 1))
					print_error("CPXdelrows() error");

			continue;
		}

		if (CPXgetbestobjval(env, lp, &inst->best_lb)) print_error("CPXgetbestobjval() error");
		if (CPXgetobjval(env, lp, &inst->zbest)) print_error("CPXgetobjval() error");

		if (VERBOSE >= 60)
			printf("--> iteration %d local_branching with best_lb: %f, zbest: %f\n", t + 1, inst->best_lb, inst->zbest);

		ratio = inst->best_lb / inst->zbest;
		if (t == 0 && ratio > 1.0 - EPS && ratio < 1.0 + EPS)
		{
			if (VERBOSE >= 60)
				printf(" !!! Heuristic found optimal solution !!!\n");
			break;
		}

		if (inst->zbest < last)
			count = 0;
		else
			count++;

		last = inst->zbest;

		if (count == 5)
		{
			if (k < 20)
				k++;
			count = 0;
		}

		// remove local branching
		if (t > 0)
			if (CPXdelrows(env, lp, CPXgetnumrows(env, lp) - 1, CPXgetnumrows(env, lp) - 1))
				print_error("CPXdelrows() error");

		// add local branching
		nnz = 0;
		for (int i = 0; i < inst->ncols; i++)
		{
			if (inst->best_sol[i] > 0.5)
			{
				index[nnz] = i;
				value[nnz++] = 1.0;
			}
		}
		rhs = inst->nnodes - k;
		sprintf(rname[0], "localbranching(%d)", t);
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, rname)) print_error("wrong CPXaddrows() for local_branching()");

		if (VERBOSE >= 60) CPXwriteprob(env, lp, "model.lp", NULL);

		// last iteration
		if (remaining_time < tl)
		{
			if (remaining_time > 180)
			{
				tl = remaining_time;
				CPXsetdblparam(env, CPX_PARAM_TILIM, tl);
			}
			else
				remaining_time = 0;
		}

		t++;
	}

	if (VERBOSE >= 60)
		printf("--> exiting local_branching with zbest: %f\n", inst->zbest);

	free(index);
	free(value);
	free(rname[0]);
	free(rname);

	return 0;
}


/**************************************************************************************************************************/
int nearest_neighbour_manager(instance* inst)
/**************************************************************************************************************************/
{
	double* last_sol = (double*)calloc(inst->ncols, sizeof(double));
	double zlast = CPX_INFBOUND;

	for (int i = 0; i < inst->nnodes; i++)
	{
		nearest_neighbour(inst, i, last_sol);
		zlast = cost(inst, last_sol);
		if (VERBOSE >= 60)
			printf("!!! start: %d, zlast: %f\n", i, zlast);

		if (zlast < inst->zbest)
		{
			inst->zbest = zlast;
			for (int i = 0; i < inst->ncols; i++)
				inst->best_sol[i] = last_sol[i];
		}

		// clean last_sol
		for (int i = 0; i < inst->ncols; i++)
			last_sol[i] = 0;
	}

	if (VERBOSE >= 60)
		printf("--> exiting nearest_neighbour_manager with zbest: %f\n", inst->zbest);

	free(last_sol);

	return 0;
}


/**************************************************************************************************************************/
int grasp_manager(instance* inst)
/**************************************************************************************************************************/
{
	double* last_sol = (double*)calloc(inst->ncols, sizeof(double));
	double zlast = CPX_INFBOUND;
	int start;
	double remaining_time = inst->timelimit;
	double elapsed_time = 0;
	double t1 = second();
	int t = 0;
	while (remaining_time > 0)
	{
		start = rand() % inst->nnodes;

		grasp(inst, start, last_sol);

		if (inst->heuristic == 8) {
			zlast = cost(inst, last_sol);
			if (VERBOSE >= 60) printf("!!! BEFORE 2-OPT MOVE t: %d, zlast: %f\n", t, zlast);
			two_opt_move(inst, last_sol);
		}

		zlast = cost(inst, last_sol);
		if (VERBOSE >= 60)
			printf("!!! t: %d, zlast: %f\n", t++, zlast);

		if (zlast < inst->zbest)
		{
			inst->zbest = zlast;
			for (int i = 0; i < inst->ncols; i++)
				inst->best_sol[i] = last_sol[i];
		}

		// clean last_sol
		for (int i = 0; i < inst->ncols; i++)
			last_sol[i] = 0;

		elapsed_time = second() - t1;
		remaining_time -= elapsed_time;
		t1 = second();
	}

	if (VERBOSE >= 60)
		printf("--> exiting grasp_manager with zbest: %f\n", inst->zbest);

	free(last_sol);

	return 0;
}


/**************************************************************************************************************************/
int nearest_neighbour(instance* inst, int start, double* xstar)
/**************************************************************************************************************************/
{
	int* visited = (int*)malloc(inst->nnodes * sizeof(int));
	for (int i = 0; i < inst->nnodes; i++)
		visited[i] = -1;
	int nvisited = 1;
	visited[start] = 1;
	int prev = start;
	int succ;

	while (nvisited < inst->nnodes)
	{
		succ = closest_node(inst, prev, visited);
		visited[succ] = 1;
		nvisited++;
		xstar[xpos(prev, succ, inst)] = 1;
		prev = succ;
	}
	xstar[xpos(prev, start, inst)] = 1;

	free(visited);

	return 0;
}


/**************************************************************************************************************************/
int grasp(instance* inst, int start, double* xstar)
/**************************************************************************************************************************/

{
	int* visited = (int*)malloc(inst->nnodes * sizeof(int));
	for (int i = 0; i < inst->nnodes; i++)
		visited[i] = -1;
	int nvisited = 1;
	visited[start] = 1;
	int prev = start;
	int succ;

	while (nvisited < inst->nnodes)
	{
		succ = rand_closest_node(inst, prev, visited);
		visited[succ] = 1;
		nvisited++;
		xstar[xpos(prev, succ, inst)] = 1;
		prev = succ;
	}
	xstar[xpos(prev, start, inst)] = 1;

	free(visited);

	return 0;
}


/**************************************************************************************************************************/
int closest_node(instance* inst, int node, int* visited)
/**************************************************************************************************************************/
{
	int min_idx = -1;
	double min_dist = CPX_INFBOUND;
	double d;
	for (int i = 0; i < inst->nnodes; i++)
	{
		if (visited[i] == 1 || i == node)
			continue;
		d = dist(i, node, inst);
		if (d < min_dist)
		{
			min_dist = d;
			min_idx = i;
		}
	}

	return min_idx;
}


/**************************************************************************************************************************/
int rand_closest_node(instance* inst, int node, int* visited)
/**************************************************************************************************************************/
{
	int min_idx[3] = { -1, -1, -1 };
	int k;

	for (k = 0; k < 3; k++)
	{
		min_idx[k] = closest_node(inst, node, visited);
		if (min_idx[k] == -1)
			break;
		visited[min_idx[k]] = 1;
	}

	int rnd_idx = rand() % k;

	for (k = 0; k < 3; k++)
	{
		if (min_idx[k] != -1)
			visited[min_idx[k]] = -1;
	}

	return min_idx[rnd_idx];
}


/**************************************************************************************************************************/
int insertion(instance* inst, int node1, int node2)
/**************************************************************************************************************************/
{
	double* last_sol = (double*)calloc(inst->ncols, sizeof(double));
	int* visited = (int*)malloc(inst->nnodes * sizeof(int));
	for (int i = 0; i < inst->nnodes; i++)
		visited[i] = -1;
	visited[node1] = 1;
	visited[node2] = 1;
	last_sol[xpos(node1, node2, inst)] = 1;
	int min_hab[3];
	get_hab(inst, visited, last_sol, min_hab);
	last_sol[xpos(node1, min_hab[0], inst)] = 1;
	last_sol[xpos(min_hab[0], node2, inst)] = 1;
	visited[min_hab[0]] = 1;
	int nvisited = 3;

	while (nvisited < inst->nnodes)
	{
		get_hab(inst, visited, last_sol, min_hab);
		last_sol[xpos(min_hab[1], min_hab[2], inst)] = 0;
		last_sol[xpos(min_hab[1], min_hab[0], inst)] = 1;
		last_sol[xpos(min_hab[0], min_hab[2], inst)] = 1;
		visited[min_hab[0]] = 1;
		nvisited++;
	}

	for (int i = 0; i < inst->ncols; i++)
		inst->best_sol[i] = last_sol[i];
	inst->zbest = cost(inst, last_sol);

	if (VERBOSE >= 60)
		printf("--> exiting insertion with zbest: %f\n", inst->zbest);

	free(visited);
	free(last_sol);

	return 0;
}


/**************************************************************************************************************************/
void get_hab(instance* inst, int* visited, double* xstar, int* min_hab)
/**************************************************************************************************************************/
{
	min_hab[0] = -1;
	min_hab[1] = -1;
	min_hab[2] = -1;
	double min_delta = CPX_INFBOUND;
	double delta;

	for (int h = 0; h < inst->nnodes; h++)
	{
		// for each h not visited yet
		if (visited[h] == 1) continue;

		for (int a = 0; a < inst->nnodes; a++)
		{
			for (int b = a + 1; b < inst->nnodes; b++)
			{
				// for each edge (a,b) within existing cycle
				if (xstar[xpos(a, b, inst)] > 0.5)
				{
					delta = dist(a, h, inst) + dist(h, b, inst) - dist(a, b, inst);
					if (delta < min_delta)
					{
						min_delta = delta;
						min_hab[0] = h;
						min_hab[1] = a;
						min_hab[2] = b;
					}
				}
			}
		}
	}
}


/**************************************************************************************************************************/
int rand_insertion_manager(instance* inst)
/**************************************************************************************************************************/
{
	double* last_sol = (double*)calloc(inst->ncols, sizeof(double)); // devo usare calloc perchè nearest neighbour mette solo gli 1 e non gli zeri!
	double zlast = CPX_INFBOUND;
	double remaining_time = inst->timelimit;
	double elapsed_time = 0;
	double t1 = second();
	int t = 0;

	while (remaining_time > 0)
	{
		rand_insertion(inst, 0, 1, last_sol);
		zlast = cost(inst, last_sol);
		if (VERBOSE >= 60)
			printf("!!! t: %d, zlast: %f\n", t++, zlast);

		if (zlast < inst->zbest)
		{
			inst->zbest = zlast;
			for (int i = 0; i < inst->ncols; i++)
				inst->best_sol[i] = last_sol[i];
		}

		// clean last_sol
		for (int i = 0; i < inst->ncols; i++)
			last_sol[i] = 0;

		elapsed_time = second() - t1;
		remaining_time -= elapsed_time;
		t1 = second();
	}

	if (VERBOSE >= 60)
		printf("--> exiting rand_insertion_manager with zbest: %f\n", inst->zbest);

	free(last_sol);

	return 0;
}


/**************************************************************************************************************************/
int rand_insertion(instance* inst, int node1, int node2, double* xstar)
/**************************************************************************************************************************/
{
	int* visited = (int*)malloc(inst->nnodes * sizeof(int));
	for (int i = 0; i < inst->nnodes; i++)
		visited[i] = -1;
	visited[node1] = 1;
	visited[node2] = 1;
	xstar[xpos(node1, node2, inst)] = 1;
	int min_hab[3];

	rand_get_hab(inst, visited, xstar, min_hab);
	xstar[xpos(node1, min_hab[0], inst)] = 1;
	xstar[xpos(min_hab[0], node2, inst)] = 1;
	visited[min_hab[0]] = 1;
	int nvisited = 3;

	while (nvisited < inst->nnodes)
	{
		rand_get_hab(inst, visited, xstar, min_hab);
		xstar[xpos(min_hab[1], min_hab[2], inst)] = 0;
		xstar[xpos(min_hab[1], min_hab[0], inst)] = 1;
		xstar[xpos(min_hab[0], min_hab[2], inst)] = 1;
		visited[min_hab[0]] = 1;
		nvisited++;
	}

	free(visited);

	return 0;
}


/**************************************************************************************************************************/
void rand_get_hab(instance* inst, int* visited, double* xstar, int* min_hab)
/**************************************************************************************************************************/
{
	min_hab[0] = -1;
	min_hab[1] = -1;
	min_hab[2] = -1;

	int** hab = (int**)malloc(3 * sizeof(int*));
	hab[0] = (int*)malloc(3 * sizeof(int));
	hab[1] = (int*)malloc(3 * sizeof(int));
	hab[2] = (int*)malloc(3 * sizeof(int));
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			hab[i][j] = -1;

	int k;

	for (k = 0; k < 3; k++)
	{
		get_hab(inst, visited, xstar, hab[k]);
		if (hab[k][0] == -1)
			break;
		visited[hab[k][0]] = 1;
	}

	int rnd_idx = rand() % k;

	for (k = 0; k < 3; k++)
	{
		if (hab[k][0] != -1)
			visited[hab[k][0]] = -1;
	}

	min_hab[0] = hab[rnd_idx][0];
	min_hab[1] = hab[rnd_idx][1];
	min_hab[2] = hab[rnd_idx][2];

	for (int i = 0; i < 3; i++)
		free(hab[i]);
	free(hab);
}


/**************************************************************************************************************************/
int two_opt_move(instance * inst, double *xstar)
/**************************************************************************************************************************/
{
	int a = -1, b = -1, c = -1, d = -1;
	int A = -1, B = -1, C = -1, D = -1;
	int e1 = -1, e2 = -1, e3 = -1, e4 = -1;
	double min_delta, delta;

	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	int* comp = (int*)malloc(inst->nnodes * sizeof(int));
	int ncomp = 0;
	build_sol(xstar, inst, succ, comp, &ncomp);
	if (ncomp != 1) print_error("solution in two_opt_move has subtour");
	
	int t = 0;
	while (1)
	{
		min_delta = CPX_INFBOUND;
		if (VERBOSE >= 60) printf(" ... iteration t=%d  ... \n", t++);
		for ( a = 0; a < inst->nnodes; a++) {
			b = succ[a];
			for (c = a + 1; c < inst->nnodes; c++) {
				if (c == b) continue;

				d = succ[c];
				if (d == a) continue;

				delta = (dist(a, c, inst) + dist(b, d, inst)) - (dist(a, b, inst) + dist(c, d, inst));

				if (delta < 0 && delta < min_delta)
				{
					if (VERBOSE >= 60) printf(" ... delta=%lf  min_delta=%lf  ... \n", delta, min_delta);
					min_delta = delta;
					e1 = xpos(a, b, inst);
					e2 = xpos(c, d, inst);
					e3 = xpos(a, c, inst);
					e4 = xpos(b, d, inst);
					A = a;
					B = b;
					C = c;
					D = d;
				}
			}
		}

		if (min_delta < 0) {
			xstar[e1] = 0;
			xstar[e2] = 0;
			xstar[e3] = 1;
			xstar[e4] = 1;
		}
		else
			break;

		if (VERBOSE >= 60) printf(" ... min_delta=%lf  ... \n", min_delta);

		int prev = B;
		int curr = succ[B];
		int next;

		while (curr != C) 
		{
			next = succ[curr];
			succ[curr] = prev;
			prev = curr;
			curr = next;
		}

		succ[A] = C;
		succ[C] = prev;
		succ[B] = D;

	}

	free(succ);
	free(comp);

	return 0;
}


/**************************************************************************************************************************/
int vns(instance* inst)
/**************************************************************************************************************************/
{
	double* xstar = (double*)calloc(inst->ncols, sizeof(double));
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	int* comp = (int*)malloc(inst->nnodes * sizeof(int));
	int ncomp = 0;
	double zlast = CPX_INFBOUND;
	inst->zbest = CPX_INFBOUND;
	int k = 5;

	int* node1 = (int*)malloc(k * sizeof(int)); //prevs
	int* node2 = (int*)malloc(k * sizeof(int)); //succs
	int* linkedList = (int*)malloc(2 * k * sizeof(int)); //sorted list of the nodes
	//int node1[k]; // prevs
	//int node2[k]; // succs


	double remaining_time = inst->timelimit;
	double elapsed_time = 0;
	double t1 = second();

	nearest_neighbour(inst, 0, xstar);
	//inst->nplot = 0;

	while (remaining_time > 0)
	{
		// intensification
		two_opt_move(inst, xstar);
		zlast = cost(inst, xstar);

		// save best solution
		if (zlast < inst->zbest)
		{
			inst->zbest = zlast;
			for (int i = 0; i < inst->ncols; i++)
				inst->best_sol[i] = xstar[i];

			//write_sol_name(inst, xstar, inst->nplot++);
		}

		//write_sol_name(inst, xstar, inst->nplot++);

		// diversification
		// randomly choose k edges
		build_sol(xstar, inst, succ, comp, &ncomp);
		int done = 0;
		int extraction;

		while (done < k)
		{
			extraction = rand() % inst->nnodes;
			node1[done] = extraction;
			node2[done] = succ[extraction];
			done++;
			for (int i = 0; i < done; i++)
			{
				if (((done - 1) > i) &&
					(node1[done - 1] == node1[i] || node1[done - 1] == node2[i] || node2[done - 1] == node1[i] || node2[done - 1] == node2[i]))
				{
					done--;
					continue;
				}
			}
		}

		if (VERBOSE >= 90) {
			printf("node1: ");
			for (int i = 0; i < k; i++) {
				printf("%d, ", node1[i]);
			}
			printf("\nnode2: ");
			for (int i = 0; i < k; i++) {
				printf("%d, ", node2[i]);
			}
			printf("\n");
		}

		linkedList[0] = node1[0];
		for (int i = 1; i < k; i++) {
			linkedList[2 * i - 1] = succ[linkedList[2 * i - 2]];
			linkedList[2 * i] = findNext(linkedList[2 * i - 1], succ, node1, k);
			if(VERBOSE >= 90) printf(" ...... primo: %d   secondo: %d \n", linkedList[2 * i - 1], linkedList[2 * i]);
		}
		linkedList[2 * k - 1] = succ[linkedList[2 * k - 2]];

		if (VERBOSE >= 90) printf("\nLinked List: ");
		for (int i = 0; i < 2*k; i++) {
			if (VERBOSE >= 90) printf("%d, ", linkedList[i]);
		}
		if (VERBOSE >= 90) printf("\n");

		for (int e = 1; e < k-1; e++)
		{
			if (VERBOSE >= 90)  printf("inserisco {%d,%d}\n", linkedList[2 * e - 1], linkedList[2 * (e + 1)]);
			xstar[xpos(linkedList[2 * e - 1], linkedList[2 * (e + 1)], inst)] = 1;
		}
		xstar[xpos(linkedList[0], linkedList[2 * k - 3], inst)] = 1;
		xstar[xpos(linkedList[2 * k - 1], linkedList[2], inst)] = 1;
		if (VERBOSE >= 90) printf("inserisco {%d,%d}\n", linkedList[0], linkedList[2 * k - 3]);
		if (VERBOSE >= 90) printf("inserisco {%d,%d}\n", linkedList[2 * k - 1], linkedList[2]);

		for (int e = 0; e < k; e++) {
			xstar[xpos(node1[e], node2[e], inst)] = 0;
			if (VERBOSE >= 90) printf("elimino {%d,%d}\n", node1[e], node2[e]);
		}
		
		//write_sol_name(inst, xstar, inst->nplot++);

		elapsed_time = second() - t1;
		remaining_time -= elapsed_time;
		t1 = second();
		if (VERBOSE >= 60) printf("elapsed time: %f, rem time: %f\n", elapsed_time, remaining_time);
	}

	if (VERBOSE >= 60)
		printf("--> exiting vns with zbest: %f\n", inst->zbest);

	free(linkedList);
	free(node1);
	free(node2);
	free(xstar);

	return 0;
}


/**************************************************************************************************************************/
int tabu(instance* inst)
/**************************************************************************************************************************/
{
	double* xstar = (double*)calloc(inst->ncols, sizeof(double));
	double zlast = CPX_INFBOUND;
	inst->zbest = CPX_INFBOUND;

	int tenure = inst->nnodes/20;
	int* tabu_list = (int*)malloc(tenure * sizeof(int));
	for (int i = 0; i < tenure; i++)
		tabu_list[i] = -1;
	int rear = 0;

	int a = -1, b = -1, c = -1, d = -1;
	int A = -1, B = -1, C = -1, D = -1;
	int e1 = -1, e2 = -1, e3 = -1, e4 = -1;
	double min_delta, delta;
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	int* comp = (int*)malloc(inst->nnodes * sizeof(int));
	int ncomp = 0;

	double remaining_time = inst->timelimit;
	double elapsed_time = 0;
	double t1 = second();

	nearest_neighbour(inst, 0, xstar);
	build_sol(xstar, inst, succ, comp, &ncomp);
	if (ncomp != 1) print_error("tabu: starting solution has subtour");
	int t = 0;
	while (remaining_time > 0)
	{
		min_delta = CPX_INFBOUND;
		if (VERBOSE >= 60) printf(" ... iteration t=%d  ... \n", t++);
		for (a = 0; a < inst->nnodes; a++) {
			b = succ[a];
			for (c = a + 1; c < inst->nnodes; c++) {
				if (c == b) continue;

				d = succ[c];
				if (d == a) continue;

				delta = (dist(a, c, inst) + dist(b, d, inst)) - (dist(a, b, inst) + dist(c, d, inst));

				if (delta < min_delta)
				{
					if (is_in_tabu_list(xpos(a, b, inst), xpos(c, d, inst), tabu_list, tenure) > 0)
						continue;

					if (VERBOSE >= 60) printf(" ... delta=%lf  min_delta=%lf  ... \n", delta, min_delta);
					min_delta = delta;
					e1 = xpos(a, b, inst);
					e2 = xpos(c, d, inst);
					e3 = xpos(a, c, inst);
					e4 = xpos(b, d, inst);
					A = a;
					B = b;
					C = c;
					D = d;
				}
			}
		}

		xstar[e1] = 0;
		xstar[e2] = 0;
		xstar[e3] = 1;
		xstar[e4] = 1;

		if (min_delta < 0) 
		{
			zlast = cost(inst, xstar);
			if (zlast < inst->zbest)
			{
				inst->zbest = zlast;
				for (int i = 0; i < inst->ncols; i++)
					inst->best_sol[i] = xstar[i];
			}
		}
		else 
		{
			if ((double) rand()/RAND_MAX > 0.5)
				tabu_list[rear] = e3;
			else
				tabu_list[rear] = e4;
			rear = (rear + 1)%tenure;
		}

		if (VERBOSE >= 60) printf(" ... min_delta=%lf  ... \n", min_delta);

		int prev = B;
		int curr = succ[B];
		int next;

		while (curr != C) 
		{
			next = succ[curr];
			succ[curr] = prev;
			prev = curr;
			curr = next;
		}

		succ[A] = C;
		succ[C] = prev;
		succ[B] = D;


		elapsed_time = second() - t1;
		remaining_time -= elapsed_time;
		t1 = second();
		if (VERBOSE >= 60) printf("elapsed time: %f, rem time: %f\n", elapsed_time, remaining_time);
	}

	if (VERBOSE >= 60)
		printf("--> exiting tabu with zbest: %f\n", inst->zbest);

	free(xstar);
	free(succ);
	free(comp);
	free(tabu_list);

	return 0;
}


/**************************************************************************************************************************/
int simulated_annealing(instance* inst)
/**************************************************************************************************************************/
{
	double* xstar = (double*)calloc(inst->ncols, sizeof(double));
	double zlast = CPX_INFBOUND;
	inst->zbest = CPX_INFBOUND;

	int N = inst->timelimit; // outer loop
	int M = 50; // inner loop

	double p_start = 0.7;
	double p_end = 0.0001;

	double t_start = -1 / (log(p_start));
	double t_end = -1 / (100*log(p_end));
	if (VERBOSE >= 60)
		printf("!!!!!!!!!!!!!!!!!!!!!!tstart: %lf    tend: %lf\n", t_start, t_end);
	double frac = pow(t_end / t_start, 1 / (N - 1.0));
	if (VERBOSE >= 60)
		printf("!!!!!!!!!!!!!!!!!!!!!!frac: %lf\n", frac);

	nearest_neighbour(inst, 0, xstar);

	int* best_succ = (int*)malloc(inst->nnodes * sizeof(int));
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	int* comp = (int*)malloc(inst->nnodes * sizeof(int));
	int ncomp = 0;
	build_sol(xstar, inst, best_succ, comp, &ncomp);
	for (int i = 0; i < inst->nnodes; i++)
		succ[i] = best_succ[i];
	for (int i = 0; i < inst->ncols; i++)
		inst->best_sol[i] = xstar[i];
	zlast = cost(inst, xstar);
	inst->zbest = zlast;

	double t = t_start;
	double delta;
	double delta_avg = 0;
	double p;
	double coin;

	int a = -1, b = -1, c = -1, d = -1;
	int accept = 0;

	for (int n = 0; n < N; n++)
	{
		for (int m = 0; m < M; m++)
		{
			a = rand() % inst->nnodes;
			b = succ[a];
			c = rand() % inst->nnodes;
			d = succ[c];
			while (c == a || c == b || d == a || d == b)
			{
				c = rand() % inst->nnodes;
				d = succ[c];
			}

			delta = (dist(a, c, inst) + dist(b, d, inst)) - (dist(a, b, inst) + dist(c, d, inst));
			
			delta_avg = (delta_avg*(n+m) + delta)/(n+m+1);

			if (VERBOSE >= 61)
				printf(" \t \t delta: %10.2lf   delta_avg: %10.2lf\n", delta, delta_avg);

			if (delta < 0)
				accept = 1;
			else
			{
				p = exp(-delta / (delta_avg * t));
				coin = (double)rand() / (RAND_MAX+0.0);
				if (coin <= p)
					accept = 1;
				else
					accept = 0;
			}

			if (accept == 1)
			{
				xstar[xpos(a, b, inst)] = 0;
				xstar[xpos(c, d, inst)] = 0;
				xstar[xpos(a, c, inst)] = 1;
				xstar[xpos(b, d, inst)] = 1;

				int curr, next, prev;
				curr = succ[b];
				next = succ[curr];
				prev = b;

				while (curr != c) {
					succ[curr] = prev;
					prev = curr;
					curr = next;
					next = succ[next];

				}
				succ[c] = prev;
				succ[a] = c;
				succ[b] = d;

				zlast = zlast + delta;
				if (VERBOSE >= 61)
					printf(" \t \t \t zlast: %10.2lf   p: %1.6lf\n", zlast,p);
				if (zlast < inst->zbest)
				{
					for (int i = 0; i < inst->nnodes; i++) 
						best_succ[i] = succ[i];

					inst->zbest = zlast;
				}
			}
		}

		if (VERBOSE >= 61)
			printf(" !!! n: %5d  t: %10.6lf    zbest: %lf   frac: %lf\n", n, t, inst->zbest, frac);
		// decrease temperature
		t = frac * t;
	}

	if (VERBOSE >= 60)
		printf("--> exiting simulated annealing with zbest: %f\n", inst->zbest);

	for (int i = 0; i < inst->ncols; i++)
		inst->best_sol[i] = 0;
	int next = succ[0];
	int prev = 0;
	while (next != 0) {
		inst->best_sol[xpos(prev, next, inst)] = 1;
		prev = next;
		next = succ[next];
	}
	inst->best_sol[xpos(prev, next, inst)] = 1;

	free(succ);
	free(comp);
	free(xstar);

	return 0;
}














