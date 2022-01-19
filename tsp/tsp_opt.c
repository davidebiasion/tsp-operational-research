#include "tsp.h"
#include <time.h>


static int CPXPUBLIC mylazycallback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p);
static int CPXPUBLIC generic_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle);
static int CPXPUBLIC myheuristicfunc(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, double* objval_p, double* x, int* checkfeas_p, int* useraction_p);


/**************************************************************************************************************************/
int TSPopt(instance *inst)
/**************************************************************************************************************************/
{ 
	// open cplex model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP"); 
	if ( VERBOSE >= 60 ) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // cplex output on screen

	// cplex parameters
	CPXsetintparam(env, CPXPARAM_RandomSeed, inst->randomseed); // Cplex random seed
	if ( inst->timelimit < CPX_INFBOUND ) CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit);

	// build model
	build_model(inst, env, lp);

	inst->ncols = CPXgetnumcols(env, lp);
	inst->best_sol = (double*)calloc(inst->ncols, sizeof(double)); 	// all entries to zero  
	inst->zbest = CPX_INFBOUND;
	
	// compute solution
	if ( CPXmipopt(env,lp) ) print_error("CPXmipopt");

	// use solution
	if ( CPXgetx(env, lp, inst->best_sol, 0, inst->ncols-1) ) print_error("CPXgetx");

	// find succ, comp and ncomp
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	int* comp = (int*)malloc(inst->nnodes * sizeof(int));
	int ncomp;
	build_sol(inst->best_sol, inst, succ, comp, &ncomp);

	if (ncomp > 1)
	{
		if ( inst->model_type != 0 ) print_error("calling addSEC with model_type != 0");
		if ( addSEC(inst, env, lp, succ, comp, &ncomp) ) print_error("addSEC");
	}

	// write .dat file for gnuplot
	write_sol(inst);

	// print solution on screen
	if ( VERBOSE >= 70 )
	{
		int x = 0;
		for (int i = 0; i < inst->nnodes; i++)
		{
			printf("\nEdges for node %d: ", x + 1);
			printf("[%d, %d] ", x + 1, succ[x] + 1);
			x = succ[x];
		}
		printf("\nThe number of connected components is: %2d\n", ncomp);
	}
	
	free(succ);
	free(comp);

	// free pools and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 	
	
	return 0;
}


/***************************************************************************************************************************/
void build_model(instance* inst, CPXENVptr env, CPXLPptr lp) 
/***************************************************************************************************************************/
{
	//initialization of the model
	inst->xstart = -1;
	inst->ustart = -1;
	inst->bigqstart = -1;
	inst->sstart = -1;
	inst->bigsstart = -1;
	inst->ystart = -1;
	inst->fstart = -1;
	inst->zstart = -1;

	int ncores = 1;
	int ncols;

	//build the model depending on the model desired
	switch (inst->model_type) 
	{
		case 0: // Loop
			build_model_0(inst, env, lp);
			break;

		case 1: // MTZ
			build_model_1(inst, env, lp);
			break;

		case 2: // F1
			build_model_2(inst, env, lp);
			break;

		case 3:	//callback
			build_model_0(inst, env, lp);

			//install the lazy callback
			CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);	// let MIP callbacks work on the original model
			if (CPXsetlazyconstraintcallbackfunc(env, mylazycallback, inst)) print_error("CPXsetlazyconstraintcallbackfunc");
			CPXgetnumcores(env, &ncores);
			CPXsetintparam(env, CPX_PARAM_THREADS, ncores); 	// it was reset after callback
			inst->nthread = ncores;
		
			break;

		case 30:	// lazy and heuristic callbacks
			build_model_0(inst, env, lp);

			//install the lazy callback
			CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);	// let MIP callbacks work on the original model
			if (CPXsetlazyconstraintcallbackfunc(env, mylazycallback, inst)) print_error("CPXsetlazyconstraintcallbackfunc");
			if (CPXsetheuristiccallbackfunc(env, myheuristicfunc, inst)) print_error("CPXsetheuristiccallbackfunc");
				CPXgetnumcores(env, &ncores);
				CPXsetintparam(env, CPX_PARAM_THREADS, ncores); 	// it was reset after callback
				inst->nthread = ncores;

				//static double patch_vec[nthread][inst->ncols];
				inst->patch_vec = (double**)malloc(inst->nthread * sizeof(double*));
				inst->patch_cost = (double*)malloc(inst->nthread * sizeof(double));
				ncols = CPXgetnumcols(env, lp);
				for (int i = 0; i < inst->nthread; i++) {
					inst->patch_vec[i] = (double*)malloc(ncols * sizeof(double));
						inst->patch_cost[i] = -1;
				}

			break;

		case 31:	// generic lazy callback
			build_model_0(inst, env, lp);

			// install the lazy callback
			CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);	// let MIP callbacks work on the original model
			if (CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE, generic_callback, inst)) print_error("CPXcallbacksetfunc");
			/*int ncores = 1; */ CPXgetnumcores(env, &ncores);
			CPXsetintparam(env, CPX_PARAM_THREADS, ncores); 	// it was reset after callback

			break;

		case 310:	// generic lazy+heuristic callback
			build_model_0(inst, env, lp);

			// install the lazy+heuristic callback
			CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);	// let MIP callbacks work on the original model
			if (CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_GLOBAL_PROGRESS, generic_callback, inst)) print_error("CPXcallbacksetfunc");
			/*int ncores = 1; */ CPXgetnumcores(env, &ncores);
			CPXsetintparam(env, CPX_PARAM_THREADS, ncores); 	// it was reset after callback

			inst->nthread = ncores;

			//static double patch_vec[nthread][inst->ncols];
			inst->patch_vec = (double**)malloc(inst->nthread * sizeof(double*));
			inst->patch_cost = (double*)malloc(inst->nthread * sizeof(double));
			ncols = CPXgetnumcols(env, lp);
			for (int i = 0; i < inst->nthread; i++) {
				inst->patch_vec[i] = (double*)malloc(ncols * sizeof(double));
				inst->patch_cost[i] = -1;
			}

			break;

		case 11: // MTZ with lazy constraints
			build_model_11(inst, env, lp);
			break;

		case 21: // F1 with n-2
			build_model_21(inst, env, lp);
			break;

		case 22: // F1 with n-2 and lazy constraints
			build_model_22(inst, env, lp);
			break;

		default:
			print_error("model type unknown");
	}
}


/***************************************************************************************************************************/
void build_model_0(instance *inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{   
	char binary = 'B'; 
	char **cname = (char **) calloc(1, sizeof(char *));		// (char **) required by cplex...
	cname[0] = (char *) calloc(100, sizeof(char));
    int izero = 0;
    int* index = (int*)malloc(inst->nnodes * sizeof(int));			
	double* value = (double*)malloc(inst->nnodes * sizeof(double));
    int nnz = 0;
	double rhs = 2.0; 
	char sense = 'E';
	double lb = 0.0;
	double ub = 1.0;
	double obj;

	// add binary var.s x(i,j) for i < j  
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = i+1; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);
			obj = dist(i, j, inst); // cost == distance   
			if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error("wrong CPXnewcols on x var.s");
    		if ( CPXgetnumcols(env,lp)-1 != xpos(i, j, inst) ) print_error("wrong position for x var.s");
		}
	} 


	// add the degree constraints
	for ( int h = 0; h < inst->nnodes; h++ )  // degree constraints
	{
		nnz = 0;
		sprintf(cname[0], "degree(%d)", h+1);   
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			if ( i == h ) continue;
			index[nnz] = xpos(i, h, inst);
			value[nnz++] = 1.0;
		}
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
			print_error("wrong CPXaddrows [degree]");
	}

	if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model.lp", NULL);   

	free(index);
	free(value);
	free(cname[0]);
	free(cname);
}  


// MTZ
/***************************************************************************************************************************/
void build_model_1(instance *inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{ 
	double zero = 0.0;  
	char binary = 'B';
	char continuous = 'C';
	//char integer = 'I';
	char **cname = (char **) calloc(1, sizeof(char *));
	cname[0] = (char *) calloc(100, sizeof(char));
    int izero = 0;
    int* index = (int*)malloc(inst->nnodes * sizeof(int));			
	double* value = (double*)malloc(inst->nnodes * sizeof(double));
    int nnz = 0;
	double rhs = 1.0; 
	char sense = 'E';
	double obj;
	double lb;
	double ub;

	// add binary var.s x(i,j)   
	if ( inst->xstart != -1 ) print_error("build_model: var. x cannot be redefined!");
	inst->xstart = CPXgetnumcols(env,lp); 		// position of the first x(,) variable  
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = 0; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);
			obj = dist(i, j, inst);   
			ub = ( i == j ) ? 0.0 : 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &zero, &ub, &binary, cname) ) print_error("wrong CPXnewcols on x var.s");
    		if ( CPXgetnumcols(env,lp)-1 != xpos_1(i, j, inst) ) print_error("wrong position for x var.s");
		}
	}

	// add continuous var.s u(i) 
	obj = 0.0;
	lb = 0.0;  
	if ( inst->ustart != -1 ) print_error("build_model: var. u cannot be redefined!");
	inst->ustart = CPXgetnumcols(env,lp); 		// position of the first u() variable   
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		sprintf(cname[0], "u(%d)", i+1);
		ub = ( i == 0 ) ? 0.0 : inst->nnodes-2;
		if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &continuous, cname) ) print_error("wrong CPXnewcols on u var.s");
    	if ( CPXgetnumcols(env,lp)-1 != upos(i, inst) ) print_error("wrong position for u var.s");
	}

	// x-constraints (out- and in-degree at nodes)
	for ( int h = 0; h < inst->nnodes; h++ )  // out-degree
	{
		nnz = 0;
		sprintf(cname[0], "outdeg(%d)", h+1);   
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			index[nnz] = xpos_1(h, i, inst);
			value[nnz++] = 1.0;
		}
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
			print_error("wrong CPXaddrows [x1]");
	} 
	for ( int h = 0; h < inst->nnodes; h++ ) // in-degree
	{
		nnz = 0;
		sprintf(cname[0], "indeg(%d)", h+1);
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			index[nnz] = xpos_1(i, h, inst);
			value[nnz++] = 1.0;
		}
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
			print_error("wrong CPXaddrows [x2]");
	} 

	// u-x linking 
	nnz = 3;
	double M = inst->nnodes;
	rhs = 1-M;
	sense = 'G';
	for ( int i = 1; i < inst->nnodes; i++ )
	{
		for ( int j = 1; j < inst->nnodes; j++ )    
		{
			if ( i == j ) continue;
			sprintf(cname[0], "link(%d,%d)", i+1,j+1);
			index[0] = upos(i, inst);
			value[0] = -1.0;
			index[1] = upos(j, inst);
			value[1] = 1.0;
			index[2] = xpos_1(i, j, inst);
			value[2] = -M;
			if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
				print_error("wrong CPXaddrows [link]");
		}
	}

	if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model.lp", NULL);   

	free(index);
	free(value);
	free(cname[0]);
	free(cname);
}


// MTZ (with lazy constraints)
/***************************************************************************************************************************/
void build_model_11(instance *inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{   
	double zero = 0.0;  
	char binary = 'B';
	char continuous = 'C';
	//char integer = 'I';
	char **cname = (char **) calloc(1, sizeof(char *));
	cname[0] = (char *) calloc(100, sizeof(char));
    int izero = 0;
    int* index = (int*)malloc(inst->nnodes * sizeof(int));			
	double* value = (double*)malloc(inst->nnodes * sizeof(double));
    int nnz = 0;
	double rhs = 1.0; 
	char sense = 'E';
	double obj;
	double ub;
	double lb;

	// add binary var.s x(i,j)   
	if ( inst->xstart != -1 ) print_error("build_model: var. x cannot be redefined!");
	inst->xstart = CPXgetnumcols(env,lp); 		// position of the first x(,) variable  
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = 0; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);
			obj = dist(i,j,inst);   
			ub = ( i == j ) ? 0.0 : 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &zero, &ub, &binary, cname) ) print_error("wrong CPXnewcols on x var.s");
    		if ( CPXgetnumcols(env,lp)-1 != xpos_1(i, j, inst) ) print_error("wrong position for x var.s");
		}
	}

	// add continuous var.s u(i)  
	obj = 0.0;
	lb = 0.0; 
	if ( inst->ustart != -1 ) print_error("build_model: var. u cannot be redefined!");
	inst->ustart = CPXgetnumcols(env,lp); 		// position of the first u() variable   
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		sprintf(cname[0], "u(%d)", i+1);
		ub = ( i == 0 ) ? 0.0 : inst->nnodes-2;
		if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &continuous, cname) ) print_error("wrong CPXnewcols on u var.s");
    	if ( CPXgetnumcols(env,lp)-1 != upos(i, inst) ) print_error("wrong position for u var.s");
	}

	// x-constraints (out- and in-degree at nodes)
	for ( int h = 0; h < inst->nnodes; h++ )  // out-degree
	{
		nnz = 0;
		sprintf(cname[0], "outdeg(%d)", h+1);   
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			index[nnz] = xpos_1(h, i, inst);
			value[nnz++] = 1.0;
		}
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
			print_error("wrong CPXaddrows [x1]");
	} 
	for ( int h = 0; h < inst->nnodes; h++ ) // in-degree
	{
		nnz = 0;
		sprintf(cname[0], "indeg(%d)", h+1);
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			index[nnz] = xpos_1(i, h, inst);
			value[nnz++] = 1.0;
		}
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
			print_error("wrong CPXaddrows [x2]");
	}

	// u-x linking with lazy constraints
	// add lazy constraints  u_i - u_j + Mx_ij <= M - 1, for each arc (i,j) not touching node 0	
	nnz = 3;
	double M = inst->nnodes-1;
	rhs = M-1;
	sense = 'L';
	for ( int i = 1; i < inst->nnodes; i++ ) 
	{
		for ( int j = 1; j < inst->nnodes; j++ )
		{
			if ( i == j ) continue;
			sprintf(cname[0], "link(%d,%d)", i+1, j+1);
			index[0] = upos(i, inst);	
			value[0] = 1.0;	
			index[1] = upos(j, inst);
			value[1] = -1.0;
			index[2] = xpos_1(i, j, inst);
			value[2] = M;
			if ( CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname) ) 
				print_error("wrong CPXlazyconstraints for u-consistency");
		}
	}

	// other lazy constraints: x_ij + x_ji <= 1, for each arc (i,j) with i < j
	rhs = 1.0; 
	sense = 'L';
	nnz = 2;
	for ( int i = 0; i < inst->nnodes; i++ ) 
	{
		for ( int j = i+1; j < inst->nnodes; j++ ) 
		{
			sprintf(cname[0], "SEC on node pair (%d,%d)", i+1, j+1);
			index[0] = xpos_1(i, j, inst);
			value[0] = 1.0;
			index[1] = xpos_1(j, i, inst);
			value[1] = 1.0;
			if ( CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname) ) 
				print_error("wrong CPXlazyconstraints on 2-node SECs");
		}
	}

	if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model.lp", NULL);

	free(index);
	free(value);
	free(cname[0]);
	free(cname);
}


// Flow 1 (GG)
/***************************************************************************************************************************/
void build_model_2(instance *inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{  
	double zero = 0.0;  
	char binary = 'B';
	//char continuous = 'C';
	char integer = 'I';
	char **cname = (char **) calloc(1, sizeof(char *));
	cname[0] = (char *) calloc(100, sizeof(char));
    int izero = 0;
    int* index = (int*)malloc( (2*inst->nnodes - 2) * sizeof(int));			
	double* value = (double*)malloc( (2*inst->nnodes - 2) * sizeof(double));
    int nnz = 0;
	double rhs = 1.0; 
	char sense = 'E';
	double obj;
	double lb;
	double ub;

	// add binary var.s x(i,j)  
	lb = 0.0; 
	if ( inst->xstart != -1 ) print_error("build_model: var. x cannot be redefined!");
	inst->xstart = CPXgetnumcols(env,lp); 		// position of the first x(,) variable  
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = 0; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);
			obj = dist(i, j, inst);
			ub = ( i == j ) ? 0.0 : 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error("wrong CPXnewcols on x var.s");
    		if ( CPXgetnumcols(env,lp)-1 != xpos_1(i, j, inst) ) print_error("wrong position for x var.s");
		}
	}

	// add flow var.s y(i,j)
	obj = 0.0; 
	if ( inst->ystart != -1 ) print_error("build_model: var. y cannot be redefined!");
	inst->ystart = CPXgetnumcols(env,lp); 		// position of the first y(,) variable  
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = 0; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "y(%d,%d)", i+1,j+1);		  
			ub = ( (i == j) || (j == 0) ) ? 0.0 : /*CPX_INFBOUND*/ inst->nnodes-1; // modificato ub
			if ( CPXnewcols(env, lp, 1, &obj, &zero, &ub, &integer, cname) ) print_error("wrong CPXnewcols on y var.s");
    		if ( CPXgetnumcols(env,lp)-1 != ypos(i, j, inst) ) print_error("wrong position for y var.s");
    	}
    }

 	// x-constraints (out- and in-degree at nodes)
	for ( int h = 0; h < inst->nnodes; h++ )  // out-degree
	{
		nnz = 0;
		sprintf(cname[0], "outdeg(%d)", h+1);   
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			index[nnz] = xpos_1(h, i, inst);
			value[nnz++] = 1.0;
		}
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
			print_error("wrong CPXaddrows [x1]");
	}
	for ( int h = 0; h < inst->nnodes; h++ ) // in-degree
	{
		nnz = 0;
		sprintf(cname[0], "indeg(%d)", h+1);
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			index[nnz] = xpos_1(i, h, inst);
			value[nnz++] = 1.0;
		}
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
			print_error("wrong CPXaddrows [x2]");	
	}


	// y-constraints
	rhs = inst->nnodes-1;
	nnz = 0;
	sprintf(cname[0], "flow(%d)", 1);

	for ( int j = 1; j < inst->nnodes; j++ )
	{
		index[nnz] = ypos(0, j, inst);
		value[nnz++] = 1.0;
	}
	if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
		print_error("wrong CPXaddrows [y1]");

	rhs = 1.0;
	for ( int h = 1; h < inst->nnodes; h++ )
	{
		nnz = 0;
		sprintf(cname[0], "flow(%d)", h+1);
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			if ( i == h ) continue;
			index[nnz] = ypos(i, h, inst);
			value[nnz++] = 1.0;
		}
		for ( int j = 0; j < inst->nnodes; j++ )
		{
			if ( j == h ) continue;
			index[nnz] = ypos(h, j, inst);
			value[nnz++] = -1.0;
		}
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
			print_error("wrong CPXaddrows [y2]");		
	}

	// y-x linking 	
	double M = inst->nnodes-1;
	nnz = 2;
	rhs = 0.0;
	sense = 'L';

	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = 1; j < inst->nnodes; j++ )
		{
			if ( i == j ) continue;
			sprintf(cname[0], "link(%d,%d)", i+1, j+1);
			index[0] = ypos(i, j, inst);
			value[0] = 1.0;
			index[1] = xpos_1(i, j, inst);
			value[1] = -M;
			if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
				print_error("wrong CPXaddrows [link]");	
		}
	}

	if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model.lp", NULL);  

	free(index);
	free(value);
	free(cname[0]);
	free(cname);
}

// Flow 1' (with n-2)
/***************************************************************************************************************************/
void build_model_21(instance *inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{  
	double zero = 0.0;  
	char binary = 'B';
	//char continuous = 'C';
	char integer = 'I';
	char **cname = (char **) calloc(1, sizeof(char *));
	cname[0] = (char *) calloc(100, sizeof(char));
    int izero = 0;
    int* index = (int*)malloc( (2*inst->nnodes - 2) * sizeof(int));			
	double* value = (double*)malloc( (2*inst->nnodes - 2) * sizeof(double));
    int nnz = 0;
	double rhs = 1.0; 
	char sense = 'E';
	double obj;
	double lb = 0.0;
	double ub;

	// add binary var.s x(i,j)   
	if ( inst->xstart != -1 ) print_error("build_model: var. x cannot be redefined!");
	inst->xstart = CPXgetnumcols(env,lp); 		// position of the first x(,) variable  
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = 0; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);
			obj = dist(i,j,inst);
			ub = ( i == j ) ? 0.0 : 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error("wrong CPXnewcols on x var.s");
    		if ( CPXgetnumcols(env,lp)-1 != xpos_1(i, j, inst) ) print_error("wrong position for x var.s");
		}
	}

	// add flow var.s y(i,j)
	if ( inst->ystart != -1 ) print_error("build_model: var. y cannot be redefined!");
	inst->ystart = CPXgetnumcols(env,lp); 		// position of the first y(,) variable  
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = 0; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "y(%d,%d)", i+1,j+1);
			obj = 0.0;   
			ub = ( (i == j) || (j == 0) ) ? 0.0 : /*CPX_INFBOUND*/ inst->nnodes-1; // modificato ub
			if ( CPXnewcols(env, lp, 1, &obj, &zero, &ub, &integer, cname) ) print_error("wrong CPXnewcols on y var.s");
    		if ( CPXgetnumcols(env,lp)-1 != ypos(i, j, inst) ) print_error("wrong position for y var.s");
    	}
    }

 	// x-constraints (out- and in-degree at nodes)
	for ( int h = 0; h < inst->nnodes; h++ )  // out-degree
	{
		nnz = 0;
		sprintf(cname[0], "outdeg(%d)", h+1);   
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			index[nnz] = xpos_1(h, i, inst);
			value[nnz++] = 1.0;
		}
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
			print_error("wrong CPXaddrows [x1]");
	}
	for ( int h = 0; h < inst->nnodes; h++ ) // in-degree
	{
		nnz = 0;
		sprintf(cname[0], "indeg(%d)", h+1);
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			index[nnz] = xpos_1(i, h, inst);
			value[nnz++] = 1.0;
		}
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
			print_error("wrong CPXaddrows [x2]");	
	}


	// y-constraints
	rhs = inst->nnodes-1;
	nnz = 0;
	sprintf(cname[0], "flow(%d)", 1);

	for ( int j = 1; j < inst->nnodes; j++ )
	{
		index[nnz] = ypos(0, j, inst);
		value[nnz++] = 1.0;
	}
	if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
		print_error("wrong CPXaddrows [y1]");

	rhs = 1.0;
	for ( int h = 1; h < inst->nnodes; h++ )
	{
		nnz = 0;
		sprintf(cname[0], "flow(%d)", h+1);
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			if ( i == h ) continue;
			index[nnz] = ypos(i, h, inst);
			value[nnz++] = 1.0;
		}
		for ( int j = 0; j < inst->nnodes; j++ )
		{
			if ( j == h ) continue;
			index[nnz] = ypos(h, j, inst);
			value[nnz++] = -1.0;
		}
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
			print_error("wrong CPXaddrows [y2]");		
	}

	// y-x linking 	
	double M;
	nnz = 2;
	rhs = 0.0;
	sense = 'L';

	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = 1; j < inst->nnodes; j++ )
		{
			if ( i == j ) continue;
			M = ( i == 0 ) ? inst->nnodes-1 : inst->nnodes-2;
			sprintf(cname[0], "link(%d,%d)", i+1, j+1);
			index[0] = ypos(i, j, inst);
			value[0] = 1.0;
			index[1] = xpos_1(i, j, inst);
			value[1] = -M;
			if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
				print_error("wrong CPXaddrows [link]");	
		}
	}

	if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model.lp", NULL);  

	free(index);
	free(value);
	free(cname[0]);
	free(cname);
}


// Flow 1'' (with n-2 and lazy constraints)
/***************************************************************************************************************************/
void build_model_22(instance *inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{  
	double zero = 0.0;  
	char binary = 'B';
	//char continuous = 'C';
	char integer = 'I';
	char **cname = (char **) calloc(1, sizeof(char *));
	cname[0] = (char *) calloc(100, sizeof(char));
    int izero = 0;
    int* index = (int*)malloc( (2*inst->nnodes - 2) * sizeof(int));			
	double* value = (double*)malloc( (2*inst->nnodes - 2) * sizeof(double));
    int nnz = 0;
	double rhs = 1.0; 
	char sense = 'E';
	double obj;
	double lb = 0.0;
	double ub;

	// add binary var.s x(i,j)   
	if ( inst->xstart != -1 ) print_error("build_model: var. x cannot be redefined!");
	inst->xstart = CPXgetnumcols(env,lp); 		// position of the first x(,) variable  
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = 0; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);
			obj = dist(i,j,inst);   
			ub = ( i == j ) ? 0.0 : 1.0;
			if ( CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname) ) print_error("wrong CPXnewcols on x var.s");
    		if ( CPXgetnumcols(env,lp)-1 != xpos_1(i, j, inst) ) print_error("wrong position for x var.s");
		}
	}

	// add flow var.s y(i,j)
	if ( inst->ystart != -1 ) print_error("build_model: var. y cannot be redefined!");
	inst->ystart = CPXgetnumcols(env,lp); 		// position of the first y(,) variable  
	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = 0; j < inst->nnodes; j++ )
		{
			sprintf(cname[0], "y(%d,%d)", i+1,j+1);
			obj = 0.0;   
			ub = ( (i == j) || (j == 0) ) ? 0.0 : /*CPX_INFBOUND*/ ( i==0 & j!=0 ) ? inst->nnodes-1 : inst->nnodes-2; // modificato ub
			if ( CPXnewcols(env, lp, 1, &obj, &zero, &ub, &integer, cname) ) print_error("wrong CPXnewcols on y var.s");
    		if ( CPXgetnumcols(env,lp)-1 != ypos(i, j, inst) ) print_error("wrong position for y var.s");
    	}
    }

 	// x-constraints (out- and in-degree at nodes)
	for ( int h = 0; h < inst->nnodes; h++ )  // out-degree
	{
		nnz = 0;
		sprintf(cname[0], "outdeg(%d)", h+1);   
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			index[nnz] = xpos_1(h, i, inst);
			value[nnz++] = 1.0;
		}
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
			print_error("wrong CPXaddrows [x1]");
	}
	for ( int h = 0; h < inst->nnodes; h++ ) // in-degree
	{
		nnz = 0;
		sprintf(cname[0], "indeg(%d)", h+1);
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			index[nnz] = xpos_1(i, h, inst);
			value[nnz++] = 1.0;
		}
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
			print_error("wrong CPXaddrows [x2]");	
	}


	// y-constraints
	rhs = inst->nnodes-1;
	nnz = 0;
	sprintf(cname[0], "flow(%d)", 1);

	for ( int j = 1; j < inst->nnodes; j++ )
	{
		index[nnz] = ypos(0, j, inst);
		value[nnz++] = 1.0;
	}
	if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
		print_error("wrong CPXaddrows [y1]");

	rhs = 1.0;
	for ( int h = 1; h < inst->nnodes; h++ )
	{
		nnz = 0;
		sprintf(cname[0], "flow(%d)", h+1);
		for ( int i = 0; i < inst->nnodes; i++ )
		{
			if ( i == h ) continue;
			index[nnz] = ypos(i, h, inst);
			value[nnz++] = 1.0;
		}
		for ( int j = 0; j < inst->nnodes; j++ )
		{
			if ( j == h ) continue;
			index[nnz] = ypos(h, j, inst);
			value[nnz++] = -1.0;
		}
		if ( CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, cname) )
			print_error("wrong CPXaddrows [y2]");		
	}

	// y-x linking 	
	double M;
	nnz = 2;
	rhs = 0.0;
	sense = 'L';

	for ( int i = 0; i < inst->nnodes; i++ )
	{
		for ( int j = 1; j < inst->nnodes; j++ )
		{
			if ( i == j ) continue;
			M = ( i == 0 ) ? inst->nnodes-1 : inst->nnodes-2;
			sprintf(cname[0], "link(%d,%d)", i+1, j+1);
			index[0] = ypos(i, j, inst);
			value[0] = 1.0;
			index[1] = xpos_1(i, j, inst);
			value[1] = -M;
			if ( CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &izero, index, value, cname) )
				print_error("wrong CPXaddrows [link]");	
		}
	}

	if ( VERBOSE >= -100 ) CPXwriteprob(env, lp, "model.lp", NULL);  

	free(index);
	free(value);
	free(cname[0]);
	free(cname);
}


/*********************************************************************************************************************************/
void build_sol(const double* xstar, instance* inst, int* succ, int* comp, int* ncomp) 
/*********************************************************************************************************************************/
{
	switch (inst->model_type) 
	{
		case 0: build_sol_0(xstar, inst, succ, comp, ncomp);
			break;

		case 3: build_sol_0(xstar, inst, succ, comp, ncomp);
			break;

		case 30: build_sol_0(xstar, inst, succ, comp, ncomp);
			break;

		case 31: build_sol_0(xstar, inst, succ, comp, ncomp);
			break;

		case 310: build_sol_0(xstar, inst, succ, comp, ncomp);
			break;

		default: build_sol_1(xstar, inst, succ, comp, ncomp);
			break;
	}
}


#define DEBUG   //da commentare se non si vuole il debug
/*********************************************************************************************************************************/
void build_sol_0(const double* xstar, instance* inst, int* succ, int* comp, int* ncomp) // for each node finds the successor and the conected component
																					  // those are saved in succ[] and comp[]. ncomp is the number of components
/*********************************************************************************************************************************/
{

#ifdef DEBUG
	int* degree = (int*)calloc(inst->nnodes, sizeof(int));
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = i + 1; j < inst->nnodes; j++)
		{
			int k = xpos(i, j, inst);
			if (fabs(xstar[k]) > EPS && fabs(xstar[k] - 1.0) > EPS ) print_error("wrong xstar in build_sol_0");
			if (xstar[k] > 0.5)
			{
				++degree[i];
				++degree[j];
			}
		}
	}
	for (int i = 0; i < inst->nnodes; i++)
	{
		if (degree[i] != 2) print_error("wrong degree in build_sol_0");
	}
	free(degree);
#endif

	* ncomp = 0;
	for (int i = 0; i < inst->nnodes; i++)
	{
		succ[i] = -1;
		comp[i] = -1;
	}

	for (int start = 0; start < inst->nnodes; start++)
	{
		if (comp[start] >= 0) continue;  // node "start" was already visited, just skip it

		// a new component is found
		(*ncomp)++;
		int i = start;				 //i is the current node in which we define the successor
		int prev = -1;				 //prev is the previous node
		while ( comp[start] == -1 )  // go and visit the current component ... the component of the start node is defined only at the end of the cycle
		{
			for (int j = 0; j < inst->nnodes; j++)
			{
				if (i != j && xstar[xpos(i, j, inst)] > 0.5 && j != prev) // the edge [i,j] is selected in xstar and j was not visited before 
				{
					succ[i] = j;	 // during each step I define the successor of i
					comp[j] = *ncomp - 1;// and the component of the successor j
					prev = i;
					i = j;
					break;
				}
			}
		}

		// go to the next component...
	}
}


/*********************************************************************************************************************************/
void build_sol_1(const double* xstar, instance* inst, int* succ, int* comp, int* ncomp) // for each node finds the successor and the conected component
																					  // those are saved in succ[] and comp[]. ncomp is the number of components
/*********************************************************************************************************************************/
{

#ifdef DEBUG
	int* indegree = (int*)calloc(inst->nnodes, sizeof(int));
	int* outdegree = (int*)calloc(inst->nnodes, sizeof(int));
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = 0; j < inst->nnodes; j++)
		{
			int k = xpos_1(i, j, inst);
			if (fabs(xstar[k]) > EPS && fabs(xstar[k] - 1.0) > EPS) print_error("wrong xstar in build_sol_1");
			if (xstar[k] > 0.5)
			{
				++outdegree[i];
				++indegree[j];
			}
		}
	}
	for (int i = 0; i < inst->nnodes; i++)
	{
		if (indegree[i] != 1 || outdegree[i] != 1) print_error("wrong degree in build_sol_1");
	}
	free(indegree);
	free(outdegree);
#endif

	* ncomp = 0;
	for (int i = 0; i < inst->nnodes; i++)
	{
		succ[i] = -1;
		comp[i] = -1;
	}

	for (int start = 0; start < inst->nnodes; start++)
	{
		if (comp[start] >= 0) continue;  // node "start" was already visited, just skip it

		// a new component is found
		(*ncomp)++;
		int i = start;				 //i is the current node in which we define the successor
		while (comp[start] == -1)  // go and visit the current component ... the component of the start node is defined only at the end of the cycle
		{
			for (int j = 0; j < inst->nnodes; j++)
			{
				if ( xstar[xpos_1(i, j, inst)] > 0.5 ) // the edge [i,j] is selected in xstar and j was not visited before 
				{
					succ[i] = j;	 // during each step I define the successor of i
					comp[j] = *ncomp - 1;// and the component of the successor j
					i = j;
					break;
				}
			}
		}

		// go to the next component...
	}
}


/**************************************************************************************************************************/
/* This function add the SEC constraints until the the number of connected components is 1*/
int addSEC(instance* inst, CPXENVptr env, CPXLPptr lp, int* succ, int* comp, int* ncomp) 
/**************************************************************************************************************************/
{
	//array to contain the constraints'name
	char** rname = (char**)calloc(1, sizeof(char*));
	rname[0] = (char*)calloc(100, sizeof(char));
	int* index = (int*)malloc(inst->nnodes*(inst->nnodes-1)/2  * sizeof(int));			
	double* value = (double*)malloc(inst->nnodes*(inst->nnodes-1)/2  * sizeof(double));
	int comp_nnodes;
	char sense = 'L';
	double rhs;
	int izero = 0;
	int nnz;

	while (*ncomp > 1)
	{
		for (int c = 0; c < *ncomp; c++)
		{
			// add SEC with sparse representation
			comp_nnodes = 0;

			for (int i = 0; i < inst->nnodes; i++) 
			{
				if (comp[i] != c) continue;
				comp_nnodes++;
			}

			nnz = 0;

			for (int i = 0; i < inst->nnodes; i++) 
			{
				if (comp[i] != c) continue;
				for (int j = i+1; j < inst->nnodes; j++)
				{
					if (comp[j] != c) continue;
					//one more non zero component
					//position and value of the coefficient
					index[nnz] = xpos(i, j, inst);
					value[nnz] = 1.0;
					nnz++;
				}
			}

			//printf("adding the SEC(%d) with %d nodes\n", lastrow - inst->nnodes + 1, comp_nnodes);
			int lastrow = CPXgetnumrows(env, lp);
			rhs = comp_nnodes - 1;
			sprintf(rname[0], "SEC(%d)", lastrow - inst->nnodes + 1);
			if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, rname)) print_error("CPXaddrows");
		}

		if (VERBOSE >= -100) CPXwriteprob(env, lp, "model.lp", NULL);

		if (CPXmipopt(env, lp)) print_error("CPXmipopt");

		int ncols = CPXgetnumcols(env, lp);
		if (CPXgetx(env, lp, inst->best_sol, 0, ncols - 1)) print_error("CPXgetx");
		build_sol(inst->best_sol, inst, succ, comp, ncomp);

		if ( VERBOSE >= 60 ) printf("\n &&& NCOMP: %d\n", *ncomp);
	}

	free(index);
	free(value);
	free(rname[0]);
	free(rname);

	return 0;
}


/**************************************************************************************************************************/
static int CPXPUBLIC mylazycallback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
/**************************************************************************************************************************/
{
	*useraction_p = CPX_CALLBACK_DEFAULT;
	instance* inst = (instance*)cbhandle; 			// casting of cbhandle

	// get solution xstar
	double* xstar = (double*)malloc(inst->ncols * sizeof(double));
	if (CPXgetcallbacknodex(env, cbdata, wherefrom, xstar, 0, inst->ncols - 1)) return 1; // xstar = current x from CPLEX-- xstar starts from position 0

	// get some random information at the node (as an example)
	double objval = CPX_INFBOUND; CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);
	int mythread = -1; CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &mythread);
	inst->nthread = -1; CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_USER_THREADS, &inst->nthread);
	double zbest; CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &zbest);

	//apply cut separator and possibly add violated cuts
	int ncuts = myseparation(inst, xstar, env, cbdata, wherefrom);

	if (inst->model_type == 30 && ncuts >= 1) {
		patching(inst, xstar, mythread);
	}

	free(xstar);

	if (ncuts >= 1) *useraction_p = CPX_CALLBACK_SET; 		// tell CPLEX that cuts have been created
	return 0; 												// return 1 would mean error --> abort Cplex's execution
}


/**************************************************************************************************************************/
void patching(instance* inst, double* xstar, int thread_idx) {
/**************************************************************************************************************************/
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	int* comp = (int*)malloc(inst->nnodes * sizeof(int));
	int ncomp;
	double min_delta;
	//int patch_idx[4] = {-1,-1,-1,-1};
	int A = -1, B = -1, C = -1, D = -1;
	int flag = 0;

	build_sol_0(xstar, inst, succ, comp, &ncomp);

	while (ncomp > 1) {

		min_delta = CPX_INFBOUND;

		for (int c1 = 0; c1 < ncomp; c1++) {
			for (int c2 = c1 + 1; c2 < ncomp; c2++) {

				for (int a = 0; a < inst->nnodes; a++) {
					if (comp[a] != c1) continue;
					int b = succ[a];

					for (int c = 0; c < inst->nnodes; c++) {
						if (comp[c] != c2) continue;
						int d = succ[c];

						double delta1 = (dist(a, c, inst) + dist(b, d, inst)) - (dist(a, b, inst) + dist(c, d, inst));
						double delta2 = (dist(a, d, inst) + dist(b, c, inst)) - (dist(a, b, inst) + dist(c, d, inst));

						if (delta1 < delta2) {
							if (delta1 < min_delta) {
								flag = 0;
								min_delta = delta1;
								/*patch_idx[0] = xpos(a, b, inst);
								patch_idx[1] = xpos(c, d, inst);
								patch_idx[2] = xpos(a, c, inst);
								patch_idx[3] = xpos(b, d, inst);*/
								A = a; B = b; C = c; D = d;
							}
						}
						else {
							if (delta2 < min_delta) {
								flag = 1;
								min_delta = delta2;
								/*patch_idx[0] = xpos(a, b, inst);
								patch_idx[1] = xpos(c, d, inst);
								patch_idx[2] = xpos(a, c, inst);
								patch_idx[3] = xpos(b, d, inst);*/
								A = a; B = b; C = c; D = d;
							}
						}
					}

				}
			}
		}

		/*xstar[patch_idx[0]] = 0;
		xstar[patch_idx[1]] = 0;
		xstar[patch_idx[2]] = 1;
		xstar[patch_idx[3]] = 1;*/
		xstar[xpos(A, B, inst)] = 0;
		xstar[xpos(C, D, inst)] = 0;

		//printf("A:%d, B:%d, C:%d, D:%d\n", A, B, C, D);

		if (flag == 0) {
			xstar[xpos(A, C, inst)] = 1;
			xstar[xpos(B, D, inst)] = 1;

			int prev = D;
			int curr = succ[D];
			int next;

			while (curr != C)
			{
				next = succ[curr];
				succ[curr] = prev;
				prev = curr;
				curr = next;
			}
			succ[C] = prev;
			succ[D] = B;
			succ[A] = C;

		}
		else {
			xstar[xpos(A, D, inst)] = 1;
			xstar[xpos(B, C, inst)] = 1;
			succ[A] = D;
			succ[C] = B;
		}

		int temp;
		//updating the comp[] of the components > c2
		for (int i = comp[C] + 1; i < ncomp; i++) {
			for (int n = 0; n < inst->nnodes; n++) {
				if (comp[n] == i) {
					temp = succ[n];
					while (temp != n) {
						comp[temp] = i - 1;
						temp = succ[temp];
					}
					comp[n] = i - 1;
				}
			}
		}

		//merging the comp[] of the components of A-B and C-D
		temp = succ[A];
		//printf("A = %d    ", A);
		while (temp != B)
		{
			comp[temp] = comp[A];
			//printf(", %d ", temp);
			temp = succ[temp];
		}
		//printf(", %d \n\n", temp);

		/*for (int i = 0; i < ncomp; i++) {
			printf("%d ", i);
		}
		printf("\n\n");*/

		ncomp--;

	}
	double costo = cost(inst, xstar);
	
	if (inst->patch_cost[thread_idx] == -1 || (costo < inst->patch_cost[thread_idx])) {
		if (VERBOSE >= 60) printf(" ... Patching heuristic added in patch_vec[%d] ...\n", thread_idx);
		for (int i = 0; i < inst->ncols; i++) {
			inst->patch_vec[thread_idx][i] = xstar[i];
		}
		inst->patch_cost[thread_idx] = costo;
	}

	free(succ);
	free(comp);

}


/**************************************************************************************************************************/
static int CPXPUBLIC myheuristicfunc(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, double* objval_p, double* x, int* checkfeas_p, int* useraction_p)
/**************************************************************************************************************************/
{
	//if (VERBOSE >= 60) printf(" ... Inside myheuristicfunc ...\n");
	*useraction_p = CPX_CALLBACK_DEFAULT;
	instance* inst = (instance*)cbhandle; 			// casting of cbhandle

	int found = 0;
	double min; if( CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_BEST_INTEGER, &min)) min = CPX_INFBOUND;
	int idx;

	//int mythread = -1; CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &mythread);

	for (int i = 0; i < inst->nthread; i++) {
		if (inst->patch_cost[i] != -1 && inst->patch_cost[i] < min-EPS) {
			min = inst->patch_cost[i];
			inst->patch_cost[i] = -1;
			idx = i;
			found = 1;
		}	
	}

	if (found == 1) {
		//tell cplex what is the heuristic integer solution
		for (int i = 0; i < inst->ncols; i++) {
			x[i] = inst->patch_vec[idx][i];
		}

		//update the objective cost
		*objval_p = min;

		//tell cplex to check if x violate any constraint
		*checkfeas_p = 1;

		//tell cplex to check the vector x
		*useraction_p = CPX_CALLBACK_SET;
		if (VERBOSE >= 60) printf(" ... Found something in patch_vec[%d], cost = %lf ...\n", idx, *objval_p);
	}

	return 0;
}


/**************************************************************************************************************************/
int myseparation(instance* inst, double* xstar, CPXCENVptr env, void* cbdata, int wherefrom) {
/**************************************************************************************************************************/

	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	int* comp = (int*)malloc(inst->nnodes * sizeof(int));
	int ncomp = 0;
	int* index = (int*)malloc(inst->nnodes * (inst->nnodes - 1) / 2 * sizeof(int));
	double* value = (double*)malloc(inst->nnodes * (inst->nnodes - 1) / 2 * sizeof(double));
	int comp_nnodes;
	int nnz;
	char sense = 'L';

	build_sol_0(xstar, inst, succ, comp, &ncomp);

	if (ncomp > 1)
	{

		for (int c = 0; c < ncomp; c++)
		{
			// add SEC with sparse representation
			comp_nnodes = 0;
			nnz = 0; //the constraint is a complete graph

			for (int i = 0; i < inst->nnodes; i++)
			{
				if (comp[i] != c) continue;
				comp_nnodes++;
				for (int j = i + 1; j < inst->nnodes; j++)
				{
					if (comp[j] != c) continue;

					//position and value of the coefficient
					index[nnz] = xpos(i, j, inst);
					value[nnz++] = 1.0;
				}
			}
			double rhs = comp_nnodes - 1;
			if (CPXcutcallbackadd(env, cbdata, wherefrom, nnz, rhs, sense, index, value, 0)) print_error("CPXcutcallbackadd");
		}
	}

	free(index);
	free(value);
	free(succ);
	free(comp);

	return (ncomp > 1) ? ncomp : 0;
}


/**************************************************************************************************************************/
static int CPXPUBLIC generic_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle)
/**************************************************************************************************************************/
{
	//if (VERBOSE >= 60) printf(" ... Inside Generic_callback switch ... \n");
	instance* inst = (instance*)userhandle; 			// casting of userhandle

	switch (contextid) {

	case CPX_CALLBACKCONTEXT_CANDIDATE:
		return callback_candidate(inst, context);
		break;

	case CPX_CALLBACKCONTEXT_GLOBAL_PROGRESS:
		return callback_heuristic(inst, context);

		break;
	}

	return 0;
}


/**************************************************************************************************************************/
int callback_candidate(instance* inst, CPXCALLBACKCONTEXTptr context)
/**************************************************************************************************************************/
{
	//if (VERBOSE >= 60) printf(" ... Inside Generic_callback_candidate ... \n");
	// get solution xstar
	double* xstar = (double*)malloc(inst->ncols * sizeof(double));
	double objval = CPX_INFBOUND;
	if (CPXcallbackgetcandidatepoint(context, xstar, 0, inst->ncols - 1, &objval)) return 1;
	// get some random information at the node (as an example)
	int mythread = -1; CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread);
	double zbest; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_BND, &zbest);

	//apply cut separator and possibly add violated cuts
	if (salvagnin_separation(inst, xstar, context, mythread)) print_error("salvagnin_separation");

	free(xstar);

	return 0; 											// return 1 would mean error --> abort Cplex's execution
}


/**************************************************************************************************************************/
int callback_heuristic(instance* inst, CPXCALLBACKCONTEXTptr context)
/**************************************************************************************************************************/
{
	//if (VERBOSE >= 60) printf(" ... Inside Generic_callback_heuristic ... \n");
	int* index = (int*)malloc(inst->ncols *  sizeof(int));
	double* value = (double*)malloc(inst->ncols * sizeof(double));
	int found = 0;
	int idx;

	//inject the heuristic solution only if it improves the best solution found
	double min; if (CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &min)) { printf("update min variable!!!!\n"); min = CPX_INFBOUND; }
	//double min = CPX_INFBOUND;

	for (int i = 0; i < inst->nthread; i++) {
		if (inst->patch_cost[i] != -1 && inst->patch_cost[i] < min-EPS) {
			min = inst->patch_cost[i];
			inst->patch_cost[i] = -1;
			idx = i;
			found = 1;
		}
	}

	if (found == 1) {
		//tell cplex what is the heuristic integer solution
		for (int i = 0; i < inst->ncols; i++) {
			index[i] = i;
			value[i] = inst->patch_vec[idx][i];
		}


		if (VERBOSE >= 60) printf(" ... Found something in patch_vec[%d], cost = %lf ...\n", idx, min);
	
		CPXCALLBACKSOLUTIONSTRATEGY strat = CPXCALLBACKSOLUTION_NOCHECK; //could also be CPXCALLBACKSOLUTION_CHECKFEAS, see documentation for other options

		CPXcallbackpostheursoln(context, inst->ncols, index, value, inst->patch_cost[idx], strat);
	}

	free(index);
	free(value);
	return 0;
}


/**************************************************************************************************************************/
int salvagnin_separation(instance* inst, double* xstar, CPXCALLBACKCONTEXTptr context, int thread_idx)
/**************************************************************************************************************************/
{
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	int* comp = (int*)malloc(inst->nnodes * sizeof(int));
	int ncomp = 0;

	build_sol_0(xstar, inst, succ, comp, &ncomp);

	int* index = (int*)malloc(inst->nnodes * (inst->nnodes - 1) / 2 * sizeof(int));
	double* value = (double*)malloc(inst->nnodes * (inst->nnodes - 1) / 2 * sizeof(double));
	int zero = 0;
	int comp_nnodes;
	int nnz;


	if (ncomp > 1)
	{
		for (int c = 0; c < ncomp; c++)
		{
			comp_nnodes = 0;
			nnz = 0;

			// add SEC with sparse representation
			for (int i = 0; i < inst->nnodes; i++)
			{
				if (comp[i] != c) continue;

				comp_nnodes++;

				for (int j = i + 1; j < inst->nnodes; j++)
				{
					if (comp[j] != c) continue;

					//position and value of the coefficient
					index[nnz] = xpos(i, j, inst);
					value[nnz] = 1.0;
					nnz++;
				}
			}

			char sense = 'L';
			double rhs = comp_nnodes - 1;

			if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &zero, index, value)) print_error("CPXcallbackrejectcandidate");
		}

		if (inst->model_type == 310)
			patching(inst, xstar, thread_idx);

	}

	free(index);
	free(value);
	free(succ);
	free(comp);

	return 0;
}

