#include "tsp.h"           


int main(int argc, char **argv) 
{ 
	if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }       
	if ( VERBOSE >= 2 ) { for (int a = 0; a < argc; a++) printf("%s ", argv[a]); printf("\n"); }

	double t1 = second(); 
	instance inst;

	parse_command_line(argc,argv, &inst);     
		  
	read_input(&inst); 
	if ( inst.heuristic != -1 )
	{
		if ( TSPapprox(&inst) ) print_error(" error within TSPapprox()");
	}
	else 
	{
		if ( TSPopt(&inst) ) print_error(" error within TSPopt()");	
	}
	double t2 = second(); 
    
    if ( VERBOSE == 0 ) // for csv file
    {
		printf("%f\n", inst.zbest);  
		//printf("%f\n", t2-t1);  
    }
	if ( VERBOSE >= 1 )   
	{
		printf("... TSP problem solved in %lf sec.s\n", t2-t1);  
	}

	//plot_sol(&inst);
	
	free_instance(&inst);
	return 0; 
}         
