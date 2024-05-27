
static unsigned long int num_steps = (unsigned long int)((1<<30)-1) ;
double step;

int main (int argc, char** argv) { 
    unsigned long int i;
    double pi;
    double PI25DT = 3.141592653589793238462643;
    double x;
    double sum=0.0;
    int factor=1;

    struct timeval t1, t2;
    double segundos;

	int __taskid = -1, __numprocs = -1;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&__numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&__taskid);
	MPI_Get_processor_name(processor_name,&namelen);

	if (__taskid == 0) {

	    if (argc == 1) {
	        printf("num_steps %ld, ¿Factor de escala (1..4)?\n", num_steps);
	        if (scanf("%d", &factor) <= 0)
		    printf("scanf error, factor sin cambio %d\n", factor);
	    }
	    else factor = atoi(argv[1]);

	    num_steps = num_steps * factor;

	    printf("%ld num_steps, %25.23f step size i: %ld size num_steps: %ld\n",
			num_steps, (double)1.0/(double) num_steps, sizeof(i), sizeof(num_steps));
			
	    gettimeofday(&t1, NULL);

	    step = 1.0/(double) num_steps;

	}
	MPI_Bcast(&num_steps, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(step, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

                for (i=0;i< num_steps; i++) {
        x = (i+0.5)*step;
        sum += 4.0/(1.0+x*x);
    }
    pi = step * sum;

    gettimeofday(&t2, NULL);
    segundos = (((t2.tv_usec - t1.tv_usec)/1000000.0f)  + (t2.tv_sec - t1.tv_sec));

    printf("Pi es %25.23f, calc con %ld pasos en %f segundos\n", pi,num_steps,segundos);
    printf("Pi es %25.23f, Error relativo %10.8e\n", PI25DT, (double) (pi - PI25DT)/PI25DT);

	MPI_Finalize();

    return(0);
}

