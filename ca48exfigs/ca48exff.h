#ifndef __CA48EXFF_H
#define __CA48EXFF_H
#define NBIN 200
#define REACH 3

#define __CA48NFILE 10
int    __ca48ex_n[__CA48NFILE];
double __ca48ex_q[__CA48NFILE][5000];
double __ca48ex_f[__CA48NFILE][5000];
double __ca48ex_min[__CA48NFILE];
double __ca48ex_max[__CA48NFILE];

void LoadCa48ExState(){
    printf("Loading Ca48 excited state data\n");
    int i, nscan;
    FILE *file;

/*
    char fn[__CA48NFILE][25] = {
	"2p", "3m1", "3p1", "3m2", "5m","5mC", "4p1", "4p2", "5p", "3m2"
    };
*/

//My Changes
    char fn[__CA48NFILE][25] = {
	"2+", "4+", "6+", "8+"
    };
//


    int n;
    for( i = 0; i < __CA48NFILE; i++ ){
	file = fopen(Form("ca48exfigs/%s_sm.txt", fn[i]), "r");
	if( !file ){ printf("Couldn'd load ca48exfigs/%s_sm.txt\n", fn[i] ); exit(1);}

	n = 0;
	__ca48ex_max[i] = -1e9;
	__ca48ex_min[i] =  1e9;

	nscan = 2;
	while( nscan==2 && !feof(file)){
	    nscan = fscanf(file, "%lf%lf", &__ca48ex_q[i][n], &__ca48ex_f[i][n] );
	    if( nscan == 2 ){ 
		if( __ca48ex_q[i][n] > __ca48ex_max[i] ){ __ca48ex_max[i] = __ca48ex_q[i][n]; }
		if( __ca48ex_q[i][n] < __ca48ex_min[i] ){ __ca48ex_min[i] = __ca48ex_q[i][n]; }
		n++; 
	    }
	}

	__ca48ex_n[i] = n;
	fclose(file);
    }

    printf("Loaded\n");
    return;
}

double GetCa48ExFF( int i, double q2 ){
    double qf = sqrt(q2)/0.197; // but into inverse fermi


    double scale = (qf - __ca48ex_min[i])/(__ca48ex_max[i] - __ca48ex_min[i]);
//    printf("max min = %f %f\n", __ca48ex_max[i], __ca48ex_min[i] );
//    printf("Q2 = %f  qf =  %f (%f)\n", Q2, qf, scale);
//   if( scale < 0 || scale > 1 ){ printf("Ca48 excited state %d scale out of range (qf %f, scale %f)\n", i, qf, scale); }
    if( scale < 0 || scale > 1 ){ return 0.0; }
   
    int idx = (int) (scale*__ca48ex_n[i]);
    double x = scale*__ca48ex_n[i] - (double) idx;

    if( i == 1  || i == 0){
//	printf("state %d q = %f, F2= %e (%f) (%e - %e)\n", i, qf, exp(log(__ca48ex_f[i][idx])*(1.0-x) + log(__ca48ex_f[i][idx+1])*x), x, __ca48ex_f[i][idx], __ca48ex_f[i][idx+1] );
    }

    return exp(log(__ca48ex_f[i][idx])*(1.0-x) + log(__ca48ex_f[i][idx+1])*x);
}

#endif//__CA48EXFF_H

