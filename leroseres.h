#ifndef __LEROSERES_H
#define __LEROSERES_H

#define __LEROSERES_N 1000

int    __leroseres_n;
double __leroseres_p[__LEROSERES_N];
double __leroseres_c[__LEROSERES_N];

#define __LEROSE_BEAM_E 2.2

#define __LEROSE_m_TO_GeV 0.1276 //  Curve we read in is in m, which
				 //  we need to convert to GeV

void LoadLeroseRes(){
    FILE *f = fopen("./brindza_4_degrees_sampdist.txt", "r");

    int n = 0;
    int nread = 2;
    double p, c;


    while( nread == 2 && !feof(f) && n < __LEROSERES_N ){
	nread = fscanf(f, "%lf%lf", &c, &p );
	if( nread == 2 ){
	    __leroseres_p[n] = p;
	    __leroseres_c[n] = c;
	    n++;
	}
    }

    fclose(f);
    __leroseres_n = n;

    return;
}

double getlerose_ressmear(){
    int i = 0;
    double x;
    double u = gRandom->Uniform();

    while( i < __leroseres_n && __leroseres_c[i] < u ){i++;}

    x = (u - __leroseres_c[i-1])/(__leroseres_c[i] - __leroseres_c[i-1]);

    return (__leroseres_p[i-1]*(1.0-x) + __leroseres_p[i]*x)
	   *__LEROSE_m_TO_GeV/__LEROSE_BEAM_E;
}


#endif//__LEROSERES_H
