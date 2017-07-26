#ifndef __EXSTATE_H
#define __EXSTATE_H

#include "ca48exff.h"

int nState( int Z, int A ){
    // Number of excited states PLUS THE GROUND STATE
    if( Z == 20 && A == 40 ){
	// Ca40
	return 3;
    }

    if( Z == 20 && A == 48 ){
	// Ca48
	return 11;
    }

//My Changes
    if( Z == 82 && A == 126){
	// Pb208
	return 5;
    }
//

    return 1;
}

double ExEnergy( int Z, int A, int j ){
    if( j == 0 ){ return 0.0; }

    if( Z == 20 && A == 40 ){
	// Ca40
	
	switch(j){
	    case 1:
		return 0.00373;
		break;
	    case 2:
		return 0.00390;
		break;
	    default:
		return 0.0;
	}
    }

    if( Z == 20 && A == 48 ){
	// Ca48
	
	switch(j){
	    case 1:
		return 0.003837;
		break;
	    case 2:
		return 0.004507;
		break;
	    case 3:
		return 0.004608;
		break;
	    case 4:
		return 0.005372;
		break;
	    case 5:
		return 0.005726;
		break;
	    case 6:
		return 0.005726;
		break;
	    case 7:
		return 0.006340;
		break;
	    case 8:
		return 0.006647;
		break;
	    case 9:
		return 0.006893;
		break;
	    case 10:
		return 0.007657;
		break;
	    default:
		return 0.0;
	}
    }

    if( Z == 82 && A == 126){
	// Pb208
	switch(j){
	    case 1:
		return 0.004085;
		break;
	    case 2:
		return 0.004323;
		break;
	    case 3:
		return 0.004424;
		break;
	    case 4:
		return 0.004610;
		break;
	    default:
		return 0.0;
	}
    }

    return 0.0;
}

double exFF( int Z, int A, int j, double Q2, double th, double q3v2 ){
    double ff2;
    if( j == 0 ){ return 1.0; }

    if( Z == 20 && A == 40 ){
	// Ca40

	switch(j){
	    case 1:
		return -4.328e-4 + 8.72449e-2*Q2 + 1.1217*Q2*Q2;
		break;
	    case 2:
		return -1.70443-4 + 8.10625e-2*Q2 - 1.70236*Q2*Q2;
		break;
	    default:
		return 0.0;
	}
    }

    if( Z == 20 && A == 48 ){
	// Ca48

	//  These are are scaled by 4piZ^2 (for no good reason)
	ff2 = GetCa48ExFF(j-1,q3v2)*4.0*3.14159/(Z*Z);

	///    3p1     5m        5p         2m
	if( j == 3 || j == 5 || j == 9 || j == 10 ){
		// These are in transverse form factors
//	    printf("ff before (%d) %f\n",j, ff2);
		ff2 *= (0.5*Q2/q3v2 + tan(th/2.0)*tan(th/2.0));
//	    printf("ff after %f\n", ff2);
	}

	return ff2;
    }

    if( Z == 82 && A == 126){
	// Pb208
	ff2 = GetCa48ExFF(j-1,q3v2);
	return ff2;	
    }

    return 0.0;
}

#endif//__EXSTATE_H
