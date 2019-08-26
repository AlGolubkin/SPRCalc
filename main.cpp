#include <cstdlib>
#include <stdio.h>
#include <iostream>

using namespace std;

/*
 * 
 */
double intencity(double REFRACTIVE_INDEX, double ENERGY);

int main(int argc, char** argv) {
    
    double I = intencity(1,2);
    
    cout << I ;

    return 0;
}

double intencity(double REFRACTIVE_INDEX, double ENERGY){
    
    double n = REFRACTIVE_INDEX;
    double En = ENERGY;
    
    double I = n * En;
    
    return I;
    
}
