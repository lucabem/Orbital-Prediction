#include <stdio.h>
#include <stdlib.h>
#include <assert.h>



double Norma2(double vector[]);
double* str2num(char *palabra);

void Test_Norma2();
void Test_Str2num();

int main()
{

    Test_Norma2();
    printf(">>>>>>>>>>> Test_Norma() superado \n");

    Test_Str2num();
    printf(">>>>>>>>>>> Test_Str2num() superado \n");


    return 0;
}


double Norma2(double vector[]){
    return sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
}

double* str2num(char *palabra){
    //return atoi(palabra);
    int i = 0;
    double resultado[10];
    printf("%s", palabra[0]);
    while (palabra[i] != NULL){
        resultado[i] = atof(palabra[i]);

        i++;
    }

    return resultado;
}

void Test_Norma2(){
    double v1[3] = {0.0, 0.0, 0.0};
    assert(Norma(v1) == 0);

    double v2[3] = {1.0, 1.0, 1.0};
    assert(Norma(v2) == sqrt(3));

    double v3[3] = {4.0, 3.0, 0.0};
    assert(Norma(v3) == 5);
}

void Test_Str2num()
{
    char *a[3]= {"123", "0.0", "-0.5"};
    double *p = str2num(a);

    assert(p[0] == 123);
    assert(p[1] == 0.5);
    assert(p[2] == -0.5);
}
