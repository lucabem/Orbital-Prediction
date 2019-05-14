
#include <stdio.h>
#include <stdlib.h>

int main()
{
    int opcion;

    do{
    printf("\n|***************************|\n");
    printf("|    Menu de seleccion      |\n");
    printf("|***************************|\n");
    printf("| 1 - Test MatlabUtilities  |\n");
    printf("| 2 - Test Ficheros .m      |\n");
    printf("| 3 - Ejemplo 1             |\n");
    printf("| 4 - Ejemplo 2             |\n");
    printf("| 5 - Ejemplo 3             |\n");
    printf("| 6 - Ejemplo 5             |\n");
    printf("| 7 - Ejemplo 6             |\n");
    printf("| 8 - Ejemplo 7             |\n");
    printf("| 9 - Salir                 |\n");
    printf("*****************************\n");
    printf("Introduce la opcion: ");
    scanf("%d", &opcion);


    switch(opcion)
    {
        case 1:
            Test_Matlab();
            printf("\n Test MatlabUtilities superados\n\n");
            break;

        case 2:
            Test_Funciones();
            printf("\n Test Funciones .m superados \n\n");
            break;

        case 3:
            Example1();
            break;

        case 4:
            Example2();
            break;

        case 5:
            Example3();
            break;

        case 6:
            Example5();
            break;

        case 7:
            Example6();
            break;

        case 8:
            Example7();
            break;
    }


    }while (opcion != 9);


    return 0;
}


