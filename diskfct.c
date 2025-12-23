#include <stdio.h>
#include <stdlib.h>

int Write(double *X, double *Y, char *namefile, int size)
{
	FILE *file = NULL;
	file = fopen(namefile, "w");
	for (int i = 0; i < size; i++)
	{
		fprintf(file, "%f ; %f\n", X[i], Y[i]);
	}
	return 0;
}
