#include <iostream>
#include <fstream>
#include "matrix.h"

using namespace std;

double EuclideanNorm(double *vector, int lenght) {
	double sum = 0;
	for (int i = 0; i < lenght; i++) {
		sum += powl(vector[i], 2);
	}
	sum = sqrt(sum);
	return sum;
}

double* Check(double *X, double *A, int numOfVariables, int numOfEqiations) {
	double *res = new double[numOfEqiations];

	for (int i = 0; i < numOfEqiations; i++) {
		res[i] = 0;
	}

	for (int i = 0; i < numOfEqiations; i++) {
		for (int j = 0; j < numOfVariables; j++) {
			res[i] += X[j] * A[i * numOfVariables + j];
		}
	}

	return res;
}

int main(void) {

	ifstream file;
	file.open("matr1.txt");

	double *dataA, *dataB;
	int numOfVariables, numOfEquations;

	file >> numOfEquations >> numOfVariables;

	dataA = new double[numOfVariables * numOfEquations];
	dataB = new double[1 * numOfEquations];

	//считыванием матрицу коэффициентов
	for (int i = 0; i < numOfEquations; i++) {
		for (int j = 0; j < numOfVariables; j++) {
			file >> dataA[i * numOfVariables + j];
		}
	}

	//считываем вектор правых частей
	for (int i = 0; i < numOfEquations; i++) {
		for (int j = 0; j < 1; j++) {
			file >> dataB[i * 1 + j];
		}
	}

	Matrix matrix = Matrix(dataA, dataB, numOfVariables, numOfEquations);

	matrix.CalculateLU();
	matrix.ForvardElimination(matrix.L, matrix.constTermsVec, numOfVariables, numOfEquations);
	matrix.BackElimination(matrix.U, matrix.Y, numOfVariables, numOfEquations);
	//matrix.PrintMatrix(matrix.X, numOfEquations, 1);
	//cout << "----------------------------------------------------------" << endl;

	double *Diff = Check(matrix.X, dataA, numOfVariables, numOfEquations);
	matrix.PrintMatrix(Diff, numOfEquations, 1);

	return 0;

}