#pragma once
#include <iostream>
#include <cstdlib>
#include <fstream>
using namespace std;

class Matrix {

public:
	double *coefficients;
	double *constTermsVec;
	double *L;
	double *U;
	double *Y;
	double *X;

	Matrix(double *coeffData, double *constTermsData, int numOfColumns, int numOfRows) {
		numOfEquations = numOfRows;
		numOfVariables = numOfColumns;

		coefficients = new double[numOfEquations * numOfVariables];
		constTermsVec = new double[numOfEquations * 1];

		CopyMatrix(coeffData, coefficients, numOfVariables, numOfEquations);
		CopyMatrix(constTermsData, constTermsVec, 1, numOfEquations);
	}
	~Matrix() {
		delete[] coefficients;
		coefficients = NULL;

		delete[] constTermsVec;
		constTermsVec = NULL;

		delete[] U;
		U = NULL;

		delete[] L;
		L = NULL;
	}

	void PrintMatrix(double *matrix, int numOfEquations, int numOfVariables) {
		ofstream file;
		file.open("result.txt");

		file.precision(10);

		for (int i = 0; i < numOfEquations; i++) {
			for (int j = 0; j < numOfVariables; j++) {
				file << matrix[i * numOfVariables + j] << "     ";
			}
			file << endl;
		}
		file.close();
	}
	void CalculateLU() {
		U = new double[numOfEquations * numOfVariables];
		L = new double[numOfEquations * numOfVariables];

		InitMatrixByZero(U, numOfVariables, numOfEquations);
		InitMatrixByZero(L, numOfVariables, numOfEquations);

        CopyRow(coefficients, U, numOfVariables, 0);
		
		for (int i = 0; i < numOfEquations; i++) {
			for (int j = 0; j < numOfVariables; j++) {
				if (i == j) {
					L[i * numOfVariables + j] = 1;
				}
			}
		}

		
		for (int i = 1; i < numOfEquations; i++) {
			for (int j = 0; j < 1; j++) {
				L[i * numOfVariables + j] = coefficients[i * numOfVariables + j] / U[j];
			}
		}

		for (int i = 1; i < numOfEquations; i++) {
			
			for (int j = i; j < numOfVariables; j++) {
                U[i * numOfVariables + j] = coefficients[i * numOfVariables + j];

				for (int k = 0; k < i; k++) {
					U[i * numOfVariables + j] -= L[i * numOfVariables + k] * U[k * numOfVariables + j];
				}
			}

			for (int j = i + 1; j < numOfVariables; j++) {
				L[j * numOfVariables + i] = coefficients[j * numOfVariables + i];
				
				for (int k = 0; k < i; k++) {
					L[j * numOfVariables + i] -= L[j * numOfVariables + k] * U[k * numOfVariables + i];
				}

				L[j * numOfVariables + i] /= U[i * numOfVariables + i];
			}
        }
	}

	void BackElimination(double *matrix, double *constTermsVec, int numOfVariables, int numOfEquations) {
		double *solutions = new double[numOfVariables];
		int numOfSolutions = numOfVariables - 1;
		int curRow = numOfEquations - 1;
		int curColumn = numOfVariables - 1;

		while (curRow >= 0) {
			double sum = 0;
			for (int i = curColumn + 1; i < numOfVariables; i++) {
				sum += matrix[curRow * numOfVariables + i] * solutions[i];
			}
			solutions[numOfSolutions] = (constTermsVec[curRow] - sum) / matrix[curRow * numOfVariables + curColumn];
			curRow--;
			curColumn--;
			numOfSolutions--;
		}
		X = solutions;
	}

	void ForvardElimination(double *matrix, double *constTermsVec, int numOfVariables, int numOfEquations) {
		double *solutions = new double[numOfVariables];
		int numOfSolutions = 0;
		int curRow = 0;
		int curColumn = 0;

		while (curRow < numOfEquations) {
			double sum = 0;
			for (int i = 0; i < curColumn; i++) {
				sum += matrix[curRow * numOfVariables + i] * solutions[i];
			}
			solutions[numOfSolutions] = (constTermsVec[curRow] - sum) / matrix[curRow * numOfVariables + curColumn];
			curRow++;
			curColumn++;
			numOfSolutions++; 
		}
		Y = solutions;
	}
private:
	int numOfVariables;
	int numOfEquations;
	

	void SwapRows(double *matrix, int from, int to) {
		for (int j = 0; j < numOfVariables; j++) {
			Swap(&matrix[from + j], &matrix[to + j]);
		}
	}
	void Swap(double *a, double *b) {
		double tmp = *a;
		*a = *b;
		*b = tmp;
	}
	void CopyMatrix(double *from, double *to, int numOfVariables, int numOfEquations) {
		for (int i = 0; i < numOfEquations; i++) {
			for (int j = 0; j < numOfVariables; j++) {
				to[i * numOfVariables + j] = from[i * numOfVariables + j];
			}
		}
	}
	void CopyRow(double *from, double *to, int numOfVariables, int row) {
		for (int j = 0; j < numOfVariables; j++) {
			to[row * numOfVariables + j] = from[row * numOfVariables + j];
		}
	}
	void InitMatrixByZero(double *matrix, int numOfVariables, int numOfEquations) {
		for (int i = 0; i < numOfEquations; i++) {
			for (int j = 0; j < numOfVariables; j++) {
				matrix[i * numOfVariables + j] = 0;
			}
		}
	}
};