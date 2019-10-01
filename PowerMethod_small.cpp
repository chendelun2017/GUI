#include <iostream>

struct Matrix {
	double* M;
	int n_rows, n_columns;
};

struct Vector {
	double* v;
	int n_rows;
};

void SetUpMatrix(Matrix &A, int rows, int columns) {
	A.n_rows = rows, A.n_columns = columns;
	A.M = (double*)calloc(rows * columns, sizeof(double));
}

void SetUpVector(Vector& a, int rows) {
	a.n_rows = rows;
	a.v = (double*)calloc(rows, sizeof(double));
}

void TransposeMatrix(Matrix A,Matrix &B) {
	B.n_rows = A.n_columns, B.n_columns = A.n_rows;
	B.M = (double*)calloc(A.n_columns * A.n_rows, sizeof(double));
	for (int i = 0; i < A.n_columns; i++) {
		for (int j = 0; j < A.n_rows; j++) B.M[j * A.n_rows + i] = 0.;
	}
	for (int i = 0; i < A.n_columns; i++) {
		for (int j = 0; j < A.n_rows; j++) B.M[j * A.n_columns + i] = A.M[i * A.n_rows + j]; 
	}
}

void MatrixProduct_Standard(Matrix A, Matrix B, Matrix &C) {
	C.n_rows = A.n_rows, C.n_columns = B.n_columns;
	C.M = (double*)calloc(C.n_rows * C.n_columns, sizeof(double));
	for (int i = 0; i < C.n_columns; i++) {
		for (int j = 0; j < C.n_rows; j++) C.M[i * C.n_rows + j] = 0.;
	}
	if (A.n_columns == B.n_rows) {
		for (int i = 0; i < C.n_columns; i++) {
			for (int j = 0; j < C.n_rows; j++) {
				for (int k = 0; k < A.n_columns; k++) C.M[i * C.n_rows + j] += A.M[k * A.n_rows + j] * B.M[i * B.n_rows + k];
			}
		}
	}
	else std::cout << "The dimensions of the matrices are not compatible for multiplying them!!!!  A matrix of zeros is going to be returned\n";
}

void MatrixProduct_XtransX(Matrix A, Matrix &B) {
	B.n_rows = A.n_columns, B.n_columns = A.n_columns;
	B.M = (double*)calloc(B.n_rows * B.n_columns, sizeof(double));
	for (int i = 0; i < B.n_columns; i++) {
		for (int j = 0; j < B.n_rows; j++) B.M[i * B.n_rows + j] = 0.;
	}
	for (int i = 0; i < B.n_columns; i++) {
		for (int j = 0; j < B.n_rows; j++) {
			for (int k = 0; k < A.n_rows; k++) B.M[i * B.n_rows + j] += A.M[i * A.n_rows + k] * A.M[j * A.n_rows + k];
		}
	}
}

void MatrixProduct_XXtrans(Matrix A, Matrix &B) {
	B.n_rows = A.n_rows, B.n_columns = A.n_rows;
	B.M = (double*)calloc(B.n_rows * B.n_columns, sizeof(double));
	for (int i = 0; i < B.n_columns; i++) {
		for (int j = 0; j < B.n_rows; j++) {
			for (int k = 0; k < A.n_columns; k++) B.M[i * B.n_rows + j] += A.M[k * A.n_rows+j] * A.M[k * A.n_rows+i];
		}
	}
}

void MatrixVectorProduct(Matrix A, Vector x, Vector& y) {
	y.n_rows = A.n_rows;
	y.v = (double*)calloc(x.n_rows, sizeof(double));
	for (int i = 0; i < x.n_rows; i++) y.v[i] = 0.;
	if (A.n_columns == x.n_rows) {
		for (int i = 0; i < A.n_rows; i++) {
			for (int j = 0; j < A.n_columns; j++) {
				y.v[i] += A.M[j * A.n_rows + i] * x.v[j];
			}
		}
	}
	else std::cout << "The dimensions of the matrix and the vector are not compatible!!!!  A vector of zeros is going to be returned\n";
}

void DisplayMatrix(Matrix A) {
	std::cout << "The matrix has " << A.n_rows << " rows and " << A.n_columns << " columns and the elements:\n";
	for (int i = 0; i < A.n_rows; i++) {
		for (int j = 0; j < A.n_columns; j++) { std::cout << A.M[j * A.n_rows + i] << "\t"; }
		std::cout << "\n";
	}
	std::cout << "\n";
}

void DisplayVector(Vector a){
	std::cout << "The vector has " << a.n_rows << " and the elements: \n";
	for (int i = 0; i < a.n_rows; i++) std::cout << a.v[i] << "\n";
	std::cout << "\n";
}

int main(void){
	Matrix X,Y,Z1,Z2;
	Vector a,b;

	SetUpMatrix(X, 2, 3);
	X.M[0] = 1, X.M[1] = 2, X.M[2] = 3;
	X.M[3] = 4, X.M[4] = 5, X.M[5] = 6;
	SetUpVector(a, 3);
	a.v[0] = 3, a.v[1] = 0, a.v[2] = 2;

	DisplayMatrix(X);
	TransposeMatrix(X,Y);
	DisplayMatrix(Y);
	MatrixProduct_Standard(X, Y, Z1);
	DisplayMatrix(Z1); 
	MatrixProduct_XtransX(Y, Z1);
	DisplayMatrix(Z1);
	MatrixProduct_XXtrans(Y, Z2);
	DisplayMatrix(Z2);
	MatrixVectorProduct(X, a, b);
	DisplayVector(b);

	free(X.M);
	free(Y.M);
	free(Z1.M);
	free(Z2.M);
}