#include <iostream>
#include <iomanip> 
#include <stdlib.h>


struct Matrix {
	double* M;
	int n_rows, n_columns;
};

struct Vector {
	double* v;
	int n_rows;
};

void DisplayMatrix(Matrix A) {			// setprecision:控制输出流显示浮点数的数字个数，setprecision(n)就是输出的n个数，会有四舍五入
	std::cout << std::fixed;
	std::cout << "The matrix has " << A.n_rows << " rows and " << A.n_columns << " columns and the elements:\n";
	for (int i = 0; i < A.n_rows; i++) {
		for (int j = 0; j < A.n_columns; j++) { std::cout << std::setprecision(6) << A.M[j * A.n_rows + i] << "\t"; }
		std::cout << "\n";
	}
	std::cout << "\n";
}

void DisplayVector(Vector a) {
	std::cout << std::fixed;
	std::cout << "The vector has " << a.n_rows << " and the elements: \n";
	for (int i = 0; i < a.n_rows; i++) std::cout << std::setprecision(4) << a.v[i] << "\n";
	std::cout << "\n";
}

void GetRandomMatrix(Matrix& A) {		// 生成矩阵元素值为(0,1)之间的随机矩阵
	for (int i = 0; i < A.n_columns; i++) {
 		for (int j = 0; j < A.n_rows; j++) A.M[i * A.n_rows + j] = ((double)rand() - (double)(RAND_MAX / 2)) / (double)RAND_MAX;;
	}
}

void GetRandomVector(Vector& a) {
	for (int i = 0; i < a.n_rows; i++) a.v[i] = ((double)rand() - (double)(RAND_MAX / 2)) / (double)RAND_MAX;
}

void SetUpMatrix(Matrix& A, int rows, int columns) {
	A.n_rows = rows, A.n_columns = columns;
	A.M = (double*)calloc(rows * columns, sizeof(double));
}

void SetUpVector(Vector& a, int rows) {
	a.n_rows = rows;
	a.v = (double*)calloc(rows, sizeof(double));
}

void TransposeMatrix(Matrix A,Matrix& B) {
	for (int i = 0; i < A.n_columns; i++) {
		for (int j = 0; j < A.n_rows; j++) B.M[j * A.n_columns + i] = A.M[i * A.n_rows + j]; 
	}
}

void MatrixProduct_Standard(Matrix A, Matrix B, Matrix& C) {
	if (A.n_columns == B.n_rows) {
		for (int i = 0; i < C.n_columns; i++) {
			for (int j = 0; j < C.n_rows; j++) {
				C.M[i * C.n_rows + j] = 0;
				for (int k = 0; k < A.n_columns; k++)  C.M[i * C.n_rows + j] += A.M[k * A.n_rows + j] * B.M[i * B.n_rows + k];
			}
		}
	}
	else std::cout << "The dimensions of the matrices are not compatible for multiplying them!!!!  A matrix of zeros is going to be returned\n";
}

void MatrixProduct_Partial(Matrix A, Matrix B, Matrix& C, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < A.n_rows; j++) {
			C.M[i*C.n_rows+j] = 0;
			for (int k = 0; k < A.n_columns; k++) C.M[i * C.n_rows + j] += A.M[k * A.n_rows + j] * B.M[i*B.n_rows+k];
		}
	}
}

void MatrixProduct_ABtrans(Matrix A, Matrix B, Matrix& C, int n) {
	for (int i = 0; i < A.n_rows; i++) {
		for (int j = 0; j < B.n_rows; j++) {
			C.M[j * C.n_rows + i] = 0;
			for (int k = 0; k < n; k++) {
				C.M[j * C.n_rows + i] += A.M[k * A.n_rows + i] * B.M[k * B.n_rows + j];
			}
		}
	}
}

void MatrixProduct_XtransX(Matrix A, Matrix& B) {
	for (int i = 0; i < B.n_columns; i++) {
		for (int j = 0; j < B.n_rows; j++) {
			if (j >= i) {
				B.M[i * B.n_rows + j] = 0.;
				for (int k = 0; k < A.n_rows; k++) B.M[i * B.n_rows + j] += A.M[i * A.n_rows + k] * A.M[j * A.n_rows + k];
			}
			else {
				B.M[i * B.n_rows + j] = B.M[j * B.n_rows + i];
			}
		}
	}
}

void MatrixProduct_XXtrans(Matrix A, Matrix& B) {
	for (int i = 0; i < B.n_columns; i++) {
		for (int j = 0; j < B.n_rows; j++) {
			B.M[i * B.n_rows + j] = 0.;
			for (int k = 0; k < A.n_columns; k++) B.M[i * B.n_rows + j] += A.M[k * A.n_rows+j] * A.M[k * A.n_rows+i];
		}
	}
}

void MatrixVectorProduct(Matrix A, Vector x, Vector& y) {
	if (A.n_columns == x.n_rows) {
		for (int i = 0; i < A.n_rows; i++) {
			y.v[i] = 0.;
			for (int j = 0; j < A.n_columns; j++) {
				y.v[i] += A.M[j * A.n_rows + i] * x.v[j];
			}
		}
	}
	else std::cout << "The dimensions of the matrix and the vector are not compatible!!!!  A vector of zeros is going to be returned\n";
}
//sqrt(a*a + b*b + c*c)
void GetVectorLength(Vector a, double& b) {
	b = 0.;
	for (int i = 0; i < a.n_rows; i++) b += a.v[i] * a.v[i];
	b = sqrt(b);
}

void ScaleVector(Vector& a, double b) {
	for (int i = 0; i < a.n_rows; i++) a.v[i] /= b;
}

void GetVectorDifferenceNorm(Vector a, Vector b, double& c) {
	c = 0.;
	for (int i = 0; i < a.n_rows; i++) c += (a.v[i] - b.v[i]) * (a.v[i] - b.v[i]);
	c = sqrt(c);
}

void ChangeVectorOrder(Vector& a, Vector b){
	for (int i = 0; i < a.n_rows; i++) a.v[i] = b.v[i];
}

void StoreVectorInMatrix(Matrix& A, Vector a, int k) {
	for (int i = 0; i < a.n_rows; i++) A.M[k * a.n_rows + i] = a.v[i];
}

void VectorDifference(Vector a, Vector b, Vector& c) {
	for (int i = 0; i < a.n_rows; i++) c.v[i] = a.v[i] - b.v[i];
}

void DeflateMatrix(Matrix& A, Vector a, double b) {
	for (int i = 0; i < A.n_columns; i++) {
		for (int j = 0; j < A.n_rows; j++) A.M[i * A.n_rows + j] -= b * a.v[j] * a.v[i];
	}
}

//A 要特征分解的矩阵，s 前n个最大的奇异值，P 是载荷矩阵
void GetEigenDecomposition(Matrix A,Vector& s,Matrix& P, int n){
	Vector x0, x1;
	double error = 100.;
	x0.n_rows = A.n_rows, x1.n_rows = A.n_rows;
	x0.v = (double*)calloc(x0.n_rows, sizeof(double)), x1.v = (double*)calloc(x1.n_rows, sizeof(double));
	for (int k = 0; k < n; k++) {
		GetRandomVector(x0);
		GetVectorLength(x0, s.v[k]);
		ScaleVector(x0, s.v[k]);
		do {
			MatrixVectorProduct(A, x0, x1);			//矩阵向量乘法，x1 = A*x0
			GetVectorLength(x1, s.v[k]);			//得到 x1 的向量长度
			ScaleVector(x1, s.v[k]);				//将 x1 化为单位向量
			GetVectorDifferenceNorm(x1, x0, error);	//判断 x1,x0 之间的误差
			ChangeVectorOrder(x0, x1);				//把 x1 的值赋予 x0
		} while (error > 1e-8);

		StoreVectorInMatrix(P, x1, k);				//将第 k 列载荷向量存入矩阵 P 中
		DeflateMatrix(A, x1, s.v[k]);				//缩小矩阵
	}
	free(x0.v), free(x1.v);
	x0.v= NULL , x1.v = NULL;
}

void GetRowVector(Matrix A, Vector& a, int j) {
	for (int i = 0; i < A.n_columns; i++) a.v[i] = A.M[i * A.n_rows + j];
}

void GetColumnVector(Matrix A, Vector& a, int i) {
	for (int j = 0; j < A.n_rows; j++) a.v[j] = A.M[i*A.n_rows+j];
}

//i 是 Z.n_colunms
//j 是 第 m 份数据
void SegmentDataMatrix(Matrix Z, Matrix& X_i_j, Vector& Xi_j, Matrix& X_ij, Vector Xij, int i, int j) {
	int ell = Xij.n_rows, kj_counter, ki_counter = 0, m = (int)(Z.n_rows/ell);	//ell 是每份训练数据的大小，一共 m 份
	for (int ki = 0; ki < Z.n_columns; ki++) {// ki 为z的列数据中做循环
		if (ki == i) {			//判断 ki列 是否是 y 的数据(ytrain,ytest)
			kj_counter = 0;
			for (int kj = 0; kj < m; kj++) {// kj 为在份数 m 中做循环
				if (kj == j) {	//在 ki列代表 y 的情况下 判断 kj份 是否是第 m 份数据，是则代表是 ytest
					for (int k = 0; k < ell; k++) Xij.v[k] = Z.M[ki * Z.n_rows + kj * ell + k];
				}
				else {			//如果 kj份 不是第 m 份数据，则代表为 ytrain
					for (int k = 0; k < ell; k++) Xi_j.v[kj_counter * ell + k] = Z.M[ki * Z.n_rows + kj * ell + k];
					kj_counter++;
				}
			}
		}
		else {					//否则 ki列 是 x 的数据(xtrain, xtest)
			kj_counter = 0;
			for (int kj = 0; kj < m; kj++) {
				if (kj == j) {	//在 ki列代表 x 的情况下 判断 kj份 是否是第 m 份数据，是则代表是 xtest
					for (int k = 0; k < ell; k++) X_ij.M[ki_counter * ell + k] = Z.M[ki * Z.n_rows + kj * ell + k];
				}
				else {			//如果 kj份 不是第 m 份数据，则代表为 xtrain
					for (int k = 0; k < ell; k++) {
						X_i_j.M[ki_counter * X_i_j.n_rows + kj_counter * ell + k] = Z.M[ki * Z.n_rows + kj * ell + k];
					}
					kj_counter++;
				}
			}
		ki_counter++;
		}
	}
}

double GetKernelMatrixElement(Matrix A, int i, int j, double sigma) {
	Vector xi, xj;
	double kij = 0;
	
	SetUpVector(xi, A.n_columns);
	SetUpVector(xj, A.n_columns);
	GetRowVector(A, xi, i);
	GetRowVector(A, xj, j);

	for (int i = 0; i < A.n_columns; i++) {
		kij += (xi.v[i] - xj.v[i]) * (xi.v[i] - xj.v[i]);	//这里 kij 就是 omega 中 第 ij 位置元素
	}
	free(xi.v), free(xj.v);
	kij = kij / sigma;		
	kij = exp(-kij);										//这里 kij 就是 Ktrain或ktest 中 第 ij 位置元素

	return kij;
}
//A 是 xtrain
void GetGramMatrixForTraining(Matrix A, Vector& Mean, Matrix& G, double& mean, double Sigma) {
	Matrix K;
	int elements = A.n_rows;		//elements 是 xtrain 的 n_rows
	SetUpMatrix(K, elements, elements);

	// Constructing the kernel matrix kij = exp ( - || xi - xj || / Sigma )
	for (int i = 0; i < elements; i++) {
		for (int j = 0; j < elements; j++) {
			if (j >= i) {	//因为 矩阵K 也是对称的，所以只需计算 j >= i 的这部分即可
				K.M[i * elements + j] = GetKernelMatrixElement(A, i, j, Sigma);
			}
			else {
				K.M[i * elements + j] = K.M[j * elements + i];
			}
			G.M[i * elements + j] = K.M[i * elements + j];
		}
	}
	// Mean centering the kernel matrix to product the Gram matrix
	mean = 0;	//mean是Ktrain中所有元素值的平均值，对应matlab中求 G 的OneTest * KTrain *OneTrain
	for (int i = 0; i < elements; i++) {
		for (int j = i; j < elements; j++) {
			if (i == j) mean += K.M[i * elements + j];
			else mean += 2 * K.M[i * elements + j];
		}
	}
	mean /= ( (double)elements * (double)elements );
	for (int i = 0; i < elements; i++) {
		Mean.v[i] = 0;
		for (int j = 0; j < elements; j++) Mean.v[i] += K.M[i * elements + j];
		Mean.v[i] /= elements;//Mean 是Ktrain每行的平均值,对应matlab中求 G 的 OneTest * KTrain 或 KTest * OneTrain
	}
	for (int i = 0; i < elements; i++) {
		for (int j = 0; j < elements; j++) {
			G.M[i * elements + j] -= Mean.v[i];		//Kruger NonlinearPCA 1.62
			G.M[i * elements + j] -= Mean.v[j];
			G.M[i * elements + j] += mean;
		}
	}
	free(K.M);
}
//Atrain 是 Xtrain, Atest 是 Xtest, G 是 Gtest
void GetGramMatrixForTesting(Matrix Atrain, Matrix Atest, Vector Mean, double mean, Matrix& G, double Sigma) {
	Matrix Ktest;
	Vector xi,xj,Temp;
	int rows = Atest.n_rows, columns = Atrain.n_rows, jIndex = 0;
	SetUpMatrix(Ktest,rows,columns);	//Ktest
	SetUpVector(xi, Atrain.n_columns);	//xi 是Xtrain的列数
	SetUpVector(xj, Atrain.n_columns);	//xj 是Xtrain的列数
	SetUpVector(Temp, rows);			//Temp是Xtest的行数

	// Constructing kernel matrix ki = exp ( || xi_test - xj_train || / Sigma )
	for (int i = 0; i < rows; i++) {
		GetRowVector(Atest, xi, i);
		for (int j = 0; j < columns; j++) {
			GetRowVector(Atrain, xj, j);
			Ktest.M[j * rows + i] = 0;
			for (int k = 0; k < Atrain.n_columns; k++) Ktest.M[j * rows + i] += (xi.v[k] - xj.v[k]) * (xi.v[k] - xj.v[k]);
			Ktest.M[j * rows + i] /= Sigma;
			Ktest.M[j * rows + i] = exp(-Ktest.M[j * rows + i]);
			G.M[j * rows + i] = Ktest.M[j * rows + i];
		}
	}
	// Mean centering the kernel matrix to produce the Gram matrix
	for (int i = 0; i < rows; i++) {
		Temp.v[i] = 0;
		for (int j = 0; j < columns; j++) Temp.v[i] += Ktest.M[j * rows + i];	
		Temp.v[i] /= columns;		//Temp存储Ktest每行的平均值
	}
	for (int i = 0; i < columns; i++) {
		for (int j = 0; j < rows; j++) {
			G.M[i * rows + j] -= Mean.v[i];
			G.M[i * rows + j] -= Temp.v[j];
			G.M[i * rows + j] += mean;
		}
	}
	free(Ktest.M), free(xi.v), free(xj.v), free(Temp.v);
}

void StoreSymmetricMatrix(Matrix& G, Vector &g) {
	int rows = G.n_rows, flag = 0;
	for (int col = 0; col < rows; col++) {
		for (int row = col; row < rows; row++) {
			g.v[flag] = G.M[col * rows + row];
			flag++;
		}
	}
	free(G.M);
	G.M = NULL;
}
Matrix GetFullSymmetricMatrix(Vector& g) {
	int n = g.n_rows;
	int l = (sqrt(8 * n + 1) - 1) / 2;
	Matrix G;
	SetUpMatrix(G, l, l);

	for (int row = 0; row < l; row++) {
		for (int col = 0; col < l; col++) {
			if (row <= col) {
				G.M[row * l + col] = g.v[row * (l + l - (row - 1)) / 2 + col - row];
			}
			else{
				G.M[row * l + col] = g.v[col * (l + l - (col - 1)) / 2 + row - col];
			}
		}
	}
	return G;
}

int ReadSymmetricElement(Vector &g ,int i, int j) {
	int temp = g.n_rows;
	int N = (sqrt(8 * temp + 1) - 1) / 2;
	if (i <= j) return g.v[i * N - (i - 1) * i / 2 + j - i];
	else return g.v[j * N - (j - 1) * j / 2 + i - j];
}

Vector ReadSymmetricColumns(Vector& g, int col) {
	int temp = g.n_rows;
	int N = (sqrt(8 * temp + 1) - 1) / 2;
	Vector g_col;
	SetUpVector(g_col, N);
	for (int row = 0; row < N; row++) {
		g_col.v[row] = ReadSymmetricElement(g, row, col);
	}
	return g_col;
}

Vector ReadSymmetricRows(Vector& g, int row) {
	int temp = g.n_rows;
	int N = (sqrt(8 * temp + 1) - 1) / 2;
	Vector g_row;
	SetUpVector(g_row, N);
	for (int col = 0; col < N; col++) {
		g_row.v[col] = ReadSymmetricElement(g, row, col);
	}
	return g_row;
}

Matrix SymmetricAtimesB(Vector& g, Matrix& B) {
	int temp = g.n_rows;
	int N = (sqrt(8 * temp + 1) - 1) / 2;
	
	if (N == B.n_rows) {
		Matrix A = GetFullSymmetricMatrix(g);
		Matrix C;
		SetUpMatrix(C, N, B.n_columns);
		MatrixProduct_Standard(A, B, C);
		return C;
	}
	else std::cout << "The dimensions of the matrices are not compatible for multiplying them!!!!  A matrix of zeros is going to be returned\n";
}

Matrix AtimesSymmetricB(Matrix& A, Vector& g) {
	int temp = g.n_rows;
	int N = (sqrt(8 * temp + 1) - 1) / 2;

	if (A.n_columns == N) {
		Matrix B = GetFullSymmetricMatrix(g);
		Matrix C;
		SetUpMatrix(C, A.n_rows, N);
		MatrixProduct_Standard(A, B, C);
		return C;
	}
	else std::cout << "The dimensions of the matrices are not compatible for multiplying them!!!!  A matrix of zeros is going to be returned\n";
}

Vector SymmetricAtimesb(Vector& g, Vector& b) {
	int temp = g.n_rows;
	int N = (sqrt(8 * temp + 1) - 1) / 2;
	Matrix A = GetFullSymmetricMatrix(g);
	Vector y;
	SetUpVector(y, N);
	MatrixVectorProduct(A, b, y);
	return y;
}

void KernelPrincipalComponentAnalysis(Matrix Z, Matrix& P, Matrix& T, Vector& s, int n) {
	Matrix X_ij, X_i_j, G_i_j, Gi_j, V, Temp_Matrix_1, Temp_Matrix_2, Temp_Matrix_3, Ree;
	Vector xij, xi_j, xij_pred, e, Mean, lambda, see;
	double sigma = 0.5, ds = 0.1, sigma_end = 4, Sigma, mean, length_e;
	int m = 10, k = 0, ell = (int)(Z.n_rows/m), counter=0;		//ell 是每份训练数据的大小，一共 m 份
	
	// Note that Sigma = sigma * sigma to avoid unneccessary multiplications when setting up Kernel matrices
	SetUpMatrix(X_ij, ell, Z.n_columns-1);						//X_ij 是 Xte_st
	SetUpMatrix(X_i_j, Z.n_rows - ell, Z.n_columns - 1);		//X_i_j 是 Xtrain
	SetUpMatrix(G_i_j, Z.n_rows - ell, Z.n_rows - ell);			//G_i_j 是 G_train
	SetUpMatrix(Gi_j, ell, Z.n_rows - ell);						//Gi_j 是 G_te_st
	SetUpMatrix(V, Z.n_rows - ell, n);							//V 是 V_train
	SetUpMatrix(Ree, n, 1 + (int)((sigma_end - sigma) / ds));	//Ree 是不同sigma下的前n个主成分的预测误差
	SetUpVector(xij, ell);										//xij 是 yte_st
	SetUpVector(xij_pred, ell);									//xij_pred 是 ypred
	SetUpVector(e, ell);										//e = xij - xij_pred
	SetUpVector(xi_j, Z.n_rows - ell);							//xi_j 是 ytrain
	SetUpVector(Mean, Z.n_rows - ell);
	SetUpVector(lambda, n);										//n是主成分，lambda存取前n个最大特征值
	SetUpVector(see, n);										//see 是 Ree 的某一列

	// Loop to find optimal Sigma based on a simple grid search
	for (int i = 0; i < n; i++) see.v[i] = 1000;
	do {
		Sigma = sigma * sigma;
		for (int j = 0; j < m; j++) {
			for (int i = 0; i < Z.n_columns; i++) {
				SegmentDataMatrix(Z, X_i_j, xi_j, X_ij, xij, i, j);
				GetGramMatrixForTraining(X_i_j, Mean, G_i_j, mean, Sigma);
				GetGramMatrixForTesting(X_i_j, X_ij, Mean, mean, Gi_j, Sigma);
				GetEigenDecomposition(G_i_j, lambda, V, n);				//V 是loading载荷矩阵
				for (int k = 0; k < n; k++) {
					SetUpMatrix(Temp_Matrix_1, Gi_j.n_rows, k+1);		//Temp_Matrix_1 是Gtest的行，k+1主成分的列
					MatrixProduct_Partial(Gi_j, V, Temp_Matrix_1, k+1);	//前 k+1 行列相乘
					for (int li = 0; li < k; li++) {
						for (int lj = 0; lj < Temp_Matrix_1.n_rows; lj++) Temp_Matrix_1.M[li * Temp_Matrix_1.n_rows + lj] /= ( lambda.v[li] * lambda.v[li] );
					}
					SetUpMatrix(Temp_Matrix_2, ell, Z.n_rows - ell);	//Temp_Matrix_2  （10*90）
					MatrixProduct_ABtrans(Temp_Matrix_1, V, Temp_Matrix_2, k+1);
					SetUpMatrix(Temp_Matrix_3, ell, Z.n_rows - ell);	//Temp_Matrix_3   (10*90)
					MatrixProduct_Standard(Temp_Matrix_2, G_i_j, Temp_Matrix_3);
					MatrixVectorProduct(Temp_Matrix_3, xi_j, xij_pred);	//xij_pred        (10*1)
					VectorDifference(xij, xij_pred, e);
					GetVectorLength(e,length_e);
					Ree.M[counter*Ree.n_rows+k] += length_e * length_e;
				}
			}
		}
		sigma += ds;
		for (int i = 0; i < n; i++) {
			if ((Ree.M[counter * Ree.n_rows + i] / ((double)Z.n_rows * (double)Z.n_columns)) < see.v[i]) see.v[i] = ( Ree.M[counter * Ree.n_rows + i] / ((double)Z.n_rows * (double)Z.n_columns));
		}
		counter++;
		std::cout << std::setprecision(5) << "I have just checked : Sigma = " << Sigma << "\t see[1] = " << see.v[1] << "\t see[2] = " << see.v[2] << "\t see[3] = " << see.v[3] << "\n";
	} while (sigma <= sigma_end);
}

int main(void){

	Matrix G_;
	Vector g;
	int rows = 4;  										// rows 假设生成对称矩阵的行, i 是恢复的行，j是恢复的列

	SetUpMatrix(G_, rows, rows);
	SetUpVector(g, rows * (rows + 1) / 2);
	for (int i = 0; i < rows * rows; i++) G_.M[i] = i;	// 生成一个对称矩阵
	StoreSymmetricMatrix(G_, g);						// Already free(G.M)
	Matrix G = GetFullSymmetricMatrix(g);				// Here we have a symmetrci matrix from vector g as an example
	std::cout << "Assuming we have a symmetric matrix G" << std::endl;
	DisplayMatrix(G);									// Show symmetric matrix

	int i = 2, j = 3;									// i represent the row, j represent the columns of symmetric matrix G
	std::cout << "The G's element of (" << i << ',' << j << ") is " << ReadSymmetricElement(g, i, j) << '\n' << std::endl; //read element of G
	std::cout << "Now we see the (column "<< j << ") of G " << std::endl;				// read column of G	
	DisplayVector(ReadSymmetricColumns(g, j));

	std::cout << "Now we see the (row " << i << ") of G " << std::endl;				// read row of G	
	DisplayVector(ReadSymmetricRows(g, i));

	std::cout << "Now we see the matrix B" << std::endl;
	Matrix B;											// Matrix B is used to: G times B
	int cols = 3;
	SetUpMatrix(B, rows, cols);
	for (int i = 0; i < rows * cols; i++) B.M[i] = i;
	DisplayMatrix(B);
	std::cout << "Now we see the result of symmetric matrix (vector g) times matrix B, we can check the rusult on Matlab" << std::endl;
	Matrix T = SymmetricAtimesB(g, B);
	DisplayMatrix(T);

	std::cout << "Now we see the matrix A" << std::endl;
	Matrix A;											// Matrix B is used to: G times B
	cols = 5;
	SetUpMatrix(A, cols, rows);
	for (int i = 0; i < rows * cols; i++) A.M[i] = i;
	DisplayMatrix(A);
	std::cout << "Now we see the result of matrix A times symmetric matrix (vector g), we can check the rusult on Matlab too" << std::endl;
	T = AtimesSymmetricB(A, g);
	DisplayMatrix(T);

	std::cout << "How about G times vector b" << std::endl;
	Vector b;
	SetUpVector(b, rows);
	for (int i = 0; i < rows; i++) b.v[i] = i;
	DisplayVector(b);
	Vector a = SymmetricAtimesb(g, b);
	std::cout << "Now we see the result of symmetric matrix (vector g) times vextor b, we can check the rusult on Matlab too" << std::endl;
	DisplayVector(a);
	
	//Matrix Z, P, T;
	//Vector s;
	//int rows = 100, columns = 3, n = 3;	//此时的rows就是matlab中的npts，n是主成分

	//srand(1); //随机数种子
	//SetUpMatrix(Z, rows, columns);
	//SetUpMatrix(P, columns, n);
	//SetUpVector(s, n);
	//GetRandomMatrix(Z);
	//DisplayMatrix(Z);
	//KernelPrincipalComponentAnalysis(Z,P,T,s,n);
}