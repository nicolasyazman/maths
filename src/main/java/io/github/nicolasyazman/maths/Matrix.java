package io.github.nicolasyazman.maths;

public class Matrix {

	/**
	 * The number of rows in the Matrix.
	 */
	private int rows;
	
	/**
	 * The number of cols in the Matrix.
	 */
	private int cols;
	
	/**
	 * The values contained inside the matrix.
	 */
	private double[][] elements;
	
	/**
	 * 
	 * @param numberOfRows
	 * @param numberOfColumns
	 * @param values
	 */
	public Matrix(double[][] values) {
		if (values == null) {
			throw new IllegalArgumentException("Trying to create a Matrix with null values.");
		}	
	}
	
	/**
	 * Matrix multiplication. Performs the multiplication A * B.
	 * It is assumed the sizes of number of columns for each row is the same.
	 * @param a The first Matrix
	 * @param b The second Matrix
	 * @return
	 */
	public static double[][] MatMul(double[][] a, double[][] b) {
		
		// Sanity check that verifies 
		if ((a == null) || (b == null)) {
			throw new IllegalArgumentException("In the matrix multiplication, A or B is null ");
		}
		int m = a.length;
		int n1 = a[0].length;
		int n2 = b.length;
		int p = b[0].length;
		
		// A is a matrix [m, n1] with m the number of rows, and n1 the number of cols.
		// B is a matrix [n2, p] with n2 the number of rows and p the number of cols.
		// To be able to multiply matrices, the number of columns of the first matrix should be equal
		// to the number of rows of the second.
		if (n1 != n2) {
			throw new IllegalArgumentException("The number of cols of the first matrix must be equal to the number of rows of the second matrix.");
		}
		
		double res[][] = new double[m][];
		for (int i = 0; i < m; i++) {
			res[i] = new double[p];
		}
		
		// For each resulting matrix row.
		for (int i = 0; i < m; i++) {
			// For each resulting matrix col.
			for (int j = 0; j < p;j++) {
				double val = 0.0;
				for (int k = 0; k < n1; k++) {
					val += a[i][k] * b[k][j];
				}
				res[i][j] = val;
			}
		}
		
		return res;
	}
}
