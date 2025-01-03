package io.github.nicolasyazman.maths;

import java.lang.Math;

public class Approximation {

	/**
	 * Since we know that 1 / (2^n (n+1)! ) * (pi / 4)n+1 <= 10⁻k
	 * With n the degree of the polynomial and k the number of accuracy digits we wish for.
	 */
	public static double[] ChebyshevTableSine = {1, 2, 3, 4, 5, 7, 9, 11, 12, 14, 15, 17, 18, 20, 22, 23, 25, 27, 28, 30, 32, 34};
	
	/**
	 * This function returns the Lambert W function using Taylor approximation.
	 * It can be used, for example, to k
	 * @param x
	 * @return
	 */
	public static double Lambert_W(double x) {
		// Let's return a Taylor expansion
		return x - x*x + 3/2 *x*x*x - 16/6 * x*x*x*x + 125/24 * x*x*x*x*x;	
	}
	
	/**
	 * Returns the number of nodes needed to get a certain accuracy.
	 * @param numberOfDigits Number of digits of accuracy.
	 * @return The number of Chebyshev nodes needed to reach a certain accuracy.
	 */
	public static int ChebyshevNumberOfNodes(int numberOfDigits) {
		for (int i = 0; i < ChebyshevTableSine.length; i++) {
			if (ChebyshevTableSine[i] > numberOfDigits) {
				return i;
			}
		}
		return ChebyshevTableSine.length;
	}
	
	/**
	 * Computes the degree of the Chebyshev Polynomial needed to reach a precision of 10^-k.
	 * @param k The power of the accuracy one wants to calculate the Chebyshev Polynomial for.
	 * @return The degree of the polynomial.
	 */
	public static double Chebyshev_Degree(int k) {
		double t = 4 / (Math.E * Math.PI) * Math.log(8 * Math.pow(10,2*k) / (Math.PI * Math.PI));
		System.out.println(t);
		System.out.println(Lambert_W(t));
		double n = Math.PI / 8 * Math.exp(1 + Lambert_W(t));
		return n;
	}
	
	public static int factorial(int n) {
		int res = 1;
		for (int i = 1; i <= n; i++) {
			res *= i;
		}
		return res;
	}
	/**
	 * Computes the error of the Chebyshev polynomial.
	 * @param The number of nodes.
	 * @return
	 */
	public static double Chebyshev_Error(int n) {
		return 1.0 / (Math.pow(2, n) * factorial(n + 1)) * Math.pow(Math.PI / 4, n+1);
	}
	
	 // Function to compute Chebyshev polynomial T_k(x)
    public static double chebyshevPolynomial(int k, double x) {
        if (k == 0) {
            return 1;
        } else if (k == 1) {
            return x;
        } else {
            double T0 = 1;
            double T1 = x;
            double T2 = 0;
            for (int i = 2; i <= k; i++) {
                T2 = 2 * x * T1 - T0;
                T0 = T1;
                T1 = T2;
            }
            return T2;
        }
    }
    
    public static void ComputeChebyshevCoeffsCos(int N) {
    	 // Number of Chebyshev nodes
        int n = 22;

        // Array to store Chebyshev coefficients
        double[] chebyshevCoeffs = new double[n];

        // Function f(x) = sin(x)
        // Chebyshev nodes
        double[] xNodes = new double[n];
        for (int i = 0; i < n; i++) {
            xNodes[i] = Math.cos((2 * i + 1) * Math.PI / (2 * n));
        }

        // Compute Chebyshev coefficients
        for (int k = 0; k < n; k++) {
            double sum = 0;
            for (int i = 0; i < n; i++) {
                double x = xNodes[i];
                double T_k_x = chebyshevPolynomial(k, x);
                sum += Math.cos(x) * T_k_x;
            }
            chebyshevCoeffs[k] = (2.0 / n) * sum;
        }

        // Print the Chebyshev coefficients
        System.out.println("Chebyshev Coefficients for cos(x) approximation:");
        for (int k = 0; k < n; k++) {
            System.out.println(chebyshevCoeffs[k]);
        }
        
    }
    
    
    // Method to calculate Chebyshev nodes in the interval [1, 2]
    public static double[] chebyshevNodes(int n, double a, double b) {
        double[] nodes = new double[n];
        for (int i = 0; i < n; i++) {
            // Chebyshev nodes in the interval [-1, 1]
            double x_i = Math.cos(Math.PI * (2 * (i + 1) - 1) / (2 * n));
            // Map nodes from [-1, 1] to [a, b]
            nodes[i] = 0.5 * (x_i + 1) * (b - a) + a;
        }
        return nodes;
    }

    // Method to calculate the natural logarithm of x
    public static double ln(double x) {
        return Math.log(x);
    }

    // Method to compute the Lagrange interpolation polynomial at point x
    public static double lagrangeInterpolation(double x, double[] nodes, double[] values) {
        double result = 0.0;
        int n = nodes.length;
        for (int i = 0; i < n; i++) {
            double term = values[i];
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    term *= (x - nodes[j]) / (nodes[i] - nodes[j]);
                }
            }
            result += term;
        }
        return result;
    }

    
    
	public static void ComputeChebyshevCoeffsSin(int numberOfNodes) {
		 // Number of Chebyshev nodes
        int n = 22;

        // Array to store Chebyshev coefficients
        double[] chebyshevCoeffs = new double[n];

        // Function f(x) = sin(x)
        // Chebyshev nodes
        double[] xNodes = new double[n];
        for (int i = 0; i < n; i++) {
            xNodes[i] = Math.cos((2 * i + 1) * Math.PI / (2 * n));
        }

        // Compute Chebyshev coefficients
        for (int k = 0; k < n; k++) {
            double sum = 0;
            for (int i = 0; i < n; i++) {
                double x = xNodes[i];
                double T_k_x = chebyshevPolynomial(k, x);
                sum += Math.sin(x) * T_k_x;
            }
            chebyshevCoeffs[k] = (2.0 / n) * sum;
        }

        // Print the Chebyshev coefficients
        System.out.println("Chebyshev Coefficients for sin(x) approximation:");
        for (int k = 0; k < n; k++) {
            System.out.println("c_" + k + " = " + chebyshevCoeffs[k]);
        }
        
        
	}
}
