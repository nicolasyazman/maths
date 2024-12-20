package io.github.nicolasyazman;

public class Trigonometry {

	private static double SIN_PI_4 = 0.70710678118654752440084436210484903928483593768847403658833986899536623923;
	private static double COS_PI_4 = 0.70710678118654752440084436210484903928483593768847403658833986899536623923;
	/**
	 * These factorial are useful for the calculation of sine function. No need to compute them each time.
	 * Since we only go to the sixth order of the Taylor series no need for more.
	 */
	private static double FACT3 = 1*2*3;
	private static double FACT5 = FACT3*4*5;
	private static double FACT7 = FACT5*6*7;
	private static double FACT9 = FACT7*8*9;
	private static double FACT11 = FACT9*10*11;
	private static double FACT13 = FACT11*12*13;
	
	/**
	 * These factorial are useful for the calculation of the cosine function. No need to compute them each time.
	 * Since we only go to the sixth order of the Taylor series, no need for more.
	 */
	private static double FACT2 = 1*2;
	private static double FACT4 = 1*2*3*4;
	private static double FACT6 = 1*2*3*4*5*6;
	private static double FACT8 = 1*2*3*4*5*6*7*8;
	private static double FACT10 = 1*2*3*4*5*6*7*8*9*10;
	private static double FACT12 = 1*2*3*4*5*6*7*8*9*10*11*12;
	
	public static double PI = 3.141592653589793238462643383279502884197;
	
	
	private static double CHEBYSHEV_COEFFS_SIN[] = {
			5.046468293750712E-17,
			0.880101171489867,
			7.065055611250996E-17,
			-0.039126707965337015,
			2.321375415125327E-16,
			4.995154604219994E-4,
			2.9269516103754126E-16,
			-3.0046516354560864E-6,
			2.220446049250313E-16,
			1.0498499584828251E-8,
			2.0185873175002848E-16,
			-2.3960883782143067E-11,
			9.083642928751281E-17,
			3.805037093488037E-14,
			-3.5325278056254984E-16,
			-3.482063122687991E-16,
			-5.702509171938304E-16,
			-6.863196879500968E-16,
			-3.3054367324067163E-16,
			-1.0875139173032783E-15,
			-1.6779507076721117E-16,
			-2.115731832154986E-15
	};
	
	private static double[] CHEBYSHEV_COEFFS_COS = {
			1.530395373115933,
			-2.0185873175002847E-17,
			-0.22980696986380092,
			1.312081756375185E-16,
			0.0049532779282195705,
			1.4130111222501992E-16,
			-4.187667600510673E-5,
			1.8671932686877633E-16,
			1.8844688304264778E-7,
			8.074349270001139E-17,
			-5.261236042790204E-10,
			2.0185873175002848E-16,
			9.994984637919722E-13,
			-5.046468293750712E-18,
			-1.4685222734814572E-15,
			-4.794144879063176E-16,
			-2.296143073656574E-16,
			-3.835315903250541E-16,
			-7.09028795271975E-16,
			-3.027880976250427E-17,
			-5.412337245047638E-16,
			-8.894400367735629E-17
        };
	/**
	 * A storage for the Taylor polynomials
	 */
	private static double pows[] = new double[6];
	
	/**
	 * Power function.
	 * @param x A value
	 * @param p The power of the value
	 * @return x^p
	 */
	public static double pow(double x, int p) {
		double res = 1;		
		while (p < 0) {
			res /= x;
			p++;
		}
		while (p > 0) {
			res *= x;
			p--;
		}
		return res;
	}
	
	 // Function to compute the Chebyshev polynomial T_k(x)
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

    // Function to compute the Chebyshev interpolation
    public static double chebyshevInterpolation(double x, double[] chebyshevCoeffs) {
        double result = 0;
        int n = chebyshevCoeffs.length;

        // Sum the terms of the Chebyshev series
        for (int k = 0; k < n; k++) {
            result += chebyshevCoeffs[k] * chebyshevPolynomial(k, x);
        }

        return result;
    }
    
 // Method to evaluate the Chebyshev series for cos(x)
    public static double chebyshevInterpolationcos(double x, double[] coeffs) {
        int n = coeffs.length;

        // Initialize T_0(x) and T_1(x)
        double T0 = 1.0;  // T_0(x)
        double T1 = x;    // T_1(x)

        // Start with the first coefficient (c0 * T0(x))
        double result = coeffs[0] * T0 + coeffs[1] * T1;

        // Compute higher order terms using the recurrence relation
        for (int i = 2; i < n; i++) {
            double Tn = 2 * x * T1 - T0; // Recurrence: T_n(x) = 2 * x * T_(n-1)(x) - T_(n-2)(x)
            result += coeffs[i] * Tn;
            T0 = T1;  // Update for next iteration
            T1 = Tn;  // Update for next iteration
        }

        return result;
    }
    
 // Payne-Hanek range reduction algorithm for trigonometric functions
    public static double[] payneHanekReduction(double x) {
        double pi = PI;

        // Determine the integer n such that x is in the range [-pi/4, pi/4]
        double n = Math.floor((x + pi) / (2 * pi));  // Find the closest n

        // Compute the reduced argument x' = x - 2n * pi
        double reducedX = x - 2 * n * pi;

        // If reducedX is outside [-pi/4, pi/4], adjust it further
        if (reducedX > pi / 4) {
            reducedX -= pi;
            n += 1;
        } else if (reducedX < -pi / 4) {
            reducedX += pi;
            n -= 1;
        }

        // Return the reduced x and the integer n
        return new double[] {reducedX, n};
    }
    
	/**
	 * Computes the sine of a number.
	 * Uses Chebyshev Interpolation with 22 nodes to compute the sine.
	 * @param x The angle (in radians)
	 */
	public static double sine(double x) {
		/*
		 * Algorithm:
		 * 1) Calculate the Chebyshev coefficients. These coefficients have been computed with the Approximation class.
		 * 2) Reduce the interpolation to the interval [-PI/4 ; PI/4]. This is because the quality of the interpolation
		 * is lower the greater the angle x.
		 * 3) 
		 */
		if (x >= 2 * PI) {
			x = x % (2*PI);
		}
		double resReduced = x;
		if ((x < 0.000000001 && x > -0.000000001) || (x > PI - 0.000000001 && x < PI + 0.000000001)) {
			return 0.0;
		}
		if ((x > PI/2 - 0.000000001) && (x < PI / 2 + 0.00000001)) {
			return 1;
		}
		if ((x > 3*PI / 2 - 0.00000001) && (x < 3*PI / 2 + 0.000000001)) {
			return -1;
		}
		
		if (x >= PI) {
			return -Trigonometry.sine(x - PI);
		}
		else if (x >= PI / 2) {
	
			// We can use the identity sin(A+B)=sin(A)cos(B)+cos(A)sin(B)
			// With here A is PI / 4 and B the remainder over PI/2
			double A = PI / 2;
			double B = x - PI / 2;
			resReduced = cos(B);
		}
		else if (x >= (PI / 4))
		{
			// We can use the identity sin(A+B)=sin(A)cos(B)+cos(A)sin(B)
			// With here A is PI / 4 and B the remainder over PI/4
			double A = PI / 4;
			double B = x - PI / 4;
			resReduced = SIN_PI_4*cos(B)+COS_PI_4*chebyshevInterpolation(B, CHEBYSHEV_COEFFS_SIN);
		}
		else
		{
			resReduced = chebyshevInterpolation(x, CHEBYSHEV_COEFFS_SIN);
		}
		return resReduced;
	}
	
	// Method to calculate the factorial of a number
    public static long factorial(int n) {
        long fact = 1;
        for (int i = 1; i <= n; i++) {
            fact *= i;
        }
        return fact;
    }
    
	/**
	 * Computes the cosine of an angle.
	 * Uses a Taylor approximation up to x^12 to compute the cosine.
	 * @param x The angle in radians.
	 * @return The cosine of an angle.
	 */
	public static double cos(double x) {
			
			
			if (x < PI/2+0.000000001 && x > PI/2-0.000000001 ) {
				return 0;
			}
			if (x > 3*PI/2-0.000000001 && x < 3*PI/2 + 0.00000000001) {
				return 0;
			}
			
			if (x < 0.0000000001 && x >-0.0000000001) {
				return 1;
			}
			// cos(PI) = -1
			if (x >= PI - 0.000000001 && x <= PI + 0.00000001) {
				return -1;
			}
			
			if (x > PI) {
				return -Trigonometry.cos(x - PI);
			}
			else if (x > PI/2) {
				return -Trigonometry.cos(PI - x);
			}
			else if (x > PI / 4 && x < PI/2) {
				// Cos (a + b) = cos a * cos b - sin a * sin b
				double A = PI / 4;
				double B = x - PI / 4;
				return COS_PI_4 * Trigonometry.cos(B) - SIN_PI_4 * Trigonometry.sine(B);
			}
		    double result = 1.0;  // The first term of the series: 1
	        double term;
	        
	        // Loop for each term in the series, up to x^12
	        for (int n = 1; n <= 6; n++) {
	            // (-1)^n term for alternating signs
	            term = Math.pow(-1, n) * Math.pow(x, 2 * n) / factorial(2 * n);
	            result += term;
	        }
	        
	        return result;
	}
}
