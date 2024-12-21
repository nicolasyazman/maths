package io.github.nicolasyazman;

import java.math.BigDecimal;

public class Trigonometry {

	/**
	 * The value of exp¹
	 */
	public static double E = 2.71828182845904523536028747135266249775724709369995957496696762772407663035;
	
	/**
	 * The value of natural logarithm of 2. y = ln(2)
	 */
	public static double Ln2 = 0.69314718055994530941723212145817656807550013436025525412068000949339362196;
	/**
	 * The value of the irrational number Pi.
	 */
	public static double PI = 3.141592653589793238462643383279502884197;
	
	/**
	 * The value of the cosine and sine of the angle Pi/4 radians (45 degrees)
	 */
	private static double SIN_PI_4 = 0.70710678118654752440084436210484903928483593768847403658833986899536623923;
	private static double COS_PI_4 = 0.70710678118654752440084436210484903928483593768847403658833986899536623923;
	
	/**
	 * A helper function that is used to compute the factorial of a number.
	 * As a reminder fact(N) = N! = 1*...*N
	 * @param n Number to compute the factorial for, in N (natural numbers group)
	 * @return The factorial of n (N!)
	 */
    private static long factorial(int n) {
        long fact = 1;
        for (int i = 1; i <= n; i++) {
            fact *= i;
        }
        return fact;
    }
    
	/**
	 * These factorial are useful for the calculation of the cosine function. No need to compute them each time.
	 * Since we only go to the sixth order of the Taylor series, no need for more.
	 */
	private static long TAYLOR_COEFFS_SIN[] = {
			1, // x¹
			-factorial(3), // x³
			factorial(5), // x⁵
			-factorial(7), // x⁷
			factorial(9), // x⁹
			-factorial(11), // x^11
			factorial(13), // x^13
			-factorial(15), // x^15
			factorial(17), // x^17
	};

	/**
	 * These factorial are useful for the calculation of the cosine function. No need to compute them each time.
	 * Since we only go to the sixth order of the Taylor series, no need for more.
	 */
	private static long TAYLOR_COEFFS_COS[] = {
			1, // x¹
			-factorial(2), // x²
			factorial(4), // x⁴
			-factorial(6), // x⁶
			factorial(8), // x⁸
			-factorial(10), // x^10
			factorial(12), // x^12
			-factorial(14), // x^14
			factorial(16), // x^16
	};

	private static long TAYLOR_COEFFS_EXP[] = {
			1, // x¹
			factorial(2), // x²
			factorial(3), // x⁴
			factorial(4), // x⁶
			factorial(5), // x⁸
			factorial(6), // x^10
			factorial(7), // x^12
			factorial(8), // x^14
			factorial(9), // x^16
			factorial(10), // x^16
			factorial(11), // x^16
			factorial(12), // x^16
			factorial(13), // x^16
			factorial(14), // x^16
			factorial(15), // x^16
			factorial(16), // x^16
			factorial(17), // x^16
	};
	
	
	private static int EXPONENTIAL_TABLE_IDX[] = {
			1,
			2,
			4,
			16,
			256
	};
	
	private static double EXPONENTIAL_TABLE[] = {
			1,
			2.71828182845904523536028747135266249775724709369995957496696762772407663035,     // e(1)
			7.38905609893065022723042746057500781318031557055184732408712782252257379607,     // e(2)
			20.0855369231876677409285296545817178969879078385541501443789342296988458780,     // e(3)
			54.5981500331442390781102612028608784027907370386140687258265939585536620999,     // e(4)
			148.413159102576603421115580040552279623487667593878989046752845110912064820,     // e(5)
			403.428793492735122608387180543388279605899897357129202613967188325151180633,     // e(6)
			1096.63315842845859926372023828812143244221913483361314378273924077612176933,     // e(7)
			2980.95798704172827474359209945288867375596793913283570220896353038773072517,     // e(8)
			8103.08392757538400770999668943275996501147608783161346250015905217827251569,     // e(9)
			22026.4657948067165169579006452842443663535126185567810742354263552252028185,     // e(10)
			59874.1417151978184553264857922577816142610796957409686527715143799324211159,     // e(11)
			162754.791419003920808005204898486783170209284478720770443556248138596770835,     // e(12)
			442413.392008920503326102775949088281784391306060589715572359330709021862336,     // e(13)
			1.202604284164776777749236770767859449412486543376102240313290633197462E6,
			
	};
	
	/**
	 * The Chebyshev coeffs for sine function.
	 * TO BE IMPLEMENTED.
	 */
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
	
	/**
	 * The Chebyshev coeffs for cosine function.
	 * TO BE IMPLEMENTED.
	 * Note: There seems to be a problem with these coefficients.
	 */
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
	 * These are the Chebyshev nodes corresponding to the X in the interval
	 * [1,2]. We use 20 nodes for accuracy.
	 * We do not use a Lagrange Polynomial with points spaced out evenly to avoid Runge phenomenon.
	 * See: https://en.wikipedia.org/wiki/Chebyshev_nodes
	 * Note: If you wish to change the number of nodes, just call Approximation.chebyshevNodes(n,a,b)
	 * Example: Generating 12 nodes between values 0 and 1: Approximation.chebyshevNodes(12, 0, 1)
	 */
	private static double CHEBYSHEV_NODES_LOG[] = {
			1.998458666866564,
			1.9861849601988384,
			1.9619397662556435,
			1.926320082177046,
			1.8802029828000155,
			1.824724024165092,
			1.7612492823579744,
			1.6913417161825448,
			1.6167226819279528,
			1.5392295478639224,
			1.4607704521360776,
			1.3832773180720475,
			1.3086582838174552,
			1.2387507176420256,
			1.1752759758349083,
			1.1197970171999845,
			1.073679917822954,
			1.0380602337443565,
			1.0138150398011616,
			1.001541333133436
	};
	
	/**
	 * These are the values of Ln(x) evaluated at the point x = CHEBYSHEV_NODES_LOG.
	 * To evaluate precisely Ln(x) you can use a Taylor series, which has slower convergence.
	 */
	private static double CHEBYSHEV_VALUES_LOG[] = {
			0.6923762168770874,
			0.6862156933100475,
			0.6739336604963982,
			0.6556114896909621,
			0.6313797405880391,
			0.6014287559640469,
			0.5660233767689095,
			0.5255221288521638,
			0.4804010642962612,
			0.4312819976420859,
			0.37896400347442033,
			0.3244555517946695,
			0.26900240152437127,
			0.2141033859894615,
			0.1615029930687512,
			0.11314743423914973,
			0.07109192356966546,
			0.037353811715549044,
			0.013720482023170928,
			0.0015401464986983199

	};
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
	
    /**
     * Computes the sine function using a 13th Degree Taylor approximation.
     * Only valid in the range [0, PI / 4]
     * If it is outside this range, a combination of sine and cosine will be used.
     * Approximates using 13th degree polynomial.
     * Note: We chose not to use Horner's algorithm to avoid catastrophic cancellation.
     * @param x An angle (in radians)
     * @param chebyshevCoeffs
     * @return
     */
    public static double TaylorInterpolationSin(double x) {
    	if (x > PI / 4 || x < -PI / 4) {
    		return Trigonometry.sin(x);
    	}
    	double result = x;
    	  // Loop for each term in the series, up to x^16
    	double x2 = x*x;
    	double powx = x2*x;
        for (int n = 1; n <= 8; n++) {
            // (-1)^n term for alternating signs
            double term = powx / TAYLOR_COEFFS_SIN[n];
            powx = powx*x2;
            result += term;
        }   
        return result;
    }
    
    /**
     * Computes the cosine function using a 13th Degree Taylor approximation.
     * Only valid in the range [0, PI / 4]
     * If it is outside this range, a combination of sine and cosine will be used.
     * Approximates using 13th degree polynomial.
     * Note: We chose not to use Horner's algorithm to avoid catastrophic cancellation.
     * @param x An angle (in radians)
     * @return
     */
    public static double TaylorInterpolationCos(double x) {
    	if (x > PI / 4 || x < -PI / 4) {
    		return Trigonometry.cos(x);
    	}
    	double result = 1.0;
    	  // Loop for each term in the series, up to x^16
    	double x2 = x*x;
    	double powx = x2;
        for (int n = 1; n <= 8; n++) {
            // (-1)^n term for alternating signs
            double term = powx / TAYLOR_COEFFS_COS[n];
            powx = powx*x2;
            result += term;
        }   
        return result;
    }
    
    public static double TaylorInterpolationExp(double x) {
    	double result = 1.0;
	  	double powx = x;
	      for (int n = 0; n <= 15; n++) {
	          double term = powx / TAYLOR_COEFFS_EXP[n];
	          powx = powx*x;
	          result += term;
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
        
	/**
	 * Computes the sine of a number, y = sin(x)
	 * By default, uses a 13th degree Taylor polynomial to compute the sine.
	 * @param x The angle (in radians)
	 * @return The sine of an angle.
	 */
	public static double sin(double x) {
		/*
		 * Algorithm:
		 * 1) Calculate the Chebyshev coefficients. These coefficients have been computed with the Approximation class.
		 * 2) Reduce the interpolation to the interval [-PI/4 ; PI/4]. This is because the quality of the interpolation
		 * is lower the greater the angle x.
		 * 3) 
		 */
		
		// We use the fact that sin(-x) = -sin(x), so that we are
		// always in the positive range.
		if (x < 0.0) {
			return -Trigonometry.sin(-x);
		}
		
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
			return -Trigonometry.sin(x - PI);
		}
		else if (x >= PI / 2) {
	
			// We can use the identity sin(A+B)=sin(A)cos(B)+cos(A)sin(B)
			// With here A is PI / 4 and B the remainder over PI/2
			double A = PI / 2;
			double B = x - PI / 2;
			resReduced = Trigonometry.cos(B);
		}
		else if (x >= (PI / 4))
		{
			// We can use the identity sin(A+B)=sin(A)cos(B)+cos(A)sin(B)
			// With here A is PI / 4 and B the remainder over PI/4
			double A = PI / 4;
			double B = x - PI / 4;
			resReduced = SIN_PI_4*Trigonometry.cos(B)+COS_PI_4*TaylorInterpolationSin(B);	
		}
		else
		{
			resReduced = TaylorInterpolationSin(x);
		}
		return resReduced;
	}
	 
	/**
	 * Computes the cosine of an angle.
	 * By default, uses a 12th degree Taylor polynomial to compute the sine.
	 * @param x The angle in radians.
	 * @return The cosine of an angle.
	 */
	public static double cos(double x) {
		
		// Remember the identity that cos(x) = cos(-x) 
		if (x < 0.0) {
			return cos(-x);
		}
		
		if (x >= 2 * PI) {
			x = x % (2*PI);
		}
		
			
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
				return COS_PI_4 * Trigonometry.cos(B) - SIN_PI_4 * Trigonometry.sin(B);
			}
		    
			return TaylorInterpolationCos(x);
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
	
    /**
     * Returns the Log2(x) of an integer. See: https://stackoverflow.com/questions/3305059/how-do-you-calculate-log-base-2-in-java-for-integers
     * The main idea is to use bit shifting, since binary is a base2.
     * If we right-shift a number 16 times, and it is not equal to zero, it means that its Log2 is at least 16
     * Then we try with 8, then 4, then 2, then 1. The resulting Log2 is then just the sum of all the shifts.
     * @param bits An integer we want to find the Log2 to.
     * @return The image of the function Log2(x)
     */
    public static int binlog( int bits ) // returns 0 for bits=0
    {
        int log = 0;
        if( ( bits & 0xffff0000 ) != 0 ) { bits >>>= 16; log = 16; }
        if( bits >= 256 ) { bits >>>= 8; log += 8; }
        if( bits >= 16  ) { bits >>>= 4; log += 4; }
        if( bits >= 4   ) { bits >>>= 2; log += 2; }
        return log + ( bits >>> 1 );
    }
    
    /**
     * The natural logarithm of a function. Uses Remez algorithm to reduce it to
     * the range [1,2].
     * For the evaluation of the Log function between [1,2], uses a Lagrange
     * Interpolation with Chebyshev Nodes, with 20 Nodes.
     * 
     * @param y Any real number.
     * @return Logarithm of y.
     */
	public static double log(double y) {
		int log2;
		double divisor,x, result;
		
		log2 = binlog((int)y);
		divisor = (double)(1 << log2);
		x = y / divisor;
		
		result = lagrangeInterpolation(x, CHEBYSHEV_NODES_LOG, CHEBYSHEV_VALUES_LOG);	
		result += ((double)log2) * Ln2;
		return result;
	}
	
	/**
	 * The logarithm of a number with base b.
	 * As a reminder, logb(x) = log(x) / log(b).
	 * 
	 * @param x Any real number.
	 * @param base The base of the logarithm.
	 * @return
	 */
	public static double log(double x, double base) {
		return Trigonometry.log(x) / Trigonometry.log(base);
	}
	
	/**
	 * Return the Log base 10 of a number.
	 * 
	 * See function log(double x, double base)
	 * Special cases: 
	 * NaN: Returns NaN
	 * 
	 * @param x Any real number.
	 * @return
	 */
	public static double log10(double x) {
		return Trigonometry.log(x, 10);
	}
	
	public static double exp(double x) {
		// Algorithm: 
		// Here we use the property of exp(a+b) = exp(a) * exp(b)
		// For example, exp(17.123) = exp(17) * exp(0.123)
		// Like this, we can reduce the exponential to the range [0,1]
		// Also, since we know that the max value of a double is 10E308, we "only" have to compute the
		// exp(1)...exp(308), otherwise it is "infinity" in the sense that it goes above the storage
		// of 
		if (x == 0) {
			return 1.0;
		}
		
		if (x < 0.0) {
			return 1.0 / exp(-x);
		}
		int nearestExponential = (int)Math.floor(x);
		if (nearestExponential > EXPONENTIAL_TABLE.length) {
			return Double.POSITIVE_INFINITY;
		}
		
		double res = EXPONENTIAL_TABLE[nearestExponential] * TaylorInterpolationExp(x - (double)nearestExponential);
		
		return res;
	}
}
