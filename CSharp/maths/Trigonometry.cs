namespace maths
{

    public class Trigonometry
    {
        /**
         * The value of exp(1)
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
        private static long Factorial(int n)
        {
            long fact = 1;
            for (int i = 1; i <= n; i++)
            {
                fact *= i;
            }
            return fact;
        }

        /**
         * These factorial are useful for the calculation of the cosine function. No need to compute them each time.
         * Since we only go to the sixth order of the Taylor series, no need for more.
         */
        private static readonly long[] TAYLOR_COEFFS_SIN = [
            1, // x¹
			-6, // x³
			120, // x⁵
			-5040, // x⁷
			362880, // x⁹
			-39916800, // x^11
			6227020800, // x^13
			-1307674368000, // x^15
			355687428096000, // x^17
	        ];

        /**
         * These factorial are useful for the calculation of the cosine function. No need to compute them each time.
         * Since we only go to the sixth order of the Taylor series, no need for more.
         */
        private static readonly long[] TAYLOR_COEFFS_COS = [
            1, // x¹
			-2, // x²
			24, // x⁴
			-720, // x⁶
			40320, // x⁸
			-3628800, // x^10
			479001600, // x^12
			-87178291200, // x^14
			20922789888000, // x^16
	       ];
        
        private static long[] TAYLOR_COEFFS_EXP = [
            1, // x¹
		    2, // x²
			6, // x^3
			24, // x^4
			120, // x^5
			720, // x^6
			5040, // x^7
			40320, // x^8
			362880, // x^9
			3628800, // x^10
			39916800, // x^11
			479001600, // x^12
			6227020800, // x^13
			87178291200, // x^14
			1307674368000, // x^15
			20922789888000, // x^16
			355687428096000, // x^17
	];

        private static readonly double[] EXPONENTIAL_TABLE = [
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

    ];

        /**
         * These are the Chebyshev nodes corresponding to the X in the interval
         * [1,2]. We use 20 nodes for accuracy.
         * We do not use a Lagrange Polynomial with points spaced out evenly to avoid Runge phenomenon.
         * See: https://en.wikipedia.org/wiki/Chebyshev_nodes
         * Note: If you wish to change the number of nodes, just call Approximation.chebyshevNodes(n,a,b)
         * Example: Generating 12 nodes between values 0 and 1: Approximation.chebyshevNodes(12, 0, 1)
         */
        private static readonly double[] CHEBYSHEV_NODES_LOG = [
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
    ];

        /**
         * These are the values of Ln(x) evaluated at the point x = CHEBYSHEV_NODES_LOG.
         * To evaluate precisely Ln(x) you can use a Taylor series, which has slower convergence.
         */
        private static double[] CHEBYSHEV_VALUES_LOG = [
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

    ];
        /**
         * Power function.
         * @param x A value
         * @param p The power of the value
         * @return x^p
         */
        public static double Pow(double x, int p)
        {
            double res = 1;
            while (p < 0)
            {
                res /= x;
                p++;
            }
            while (p > 0)
            {
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
        public static double TaylorInterpolationSin(double x)
        {
            if (x > PI / 4 || x < -PI / 4)
            {
                return Trigonometry.Sin(x);
            }
            double result = x;
            // Loop for each term in the series, up to x^16
            double x2 = x * x;
            double powx = x2 * x;
            for (int n = 1; n <= 8; n++)
            {
                // (-1)^n term for alternating signs
                double term = powx / TAYLOR_COEFFS_SIN[n];
                powx = powx * x2;
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
        public static double TaylorInterpolationCos(double x)
        {
            if (x > PI / 4 || x < -PI / 4)
            {
                return Trigonometry.Cos(x);
            }
            double result = 1.0;
            // Loop for each term in the series, up to x^16
            double x2 = x * x;
            double powx = x2;
            for (int n = 1; n <= 8; n++)
            {
                // (-1)^n term for alternating signs
                double term = powx / TAYLOR_COEFFS_COS[n];
                powx = powx * x2;
                result += term;
            }
            return result;
        }

        public static double TaylorInterpolationExp(double x)
        {
            double result = 1.0;
            double powx = x;
            for (int n = 0; n <= 15; n++)
            {
                double term = powx / Trigonometry.TAYLOR_COEFFS_EXP[n];
                powx = powx * x;
                result += term;
            }
            return result;

        }

        /**
         * Constant used by the TaylorInterpolationArctan function.
         */
        private static readonly double[] atanhi = [
              4.63647609000806093515e-01, /* atan(0.5)hi 0x3FDDAC67, 0x0561BB4F */
    		  7.85398163397448278999e-01, /* atan(1.0)hi 0x3FE921FB, 0x54442D18 */
    		  9.82793723247329054082e-01, /* atan(1.5)hi 0x3FEF730B, 0xD281F69B */
    		  1.57079632679489655800e+00, /* atan(inf)hi 0x3FF921FB, 0x54442D18 */
    	];

        /**
         * Constant used by the TaylorInterpolationArctan function.
         */
        private static readonly double[] atanlo = [
              2.26987774529616870924e-17, /* atan(0.5)lo 0x3C7A2B7F, 0x222F65E2 */
    		  3.06161699786838301793e-17, /* atan(1.0)lo 0x3C81A626, 0x33145C07 */
    		  1.39033110312309984516e-17, /* atan(1.5)lo 0x3C700788, 0x7AF0CBBD */
    		  6.12323399573676603587e-17, /* atan(inf)lo 0x3C91A626, 0x33145C07 */
    	];

        /**
         * Constant used by the TaylorInterpolationArctan function.
         */
        private static double[] aT = [
              3.33333333333329318027e-01, /* 0x3FD55555, 0x5555550D */
    		 -1.99999999998764832476e-01, /* 0xBFC99999, 0x9998EBC4 */
    		  1.42857142725034663711e-01, /* 0x3FC24924, 0x920083FF */
    		 -1.11111104054623557880e-01, /* 0xBFBC71C6, 0xFE231671 */
    		  9.09088713343650656196e-02, /* 0x3FB745CD, 0xC54C206E */
    		 -7.69187620504482999495e-02, /* 0xBFB3B0F2, 0xAF749A6D */
    		  6.66107313738753120669e-02, /* 0x3FB10D66, 0xA0D03D51 */
    		 -5.83357013379057348645e-02, /* 0xBFADDE2D, 0x52DEFD9A */
    		  4.97687799461593236017e-02, /* 0x3FA97B4B, 0x24760DEB */
    		 -3.65315727442169155270e-02, /* 0xBFA2B444, 0x2C6A6C2F */
    		  1.62858201153657823623e-02, /* 0x3F90AD3A, 0xE322DA11 */
    	];

        /**
         * Computation an approximation of the arctangent of a ratio.
         * 
         * Method
         *   1. Reduce x to positive by atan(x) = -atan(-x).
         *   2. According to the integer k=4t+0.25 chopped, t=x, the argument
         *      is further reduced to one of the following intervals and the
         *      arctangent of t is evaluated by the corresponding formula:
         *
         *      [0,7/16]      atan(x) = t-t^3*(a1+t^2*(a2+...(a10+t^2*a11)...)
         *      [7/16,11/16]  atan(x) = atan(1/2) + atan( (t-0.5)/(1+t/2) )
         *      [11/16.19/16] atan(x) = atan( 1 ) + atan( (t-1)/(1+t) )
         *      [19/16,39/16] atan(x) = atan(3/2) + atan( (t-1.5)/(1+1.5t) )
         *      [39/16,INF]   atan(x) = atan(INF) + atan( -1/t )
         *
         * Constants:
         * The hexadecimal values are the intended ones for the following
         * constants. The decimal values may be used, provided that the
         * compiler will convert from decimal to binary accurately enough
         * to produce the hexadecimal values shown.
         * 
         * For more information:
         * See https://git.musl-libc.org/cgit/musl/tree/src/math/atan.c
         * @param x A real representing a ratio of adjacent side to opposite side.
         * @return The theta angle.
         */
        public static double TaylorInterpolationArctan(double x)
        {
            int id;
            bool negative_sign;
            if (x < 0.0)
            {
                negative_sign = true;
            }
            else
            {
                negative_sign = false;
            }
            if (x < 0.0)
                x = -x;
            if (x < 1.1875)
            {  /* |x| < 1.1875 */
                if (x < 11.0 / 16.0)
                {  /*  7/16 <= |x| < 11/16 */
                    id = 0;
                    x = (2.0 * x - 1.0) / (2.0 + x);
                }
                else
                {                /* 11/16 <= |x| < 19/16 */
                    id = 1;
                    x = (x - 1.0) / (x + 1.0);
                }
            }
            else
            {
                if (x < 2.4375)
                {  /* |x| < 2.4375 */
                    id = 2;
                    x = (x - 1.5) / (1.0 + 1.5 * x);
                }
                else
                {                /* 2.4375 <= |x| < 2^66 */
                    id = 3;
                    x = -1.0 / x;
                }
            }

            double z = x * x;
            double w = z * z;
            /* break sum from i=0 to 10 aT[i]z**(i+1) into odd and even poly */
            double s1 = z * (aT[0] + w * (aT[2] + w * (aT[4] + w * (aT[6] + w * (aT[8] + w * aT[10])))));
            double s2 = w * (aT[1] + w * (aT[3] + w * (aT[5] + w * (aT[7] + w * aT[9]))));
            if (id < 0)
                return x - x * (s1 + s2);
            z = atanhi[id] - (x * (s1 + s2) - atanlo[id] - x);
            return negative_sign == true ? -z : z;

        }

        // Method to evaluate the Chebyshev series for cos(x)
        public static double chebyshevInterpolationcos(double x, double[] coeffs)
        {
            int n = coeffs.Length;

            // Initialize T_0(x) and T_1(x)
            double T0 = 1.0;  // T_0(x)
            double T1 = x;    // T_1(x)

            // Start with the first coefficient (c0 * T0(x))
            double result = coeffs[0] * T0 + coeffs[1] * T1;

            // Compute higher order terms using the recurrence relation
            for (int i = 2; i < n; i++)
            {
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
        public static double Sin(double x)
        {
            /*
             * Algorithm:
             * 1) Calculate the Chebyshev coefficients. These coefficients have been computed with the Approximation class.
             * 2) Reduce the interpolation to the interval [-PI/4 ; PI/4]. This is because the quality of the interpolation
             * is lower the greater the angle x.
             * 3) 
             */

            // We use the fact that sin(-x) = -sin(x), so that we are
            // always in the positive range.
            if (x < 0.0)
            {
                return -Trigonometry.Sin(-x);
            }

            if (x >= 2 * PI)
            {
                x = x % (2 * PI);
            }
            double resReduced = x;
            if ((x < 0.000000001 && x > -0.000000001) || (x > PI - 0.000000001 && x < PI + 0.000000001))
            {
                return 0.0;
            }
            if ((x > PI / 2 - 0.000000001) && (x < PI / 2 + 0.00000001))
            {
                return 1;
            }
            if ((x > 3 * PI / 2 - 0.00000001) && (x < 3 * PI / 2 + 0.000000001))
            {
                return -1;
            }

            if (x >= PI)
            {
                return -Trigonometry.Sin(x - PI);
            }
            else if (x >= PI / 2)
            {

                // We can use the identity sin(A+B)=sin(A)cos(B)+cos(A)sin(B)
                // With here A is PI / 4 and B the remainder over PI/2
                double B = x - PI / 2;
                resReduced = Trigonometry.Cos(B);
            }
            else if (x >= (PI / 4))
            {
                // We can use the identity sin(A+B)=sin(A)cos(B)+cos(A)sin(B)
                // With here A is PI / 4 and B the remainder over PI/4
                double B = x - PI / 4;
                resReduced = SIN_PI_4 * Trigonometry.Cos(B) + COS_PI_4 * TaylorInterpolationSin(B);
            }
            else
            {
                resReduced = TaylorInterpolationSin(x);
            }
            return resReduced;
        }

        /**
         * Computes the cosine of an angle.
         * The value y = Cos(alpha) of an angle alpha expressed in radians is expressed easily using
         * a right triangle.
         * Imagine you have a circle. Inside the circle, there are two axis, the x and the y
         * axis. You select a point [x1,y1] on the unit circle, and it forms an angle alpha between the
         * x axis and the segment [xO, yO] and [x1, y1]. xO and yO being the origin.
         * Given a right triangle with an angle theta
         * We know and can use the following formulas
         * By default, uses a 12th degree Taylor polynomial to compute the sine.
         * @param x The angle in radians.
         * @return The cosine of an angle.
         */
        public static double Cos(double x)
        {

            // Remember the identity that cos(x) = cos(-x) 
            if (x < 0.0)
            {
                return Cos(-x);
            }

            if (x >= 2 * PI)
            {
                x = x % (2 * PI);
            }


            if (x < PI / 2 + 0.000000001 && x > PI / 2 - 0.000000001)
            {
                return 0;
            }
            if (x > 3 * PI / 2 - 0.000000001 && x < 3 * PI / 2 + 0.00000000001)
            {
                return 0;
            }

            if (x < 0.0000000001 && x > -0.0000000001)
            {
                return 1;
            }
            // cos(PI) = -1
            if (x >= PI - 0.000000001 && x <= PI + 0.00000001)
            {
                return -1;
            }

            if (x > PI)
            {
                return -Trigonometry.Cos(x - PI);
            }
            else if (x > PI / 2)
            {
                return -Trigonometry.Cos(PI - x);
            }
            else if (x > PI / 4 && x < PI / 2)
            {
                // Cos (a + b) = cos a * cos b - sin a * sin b
                double B = x - PI / 4;
                return COS_PI_4 * Trigonometry.Cos(B) - SIN_PI_4 * Trigonometry.Sin(B);
            }

            return TaylorInterpolationCos(x);
        }

        // Method to compute the Lagrange interpolation polynomial at point x
        public static double LagrangeInterpolation(double x, double[] nodes, double[] values)
        {
            double result = 0.0;
            int n = nodes.Length;
            for (int i = 0; i < n; i++)
            {
                double term = values[i];
                for (int j = 0; j < n; j++)
                {
                    if (i != j)
                    {
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
        public static int Binlog(int bits) // returns 0 for bits=0
        {
            int log = 0;
            if ((bits & 0xffff0000) != 0) { bits >>>= 16; log = 16; }
            if (bits >= 256) { bits >>>= 8; log += 8; }
            if (bits >= 16) { bits >>>= 4; log += 4; }
            if (bits >= 4) { bits >>>= 2; log += 2; }
            return log + (bits >>> 1);
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
        public static double Log(double y)
        {
            int log2;
            double divisor, x, result;

            log2 = Trigonometry.Binlog((int)y);
            divisor = (double)(1 << log2);
            x = y / divisor;

            result = Trigonometry.LagrangeInterpolation(x, CHEBYSHEV_NODES_LOG, CHEBYSHEV_VALUES_LOG);
            result += ((double)log2) * Ln2;
            return result;
        }

        /**
         * An estimation of the arctangent function.
         * 
         * Important: Instead of using atan, consider using atan2.
         * @param x Ratio of the opposite side over the adjacent side (dy / dx)
         * @return Theta the angle 
         */
        public static double Atan(double x)
        {
            return TaylorInterpolationArctan(x);
        }

        /**
         * Computes an approximation of the arctangent of the ratio y / x.
         * 
         * @param y Any real. Represents the opposite side in a triangle.
         * @param x Any real. Represents the adjacent side in a triangle.
         * @return
         */
        public static double Atan2(double y, double x)
        {
            if ((x > -0.0000001 && x < 0.000001) && y > 0)
            {
                return PI / 2.0;
            }
            if ((x > -0.000001 && x < 0.000001) && y < 0)
            {
                return -PI / 2.0;
            }
            if (x > 0)
            {
                return Trigonometry.Atan(y / x);
            }
            if (x < 0 && y >= 0)
            {
                return Trigonometry.Atan(y / x) + PI;
            }
            if (x < 0 && y < 0)
            {
                return Trigonometry.Atan(y / x) - PI;
            }

            return Double.NaN;
        }

        /**
         * The logarithm of a number with base b.
         * As a reminder, logb(x) = log(x) / log(b).
         * 
         * @param x Any real number.
         * @param base The base of the logarithm.
         * @return
         */
        public static double Log(double x, double @base)
        {
            return Trigonometry.Log(x) / Trigonometry.Log(@base);
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
        public static double Log10(double x)
        {
            return Trigonometry.Log(x, 10);
        }

        /**
         * Computes x to the power y.
         * @param x Any real.
         * @param y Any real
         * @return
         */
        public static double Pow(double x, double y)
        {
            return Trigonometry.Exp(y * Trigonometry.Log(x));
        }

        /**
         * The exponential function. 
         * 
         * Method:
         * Here we use the property of exp(a+b) = exp(a) * exp(b)
         * For example, exp(17.123) = exp(17) * exp(0.123)
         * Like this, we can reduce the exponential to the range [0,1]
         * Also, since we know that the max value of a double is 10E308, we "only" have to compute the
         * exp(1)...exp(308), otherwise it is "infinity" in the sense that it goes above the storage.	
         * Special cases:
         * NaN: Returns NaN
         * 
         * @param x Any floating point number.
         * @return The exponential of any number.
         */
        public static double Exp(double x)
        {

            // of 
            if (x == 0)
            {
                return 1.0;
            }

            if (x < 0.0)
            {
                return 1.0 / Exp(-x);
            }
            int nearestExponential = (int)(x);
            if (nearestExponential > EXPONENTIAL_TABLE.Length)
            {
                return Double.PositiveInfinity;
            }

            double res = EXPONENTIAL_TABLE[nearestExponential] * TaylorInterpolationExp(x - (double)nearestExponential);

            return res;
        }

    }
}
