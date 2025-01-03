using maths;

namespace mathTests
{
 
    public class TrigonometryTests
    {
        [SetUp]
        public void Setup()
        {

        }

        [Test]
    public void Compute10ToThePower3()
        {
            Assert.AreEqual(1000.0, Trigonometry.Pow(10.0, 3.0), Math.Pow(10, -10));
        }

        [Test]
    public void Compute2ToThePower8()
        {
            
            Assert.AreEqual(256, Trigonometry.Pow(2, 8), Math.Pow(10, -10));
        }

        [Test]
    public void ComputeNegativeNumberPower()
        {
            Assert.AreEqual(1.0 / 256.0, Trigonometry.Pow(2, -8), Math.Pow(10, -10));
        }


        [Test]
    public void ComputePositiveCos()
        {
            for (int i = 0; i < 1000; i++)
            {
                Assert.AreEqual(Math.Cos(i * (Math.PI / 180.0)), Trigonometry.Cos(i * (Math.PI / 180.0)), Math.Pow(10, -12));
            }
        }

        [Test]
    public void ComputeNegativeCos()
        {
            for (int i = 0; i < 1000; i++)
            {
                Assert.AreEqual(Math.Cos(i * (-Math.PI / 180.0)), Trigonometry.Cos(i * (-Math.PI / 180.0)), Math.Pow(10, -12));
            }
        }

        [Test]
    public void ComputePositiveSin()
        {
            for (int i = 0; i < 1000; i++)
            {
                Assert.AreEqual(Math.Sin(i * (Math.PI / 180.0)), Trigonometry.Sin(i * (Math.PI / 180.0)), Math.Pow(10, -12));
            }
        }

        [Test]
    public void ComputeNegativeSin()
        {
            for (int i = 0; i < 1000; i++)
            {
                Assert.AreEqual(Math.Sin(i * (-Math.PI / 180.0)), Trigonometry.Sin(i * (-Math.PI / 180.0)), Math.Pow(10, -12));
            }
        }
        [Test]
    public void ComputeExponentialApproximationAndCompareToSuns()
        {
            for (int i = 0; i < 10; i++)
            {
                Assert.AreEqual(Math.Exp((double)(i) * 0.5), Trigonometry.Exp((double)(i) * 0.5), Math.Pow(10, -12));
            }
        }

        [Test]
    public void ComputeNegativeExponentialApproximationAndCompareToSuns()
        {
            for (int i = 0; i < 10; i++)
            {
                Assert.AreEqual(Math.Exp((double)(i) * -0.5), Trigonometry.Exp((double)(i) * -0.5), Math.Pow(10, -12));
            }
        }

        [Test]
    public void ComputeAtan2AndCompareToSuns()
        {

            for (int x = -100; x < 100; x++)
            {
                for (int y = -100; y < 100; y++)
                {
                    if (x != y)
                    {
                        double sunAtan2 = Math.Atan2((double)y, (double)x);
                        double myAtan2 = Trigonometry.Atan2((double)y, (double)x);
                        Assert.AreEqual(sunAtan2, myAtan2, Math.Pow(10, -12));
                    }
                }
            }

        }

        /*
        [Test]
        public void computeArctanAndCompareToSuns() {
            for (int i = 0; i < 3; i++) {
                Assert.AreEqual(Math.atan(i* (Math.PI /180.0)), Trigonometry.atan(i* (Math.PI /180.0)), Math.Pow(10,-12));
            }
        }*/

        [Test]
    public void computeRandomExponentialBetween0and1AndCompareToSuns()
        {
            Random ran = new Random();
            double randx = ran.NextDouble();
            Assert.AreEqual(Math.Exp(randx), Trigonometry.Exp(randx), Math.Pow(10, -10));
        }

        [Test]
    public void ComputeLogBetween1and2AndCompareItToSuns()
        {
            Random ran = new Random();
            double randx = 15 + ran.NextDouble();
            Console.WriteLine("Random value selected: " + randx);
            Console.WriteLine("My Log of this random value: " + Trigonometry.Log(randx));
            Console.WriteLine("Microsoft's Log of this random value: " + Math.Log(randx));

            //System.out.println("Random number generated for natural Logarithm test:");
            //System.out.println(randx);
            Assert.AreEqual(Math.Log(randx), Trigonometry.Log(randx), Math.Pow(10, -12));
        }

        [Test]
    public void ComputeLogBase10Of938183()
        {
            double x = 938183;
            Assert.AreEqual(Math.Log10(x), Trigonometry.Log(x, 10), Math.Pow(10, -12));
        }

        [Test]
    public void SinNinetyDegreesShouldBeAround1()
        {
            Assert.AreEqual(1.0, Trigonometry.Sin(Trigonometry.PI / 2), 0.00001);
        }

        [Test]
    public void Sin270ShouldBeAroundMinus1()
        {
            Assert.AreEqual(-1.0, Trigonometry.Sin(3 * Trigonometry.PI / 2), 0.00001);
        }


        [Test]
    public void SinZeroDegreesShouldBeAround0()
        {
            Assert.AreEqual(0.0, Trigonometry.Sin(0), 0.00001);
        }

        [Test]
    public void CosZeroShouldBeAround1()
        {
            Assert.AreEqual(1.0, Trigonometry.Cos(0), 0.00001);
        }

        [Test]
    public void CosNinetyDegreesShouldBeAround0()
        {
            Assert.AreEqual(0.0, Trigonometry.Cos(Trigonometry.PI / 2), 0.0001);
        }

        [Test]
    public void Cos42PiShouldBeAround1()
        {
            Assert.AreEqual(1, Trigonometry.Cos(Trigonometry.PI * 42), 0.0001);
        }
        [Test]
    public void Cos43PiShouldBeAroundMinus1()
        {
            Assert.AreEqual(-1, Trigonometry.Cos(Trigonometry.PI * 43), 0.0001);
        }

        [Test]
    public void XPow3shouldbe27()
        {
            Assert.AreEqual(27, Trigonometry.Pow(3, 3));
        }
    }
}