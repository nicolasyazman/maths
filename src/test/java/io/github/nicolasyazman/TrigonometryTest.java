package io.github.nicolasyazman;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import io.github.nicolasyazman.Trigonometry;
import org.junit.jupiter.api.Test;
import java.lang.Math;
import java.util.Random;

public class TrigonometryTest {

	@Test
	public void computePositiveCos() {
		for (int i = 0; i < 1000; i++) {
			assertEquals(Math.cos(i* (Math.PI /180.0)), Trigonometry.cos(i* (Math.PI /180.0)), Math.pow(10,-12));
		}
	}
	
	@Test
	public void computeNegativeCos() {
		for (int i = 0; i < 1000; i++) {
			assertEquals(Math.cos(i* (-Math.PI /180.0)), Trigonometry.cos(i* (-Math.PI /180.0)), Math.pow(10,-12));
		}
	}
	
	@Test
	public void computePositiveSin() {
		for (int i = 0; i < 1000; i++) {
			assertEquals(Math.sin(i* (Math.PI /180.0)), Trigonometry.sin(i* (Math.PI /180.0)), Math.pow(10,-12));
		}
	}
	
	@Test
	public void computeNegativeSin() {
		for (int i = 0; i < 1000; i++) {
			assertEquals(Math.sin(i* (-Math.PI /180.0)), Trigonometry.sin(i* (-Math.PI /180.0)), Math.pow(10,-12));
		}
	}
	@Test
	public void computeExponentialApproximationAndCompareToSuns() {
		double x = 0;
		for (int i = 0; i < 10; i++) {
			assertEquals(Math.exp((double)(i)*0.5), Trigonometry.exp((double)(i)*0.5), Math.pow(10,-12));
		}
	}

	@Test
	public void computeNegativeExponentialApproximationAndCompareToSuns() {
		double x = 0;
		for (int i = 0; i < 10; i++) {
			assertEquals(Math.exp((double)(i)*-0.5), Trigonometry.exp((double)(i)*-0.5), Math.pow(10,-12));
		}
	}
	
	@Test
	public void computeRandomExponentialBetween0and1AndCompareToSuns() {
		Random ran = new Random();
		double randx = ran.nextDouble();
		assertEquals(Math.exp(randx), Trigonometry.exp(randx), Math.pow(10, -10));
	}
	
	@Test
	public void computeLogBetween1and2AndCompareItToSuns() {
		Random ran = new Random();
		double randx = 15 + ran.nextDouble();
		System.out.println("Random number generated for natural logarithm test:");
		System.out.println(randx);
		assertEquals(Math.log(randx), Trigonometry.log(randx), Math.pow(10, -12) );
	}
	
	@Test
	public void computeLogBase10Of938183() {
		double x = 938183;
		assertEquals(Math.log10(x), Trigonometry.log(x, 10), Math.pow(10, -12));
	}
	
	@Test
	public void sineNinetyDegreesShouldBeAround1() {
		assertEquals(1.0, Trigonometry.sin(Trigonometry.PI/2), 0.00001);
	}

	@Test
	public void sine270ShouldBeAroundMinus1() {
		assertEquals(-1.0, Trigonometry.sin(3*Trigonometry.PI/2), 0.00001);
	}
	
	
	@Test
	public void sineZeroDegreesShouldBeAround0() {
		assertEquals(0.0, Trigonometry.sin(0), 0.00001);
	}
	
	@Test
	public void cosZeroShouldBeAround1() {
		assertEquals(1.0, Trigonometry.cos(0), 0.00001);
	}

	@Test
	public void cosNinetyDegreesShouldBeAround0() {
		assertEquals(0.0, Trigonometry.cos(Trigonometry.PI/2), 0.0001);
	}
	
	@Test
	public void cos42PiShouldBeAround1() {
		assertEquals(1, Trigonometry.cos(Trigonometry.PI * 42), 0.0001);
	}
	@Test
	public void cos43PiShouldBeAroundMinus1() {
		assertEquals(-1, Trigonometry.cos(Trigonometry.PI * 43), 0.0001);
	}
	
	@Test
	public void xpow3shouldbe27() {
		assertEquals(27, Trigonometry.pow(3, 3));
	}
	
	
}
