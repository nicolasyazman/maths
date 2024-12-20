package io.github.nicolasyazman;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import io.github.nicolasyazman.Trigonometry;
import org.junit.jupiter.api.Test;
import java.lang.Math;

public class TrigonometryTest {

	@Test
	public void computeAllCosOfAnglesOfThisPackagesTrigonometricFunctionsAndCompareToSuns() {
		double angle = 0.0;
		for (int i = 0; i < 1000; i++) {
			assertEquals(Math.cos(i* (Math.PI /180.0)), Trigonometry.cos(i* (Math.PI /180.0)), Math.pow(10,-10));
		}
	}
	
	@Test
	public void computeAllAnglesOfThisPackagesTrigonometricFunctionsAndCompareToSuns() {
		double angle = 0.0;
		for (int i = 0; i < 1000; i++) {
			assertEquals(Math.sin(i* (Math.PI /180.0)), Trigonometry.sine(i* (Math.PI /180.0)), Math.pow(10,-10));
		}
	}
	@Test
	public void sineNinetyDegreesShouldBeAround1() {
		assertEquals(1.0, Trigonometry.sine(Trigonometry.PI/2), 0.00001);
	}

	@Test
	public void sine270ShouldBeAroundMinus1() {
		assertEquals(-1.0, Trigonometry.sine(3*Trigonometry.PI/2), 0.00001);
	}
	
	
	@Test
	public void sineZeroDegreesShouldBeAround0() {
		assertEquals(0.0, Trigonometry.sine(0), 0.00001);
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
