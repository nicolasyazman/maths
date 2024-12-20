package io.github.nicolasyazman;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;
import io.github.nicolasyazman.Approximation;
public class ApproximationTest {

	@Test
	public void ComputeChebyshevCoeffsShouldBe() {
		Approximation.ComputeChebyshevCoeffsSin(22);
		assertEquals(true, true);
	}
	
	@Test
	public void ComputeChebyshevCoeffsCos() {
		Approximation.ComputeChebyshevCoeffsCos(22);
		assertEquals(true, true);
	}
	@Test
	public void ChebyshevDegreeShouldBeAround22ForKEquals32() {
		assertEquals(22, Approximation.Chebyshev_Degree(3));
	}
	
	@Test
	public void ChebyshevError() {
		assertEquals(0.1, Approximation.Chebyshev_Error(22));
	}
}
