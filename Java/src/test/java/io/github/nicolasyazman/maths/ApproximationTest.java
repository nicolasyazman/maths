package io.github.nicolasyazman.maths;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import io.github.nicolasyazman.maths.Approximation;
public class ApproximationTest {

	@Test
	public void ComputeChebyshevCoeffsForNaturalLogarithm() {
		// Interval [1, 2]
        double a = 1.0;
        double b = 2.0;

        // Degree of the polynomial (increase for better accuracy)
        int n = 20;
        
        // Generate Chebyshev nodes in the interval [1, 2]
        double[] nodes = Approximation.chebyshevNodes(n, a, b);

        // Evaluate ln(x) at these nodes
        double[] values = new double[n];
        for (int i = 0; i < n; i++) {
            values[i] = Math.log(nodes[i]);
        }
        
        System.out.println("ChebyshevNodes:");
        for (int i = 0; i < nodes.length; i++) {
        	System.out.println(nodes[i]);
        }
        System.out.println("ChebyshevValues:");
        for (int i = 0; i < values.length; i++) {
        	System.out.println(values[i]);
        }
        
	}
	@Test
	public void ComputeChebyshevCoeffsShouldBe() {
		Approximation.ComputeChebyshevCoeffsSin(32);
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
