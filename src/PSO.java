import java.util.Random;

/**
* 
* Particle Swarm Optimisation Algorithm for Continuous Optimization
* Implemented by Fevzi Tugrul Varna: f.varna@sussex.ac.uk
*/

public class PSO {

	int Tmax = 10000; //max iterations
	int n = 40; // population
	int d = 100; // problem dimension
	
	final double LB = -10;
	final double UB = 10;
	final double vMax = 0.1*(UB-LB);
	final double vMin = -vMax;
	
	//PSO parameters
	final double w = 0.4;
	final double c1 = 2;
	final double c2 = 2;
	
	//particle attributes
	double[][] V = new double[n][d]; //velocity
	double[][] X = new double[n][d];  //position
	double[][] PX = new double[n][d]; //pbest
	double[] F = new double[n]; //fitness
	double[] PF = new double[n]; //pbest fitness
	double[] GX = new double[d]; //gbest
	double GF = Double.POSITIVE_INFINITY; //gbest fitness

	Random r = new Random();

	public PSO() {
		X = getInitialPositions(n, d, LB, UB);
		PX = X;
		
		for (int i = 0; i < n; i++) { // Evaluate initial fitness
			F[i] = getFitness(X[i]);
			PF[i] = F[i];
			updatePbest(i);
			updateGbest(i);
		}
		
		//Main loop
		for (int t = 0; t < Tmax; t++) { // for each time t
			
			for (int i = 0; i < n; i++) { // for each particle i

				double[] r1 = getRandomVector(d, 0, 1);
				double[] r2 = getRandomVector(d, 0, 1);
				
				for (int k = 0; k < d; k++) {
					V[i][k] = w * V[i][k] + c1 * r1[k] * PX[i][k] - X[i][k] + c2 *r2[k] * GX[k] - X[i][k];
				}
				
				V[i] = checkBounds(V[i], vMin, vMax);
				
				//update position
				for (int k = 0; k < d; k++) {
					X[i][k] = X[i][k] + V[i][k];
				}
				
				X[i] = checkBounds(X[i], LB, UB);
				F[i] = getFitness(X[i]);
				updatePbest(i);
				updateGbest(i);
			}
			System.out.println("Iteration = " + t + " - GBest = " + GF);
		}
	}
	
	/**
	 * update memory (pbest) of the particle.
	 * @param i - index of the particle
	 * */
	void updatePbest(int i) {
		if (F[i] < PF[i]) {
			PF[i] = F[i];
			PX[i] = X[i];
		}
	}
	
	/**
	 * update the best known (gbest) solution.
	 * @param i - index of the particle
	 * */
	void updateGbest(int i) {
		if (PF[i] < GF) {
			GF = PF[i];
			GX = PX[i];
		}
	}

	/**
	 * generates and returns the initial random positions for all particles
	 * 
	 * @param n - swarm size
	 * @param d - problem dimension
	 * @param LB - lower bound
	 * @param UB - upper bound
	 * @return d-dimensional positions of particles
	 * */
	double[][] getInitialPositions(int n, int d, double LB, double UB) {
		for (int i = 0; i < n; i++) {
			X[i] = getRandomVector(d, LB, UB);
		}
		return X;
	}

	/**
	 * returns a random solution vector within upper and lower bounds
	 */
	double[] getRandomVector(int d, double LB, double UB) {
		double[] v = new double[d];
		for (int i = 0; i < d; i++) {
			v[i] = LB + (UB - LB) * r.nextDouble();
		}
		return v;
	}

	/**
	 * objective function - Sphere function.
	 * @param x - position of particle
	 * @return fitness of particle
	 */
	double getFitness(double[] x) {
		double z = 0;
		for (int i = 0; i < d; i++) {
			z = z + Math.pow(x[i], 2);
		}
		return z;
	}
	
	/**
	 * corrects and returns a vector within the given bounds.
	 * @param x - vector (position or velocity)
	 * @param LB - upper bound
	 * @param UB - lower bound
	 * @return corrected vector within the provided upper and lower bounds
	 * */
	double[] checkBounds(double[] x, double min, double max){
		
		for (int i = 0; i < x.length; i++) {
			if (x[i] < min) {
				x[i] = min;
			}

			if (x[i] > max) {
				x[i] = max;
			}
		}
		return x;
	}
}