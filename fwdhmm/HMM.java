import java.text.DecimalFormat;

public class HMM {
	int N, M, maxIters, iters;
	double oldLogProb;
	double a[][], b[][], pi[];

	public HMM(int states, int obsSym, int mx) {
		N = states;
		M = obsSym;

		a = new double[N][N];
		b = new double[N][M];
		pi = new double[N];
		double nProb = 1 / N;
		double mProb = 1 / M;
		for (int i = 0; i < N; ++i) {
			pi[i] = nProb;
			for (int j = 0; j < N; ++j) {
				a[i][j] = nProb;
			}
			for (int j = 0; j < M; ++j) {
				b[i][j] = mProb;
			}
		}

		maxIters = mx;
		iters = 0;
		oldLogProb = Double.MIN_VALUE;
	}

	public double[][] forward(int[] obsSeq) {

		int T = obsSeq.length;
		double[][] alpha = new double[T][N];
		double c[] = new double[T];
		// compute a0(i)
		c[0] = 0;
		for (int i = 0; i < N; i++) {
			alpha[0][i] = pi[i] * b[i][obsSeq[0]];
			c[0] += alpha[0][i];
		}
		// scale the a0[i]
		scale(0, alpha, c);
		// compute at(i)
		for (int t = 1; t < T; t++) {
			c[t] = 0;
			for (int i = 0; i < N; i++) {
				alpha[t][i] = 0;
				for (int j = 0; j < N; j++) {
					alpha[t][i] += alpha[t - 1][j] * a[j][i];

				}
				alpha[t][i] *= b[i][obsSeq[t]];
				c[t] += alpha[t][i];
			}
			// scale at(i)
			scale(t, alpha, c);
		}

		return alpha;
	}

	public void scale(int t, double[][] pass, double[] c) {
		if (c[t] != 0) {
			c[t] = 1 / c[t];
		}
		for (int i = 0; i < N; i++) {
			pass[t][i] = c[t] * pass[t][i];
		}
	}

	public double[][] backward(int[] obsSeq) {
		int T = obsSeq.length;

		double[][] beta = new double[T][N];

		for (int i = 0; i < N; i++)
			beta[T - 1][i] = 1;

		for (int t = T - 2; t >= 0; t--) {
			for (int i = 0; i < N; i++) {
				beta[t][i] = 0;
				for (int j = 0; j < N; j++)
					beta[t][i] += a[i][j] * b[j][obsSeq[t + 1]]
							* beta[t + 1][j];
			}
		}

		return beta;
	}

	public double gamma(int t, int i, int j, int[] obsSeq, double[][] alpha,
			double[][] beta) {
		double g = alpha[t][i] * a[i][j] * b[j][obsSeq[t + 1]] * beta[t + 1][j];
		return g;
	}

	// gamma(i,t)
	public double g_sum(int t, int i, int[] obsSeq, double[][] alpha,
			double[][] beta) {
		int T = obsSeq.length;
		double res = 0;
		if (t == T - 1)
			res = alpha[t][i];
		else
			for (int j = 0; j < N; j++)
				res += gamma(t, i, j, obsSeq, alpha, beta);
		return res;
	}

	public void reestimate(int[] obsSeq, int st) {
		int T = obsSeq.length;
		double[][] alpha;
		double[][] beta;
		double ar[][] = new double[N][N];
		double br[][] = new double[N][M];
		double pir[] = new double[N];
		for (int s = 0; s < st; s++) {
			alpha = forward(obsSeq);
			beta = backward(obsSeq);

			// re-estimation of pi (initial state probabilities)
			for (int i = 0; i < N; i++)
				pir[i] = g_sum(0, i, obsSeq, alpha, beta);

			// re-estimation of A (transition probabilities)
			double gamma[][][] = newGamma(obsSeq, alpha, beta);
			double gammaSum[][] = newGammaSum(obsSeq, alpha, beta, gamma);
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					double numer = 0;
					double denom = 0;
					for (int t = 0; t <= T - 2; t++) {
						//numer += gamma(t, i, j, obsSeq, alpha, beta);
						//denom += g_sum(t, i, obsSeq, alpha, beta);
						numer += gamma[t][i][j];
						denom += gammaSum[t][i];
					}
					try {
						if ( denom == 0.0 ) {
							ar[i][j] = 0;
						} else {
							ar[i][j] = numer / denom;
						}
					} catch (java.lang.ArithmeticException e) {
						ar[i][j] = 0;
						System.out.println(e.getMessage());
					}
				}
			}

			// re-estimation of B (emission probabilities)
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < M; j++) {
					double numer = 0;
					double denom = 0;

					for (int t = 0; t <= T - 1; t++) {
						if (obsSeq[t] == j)
							numer += gammaSum[t][i];
						denom += gammaSum[t][i];
						/*
						if (obsSeq[t] == j)
							numer += g_sum(t, i, obsSeq, alpha, beta);
						denom += g_sum(t, i, obsSeq, alpha, beta);
						*/
					}
					try {
						if ( denom == 0.0 ) {
							br[i][j] = 0;
						} else {
							br[i][j] = numer / denom;
						}
					} catch (java.lang.ArithmeticException e) {
						br[i][j] = 0;
						System.out.println(e.getMessage());
					}
				}
			}
			pi = pir;
			a = ar;
			b = br;
		}

	}

	private double[][][] newGamma(int[] obsSeq, double[][] alpha, double[][] beta) {
		int T = obsSeq.length;
		double[][][] ng = new double[T][N][N];
		for (int t = 0; t <= T - 2; ++t) {
			for (int i = 0; i <= N - 1; ++i) {
				for (int j = 0; j <= N - 1; ++j) {
					ng[t][i][j] = alpha[t][i] * a[i][j] * b[j][obsSeq[t + 1]]
							* beta[t + 1][j];
				}
			}
		}
		return ng;
	}

	private double[][] newGammaSum(int[] obsSeq, double[][] alpha, double[][] beta, double[][][] ng) {
		int T = obsSeq.length;
		double[][] ngs = new double[T][N];
		for (int t = 0; t <= T - 2; ++t) {
			for (int i = 0; i <= N - 1; ++i) {
				ngs[t][i] = 0;
				for (int j = 0; j <= N - 1; ++j) {
					ngs[t][i] += ng[t][i][j];
				}
			}
		}
		// Special case
		for (int i = 0; i <= N - 1; ++i ) {
			ngs[T-1][i] = alpha[T-1][i];
		}
		return ngs;
	}

	/** prints all the parameters of an HMM */
	public void print() {
		DecimalFormat fmt = new DecimalFormat();
		fmt.setMinimumFractionDigits(5);
		fmt.setMaximumFractionDigits(5);

		for (int i = 0; i < N; i++)
			System.out.println("pi(" + i + ") = " + fmt.format(pi[i]));
		System.out.println();

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++)
				System.out.print("a(" + i + "," + j + ") = "
						+ fmt.format(a[i][j]) + "  ");
			System.out.println();
		}

		System.out.println();
		for (int i = 0; i < N; i++) {
			for (int k = 0; k < M; k++)
				System.out.print("b(" + i + "," + k + ") = "
						+ fmt.format(b[i][k]) + "  ");
			System.out.println();
		}
	}
		
		public void logCompute( int logProb)
		{
		logProb = 0;
		int T = 0;
		for(int i = 0; i < T ; i++){
			Object[] c = null;
			logProb = logProb + log(c[i]);
		}
		logProb = -logProb;
		}

		private int log(Object object) 
		{
			return 0;
		}
		public void iteration( int iters, int oldLogProb, int mx)
		{ 
		int maxIters;
		double[][] a;
		double[][] b;
		int c;
		int pi;
		double ar[][] = new double[N][N];
		double br[][] = new double[N][M];
		a = ar;
		b = br;
		
		maxIters = mx;
		iters = 0;
		oldLogProb = Integer.MIN_VALUE;
		iters = iters + 1;
		int logProb = 0;
		int i = 0,j = 0;
		if(iters < maxIters && logProb > oldLogProb)
		{
		oldLogProb = logProb;
		}
		else 
		{
		oldLogProb =log (a[i][j]) + log ( b[i][j]);
		}
		}


		}

	
