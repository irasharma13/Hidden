import java.util.*;
import java.io.*;

public class cake {
	// Initial states
	private int pi_min[];
	private int pi_max[];
	// Transitions
	private int a_min[][];
	private int a_max[][];

	// Outputs 
	private int b_min[][];
	private int b_max[][];

	// output symbols
	private static final int tiramisu = 0;
	private static final int italian_cookie = 1;
	private static final int the_green_one = 2;

	// States
	private static final int tira_pref = 0;
	private static final int italian_pref = 1;

	private void init() {
		pi_min = new int[2];
		pi_max = new int[2];

		pi_min[tira_pref] = 1;
		pi_max[tira_pref] = 1000;
		pi_min[italian_pref] = 0;
		pi_max[italian_pref] = 0;

		a_min = new int[2][2];
		a_max = new int[2][2];

		a_min[tira_pref][tira_pref] = 1;
		a_max[tira_pref][tira_pref] = 700;
		a_min[tira_pref][italian_pref] = 701;
		a_max[tira_pref][italian_pref] = 1000;
		a_min[italian_pref][tira_pref] = 1;
		a_max[italian_pref][tira_pref] = 500;
		a_min[italian_pref][italian_pref] = 501;
		a_max[italian_pref][italian_pref] = 1000;

		b_min = new int[2][3];
		b_max = new int[2][3];

		b_min[tira_pref][tiramisu] = 1;
		b_max[tira_pref][tiramisu] = 600;

		b_min[tira_pref][italian_cookie] = 601;
		b_max[tira_pref][italian_cookie] = 700;

		b_min[tira_pref][the_green_one] = 701;
		b_max[tira_pref][the_green_one] = 1000;

		b_min[italian_pref][tiramisu] = 1;
		b_max[italian_pref][tiramisu] = 100;

		b_min[italian_pref][italian_cookie] = 101;
		b_max[italian_pref][italian_cookie] = 800;

		b_min[italian_pref][the_green_one] = 801;
		b_max[italian_pref][the_green_one] = 1000;
	}

	private int[] gen_output(int size) {
		int[] data = new int[size];
		int rand = ((int) (Math.random() * 1000)) + 1;
		int state;

		for (state = 0; state < 2; state++) {
			if ((pi_min[state] <= rand) && (pi_max[state] >= rand))
				break;
		}

		for (int i = 0; i < size; i++) {
			rand = ((int) (Math.random() * 1000)) + 1;
			for (int symb = 0; symb < 3; symb++) {
				if ((b_min[state][symb] <= rand)
						&& (b_max[state][symb]) >= rand) {
					data[i] = symb;
					break;
				}
			}

			rand = ((int) (Math.random() * 1000)) + 1;
			for (int newState = 0; newState < 2; newState++) {
				if ((a_min[state][newState] <= rand)
						&& (a_max[state][newState] >= rand)) {
					state = newState;
					break;
				}
			}
		}
		return data;
	}

	/**
	 * java cake <sequence_size> and <steps>, where sequence_size is the length of the
	 * sequence to be generated, steps - re-estimation steps.
	 */
	public static void main(String args[]) {
		cake c = new cake();
		c.init();
		int[] o = c.gen_output(Integer.parseInt(args[0]));

		HMM hmm = new HMM(2, 3, 4);

		hmm.pi[0] = 1.0;
		hmm.pi[1] = 0.0;

		hmm.a[0][0] = 0.5;
		hmm.a[0][1] = 0.5;
		hmm.a[1][1] = 0.5;
		hmm.a[1][0] = 0.5;

		hmm.b[0][0] = 1.0 / 3.0;
		hmm.b[0][1] = 1.0 / 3.0;
		hmm.b[0][2] = 1.0 / 3.0;
		hmm.b[1][0] = 1.0 / 3.0;
		hmm.b[1][1] = 1.0 / 3.0;
		hmm.b[1][2] = 1.0 / 3.0;
		
		

		System.out.println("Initial Parameters:");
		hmm.print();

		hmm.reestimate(o, Integer.parseInt(args[1]));

		System.out.println();

		System.out.println("Trained Model:");
		hmm.print();
	}
}
