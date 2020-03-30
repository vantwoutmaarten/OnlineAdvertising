import java.util.ArrayList;
import java.util.Random;
import org.apache.commons.math3.optim.linear.NoFeasibleSolutionException;


public class Homework {
	private static final Random rnd = new Random(222);
	
	// Example
	public static void task0() {
		// Get CMDP model for 1 agent
		CMDP cmdp = UserGenerator.getCMDPChild();
		CMDP[] cmdps = new CMDP[]{cmdp};
		
		// construct dummy policy that always executes action 4
		ArrayList<double[][]> policies = new ArrayList<double[][]>();
		double[][] policy = new double[cmdp.getNumStates()][cmdp.getNumActions()];
		for(int s=0; s<cmdp.getNumStates(); s++) {
			policy[s][4] = 1.0;
		}
		policies.add(policy);
		
		// construct dummy solution object with dummy expectations 0.0
		Solution solution = new Solution(policies, 0.0, 0.0, new double[]{0.0}, new double[]{0.0});
		
		// use the simulator to execute on run
		Simulator sim = new Simulator(rnd);
		sim.printActions();
		sim.simulate(cmdps, solution, 1);
	}
	
	// Solve unconstrained problem for 1 agent with value iteration
	public static void task1() {		
		// Get CMDP model for 1 agent
		CMDP cmdp = UserGenerator.getCMDPChild();
		CMDP[] cmdps = new CMDP[]{cmdp};
		
		// Solve the problem without constraints
		PlanningAlgorithm alg = new PlanningAlgorithm();		
		Solution solution = alg.solveVI(cmdps, 0.95);
		System.out.println("Expected reward: "+solution.getExpectedReward());
		System.out.println("Expected cost: "+solution.getExpectedCost());
		
		// Simulate solution
		System.out.println();
		Simulator sim = new Simulator(rnd);
		sim.simulate(cmdps, solution, 1000);
		
		// Print policy of agent 0
		int agentID = 0;
		double[][] policy = solution.getPolicy(agentID);
		System.out.println();
		for(int s=0; s<cmdps[agentID].getNumStates(); s++) {
			System.out.print("State "+s+": ");
			for(int a=0; a<cmdps[agentID].getNumActions(); a++) {
				System.out.print(policy[s][a]+" ");
			}
			System.out.println();
		}
		
	}
	
	// Solve unconstrained problem for 1 agent with value iteration
		public static void task1DF() {		
			// Get CMDP model for 1 agent
			CMDP cmdp = UserGenerator.getCMDPAdult();
			CMDP[] cmdps = new CMDP[]{cmdp};
			
			
			for (double gamma = 0.00; gamma <= 1; gamma += 0.05) {
				// Solve the problem without constraints
				PlanningAlgorithm alg = new PlanningAlgorithm();		
				Solution solution = alg.solveVI(cmdps, gamma);
				gamma = gamma * 10000;
				gamma = Math.round(gamma);
				gamma = gamma / 10000;
				double sol = solution.getExpectedReward() * 10000;
				sol = Math.round(sol);
				sol = sol / 10000;
				System.out.println("(" + gamma + "," + sol + ")");
			}
			
		}
	
	// Solve unconstrained problem for 1 agent with cost
	public static void task2() {		
		// Get CMDP model for 1 agent
		CMDP cmdp = UserGenerator.getCMDPChild();
		CMDP[] cmdps = new CMDP[]{cmdp};
		
		// Assign cost
		for (int s = 0; s < cmdp.getNumStates(); s++) {
			for (int a = 0; a < cmdp.getNumActions(); a++) {
				cmdp.assignCost(s, a, 2 * (a + 1));
			}
		}
		
		// Solve the problem without constraints
		PlanningAlgorithm alg = new PlanningAlgorithm();
		Solution solution = alg.solveUnconstrained(cmdps);
		System.out.println("Expected reward: "+solution.getExpectedReward());
		System.out.println("Expected cost: "+solution.getExpectedCost());
		
		// Simulate solution
		System.out.println();
		Simulator sim = new Simulator(rnd);
		sim.simulate(cmdps, solution, 1000);
		
		// Print policy of agent 0
		int agentID = 0;
		double[][] policy = solution.getPolicy(agentID);
		System.out.println();
		for(int s=0; s<cmdps[agentID].getNumStates(); s++) {
			System.out.print("State "+s+": ");
			for(int a=0; a<cmdps[agentID].getNumActions(); a++) {
				System.out.print(policy[s][a]+" ");
			}
			System.out.println();
		}
		
	}
	
	// Solve constrained problem for 1 agent
	public static void task3() {		
		// Get CMDP model for 1 agent
		CMDP cmdp = UserGenerator.getCMDPChild();
		CMDP[] cmdps = new CMDP[]{cmdp};
		
		// Assign cost
		for (int s = 0; s < cmdp.getNumStates(); s++) {
			for (int a = 0; a < cmdp.getNumActions(); a++) {
				cmdp.assignCost(s, a, 2 * (a + 1));
			}
		}

		PlanningAlgorithm alg = new PlanningAlgorithm();
		ArrayList<Double> expectedRewards = new ArrayList<>();
		for (int i = 1; i < 101; i++) {
			try {
				Solution solution = alg.solve(cmdps, i);
				expectedRewards.add(solution.getExpectedReward());
				System.out.println("(" + i + ", " + expectedRewards.get(expectedRewards.size() - 1) + ")");
			}
			catch (NoFeasibleSolutionException e) {
				System.out.println("(" + i + ", " + 0 + ")");
			}
		}
		
	}
	
	// Solve constrained problem for 2 agents with trivial budget split
	public static void task4() {		
		// Get CMDP models
		CMDP cmdpChild = UserGenerator.getCMDPChild();
		CMDP cmdpAdult = UserGenerator.getCMDPAdult();
		
		// Assign cost to child
		for (int s = 0; s < cmdpChild.getNumStates(); s++) {
			for (int a = 0; a < cmdpChild.getNumActions(); a++) {
				cmdpChild.assignCost(s, a, 2 * (a + 1));
			}
		}
		
		// Assign cost to adult
		for (int s = 0; s < cmdpAdult.getNumStates(); s++) {
			for (int a = 0; a < cmdpAdult.getNumActions(); a++) {
				cmdpAdult.assignCost(s, a, 2 * (a + 1));
			}
		}
		
		PlanningAlgorithm alg = new PlanningAlgorithm();
		Simulator sim = new Simulator(rnd);
		
		// Solve both problems separately without constraints and print expectations
		System.out.println("=========== UNCONSTRAINED ===========");
		for(int i=0; i<2; i++) {
			try {
				CMDP cmdp = (i==0) ? cmdpChild : cmdpAdult;
				Solution sol = alg.solveUnconstrained(new CMDP[]{cmdp});
				double expectedReward0 = sol.getExpectedReward();
				double expectedCost0 = sol.getExpectedCost();
				System.out.println("Expected reward agent "+i+": "+expectedReward0);
				System.out.println("Expected cost agent "+i+": "+expectedCost0);
			}
			catch (NoFeasibleSolutionException e) {
				System.out.println("(" + i + ", " + 0 + ")");
			}
			
		}
		
		// trivial budget split: invest 10 in each agent
		double budgetPerAgent = 10;
		
		System.out.println();
		System.out.println("=========== SEPARATE PLANNING ===========");
		double expectedReward = 0.0;
		double expectedCost = 0.0;
		for(int i=0; i<2; i++) {
			try {
				CMDP cmdp = (i==0) ? cmdpChild : cmdpAdult;
				Solution sol = alg.solve(new CMDP[]{cmdp}, budgetPerAgent); // TODO replace the number with the correct limit
				double expectedReward0 = sol.getExpectedReward();
				double expectedCost0 = sol.getExpectedCost();
				System.out.println("Expected reward agent "+i+": "+expectedReward0);
				System.out.println("Expected cost agent "+i+": "+expectedCost0);
				expectedReward += expectedReward0;
				expectedCost += expectedCost0;
			}
			catch (NoFeasibleSolutionException e) {
				System.out.println("(" + i + ", " + 0 + ")");
			}
			
		}
		System.out.println("Expected reward: "+expectedReward);
		System.out.println("Expected cost: "+expectedCost);
		
		// multi-agent problem: invest 20 in total
		Solution combinedSolution = alg.solve(new CMDP[]{cmdpChild, cmdpAdult}, budgetPerAgent * 2); // TODO replace the number with the correct limit
		System.out.println();
		System.out.println("=========== MULTI-AGENT PLANNING ===========");
		System.out.println("Expected reward: "+combinedSolution.getExpectedReward());
		System.out.println("Expected reward agent 0: "+combinedSolution.getExpectedReward(0));
		System.out.println("Expected reward agent 1: "+combinedSolution.getExpectedReward(1));
		System.out.println("Expected cost total: "+combinedSolution.getExpectedCost());
		System.out.println("Expected cost agent 0: "+combinedSolution.getExpectedCost(0));
		System.out.println("Expected cost agent 1: "+combinedSolution.getExpectedCost(1));
		
		// simulate
		sim.simulate(new CMDP[]{cmdpChild, cmdpAdult}, combinedSolution, 10000);
	}
	
	// Solve constrained problem for 2 agents with trivial budget split
		public static void task5() {	
			// Number of Children
			for (int i = 1; i < 50; i++) {
				// Number of Adults
				for (int j = 1; j < 50 && i + j <= 50; j++) {
					// Get CMDP models
					CMDP[] cmdps = new CMDP[i + j];
					for (int k = 0; k < i; k++) {
						cmdps[k] = UserGenerator.getCMDPChild();
						
						// Assign cost to child
						for (int s = 0; s < cmdps[k].getNumStates(); s++) {
							for (int a = 0; a < cmdps[k].getNumActions(); a++) {
								cmdps[k].assignCost(s, a, 2 * (a + 1));
							}
						}
					}
					for (int k = i; k < i + j; k++) {
						cmdps[k] = UserGenerator.getCMDPAdult();
						
						// Assign cost to adult
						for (int s = 0; s < cmdps[k].getNumStates(); s++) {
							for (int a = 0; a < cmdps[k].getNumActions(); a++) {
								cmdps[k].assignCost(s, a, 2 * (a + 1));
							}
						}
					}
					
					long startTime = System.nanoTime();
					
					PlanningAlgorithm alg = new PlanningAlgorithm();
					Simulator sim = new Simulator(rnd);
					
					// trivial budget split: invest 10 in each agent
					double budgetPerAgent = 10;
					
					// multi-agent problem: invest 20 in total
					try {
						Solution combinedSolution = alg.solve(cmdps, budgetPerAgent * cmdps.length); // TODO replace the number with the correct limit
						System.out.println();
						
						// simulate
						sim.simulate(cmdps, combinedSolution, 1000);
						
						double runtime = (System.nanoTime() - startTime) / 1000000000;
						System.out.println("Runtime with " + i + " Children and " + j + " Adults and Budget " + budgetPerAgent * cmdps.length + " and runtime " + runtime);
					}
					catch (NoFeasibleSolutionException e) {
						System.out.println("(" + i + ", " + 0 + ")");
					}
				}
			}
		}
	
	public static void main(String[] args) {
		task3();
	}
}
