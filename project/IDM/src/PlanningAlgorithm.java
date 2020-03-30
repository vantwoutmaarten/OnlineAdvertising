import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.linear.LinearConstraintSet;
import org.apache.commons.math3.optim.linear.LinearObjectiveFunction;
import org.apache.commons.math3.optim.linear.Relationship;
import org.apache.commons.math3.optim.linear.SimplexSolver;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;


public class PlanningAlgorithm {
	
	public PlanningAlgorithm() {
		
	}
	
	public Solution solveUnconstrained(CMDP[] cmdps) {
		return solve(cmdps, Double.MAX_VALUE);
	}
	
	public Solution solve(CMDP[] cmdps, double costLimit) {
		int nAgents = cmdps.length;
		
		// assign variable IDs to state-action pairs
		int[][][] varIDs = new int[nAgents][][];
		int varCounter = 0;
		for(int i=0; i<nAgents; i++) {
			varIDs[i] = new int[cmdps[i].getNumStates()][cmdps[i].getNumActions()];
			for(int s=0; s<cmdps[i].getNumStates(); s++) {
				for(int a=0; a<cmdps[i].getNumActions(); a++) {
					varIDs[i][s][a] = varCounter;
					varCounter++;
				}
			}
		}
		int numVars = varCounter;
		
		// create objective function
		RealVector objectiveCoefficients = new ArrayRealVector(numVars);
		for(int i=0; i<nAgents; i++) {
			CMDP cmdp = cmdps[i];
			
			for(int s=0; s<cmdp.getNumStates(); s++) {
				for(int a=0; a<cmdp.getNumActions(); a++) {
					objectiveCoefficients.setEntry(varIDs[i][s][a], cmdps[i].getReward(s, a));
				}
			}
		}
		LinearObjectiveFunction objectiveFunction = new LinearObjectiveFunction(objectiveCoefficients, 0);
		
		// create constraint set
		Collection<LinearConstraint> constraints = new ArrayList<LinearConstraint>();
		
		// add flow conservation constraints
		for(int i=0; i<nAgents; i++) {
			CMDP cmdp = cmdps[i];
			
			for(int sNext=0; sNext<cmdp.getNumStates(); sNext++) {
				RealVector lhs = new ArrayRealVector(numVars);
				
				for(int aPrime=0; aPrime<cmdp.getNumActions(); aPrime++) {
					lhs.setEntry(varIDs[i][sNext][aPrime], 1.0);
				}
				
				for(int s=0; s<cmdp.getNumStates(); s++) {
					for(int a=0; a<cmdp.getNumActions(); a++) {
						lhs.addToEntry(varIDs[i][s][a], -1.0 * (cmdp.getDiscountFactor() * cmdp.getTransitionProbability(s, a, sNext)));
					}
				}
				
				double rhs = cmdp.getInitialState() == sNext ? 1.0 : 0.0;
				LinearConstraint constr = new LinearConstraint(lhs, Relationship.EQ, rhs);
				constraints.add(constr);
			}
		}
		
		// add cost constraint
		RealVector costLHS = new ArrayRealVector(numVars);
		for(int i=0; i<nAgents; i++) {
			CMDP cmdp = cmdps[i];
			
			for(int s=0; s<cmdp.getNumStates(); s++) {
				for(int a=0; a<cmdp.getNumActions(); a++) {
					costLHS.setEntry(varIDs[i][s][a], cmdp.getCost(s, a));
				}
			}
		}
		LinearConstraint costConstraint = new LinearConstraint(costLHS, Relationship.LEQ, costLimit);
		constraints.add(costConstraint);
		
		// add non-negative constraints
		for(int i=0; i<nAgents; i++) {
			CMDP cmdp = cmdps[i];
			
			for(int s=0; s<cmdp.getNumStates(); s++) {
				for(int a=0; a<cmdp.getNumActions(); a++) {
					RealVector lhs = new ArrayRealVector(numVars);
					lhs.setEntry(varIDs[i][s][a], 1.0);
					LinearConstraint constr = new LinearConstraint(lhs, Relationship.GEQ, 0.0);
					constraints.add(constr);
				}
			}	
		}
		
		// solve the problem using simplex
		SimplexSolver solver = new SimplexSolver();
		LinearConstraintSet lcs = new LinearConstraintSet(constraints);
		PointValuePair solution = solver.optimize(objectiveFunction, lcs, GoalType.MAXIMIZE);
		
		// compute expected reward and cost
		double expectedReward = 0.0;
		double expectedCost = 0.0;
		double[] expectedRewardAgent = new double[nAgents];
		double[] expectedCostAgent = new double[nAgents];
		for(int i=0; i<nAgents; i++) {
			CMDP cmdp = cmdps[i];
			
			for(int s=0; s<cmdp.getNumStates(); s++) {
				for(int a=0; a<cmdp.getNumActions(); a++) {
					int varID = varIDs[i][s][a];
					double flow = solution.getPoint()[varID];
					expectedReward += flow * cmdp.getReward(s, a);
					expectedCost += flow * cmdp.getCost(s, a);
					expectedRewardAgent[i] += flow * cmdp.getReward(s, a);
					expectedCostAgent[i] += flow * cmdp.getCost(s, a);
				}
			}
		}
		
		// get solution
		double[] solutionValues = new double[numVars];
		for(int v=0; v<numVars; v++) {
			solutionValues[v] = solution.getPoint()[v];
			
			if(Math.abs(solutionValues[v]) < 0.00000001) {
				solutionValues[v] = 0.0;
			}
		}
		
		// construct policy
		double[][][] policy = new double[nAgents][][];
		for(int i=0; i<nAgents; i++) {
			CMDP cmdp = cmdps[i];
			
			policy[i] = new double[cmdp.getNumStates()][cmdp.getNumActions()];
			
			for(int s=0; s<cmdp.getNumStates(); s++) {
				for(int a=0; a<cmdp.getNumActions(); a++) {
					double divisor = 0.0;
					for(int aPrime=0; aPrime<cmdp.getNumActions(); aPrime++) {
						int varID = varIDs[i][s][aPrime];
						divisor += solutionValues[varID];
					}
					
					int varID = varIDs[i][s][a];
					policy[i][s][a] = solutionValues[varID] / divisor;
				}
			}
		}
		
		ArrayList<double[][]> policies = new ArrayList<double[][]>();
		for(int i=0; i<nAgents; i++) {
			policies.add(policy[i]);
		}
		
		return new Solution(policies, expectedReward, expectedCost, expectedRewardAgent, expectedCostAgent);
	}
	
	public Solution solveVI(CMDP[] cmdps, double discountFactor) {
		CMDP cmdp = cmdps[0];
		
		// Array that hold maximum reward of all actions of a state for Q_{n}
		double[] V = new double[cmdp.getNumStates()];
		
		// TODO compute an optimal value function for the cmdp object
		double[][] Q = new double[cmdp.getNumStates()][cmdp.getNumActions()];
		double epsilon = 0.001;
		double delta;
		do {
			delta = 0;
			double[][] Q_next = new double[cmdp.getNumStates()][cmdp.getNumActions()];
			// Array that hold maximum reward of all actions of a state for Q_next = Q_{n+1}
			double[] V_next = new double[cmdp.getNumStates()];
			
			for (int s = 0; s < cmdp.getNumStates(); s++) {
				for (int a = 0; a < cmdp.getNumActions(); a++) {
					double sumValuesNextStep = 0;
					// Calculate sum of the returns of going to the next state.
					for (int s_next = 0; s_next < cmdp.getNumStates(); s_next++) {
						sumValuesNextStep += cmdp.getTransitionProbability(s, a, s_next) * V[s_next];
					}
					// Update expected reward for state s and action a
					Q_next[s][a] = cmdp.getReward(s, a) + discountFactor * sumValuesNextStep;
					// Update maximum reward for a state s
					if (V_next[s] < Q_next[s][a]) {
						V_next[s] = Q_next[s][a];
					}
					// Update delta
					delta = Math.max(delta, Q_next[s][a] - Q[s][a]);
				}
			}
			Q = Q_next;
			V = V_next;
			// Loop while max updated value is larger or equal to epsilon
		} while (delta >= epsilon);
		
		
		double[][] policy = new double[cmdp.getNumStates()][cmdp.getNumActions()];
		
		// TODO fill the policy array with probabilities
		// It shares probability evenly if the expected reward is equal
		for (int s = 0; s < cmdp.getNumStates(); s++) {
			// ArrayList of indices of maximum expected reward actions
			ArrayList<Integer> maxIndex = new ArrayList<Integer>();
			// maximum expected reward for any action for state s
			double maxValue = 0;
			for (int a = 0; a < cmdp.getNumActions(); a++) {
				// Update maxValue if return is higher and create new ArrayList of indices
				if (maxValue < Q[s][a]) {
					maxValue = Q[s][a];
					maxIndex = new ArrayList<Integer>();
					maxIndex.add(a);
				} else if (maxValue == Q[s][a]) {
					// Add to existing ArrayList if expected reward is equal
					maxIndex.add(a);
				}
			}
			// Assign probabilities equal for max expected reward indices
			for (Integer index : maxIndex) {
				policy[s][index] = 1.0 / (double) maxIndex.size();
			}
		}
		
		ArrayList<double[][]> policies = new ArrayList<double[][]>();
		policies.add(policy);
		
		double expectedReward = V[cmdp.getInitialState()];
		
		return new Solution(policies, expectedReward, 0.0, new double[]{expectedReward}, new double[1]);
	}
}
