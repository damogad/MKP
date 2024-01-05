package ga.ssGA;

import java.util.Arrays;

public class ProblemMKP extends Problem {
    // number of variables ('items')
    int numItems;

    // number of constraints ('knapsacks')
    int numKnapsacks;

    // optimal solution, if available
    double optimalFitness;

    // profit per Item
    double[] profitPerItem;

    // capacity (number of available resources) for each knapsack and item
    // first dimension -> knapsack index
    // second dimension -> item index
    double[][] knapsacksCoefficients;

    // total capacity (number of available resources) per knapsack
    double[] knapsacksCapacity;

    public ProblemMKP(int numItems, int numKnapsacks, double optimalFitness, double[] profitPerItem,
                      double[][] knapsacksCoefficients, double[] knapsacksCapacity) {
        this.numItems = numItems;
        this.numKnapsacks = numKnapsacks;
        this.optimalFitness = optimalFitness;
        this.profitPerItem = profitPerItem;
        this.knapsacksCoefficients = knapsacksCoefficients;
        this.knapsacksCapacity = knapsacksCapacity;
    }

    public String toString() {
        StringBuilder str = new StringBuilder("Items: " + this.numItems + ", knapsacks: " + this.numKnapsacks);
        if (this.optimalFitness > 0) {
            str.append(", optimal solution fitness: ");
            str.append(this.optimalFitness);
        } else {
            str.append(", no optimal solution fitness provided");
        }
        str.append("\nProfit per item: ");
        str.append(Arrays.toString(this.profitPerItem));
        str.append("\nKnapsacks coefficients (available resources) per item:\n");
        for (double[] knapsackCoefficients : this.knapsacksCoefficients) {
            str.append(Arrays.toString(knapsackCoefficients));
            str.append("\n");
        }
        str.append("Knapsacks total capacity: ");
        str.append(Arrays.toString(this.knapsacksCapacity));
        str.append("\n");
        return str.toString();
    }

    private boolean isFeasible(Individual indiv) {
        for (int i = 0; i < this.numKnapsacks; i++) {
            double knapsackCapacity = this.knapsacksCapacity[i];
            double indivAllocated = 0.0;
            for (int j = 0; j < this.numItems; j++) {
                if (indiv.get_allele(j) == (byte)1) {
                    indivAllocated += this.knapsacksCoefficients[i][j];
                }
            }
            if (indivAllocated > knapsackCapacity) {
                return false;
            }
        }
        return true;
    }

    private double computeFitness(Individual indiv) {
        double fitness = 0.0;
        for (int j = 0; j < this.numItems; j++) {
            if (indiv.get_allele(j) == (byte)1) {
                fitness += this.profitPerItem[j];
            }
        }
        indiv.set_fitness(fitness);
        return fitness;
    }

    @Override
    public double Evaluate(Individual Indiv) {
        if (this.isFeasible(Indiv)) {
            return this.computeFitness(Indiv);
        }
        return -1;  // In order not to be mixed up with the empty solution (all 0s)
    }
}