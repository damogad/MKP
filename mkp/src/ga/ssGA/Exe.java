///////////////////////////////////////////////////////////////////////////////
///            Steady State Genetic Algorithm v1.0                          ///
///                by Enrique Alba, July 2000                               ///
///                                                                         ///
///   Executable: set parameters, problem, and execution details here       ///
///////////////////////////////////////////////////////////////////////////////

package ga.ssGA;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

public class Exe {

    private static List<String> getLineTokens(String line) {
        List<String> lineTokens = Arrays.asList(line.replaceAll("\\s+", " ").trim().split(" "));
        if (lineTokens.isEmpty() || (lineTokens.size() == 1 && lineTokens.get(0).isEmpty())) {
            return new ArrayList<>();
        }
        return lineTokens;
    }

    private static List<List<String>> readMKPFile(String filePath) throws java.io.FileNotFoundException {
        List<List<String>> lines = new ArrayList<>();
        try (Scanner scanner = new Scanner(new File(filePath))) {
            while (scanner.hasNextLine()) {
                List<String> lineTokens = Exe.getLineTokens(scanner.nextLine());
                if (!lineTokens.isEmpty()) {
                    lines.add(lineTokens);
                }
            }
        }
        return lines;
    }

    private static List<MKPParams> getMKPInstances(int numberOfInstances, int index, List<List<String>> lines) {
        List<MKPParams> problemInstances = new ArrayList<>();

        // We will parse the complete file, even though we may have specified to execute a single problem instance
        int currentProblemIndex = 0;
        while (currentProblemIndex < numberOfInstances) {
            List<String> generalProblemInfo = lines.remove(0);
            int numItems = Integer.parseInt(generalProblemInfo.get(0));
            int numKnapsacks = Integer.parseInt(generalProblemInfo.get(1));
            double optimalFitness = Double.parseDouble(generalProblemInfo.get(2));

            List<String> profitPerItemLine = lines.remove(0);
            double[] profitPerItem = new double[numItems];
            int multipleLinesOffset = 0;  // The information may span multiple lines if too many items
            for (int j = 0; j < numItems; j++) {
                if ((j - multipleLinesOffset) >= profitPerItemLine.size()) {
                    multipleLinesOffset += profitPerItemLine.size();
                    profitPerItemLine = lines.remove(0);
                }
                profitPerItem[j] = Double.parseDouble(profitPerItemLine.get(j-multipleLinesOffset));
            }

            double[][] knapsacksCoefficientsMatrix = new double[numKnapsacks][numItems];
            for (int i = 0; i < numKnapsacks; i++) {
                List<String> knapsackCoefficientsLine = lines.remove(0);
                multipleLinesOffset = 0;  // The information may span multiple lines if too many items
                for (int j = 0; j < numItems; j++) {
                    if ((j - multipleLinesOffset) >= knapsackCoefficientsLine.size()) {
                        multipleLinesOffset += knapsackCoefficientsLine.size();
                        knapsackCoefficientsLine = lines.remove(0);
                    }
                    knapsacksCoefficientsMatrix[i][j] = Double.parseDouble(knapsackCoefficientsLine.get(j-multipleLinesOffset));
                }
            }

            List<String> knapsacksCapacityLine = lines.remove(0);
            double[] knapsacksCapacity = new double[numKnapsacks];
            multipleLinesOffset = 0;  // The information may span multiple lines if too many knapsacks
            for (int i = 0; i < numKnapsacks; i++) {
                if ((i - multipleLinesOffset) >= knapsacksCapacityLine.size()) {
                    multipleLinesOffset += knapsacksCapacityLine.size();
                    knapsacksCapacityLine = lines.remove(0);
                }
                knapsacksCapacity[i] = Double.parseDouble(knapsacksCapacityLine.get(i-multipleLinesOffset));
            }

            MKPParams problemParams = new MKPParams(numItems, numKnapsacks, optimalFitness, profitPerItem,
                    knapsacksCoefficientsMatrix, knapsacksCapacity);
            if (index == -1 || index == currentProblemIndex) {
                problemInstances.add(problemParams);
            }

            currentProblemIndex++;
        }

        return problemInstances;
    }

    private static void execute(int gn, int gl, int popsize, double pc, double pm, double tf, long max_isteps) throws Exception {
        /*
        // PARAMETERS PPEAKS
        int    gn         = 512;                           // Gene number
        int    gl         = 1;                            // Gene length
        int    popsize    = 512;                          // Population size
        double pc         = 0.8;                          // Crossover probability
        double pm  = 1.0/(double)((double)gn*(double)gl); // Mutation probability
        double tf         = (double)1 ;              // Target fitness beign sought
        long   MAX_ISTEPS = 50000;
    */
    /*
        // PARAMETERS ONEMAX
        int gn = 512;                          // Gene number
        int gl = 1;                            // Gene length
        int popsize = 512;                          // Population size
        double pc = 0.8;                          // Crossover probability
        double pm = 1.0 / (double) ((double) gn * (double) gl); // Mutation probability
        double tf = (double) gn * gl;           // Target fitness being sought
        long MAX_ISTEPS = 50000;
     */

        Problem problem;                             // The problem being solved

        // problem = new ProblemPPeaks();
        problem = new ProblemOneMax();

        problem.set_geneN(gn);
        problem.set_geneL(gl);
        problem.set_target_fitness(tf);

        Algorithm ga;          // The ssGA being used
        ga = new Algorithm(problem, popsize, gn, gl, pc, pm);

        for (int step = 0; step < max_isteps; step++) {
            ga.go_one_step();
            System.out.print(step);
            System.out.print("  ");
            System.out.println(ga.get_bestf());

            if ((problem.tf_known()) &&
                    (ga.get_solution()).get_fitness() >= problem.get_target_fitness()
            ) {
                System.out.print("Solution Found! After ");
                System.out.print(problem.get_fitness_counter());
                System.out.println(" evaluations");
                break;
            }
        }

        // Print the solution
        for (int i = 0; i < gn * gl; i++)
            System.out.print((ga.get_solution()).get_allele(i));
        System.out.println();
        System.out.println((ga.get_solution()).get_fitness());
    }

    public static void main(String[] args) throws Exception {
        if (args.length < 1) {
            throw new Exception(".csv file path expected as first argument");
        }

        String filePath = args[0];
        int instanceIndex = -1;  // By default, we assume that all problems in the file will be executed
        if (args.length > 1) {
            instanceIndex = Integer.parseInt(args[1]);
        }

        List<List<String>> fileLines = Exe.readMKPFile(filePath);
        int numProblemInstances = Integer.parseInt(fileLines.remove(0).get(0));
        if (instanceIndex >= numProblemInstances) {
            throw new IllegalArgumentException("Invalid index, make sure to specify one that exists given the " +
                    "number of instances in the file, starting at 0");
        }
        List<MKPParams> problemInstances = Exe.getMKPInstances(numProblemInstances, instanceIndex, fileLines);

        for (MKPParams instance : problemInstances) {
            System.out.println(instance);
        }
    }

}
// END OF CLASS: Exe
