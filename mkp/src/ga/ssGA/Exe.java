///////////////////////////////////////////////////////////////////////////////
///            Steady State Genetic Algorithm v1.0                          ///
///                by Enrique Alba, July 2000                               ///
///                                                                         ///
///   Executable: set parameters, problem, and execution details here       ///
///////////////////////////////////////////////////////////////////////////////

package ga.ssGA;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

public class Exe {

    public static final int NUM_EXECUTIONS_PER_COMBINATION = 30;
    public static final double[] CROSSOVER_PROBABILITIES = new double[]{0.6, 0.7, 0.8, 0.9, 1.0};
    public static final int[] NUMBER_AVG_EXPECTED_MUTATIONS = new int[]{1, 2, 3, 4, 5};
    public static final String OUTPUT_CSV_HEADER = "problem_index,population_size,crossover_probability," +
            "mutation_probability,search_optimal,max_steps,obtained_fitness,num_evaluations,execution_time\n";
    public static final long MAX_NUMBER_EVALUATIONS = 50000000L;  // Max number of evaluations before resetting when searching for optimal

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

    /**
     * Parses a list of lines (tokens) from the input file with the format specified in
     * http://people.brunel.ac.uk/~mastjjb/jeb/orlib/mknapinfo.html.
     * @param numberOfInstances
     * @param index if -1, all instances in the input file will be loaded; otherwise, the instance at the specified
     *              index (starting at 0) will be loaded.
     * @param lines list of lines from the input file, each one divided in tokens
     * @return list of MKP problem instances read from the input file
     */
    private static List<ProblemMKP> getMKPInstances(int numberOfInstances, int index, List<List<String>> lines) {
        List<ProblemMKP> problemInstances = new ArrayList<>();

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

            ProblemMKP problemParams = new ProblemMKP(numItems, numKnapsacks, optimalFitness, profitPerItem,
                    knapsacksCoefficientsMatrix, knapsacksCapacity);
            if (index == -1 || index == currentProblemIndex) {
                problemInstances.add(problemParams);
            }

            currentProblemIndex++;
        }

        return problemInstances;
    }

    public static void writeResultToCsvFile(String filePath, int problemIndex, int popSize, double pc, double pm,
                                            boolean searchForOptimal, long maxSteps, double obtainedFitness,
                                            long numEvaluations, double executionTime) {
        File f = new File(filePath);
        PrintWriter out = null;

        String header = "";
        try {
            if (f.exists() && !f.isDirectory()) {
                out = new PrintWriter(new FileOutputStream(f, true));
            } else {
                out = new PrintWriter(filePath);
                header = Exe.OUTPUT_CSV_HEADER;
            }
            if (!header.isEmpty()) {
                out.append(header);
            }
            String newLine = String.join(",", String.valueOf(problemIndex), String.valueOf(popSize),
                    String.valueOf(pc), String.valueOf(pm), String.valueOf(searchForOptimal), String.valueOf(maxSteps),
                    String.valueOf(obtainedFitness), String.valueOf(numEvaluations), String.valueOf(executionTime));
            out.append(newLine + "\n");
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        } finally {
            if (out != null) {
                out.close();
            }
        }
    }

    /**
     * Executes the ssGA algorithm on the provided problem instance
     * @param problem contains the problem instance information
     * @param gn gene number
     * @param gl gene length
     * @param popSize population size
     * @param pc crossover probability
     * @param pm mutation probability
     * @param tf target fitness value
     * @param maxSteps maximum number of steps (i.e., evaluations)
     * @param searchForOptimal true if the stopping criteria is reaching the target finess; false if the stopping
     *                         criteria is the provided maximum number of steps
     * @param outputFilePath path of the file where the result will be written
     * @throws Exception
     */
    private static void execute(Problem problem, int problemIndex, int gn, int gl, int popSize, double pc, double pm,
                                double tf, long maxSteps, boolean searchForOptimal,
                                String outputFilePath) throws Exception {
        problem.set_geneN(gn);
        problem.set_geneL(gl);
        problem.set_target_fitness(tf);

        Algorithm ga;          // The ssGA being used
        ga = new Algorithm(problem, popSize, gn, gl, pc, pm);

        long startTs = System.currentTimeMillis();  // Start timer

        boolean stop = problem.tf_known() && (ga.get_solution().get_fitness() >= problem.get_target_fitness());
        long step = 0L;
        boolean solutionFound = false;
        while (!stop) {
            ga.go_one_step();
            step++;
            solutionFound = problem.tf_known() && (ga.get_solution().get_fitness() >= problem.get_target_fitness());
            stop = (searchForOptimal && solutionFound) || (!searchForOptimal && step >= maxSteps);
            if (searchForOptimal && !stop && problem.get_fitness_counter() >= Exe.MAX_NUMBER_EVALUATIONS) {
                problem.resetFitnessCounter();
                step = 0L;
                ga = new Algorithm(problem, popSize, gn, gl, pc, pm);
                System.out.println("Resetting after not being able to find optimal solution in " +
                        Exe.MAX_NUMBER_EVALUATIONS + " evaluations");
                startTs = System.currentTimeMillis();  // Reset timer
            }
        }

        long endTs = System.currentTimeMillis();  // End timer

        long numberEvaluations = problem.get_fitness_counter();
        if (solutionFound) {
            System.out.println("Solution Found! After " + numberEvaluations + " evaluations");
        }
        double executionTime = (endTs-startTs)/1000.0;
        System.out.println("Elapsed time: " + executionTime + "s");

        Individual solution = ga.get_solution();
        double obtainedFitness = solution.get_fitness();

        // Print the solution
        for (int i = 0; i < gn * gl; i++) {
            System.out.print(solution.get_allele(i));
        }
        System.out.println();
        System.out.println(obtainedFitness);

        // Write results to the specified .csv file
        writeResultToCsvFile(outputFilePath, problemIndex, popSize, pc, pm, searchForOptimal, maxSteps, obtainedFitness,
                numberEvaluations, executionTime);
    }

    /**
     * This function will launch the execution of the specified MKP problem instances for each combination of crossover
     * and mutation 30 times.
     * @param args These are the arguments to be passed:
     *             - input file path with the problem instances
     *             - index (starting at 0) of the problem instance in the input file
     *             - population size
     *             - max number of steps
     *             - output file path (.csv)
     * @throws Exception
     */
    public static void main(String[] args) throws Exception {
        if (args.length < 5) {
            throw new Exception("Input file path, problem instance index, population size, max number of steps and " +
                    "results file path expected as arguments (in this order)");
        }

        String inputFilePath = args[0];

        int instanceIndex = Integer.parseInt(args[1]);
        if (instanceIndex < 0) {
            instanceIndex = -1;
        }

        int populationSize = Integer.parseInt(args[2]);
        if (populationSize < 10) {
            throw new IllegalArgumentException("Population size cannot be smaller than 10");
        }

        long maxSteps = Long.parseLong(args[3]);
        if (maxSteps < 0) {
            throw new IllegalArgumentException("Max number of steps cannot be smaller than 0");
        }

        String outputFilePath = args[4];

        List<List<String>> fileLines = Exe.readMKPFile(inputFilePath);
        int numProblemInstances = Integer.parseInt(fileLines.remove(0).get(0));
        if (instanceIndex >= numProblemInstances) {
            throw new IllegalArgumentException("Invalid index, make sure to specify one that exists given the " +
                    "number of instances in the file, starting at 0");
        }
        List<ProblemMKP> problemInstances = Exe.getMKPInstances(numProblemInstances, instanceIndex, fileLines);

        for (ProblemMKP instance : problemInstances) {
            System.out.println(instance);

            int geneNumber = instance.numItems;
            int geneLength = 1;
            double targetFitness = instance.optimalFitness;

            for (double crossoverProbability : Exe.CROSSOVER_PROBABILITIES) {
                for (int numAvgExpectedMutations : Exe.NUMBER_AVG_EXPECTED_MUTATIONS) {
                    double mutationProbability = numAvgExpectedMutations / ((double) geneNumber * (double) geneLength);
                    for (boolean searchForOptimal : new boolean[]{false, true}) {
                        System.out.println("\nCrossover probability: " + crossoverProbability +
                                ", mutation probability: " + mutationProbability + ", search for optimal: " +
                                searchForOptimal);
                        for (int i = 0; i < Exe.NUM_EXECUTIONS_PER_COMBINATION; i++) {
                            System.out.println("Execution " + (i+1) + "/" + Exe.NUM_EXECUTIONS_PER_COMBINATION);
                            Exe.execute(instance, instanceIndex, geneNumber, geneLength, populationSize,
                                    crossoverProbability, mutationProbability, targetFitness, maxSteps,
                                    searchForOptimal, outputFilePath);
                            instance.resetFitnessCounter();  // Set the number of evaluations back to 0 after execution
                        }
                    }
                }
            }
        }
    }

}
// END OF CLASS: Exe
