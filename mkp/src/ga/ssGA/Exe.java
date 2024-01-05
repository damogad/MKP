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

    // Number of executions for each combination of crossover probability, mutation probability and study type (search
    // with limited evaluation or search for target fitness)
    public static final int NUM_EXECUTIONS_PER_COMBINATION = 30;
    // Fixed crossover probabilities to be tested
    public static final double[] CROSSOVER_PROBABILITIES = new double[]{0.6, 0.7, 0.8, 0.9, 1.0};
    // This constant is used to compute the mutation probability for each problem instance (num / chromosome length),
    // meaning the average expected number of mutated bits (genes) by chance
    public static final int[] NUMBER_AVG_EXPECTED_MUTATIONS = new int[]{1, 2, 3, 4, 5};
    // Header of the .csv output file with the execution results
    public static final String OUTPUT_CSV_HEADER = "problem_index,population_size,crossover_probability," +
            "mutation_probability,search_optimal,max_steps,optimal_fitness,optimal_percentage,best_fitness_last_pop," +
            "avg_fitness_last_pop,shannon_entropy,num_evaluations,execution_time\n";
    // Max number of evaluations before resetting when searching for optimal
    public static final long MAX_NUMBER_EVALUATIONS = 50000000L;

    /**
     * Splits an input line into tokens.
     * @param line input file line
     * @return list of tokens
     */
    private static List<String> getLineTokens(String line) {
        List<String> lineTokens = Arrays.asList(line.replaceAll("\\s+", " ").trim().split(" "));
        if (lineTokens.isEmpty() || (lineTokens.size() == 1 && lineTokens.get(0).isEmpty())) {
            return new ArrayList<>();
        }
        return lineTokens;
    }

    /**
     * Extracts all the lines of a mknap problem instances input file.
     * @param filePath input file path (.txt)
     * @return list of lines, each one tokenized
     * @throws java.io.FileNotFoundException
     */
    private static List<List<String>> readMKPFile(String filePath) throws java.io.FileNotFoundException {
        List<List<String>> lines = new ArrayList<>();
        try (Scanner scanner = new Scanner(new File(filePath))) {
            while (scanner.hasNextLine()) {
                List<String> lineTokens = getLineTokens(scanner.nextLine());
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
     * @param numberOfInstances indicates the total number of problem instances defined in the input file
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

    /**
     * Writes the results of a ssGA execution to the provided output file in .csv format. If the file does not exist,
     * it creates it.
     * @param filePath .csv output file path
     * @param problemIndex indicates the index (starting at 0) of the problem instance in the input file
     * @param popSize size of the ssGA population
     * @param pc crossover probability (search param)
     * @param pm mutation probability (search param)
     * @param searchForOptimal whether the ssGA has been executed with a fixed number of evaluations (false) or
     *                         targeting the provided fitness (optimal or relaxed)
     * @param maxSteps number of fixed evaluations (only relevant if searchForOptimal is false)
     * @param optimalFitness the problem instance's optimal fitness (only relevant if searchForOptimal is true). The
     *                       actual target fitness may be different (if optimalPercentage is < 100)
     * @param optimalPercentage percentage of the optimalFitness that was targeted (only relevant if searchForOptimal is
     *                          true)
     * @param bestFitnessLastPop best fitness obtained with a single individual present in the last population of the
     *                           ssGA
     * @param avgFitnessLastPop average fitness obtained in the (complete) last population of the ssGA
     * @param shannonEntropy average Shannon Entropy across bits (genes) of the last population's individuals
     * @param numEvaluations number of total executed evaluations (if searchForOptimal is false, it will be equal to
     *                       maxSteps + size of the population, as the evaluations on the first, randomly generated
     *                       population is not considered for the limit of evaluations)
     * @param executionTime ssGA execution time in seconds
     */
    public static void writeResultToCsvFile(String filePath, int problemIndex, int popSize, double pc, double pm,
                                            boolean searchForOptimal, long maxSteps, double optimalFitness,
                                            double optimalPercentage, double bestFitnessLastPop,
                                            double avgFitnessLastPop, double shannonEntropy, long numEvaluations,
                                            double executionTime) {
        File f = new File(filePath);
        PrintWriter out = null;

        String header = "";
        try {
            if (f.exists() && !f.isDirectory()) {
                out = new PrintWriter(new FileOutputStream(f, true));
            } else {
                out = new PrintWriter(filePath);
                header = OUTPUT_CSV_HEADER;
            }
            if (!header.isEmpty()) {
                out.append(header);
            }
            String newLine = String.join(",", String.valueOf(problemIndex), String.valueOf(popSize),
                    String.valueOf(pc), String.valueOf(pm), String.valueOf(searchForOptimal), String.valueOf(maxSteps),
                    String.valueOf(optimalFitness), String.valueOf(optimalPercentage),
                    String.valueOf(bestFitnessLastPop), String.valueOf(avgFitnessLastPop),
                    String.valueOf(shannonEntropy), String.valueOf(numEvaluations), String.valueOf(executionTime));
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
     * Computes the average Shannon Entropy of the bits (genes) in the population's individuals.
     * For example, if there were just these two individuals in the population generated by the provided ssGA Algorithm:
     * 01100110 and 01100110, the entropy would be 0. However, if they were 10011101 and 01100010, the entropy would be
     * the maximum, 1, as it would be the average of the entropy of each bit (gene).
     * @param ga Algorithm object which contains all the logic and the state of the ssGA
     * @return average Shannon Entropy across bits (genes) of the population's individuals. It's a synonym of how
     *         diverse the genotype is in the population, the higher (the closer to 1), the more diverse.
     */
    private static double computePopulationEntropy(Algorithm ga) {
        int chromosomeLength = ga.getChromosomeLength();
        int[] sumPositiveAllelesPerChromosome = new int[chromosomeLength];
        for (Individual individual : ga.getPopulation().getIndividuals()) {
            for (int i = 0; i < chromosomeLength; i++) {
                sumPositiveAllelesPerChromosome[i] += individual.get_allele(i);
            }
        }

        // Shannon Entropy
        double totalEntropy = 0.0;
        for (int i = 0; i < chromosomeLength; i++) {
            double p = sumPositiveAllelesPerChromosome[i]/(double)ga.getPopulationSize();
            if (p > 0 && p < 1) {
                totalEntropy += (-p) * (Math.log(p) / Math.log(2)) - (1-p) * (Math.log(1-p) / Math.log(2));
            }
        }

        return totalEntropy/(double)chromosomeLength;
    }

    /**
     * Generates the metrics of the provided algorithm, writes to the standard output and to the .csv output file.
     * @param ga Algorithm object which contains all the logic and the state of the ssGA
     * @param problemIndex indicates the index (starting at 0) of the problem instance in the input file
     * @param maxSteps maximum number of steps (i.e., evaluations) if searchForOptimal is false
     * @param searchForOptimal true if the stopping criteria is reaching the target finess; false if the stopping
     *                         criteria is the provided maximum number of steps
     * @param optimalPercentage percentage of the target fitness to be reached if searchForOptimal is true
     * @param outputFilePath path of the file where the result will be written
     * @param startTs Unix timestamp corresponding to the beginning of the Algorithm execution
     * @param endTs Unix timestamp corresponding to the end of the Algorithm execution
     * @param solutionFound whether the target fitness has been reached or not
     * @throws Exception
     */
    private static void getMetricsAndWriteResult(Algorithm ga, int problemIndex, long maxSteps,
                                                 boolean searchForOptimal, double optimalPercentage,
                                                 String outputFilePath, long startTs, long endTs,
                                                 boolean solutionFound) throws Exception {
        Problem problem = ga.getProblem();

        long numberEvaluations = problem.get_fitness_counter();
        Individual solution = ga.get_solution();
        double bestFitnessLastPop = solution.get_fitness();
        double avgFitnessLastPop = ga.get_avgf();
        if (solutionFound) {
            System.out.println("Solution Found! After " + numberEvaluations + " evaluations");

            // Print the solution
            for (int i = 0; i < ga.getChromosomeLength(); i++) {
                System.out.print(solution.get_allele(i));
            }
            System.out.println();
            System.out.println("Solution fitness: " + bestFitnessLastPop);
            System.out.println("Average population fitness: " + avgFitnessLastPop);
        }
        double shannonEntropy = computePopulationEntropy(ga);
        System.out.println("Population entropy: " + shannonEntropy);
        double executionTime = (endTs-startTs)/1000.0;
        System.out.println("Elapsed time: " + executionTime + "s\n");

        writeResultToCsvFile(outputFilePath, problemIndex, ga.getPopulationSize(), ga.getCrossoverProbability(),
                ga.getMutationProbability(), searchForOptimal, maxSteps, problem.get_target_fitness(),
                optimalPercentage, bestFitnessLastPop, avgFitnessLastPop, shannonEntropy, numberEvaluations,
                executionTime);
    }

    /**
     * Executes the ssGA algorithm on the provided problem instance
     * @param problem contains the problem instance information
     * @param problemIndex indicates the index (starting at 0) of the problem instance in the input file
     * @param gn gene number
     * @param gl gene length
     * @param popSize population size
     * @param pc crossover probability
     * @param pm mutation probability
     * @param tf target fitness value (optimal)
     * @param maxSteps maximum number of steps (i.e., evaluations) if searchForOptimal is false
     * @param searchForOptimal true if the stopping criteria is reaching the target fitness; false if the stopping
     *                         criteria is the provided maximum number of steps
     * @param optimalPercentage percentage of the target fitness to be reached if searchForOptimal is true
     * @param outputFilePath path of the file where the result will be written
     * @throws Exception
     */
    private static void execute(Problem problem, int problemIndex, int gn, int gl, int popSize, double pc, double pm,
                                double tf, long maxSteps, boolean searchForOptimal, double optimalPercentage,
                                String outputFilePath) throws Exception {
        problem.set_geneN(gn);
        problem.set_geneL(gl);
        problem.set_target_fitness(tf);

        double targetFitnessRelaxation = tf * (optimalPercentage) / 100.0;

        Algorithm ga;          // The ssGA being used
        ga = new Algorithm(problem, popSize, gn, gl, pc, pm);

        long startTs = System.currentTimeMillis();  // Start timer

        boolean stop = (searchForOptimal && problem.tf_known() && (ga.get_solution().get_fitness() >= targetFitnessRelaxation));
        long step = 0L;
        boolean solutionFound = false;
        long endTs;
        while (!stop) {
            ga.go_one_step();
            step++;
            solutionFound = problem.tf_known() && (ga.get_solution().get_fitness() >= targetFitnessRelaxation);
            stop = (searchForOptimal && solutionFound) || (!searchForOptimal && step >= maxSteps);

            // if max allowed number of evaluations is reached when searching for the optimal fitness
            if (searchForOptimal && !stop && problem.get_fitness_counter() >= MAX_NUMBER_EVALUATIONS) {
                endTs = System.currentTimeMillis();  // End timer
                getMetricsAndWriteResult(ga, problemIndex, maxSteps, true, optimalPercentage,
                        outputFilePath, startTs, endTs, false);

                // reset the search
                problem.resetFitnessCounter();
                step = 0L;
                ga = new Algorithm(problem, popSize, gn, gl, pc, pm);
                System.out.println("Resetting after not being able to find optimal solution in " +
                        MAX_NUMBER_EVALUATIONS + " evaluations");
                startTs = System.currentTimeMillis();  // Reset timer
            }
        }

        endTs = System.currentTimeMillis();  // End timer

        getMetricsAndWriteResult(ga, problemIndex, maxSteps, searchForOptimal, optimalPercentage, outputFilePath,
                startTs, endTs, solutionFound);
    }

    /**
     * This function will launch the execution of the specified MKP problem instances for each combination of crossover
     * and mutation (and type of search or study, i.e., with a fixed number of evaluation or targeting the optimal
     * fitness) 30 times for each.
     * @param args These are the arguments to be passed:
     *             - input file path with the problem instances
     *             - index (starting at 0) of the problem instance in the input file. It can also be negative to
     *               indicate that we want to run all the problems (which is not recommended unless you have several
     *               hours of computation time).
     *             - population size
     *             - max number of steps
     *             - output file path (.csv)
     *             - (optional) target percentage of the optimal fitness. If this is provided and is smaller than 100,
     *               only the optimal search will be executed.
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

        // Last argument is optional, and it is the percentage of the optimal that can be considered
        // as reaching the optimal. By default, 100%, that is, if argument is not provided we will
        // try to reach the optimal when searching for the optimal.
        double optimalPercentage = 100;
        if (args.length > 5) {
            optimalPercentage = Double.parseDouble(args[5]);
            if (optimalPercentage <= 0 || optimalPercentage > 100) {
                throw new IllegalArgumentException("Invalid percentage of optimal");
            }
        }

        List<List<String>> fileLines = readMKPFile(inputFilePath);
        int numProblemInstances = Integer.parseInt(fileLines.remove(0).get(0));
        if (instanceIndex >= numProblemInstances) {
            throw new IllegalArgumentException("Invalid index, make sure to specify one that exists given the " +
                    "number of instances in the file, starting at 0");
        }
        List<ProblemMKP> problemInstances = getMKPInstances(numProblemInstances, instanceIndex, fileLines);

        for (ProblemMKP instance : problemInstances) {
            System.out.println(instance);

            int geneNumber = instance.numItems;
            int geneLength = 1;
            double targetFitness = instance.optimalFitness;

            for (double crossoverProbability : CROSSOVER_PROBABILITIES) {
                for (int numAvgExpectedMutations : NUMBER_AVG_EXPECTED_MUTATIONS) {
                    double mutationProbability = numAvgExpectedMutations / ((double) geneNumber * (double) geneLength);
                    for (boolean searchForOptimal : new boolean[]{false, true}) {
                        // If a percentage is provided, and it's not 100, we will just do the search for the optimal
                        if (args.length > 5 && !searchForOptimal && optimalPercentage < 100) {
                            continue;
                        }
                        System.out.println("\nCrossover probability: " + crossoverProbability +
                                ", mutation probability: " + mutationProbability + ", search for optimal: " +
                                searchForOptimal + ", optimal percentage: " + optimalPercentage);
                        for (int i = 0; i < NUM_EXECUTIONS_PER_COMBINATION; i++) {
                            System.out.println("Execution " + (i+1) + "/" + NUM_EXECUTIONS_PER_COMBINATION);
                            execute(instance, instanceIndex, geneNumber, geneLength, populationSize,
                                    crossoverProbability, mutationProbability, targetFitness, maxSteps,
                                    searchForOptimal, optimalPercentage, outputFilePath);
                            instance.resetFitnessCounter();  // Set the number of evaluations back to 0 after execution
                        }
                    }
                }
            }
        }
    }

}
// END OF CLASS: Exe
