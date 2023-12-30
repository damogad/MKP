# Multidimensional Knapsack Problem

Subject: _Resolución de problemas con metaheurísticas_ (Máster Universitario en Investigación en Inteligencia Artificial, UIMP).

**Author: David Mora Garrido**

This code is based on the steady-state genetic algorithm (ssGA) that can be found [here](https://neo.lcc.uma.es/software/ssga/description.php).

In order to compile the source code into JVM bytecode, ensure that you are in the root of the complete project (**not in the root of the Java project**) and execute the following command:

```
javac mkp/src/ga/ssGA/*.java -d out/
```

Then, you can directly run the java bytecode in the JVM:

```
java -cp out/ ga.ssGA.Exe
```

Or you can build a `.jar` file:

```
jar cfe mkp.jar ga.ssGA.Exe -C out/ ga/ssGA/
```

Version used in to build and tun: **java 21.0.1 2023-10-17 LTS**

This is the recommended way to execute the algorithm (using the `.jar` build, which is also included as a release):

```
java -jar mkp.jar <input_file_path> <problem_index> <population_size> <number_evaluations> <output_file_path> <optimal_percentage(optional, default=100)>
```

The executions carried out in this study correspond to the following groups of parameters:

```
"./input_files/mknap1.txt" 2 100 1000 "./results/mknap1.csv" 100
"./input_files/mknap1.txt" 4 100 10000 "./results/mknap1.csv" 100
"./input_files/mknap1.txt" 6 200 1000000 "./results/mknap1.csv" 100
"./input_files/mknap1.txt" 6 200 1000000 "./results/mknap1.csv" 99.5
"./input_files/mknap1.txt" 6 200 1000000 "./results/mknap1.csv" 99
```