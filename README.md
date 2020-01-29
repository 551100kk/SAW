# SAW
A Tool for Safety Analysis of Weakly-hard Systems

### Dependency

| Names                           | Weblinks                           |
| ------------------------------- | ---------------------------------- |
| GNU MPFR Library                | http://gmplib.org/                 |
| GNU MPFR Library                | http://www.mpfr.org/               |
| GNU Scientific Library          | http://www.gnu.org/software/gsl/   |
| GNU Linear Programming Kit      | http://www.gnu.org/software/glpk/  |
| Bison - GNU parser generator    | http://www.gnu.org/software/bison/ |
| Flex: The Fast Lexical Analyzer | http://flex.sourceforge.net/       |
| Gnuplot                         | http://www.gnuplot.info/           |
| Boost C++ Library               | https://www.boost.org/             |

### Installation

We provide three ways to install the program.

* The virtual machine where the tool is ready to run can be found at https://www.csie.ntu.edu.tw/~r08922054/SAW.ova. To reproduce the results of all the examples in the paper, please follow the instruction in ./artifact_evaluation_readme.txt

  * Account: saw, Password: saw
  * Path: /home/SAW
  * Precompiled version only. No dependency in this VM.
  * SHA1: 7d11f1fbe1ca1c65ae94f00a9cd1ad1e8502469d

* Build in your own environment.

  ```
  make
  ```

  **Warning:** We strongly suggest you using **g++-8** to build the project.

* Use the precompiled program in x86_64 linux environment.

  * **./saw_linux_x86_64**

### Usage

1. Create a model file as follow format.

   ```
   <state_dim> <input_dim> <grid_count>
   <state_var_names> <input_var_names>
   <state_ode.1>
   ...
   <state_ode.state_dim>
   <input_equa.1>
   ...
   <input_equa.input_dim>
   <period> <step_size>
   <m> <k>
   <safe_state.1>
   ...
   <safe_state.state_dim>
   <initial_state.1>
   ...
   <initial_state.state_dim>
   ```

2. You can also modify the configuration of **Flow\***.

   ```
   <order>
   <cutoff_threshold>
   <queue_size>
   <remainder_estimation>
   ```

3. Execute the verification tool with the model file.

   ```
   ./saw model.txt
   or
   ./saw_linux_x86_64 model.txt
   ```

   You can test the example models in **example/**.
   
   ```
   ./saw example/model1.txt
   or
   ./saw_linux_x86_64 example/model1.txt
   ```
   
4. The result of **example/model1** is as follow:

   ```
   [Info] Parsing model.
   [Info] Building FLOW* configuration.
   [Info] Building grids.
   [Info] Building one-step graph.
          Process: 100.00%
   [Success] Number of edges: 19354
   [Info] Building K-step graph.
   [Success] Start Region Size: 1908
             End Region: 1208
             Number of Edges: 102436
   [Info] Finding the largest closed subgraph.
   [Success] Safe Initial Region Size: 1622
   [Info] Calculating area.
          Initial state region: 4.000000
          Grids Intersection:   4.000000
   ```

    The program will plot the result region of first two dimensions to **output.svg**.

   ![output2](example/output1_(2,5).svg)




