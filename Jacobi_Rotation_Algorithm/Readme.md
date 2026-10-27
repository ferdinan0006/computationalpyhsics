## Project 2 Readme

This folder contains all code for Project 2 of the course FYS3150/FYS4150. And is developed by:

-Johannes Schäfers
-Ferdinand Nilius
-Michał Kiedrzyński

This Project discusses the Jacobi rotation algorithm, which is implemented to solve an eigenvalue problem derived from the differential equation of a buckling beam. Different solvers (including a handwritten version of this algorithm and the armadillo implementation) are compared to the analytical solution. Furthermore the scaling behaviour in dependence on the matrix size is explored.

For every task, which needed programming effort, separate files are created. A short description of all task and how to run the code is listed below.
In order to start make sure to clone the repository and navigate into the folder Project_2 (where this Readme is also listed).

# Requirements

The Python scripts require Python 3 and the packages NumPy and Matplotlib.
The C++ code requires a compiler supporting C++20.

# Task 2: Calculation of Eigenvectors using armadillo

The Code for Task 2 is contained inside of Task2.cpp, it needs to be compiled and linked. (use e.g. `g++ Task2.cpp -o Task2.exe -larmadillo -llapack -lblas` to make sure this will link properly)

The output of this code should be a sequence of 6 1´s. The program calculates the numerical and anlytical eigenvectors and then checks if they agree to a tolerance of 10^(-5). If the eigenvectors do not agree to said precision, there will be a 0 instead of a 1.
The program will also create a .txt file containing all calculated eigenvectors and eigenvalues calculated using both the analytical solution and armadillo.


# Task 3 and 4: Implementation of the Jacobi Rotation Algorithm

The code for Task 3 and 4 is contained in a single cpp file (Task3_4.cpp). This file needs to be compiled and linked (use e.g. `g++ Task3_4.cpp -o Task3_4.exe -larmadillo -llapack -lblas`).

The file prints the following information to the Console:
For Task 3: The absolute value, row index and column index of the largest element.
For Task 4: A Statement if the algorithm converged or not and how many iterations were necessary to reach the set precision. The comparison if all of 6 eigenvectors, which were found using Jacobis rotation algorithm are also contained in the analytical eigenvectors. 

Another .txt file will be created with the eigenvectors and eigenvalues calculated using the Jacobi rotation algorithm and the analytical solution.

# Task 5: Scaling of the Jacobi Rotation Algorithm

Task 5 is contained in two separate C++ codes (Task5.cpp and task5b.cpp). One code runs the algorithm created in 4 for the matrix scheme specified in the problems before for different sizes, while the other does the same for a random matrix of different sizes.
To compile and link the code run the following:

`g++ Task5.cpp -o Task5.exe -larmadillo`
`./Task5.exe`
`g++ task5b.cpp -o task5b.exe -larmadillo`
`./task5b.exe`

The results, which are output into two .txt files can be plotted using two python scripts. To run these use the following commands:

`python3 task5a.py`
`python3 task5b.py`

If path dependency issues arise, paste the absolute paths into line 18 of the Python scripts.

# Task 6: Solving the equation and plotting of the solution

The Code for Task 6 is contained inside of task6.cpp, it needs to be compiled and linked. (use e.g. `g++ task6.cpp -o task6.exe -larmadillo -llapack -lblas` to make sure this will link properly)

The program outputs four text files containing the calculated eigenvalues and eigenvectors for the numerical and analytical solutions for N=10 and N=100 (4 files in total). The associated python script `task6.py` can then parse this data and generate the requested plots (use `python3 test6.py` to run).