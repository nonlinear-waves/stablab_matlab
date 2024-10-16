This is a repository to hold the code associated with the paper "The impact of a wind switch on traveling fronts in a temperature-based model  of fire propagation." The pre-print is available on the ArXiv. 

There are two Jupyter notebooks to calculate solution trajectories and intersections for intersections of type 1 and 2, as defined in the paper. There is also a Jupyter notebook to calculate the accumulated angle, angular solution trajectories and eigenvalue minimizer for solutions corresponding to type 1 intersections.

The accumulated angle code requires several helper .py functions, also included. It takes a profile solution as input. These profile solutions are computed in matlab and the code is included in the folder "Matlab profile code." To run the accumulated angle code, you first need to generate profile solution information by running "driver_save_profiles.m". This will generate a "data" folder (that you can rename as needed) with _info.csv and _y.csv files that you will need to manually sort in to "data_profile" and "data_info." There is currently a data folder from 09/23 that you can use to run the rest of the code.

