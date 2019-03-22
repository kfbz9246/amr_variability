This repo provides the code for the Preventative Veterinary Medicine paper 'Estimation of Multidrug Resistance Variability in the National Antimicrobial Monitoring System' - doi.org/10.1016/j.prevetmed.2019.03.006

A docker container with the environment required to run this code can be found at: kfbz9246/analysis_environment:10.1016-j.prevetmed.2019.03.006

The 'amr_db' directory contains the code for building the sqlite database that the analysis code queries.
	- 'create_narms_database.py' is the program that does the database creation
	- the 'database_interface' directory is a package implementing an sqlalchemy interface to the sqlite database that the analysis code requires.
	
The 'data_analysis' directory contains the code for all the analyses in the paper.
	- 'overdispersion_factor_table.py' calculates all the values for the overdispersion tables.
	- 'overdispersion_difference_test.py' tests the differences between overdispersion values.
	- 'power_analysis.py' runs the power analyses.
	- 'q_q_plots.py' makes the qq-plots.