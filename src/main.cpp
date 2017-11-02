#include "MDSimulation.h"
#include <iostream>

using namespace std;
                                                                            
void print_usage(){                                                         
	std::cout <<  "usage: tbmd -c nc -i nd -o no" << std::endl;
 	std::cout << "     -c nc: Name of the configuration file." << std::endl;
 	std::cout << "     -i nd: Name of the input data file." << std::endl;
 	std::cout << "     -o no: Name of the output data file." << std::endl;
}


int main (int argc, char *argv[]) {
	std::string conf_fname, data_fname, output_fname;
	int option, index;
	bool conf_flag = false, data_flag = false, output_flag = false, q_flag = false;
	if (argc == 1) {
		print_usage();
		exit(EXIT_FAILURE);
	}

	while ((option = getopt(argc, argv,"-qc:i:o:")) != -1) {
		switch (option) {
			case 'i' :  
				data_flag = true;
				data_fname = std::string(optarg);
				break;
			case 'c' :
				conf_flag = true;
				conf_fname = std::string(optarg);
			 	break;
			case 'o' : 
				output_flag = true;
				output_fname = std::string(optarg);
				break;
			case 'q':
				q_flag = true;
				break;
			default: 
				print_usage(); 
				exit(EXIT_FAILURE);
		}
	}
  for (index = optind; index < argc; index++) {
		printf ("Non-option argument %s\n", argv[index]);
		print_usage();
		exit(EXIT_FAILURE);
	}

	if (!output_flag | !conf_flag | !data_flag) {
		print_usage(); 
		exit(EXIT_FAILURE);
	}

	std::cout << std::setprecision(16);
	MDSimulation mds(data_fname, conf_fname, output_fname);	
	mds.quietMode(q_flag);
	mds.run();
}
