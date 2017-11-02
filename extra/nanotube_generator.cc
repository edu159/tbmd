// Here starts file map.cc that creates .xyz file for graphene edge or carbon nanotube of chosen size and chirality (zigzag/armchair). Output is directed to file map.xyz

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;
const long double PI = 3.1415926535897, BOND_LENGTH = 1.421; // in Å
const long double SIZE_A = 11, SIZE_B = 11; // in Å; 
// Cell dimension = nanotube radius + SIZE_A

int main()
{
int n, number_of_unit_cells;
char material, type;
vector<long double> centre_of_mass(2, 0);
long double radius;

cout << "Enter g/n (graphene/nanotube): ";
cin >> material;
cout << "Enter n: ";
cin >> n;
cout << "Enter a/z (armchair/zigzag): ";
cin >> type;
cout << "Enter number of repeating units: ";
cin >> number_of_unit_cells;

// Clean the map.xyz file before writing data there
ofstream to_clean("map.xyz");
to_clean << "";
to_clean.close();

// If nanotubes is chosen
if (material == 'n') {
	// If armchair is chosen
	if (type == 'a') {
		vector<long double> shell_one(8*n, 0);
		shell_one[0] = 0, shell_one[1] = 0;
		long double inner_angle = PI*(1.0-1.0/(2.0*n));
		long double angle_step = PI - inner_angle; 
		// Calculate and output positions in the x-y plane in armchair
		cout << "Map for (" << n  << ", " << n << ") nanotube" << endl;	
		for (int i = 1; i < 4*n; i++) {
			if (i%2 != 0) {
				shell_one[2*i] = shell_one[2*i - 2] + 
BOND_LENGTH*cos(i*angle_step);
				shell_one[2*i + 1] = shell_one[2*i - 1] + 
BOND_LENGTH*sin(i*angle_step);
				centre_of_mass[0] += shell_one[2*i];
				centre_of_mass[1] += shell_one[2*i + 1];
			}
			else {
				shell_one[2*i] = shell_one[2*i - 2] + 
0.5*BOND_LENGTH*cos(i*angle_step);
				shell_one[2*i + 1] = shell_one[2*i - 1] + 
0.5*BOND_LENGTH*sin(i*angle_step);
				centre_of_mass[0] += shell_one[2*i];
				centre_of_mass[1] += shell_one[2*i + 1];
			}
		}
		centre_of_mass[0] = centre_of_mass[0]/(4*n);
		centre_of_mass[1] = centre_of_mass[1]/(4*n);
		radius = sqrt(pow(centre_of_mass[0] - shell_one[0], 2) +
 pow(centre_of_mass[1] -  shell_one[1], 2));
		radius = radius;
		// Write data to file 
	std::ofstream to_write("map.xyz", std::ios::out | std::ios::app);
	to_write << n*number_of_unit_cells*4 << endl << endl;	
	for (int j = 0; j < number_of_unit_cells; j++) {
			for (int i = 0; i < 4*n; i++) {
				if (i%4 == 0) {
				to_write << "C " << shell_one[2*i] -
 centre_of_mass[0] + radius << 
		" " << shell_one[2*i + 1] - centre_of_mass[1] + radius 
<< " " << 
		+ BOND_LENGTH*(j - 0.25)*sqrt(3) << endl;
				to_write << "C " << shell_one[2*i + 2] 
- centre_of_mass[0] + radius << 
		" " << shell_one[2*i + 3] - centre_of_mass[1] + radius 
<< " " <<
		+ BOND_LENGTH*(j - 0.25)*sqrt(3) << endl;
				}
				else if ((i+2)%4 == 0) {
				to_write << "C " << shell_one[2*i] 
- centre_of_mass[0] + radius << 
		" " << shell_one[2*i + 1] - centre_of_mass[1] + radius
 << " "  << 
		BOND_LENGTH*(j + 0.25)*sqrt(3) << endl;
				to_write << "C " << shell_one[2*i + 2] 
- centre_of_mass[0] + radius << 
		" " << shell_one[2*i + 3] - centre_of_mass[1] + radius 
<< " " << BOND_LENGTH*(j + 0.25)*sqrt(3)  << endl;
				}
			}
		}	
		to_write.close();	
	
	}
	// If zig-zag is chosen
	else if (type == 'z') {
		vector<long double> shell_one(4*n, 0);
		shell_one[0] = 0, shell_one[1] = 0;
		long double inner_angle = PI*(1.0-1.0/(1.0*n));
		long double angle_step = PI - inner_angle; 
		// Calculate and output positionsin the x-y plane in zig-zag
		cout << "Map for (" << n  << ", " << 0 << ") nanotube" << endl;	
		for (int i = 1; i < 2*n; i++) {
			shell_one[2*i] = shell_one[2*i - 2] 
+ 0.5*sqrt(3)*BOND_LENGTH*cos(i*angle_step);
			shell_one[2*i + 1] = shell_one[2*i - 1] 
+ 0.5*sqrt(3)*BOND_LENGTH*sin(i*angle_step);
			centre_of_mass[0] += shell_one[2*i];
			centre_of_mass[1] += shell_one[2*i + 1];
		}
		centre_of_mass[0] = centre_of_mass[0]/(2*n);
		centre_of_mass[1] = centre_of_mass[1]/(2*n);
		radius = sqrt(pow(centre_of_mass[0] - shell_one[0], 2) 
+ pow(centre_of_mass[1] -  shell_one[1], 2));
		radius = radius;
		// Write data to file 
	std::ofstream to_write("map.xyz", std::ios::out | std::ios::app);
	to_write << n*number_of_unit_cells*4 << endl << endl;	
		for (int j = 0; j < number_of_unit_cells; j++) {
			for (int i = 0; i < 2*n; i++)
			{
				if (i%2 == 0) {
					to_write << "C " << shell_one[2*i] 
- centre_of_mass[0] + radius <<
		 " " << shell_one[2*i + 1] - centre_of_mass[1] + radius 
<< " " << (3*j - 1.25)*BOND_LENGTH << endl;
					to_write << "C " << shell_one[2*i] 
- centre_of_mass[0] + radius << 
		" " << shell_one[2*i + 1] - centre_of_mass[1] + radius 
<< " " << (3*j - 0.25)*BOND_LENGTH << endl;
				}
				else if ((i+1)%2 == 0) {
					to_write << "C " << shell_one[2*i] 
- centre_of_mass[0] + radius <<
		 " " << shell_one[2*i + 1] - centre_of_mass[1] + radius
 << " " << (3*j + 0.25)*BOND_LENGTH << endl;
					to_write << "C " << shell_one[2*i] 
- centre_of_mass[0] + radius << 
		 " " << shell_one[2*i + 1] - centre_of_mass[1] + radius 
<< " " << (3*j + 1.25)*BOND_LENGTH << endl;
				} 
			}	
		}
		to_write.close();
		}
}
// If graphene is chosen
else if (material == 'g') {
	// If armchair is chosen
	if (type == 'a') {
		vector<long double> shell_one(8*n, 0);
		shell_one[0] = 0, shell_one[1] = 0;
		long double inner_angle = PI;
		long double angle_step = 0; 
		// Calculate and output positions in the x-y plane in armchair
		cout << "Map for " << n  << "x" << number_of_unit_cells << 
" armchair supercell" << endl;	
		for (int i = 1; i < 4*n; i++) {
			if (i%2 != 0) {
				shell_one[2*i] = shell_one[2*i - 2] 
+ BOND_LENGTH;
				shell_one[2*i + 1] = shell_one[2*i - 1];
		
			}
			else {
				shell_one[2*i] = shell_one[2*i - 2] 
+ 0.5*BOND_LENGTH;
				shell_one[2*i + 1] = shell_one[2*i - 1];	
			}
		}
		// Write data to file 
	std::ofstream to_write("map.xyz", std::ios::out | std::ios::app);
	to_write << n*number_of_unit_cells*4 << endl << endl;	
		for (int j = 0; j < number_of_unit_cells; j++) {
			for (int i = 0; i < 4*n; i++) {
				if (i%4 == 0) {
				to_write << "C " << shell_one[2*i]  << 
		" " << 
		+ BOND_LENGTH*(j)*sqrt(3) << " " << 0  << endl;
				to_write << "C " << shell_one[2*i + 2]  << 
		" "  <<
		+ BOND_LENGTH*(j)*sqrt(3) << " " << 0 << endl;
				}
				else if ((i+2)%4 == 0) {
				to_write << "C " << shell_one[2*i] << 
	" " << BOND_LENGTH*(j + 0.5)*sqrt(3) << " " << 0  << endl ;
				to_write << "C " << shell_one[2*i + 2] << 
	" " <<  BOND_LENGTH*(j + 0.5)*sqrt(3) << " " << 0 << endl;
				}
			}
		}	
		to_write.close();	
	}

// If zigzag is chosen
else if (type == 'z') {
		vector<long double> shell_one(4*n, 0);
		shell_one[0] = shell_one[1] = 0;
		// Calculate and output positions in the x-y plane in zig-zag
		cout << "Map for " << n  << "x" << number_of_unit_cells << 
" zigzag supercell" << endl;	
		for (int i = 1; i < 2*n; i++) {
			if (i%2 != 0) {
				shell_one[2*i] = shell_one[2*i - 2] + 
0.5*BOND_LENGTH ;
				shell_one[2*i+1] = shell_one[2*i - 1] +
0.5*BOND_LENGTH*sqrt(3) ;
			}
			else {
				shell_one[2*i] = shell_one[2*i - 4] 
+ 1.5*BOND_LENGTH;
				shell_one[2*i + 1] = shell_one[2*i - 3] 
- 0.5*BOND_LENGTH*sqrt(3);	
			}
		}

		// Write data to file 
	std::ofstream to_write("map.xyz", std::ios::out | std::ios::app);
	to_write << n*number_of_unit_cells*2 << endl << endl;	
		for (int j = 0; j < number_of_unit_cells; j++) {
			for (int i = 0; i < 2*n; i++)
			{
				if (i%2 == 0) {
				to_write << "C " << shell_one[2*i]  <<
		 " " << shell_one[2*i + 1] + sqrt(3)*BOND_LENGTH*j 
<< " " << 0 << endl;
			}
				else if ((i+1)%2 == 0) {
				to_write << "C " << shell_one[2*i]   <<
		 " " << shell_one[2*i + 1]  + sqrt(3)*BOND_LENGTH*j 
 << " " << 0 << endl;
				} 
			}	
		}
		to_write.close();
		}
}	

else { 
cout << "Error: wrong n or configuration!" << endl;
}
return 0;
}

