/*******************

test_harness.cpp
Authored by Alex Shows
Released under the MIT License 
(http://opensource.org/licenses/mit-license.php) 

Test harness for the tsClusters template class

********************/

#include "tsClusters.h"

#include <random>
#include <iostream>

#define TS_DIMENSIONS 5
#define TS_DATAPOINTS 1000

/*******************
Main application entry point
********************/
int main(int argc, char* argv[])
{
	srand(0xDEADBEEF);

	tsClusters<float> clusters;

	float* data_array = new float [TS_DIMENSIONS * TS_DATAPOINTS];
	
	for(int i=0; i< (TS_DIMENSIONS * TS_DATAPOINTS); i+=5)
	{
		if((i+4) < (TS_DIMENSIONS * TS_DATAPOINTS))
		{
			data_array[i] = rand()%30 + 30.f;
			data_array[i+1] = rand()%100 + 50.f;
			data_array[i+2] = rand()%50 + 100.f;
			data_array[i+3] = rand()%150 + 25.f;
			data_array[i+4] = rand()%10 + 10.f;
		}
	}

	// TODO: Built more robust tests of these functions
	unsigned int num_data_points = clusters.fill_data_array(data_array, TS_DIMENSIONS * TS_DATAPOINTS, TS_DIMENSIONS);
	clusters.initialize_clusters();

	int round_counter = 0;
	bool keep_running = true;
	while (keep_running)
	{	
		round_counter++;

		clusters.assign_clusters();
		clusters.compute_centroids();
		unsigned int moved = clusters.get_num_data_points_moved();
		std::cout << "Round " << round_counter << " complete: " << moved << " points moved." << std::endl;

		if (moved == 0)
			keep_running = false;
	}
	
	std::cout << std::endl;
	std::cout << "Convergence complete in " << round_counter << " rounds!" << std::endl;
	std::cout << "Least sum of squares found for the data set given." << std::endl;
	std::cout << std::endl;
	std::cout << "Press Enter to Exit." << std::endl;
	std::cin.get();

	return 0;
}