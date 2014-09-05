// tsClusters.h
// Authored by Alex Shows
// Released under the MIT License 
// (http://opensource.org/licenses/mit-license.php) 
//
// Definition of the tsClusters template class 
// Given a data set of N-dimensional values, find
// some M number of clusters
#ifndef _TS_CLUSTERS_H
#define _TS_CLUSTERS_H

#include <cmath>
#include <limits>
#include <thread>
#include <mutex>
#include <vector>
#include <memory>
#include <fstream>

/*******************
The idea here is to have a template class for N-dimensional arrays that
can be searched for M clusters.

Internally, each data point is comprised of some number of T values. 
Those T values are stored in a vector, which is inside a struct that pairs
the vector of T's with a cluster index to which the data point is currently
assigned.

Then all of these structs of vectors of T's are placed in a vector that 
is the entire dataset. Thus you can traverse the full data set, find a 
data point and look at its T values and cluster assigned. 
********************/

/*******************
TODO List: 
	-Add logging support [ IN PROGRESS ]
	-Migrate to a better logger (need timestamps and such)
	-Test data filling functions, bounds check, error handling
	-Add logical processor counting [ DONE ]
	-Add cluster centers data [ DONE ]
	-Add cluster starting position assignments [ DONE ]
	-Add minimum safe distance checks on cluster starting positions 
	-Add cluster assignment [ DONE ]
	-Add centroid computation [ DONE ]
	-Add checks for movement, whether a data point changed clusters
	-Add multi-threading and experiment with workload distribution 
********************/

/*******************
A template class for cluster analysis across an N-dimensional array
********************/
template <typename T> class tsClusters
{
public:
	tsClusters();
	tsClusters(const tsClusters&); // Copy constructor
	virtual ~tsClusters(); // Destructor
	tsClusters& operator=(const tsClusters&); // Assignment operator
	// TODO: What about a move operator?
	unsigned int fill_data_array(T* input, unsigned int size,  unsigned int stride);
	void set_number_of_clusters(unsigned int num_clusters);
	void initialize_clusters();
	void assign_clusters(); // For each data point, assign the closest cluster to it
	// For each cluster, recompute the position 
	// as the centroid of all associated data points
	void compute_centroids(); 
	// Return the number of data points that moved in the last round
	unsigned int get_num_data_points_moved(){ return data_points_moved; };
private:
	/* This is the internal data structure for storing each data point
	of N dimensions, where each data point has a cluster index to which
	it is assigned.	*/
	struct data_structure
	{
		std::vector<T> p; // The N-dimensional point
		unsigned int ci; // The cluster index
		T distance_squared; // Distance to the nearest cluster
	};

	/* A shared pointer to the data vector itself, of which there
	may be any number of points, each with its own assigned cluster */
	std::shared_ptr<std::vector<data_structure>> data;

	/* The clusters are referenced by index, so we just need a shared
	pointer to the vector of vectors (2 dimensional layout), thus
	clusters[0] is the first index, clusters[1] is the second, and so on */
	std::shared_ptr<std::vector<std::vector<T>>> clusters;
	
	/* Stride is number of number of dimensions to the data */
	unsigned int stride;

	/* This is not redundant, because changing this and then 
	initializing clusters will change the cluster data structure.
	Therefore, be careful not to change the cluster size _without_
	initializing clusters afterward. */
	unsigned int number_of_clusters;

	/*******************
	A lock for any thread that may need to write to the data,
	though this should only ever need to be the thread owning the
	object created with this template, it's here just in case
	********************/
	std::mutex tsLock;

	/* The log file */
	std::fstream log;

	unsigned int cpu_count = 0;

	/* Has a data point moved from one cluster assignment to another? */
	unsigned int data_points_moved = std::numeric_limits<unsigned int>::max();

	T compute_squared_distance(std::vector<T>& pointA, std::vector<T>& pointB);
};

/*
Default constructor
*/
template <typename T> tsClusters<T>::tsClusters()
{
	// By default we create this shared pointer, but we don't know the stride yet
	// until the data is filled
	data = std::make_shared<std::vector < data_structure >>(*(new std::vector< data_structure>));
	clusters = std::make_shared<std::vector<std::vector<T>>>(*(new std::vector<std::vector<T>>));

	stride = 0;
	number_of_clusters = 0;
	cpu_count = std::thread::hardware_concurrency(); // Logical processor count

#ifdef _DEBUG
	log.open("debug.log", std::fstream::out);
	log << "Initialized tsClusters debug log.";
#endif
}

template <typename T> tsClusters<T>::~tsClusters()
{
#ifdef _DEBUG
	log << std::endl;
	log << "tsClusters cleaning up and closing out the log file.";
	log.close();
#endif	
}

/*
Templated copy constructor
*/
template <typename T> tsClusters<T>::tsClusters(const tsClusters<T> &other)
{
#ifdef _DEBUG
	log << "tsClusters copy constructor called.";
#endif
}

/*
Templated assignment operator
*/
template <typename T> tsClusters<T>& tsClusters<T>::operator=(const tsClusters<T> &other)
{
#ifdef _DEBUG
	log << "tsClusters assignment operator called.";
#endif

	data.reset(new std::vector< data_structure>);
	clusters.reset(new std::vector<std::vector<T>>);
	tsLock = other.tsLock;
	return *this;
}

/*
Fill the data array with an array of variable type T of a given size, where stride is
the dimension of the array.
Obviously size%stride should be 0.
Returns the size of the internal data vector
*/
template <typename T> unsigned int tsClusters<T>::fill_data_array(T* input_data, unsigned int input_size, unsigned int input_stride)
{
	if(!input_data || !input_size || !input_stride)
		return 0;
	
	number_of_clusters = input_stride; // To begin, we assume this, but the user can change it
	stride = input_stride; // This should be internally consistent everywhere

	std::lock_guard<std::mutex> lock(tsLock);

	try
	{
		for (unsigned int i = 0; i < input_size; i++)
		{
			data_structure ds;
			ds.p.push_back(input_data[i]);
			// Fill this row to the stride
			while (ds.p.size() < stride)
			{				
				i++;
				ds.p.push_back(input_data[i]);
			}

			ds.ci = 0;
			ds.distance_squared = std::numeric_limits<T>::max();
			data->push_back(ds);
		}
	}
	catch (std::exception e)
	{
#ifdef _DEBUG
		log << "Exception in fill_data_array: " << e.what() << std::endl;
#endif
		return 0;
	}


#ifdef _DEBUG
	log << std::endl << std::endl;
	log << "Data points: " << std::endl;

	for(auto it=data->begin(); it!=data->end(); it++)
	{
		for (auto it_p = it->p.begin(); it_p != it->p.end(); it_p++)
			log << *it_p << "\t";

		log << std::endl;
	}
#endif

	return (unsigned int)(data->size() * data->begin()->p.size());
}

/*
Set the desired number of clusters.
If this isn't called, number of clusters will default to the number of dimensions. 
*/
template <typename T> void tsClusters<T>::set_number_of_clusters(unsigned int input_number)
{
#ifdef _DEBUG
	log << "Setting the number of clusters to " << input_number;
#endif

	if(input_number)
		number_of_clusters = input_number;
}

/*
Initalize the clusters to a new starting position, based on the min and max
of each dimension (of which there are N dimensions, where N is the stride)
TODO: Test this more thorougly
*/
template <typename T> void tsClusters<T>::initialize_clusters()
{
	if (!stride || !number_of_clusters)
		return;

	std::lock_guard<std::mutex> lock(tsLock);

	// Find the upper and lower bound of each dimension in the data vector
	T* ub = new T[stride];
	T* lb = new T[stride];

	for (unsigned int j = 0; j < stride; j++)
	{
		ub[j] = std::numeric_limits<T>::min();
		lb[j] = std::numeric_limits<T>::max();
	}

	std::vector<data_structure>::iterator it = data->begin();

	while (it != data->end())
	{
		std::vector<T>::iterator it_t = it->p.begin();
		unsigned int i = 0;
		while (it_t != it->p.end())
		{
			if (*it_t > ub[i])
				ub[i] = *it_t;
			if (*it_t < lb[i])
				lb[i] = *it_t;

			i++;
			it_t++;
		}

		it++;
	}

	// We should now have a lower and upper bound for every dimension in
	// the data, based on traversing all the data

	// For every cluster...
	for (unsigned int idx_c = 0; idx_c < number_of_clusters; idx_c++)
	{
		// Create a new cluster point...
		std::vector<T> cp;

		// ... and for every dimension of input (stride)...
		for (unsigned int idx_s = 0; idx_s < stride; idx_s++)
		{
			// ...put a random value into the clusters vector that is between the
			// lower bound and upper bound of this particular dimension
			// Using fmod from cmath as modulo is not defined for float
			cp.push_back((T)(std::fmod(rand(), (ub[idx_s] - lb[idx_s]))) + lb[idx_s]);

		}
		// Then add this point and move on to the next cluster...
		clusters->push_back(cp);
	}

	// TODO: Test for minimum safe distance...

#ifdef _DEBUG
	log << std::endl << std::endl;
	log << "Cluster starting positions:" << std::endl;

	unsigned int cluster_index = 0;

	// The clusters
	std::vector<std::vector<T>>::iterator it_c = clusters->begin();

	while (it_c != clusters->end())
	{
		log << std::endl;
		log << "Cluster " << cluster_index << ":" << std::endl;
		log << "     ";

		// The cluster points
		std::vector<T>::iterator it_cp = it_c->begin();

		while (it_cp != it_c->end())
		{
			log << *it_cp << " ";
			it_cp++;
		}

		cluster_index++;
		it_c++;
	}

	log << std::endl;
#endif
	
}

/*
For every data point, find the closest cluster to it, and assign that one to it.
*/
template <typename T> void tsClusters<T>::assign_clusters()
{
	data_points_moved = 0;

	T computed_distance = 0; // Accumulator for the (p1-q1)^2 part of the distance computation
	
	unsigned int current_cluster_index = 0; // For keeping track of which cluster we're checking
	unsigned int closest_cluster_index = 0; // For keeping track of which was the closest cluster so far
	T closest_cluster_distance = std::numeric_limits<T>::max();

	// For every data point in the data vector...
	for (auto dp_it = data->begin(); dp_it != data->end(); dp_it++)
	{
		current_cluster_index = 0;
		closest_cluster_distance = std::numeric_limits<T>::max();;

		// ...compare to every cluster point...
		for (auto cp_it = clusters->begin(); cp_it != clusters->end(); cp_it++)
		{
			computed_distance = 0;

			/*****************
			Basically we're doing an N-dimensional distance comparison, to find
			the index of the cluster that is the closest to the data point.
			Using a distance formula like this:
			distance(p,q) = sqrt((p1-q1)^2 + (p2-q2)^2 + ... + (pN-qN)^2)
			Except we don't need to bother with the expensive sqrt operation.
			So in this inner loop we're just accumulating the various (p1-q1)^2
			result for each cluster and data point T value.
			NOTE: The count of it_p's should match the count of it_cp's because
			of internal consistency checks for stride.
			******************/
			computed_distance = compute_squared_distance(dp_it->p, *cp_it);

			if (computed_distance < closest_cluster_distance)
			{
				closest_cluster_distance = computed_distance;
				closest_cluster_index = current_cluster_index;
			}

			current_cluster_index++;

		} // end for every cluster

		// Test and set the movement flag if the cluster changed this round
		// (this is tested and skipped if it's already true)
		if (dp_it->ci != closest_cluster_index)
			data_points_moved++;

		// Assign the cluster index to this data point
		dp_it->ci = closest_cluster_index;
		dp_it->distance_squared = closest_cluster_distance;
		
	} // end for every data point
}

/* 
Given a set of data points with clusters assigned, compute new cluster
positions as the centroid of all the points assigned to that cluster.
If a cluster has no points assigned, it needs to be moved to a new 
random location. 
*/
template <typename T> void tsClusters<T>::compute_centroids()
{
	if (!stride || !number_of_clusters)
		return;
		
	// Outside loop is each cluster by index
	for (unsigned int i = 0; i < number_of_clusters; i++)
	{		
		T* accum = new T[stride]; // accumulator for each data point
		memset(accum, 0, sizeof(T) * stride);

		unsigned int data_point_counter = 0;

		std::vector<data_structure>::iterator ds_iter = data->begin();
		while (ds_iter != data->end())
		{
			if (ds_iter->ci == i)
			{
				data_point_counter++; // Used later to compute the mean

				unsigned int Tcounter = 0;
				for (auto& iter : ds_iter->p)
				{
					accum[Tcounter] += iter;
					Tcounter++;
				}				
			}
			ds_iter++;
		} // End for each data point by iterator

		// Now update the ith cluster position as the mean
		// of the accumulator for each T value
		unsigned int j = 0;
		auto cp_iter = clusters->begin();
		// Move out to the iterator of the current index in the outside loop
		while (j != i && cp_iter != clusters->end())
		{
			j++;
			cp_iter++;
		}
		
		// Compute the mean of each accumulator and store it as the new
		// set of T values in the cluster
		j = 0;
		std::vector<T>::iterator Tval_iter = cp_iter->begin();
		while (Tval_iter != cp_iter->end())
		{
			*Tval_iter = accum[j] / data_point_counter;
			j++;
			Tval_iter++;
		}			

		delete accum;

	} // End for each cluster by index
	
}

/* Compute the squared distance, ignoring the expensive sqrt operation. 
Useful for comparing without worrying about the _actual_ distance of the two points.
If you want the _actual_ distance, just sqrt the return value of this. */
template <typename T> T tsClusters<T>::compute_squared_distance(std::vector<T>& pointA, std::vector<T>& pointB)
{
	T accum = 0;
	std::vector<T>::iterator it_a = pointA.begin();
	std::vector<T>::iterator it_b = pointB.begin();

	while (it_a != pointA.end() && it_b != pointB.end())
	{
		accum += (*it_a - *it_b) * (*it_a - *it_b);
		it_a++;
		it_b++;
	}

	return accum;
}

#endif // _TS_CLUSTERS_H