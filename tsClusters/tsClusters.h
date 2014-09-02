// tsClusters.h
// Authored by Alex Shows
// Released under the MIT License 
// (http://opensource.org/licenses/mit-license.php) 
//
// Definition of the tsClusters template class 
// Given a data set of N-dimensional values, find
// some X number of clusters
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
can be searched for X clusters.

Internally, while it would be nice to store the data as a vector of 
vectors of vectors...ad infinitum...this is intractible and would look
silly even with typedefs.

Instead, I'm trying a single vector with a stride component so that
any number of T-type values can be copied into the class and traversed
by the stride provided in order to run the cluster analysis. 
********************/

/*******************
PRIORITY TODO
=============
I just realized the fatal flaw in this design. :-(
This can't be a template T where T is in a vector, because there's not enough
descriptive information about cluster assignment. I either need to build a vector
of a custom data structure using T (probably the better path), or build and
maintain a mapping between data and clusters. Every data point needs to have a
cluster to which it is assigned.

TODO List: 
	-Add logging support [ IN PROGRESS ]
	-Migrate to a better logger (need timestamps and such)
	-Test data filling functions, bounds check, error handling
	-Add cluster centers data [ DONE ]
	-Add cluster starting position assignments [ IN PROGRESS ]
	-Add cluster assignment
	-Add centroid computation
	-Add logical processor counting [ DONE ] 
	-Add multi-threading and experiment with workload distribution 
			(probably divide the entire dataset by thread along stride boundaries?)
	-Add distance computation (test for instruction support, use intrinsics)
	-Try out sum of squares with/without the sqrt (dot product the array with itself)
	-Performance, performance, performance

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
	void tsClusters<T>::initialize_clusters();
private:
	/* This is the internal data structure for storing each data point
	of N dimensions, where each data point has a cluster index to which
	it is assigned.	*/
	struct data_structure
	{
		std::vector<T> p; // The N-dimensional point
		unsigned int ci; // The cluster index
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

	std::lock_guard<std::mutex> lock(tsLock);

	try
	{
		for (unsigned int i = 0; i < input_size; i++)
		{
			data_structure ds;
			while ((i % (input_stride + 1)) < input_stride)
			{
				ds.p.push_back(input_data[i]);
				i++;
			}

			ds.ci = 0;
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
	
	number_of_clusters = input_stride;
	stride = input_stride;

#ifdef _DEBUG
	log << std::endl << std::endl;
	log << "Data points: " << std::endl;

	unsigned int count = 0;
	unsigned int row = 0;

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
TODO: Test this
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
		while (it_t != it->p.begin())
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
			log << *it_cp << ", ";
			it_cp++;
		}

		cluster_index++;
		it_c++;
	}

	log << std::endl;
#endif

}

#endif // _TS_CLUSTERS_H