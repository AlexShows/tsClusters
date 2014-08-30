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
	/* We need a shared_ptr here so we can traverse the data across
	multiple threads */
	std::shared_ptr<std::vector<T>> data;
	/* Like the data vector, the clusters have a stride to them, to 
	retain the flatness of the vector layout */
	std::shared_ptr<std::vector<T>> clusters;

	/* Stride is number of number of dimensions to the data */
	unsigned int stride;
	unsigned int number_of_clusters;
	/*******************
	A lock for any thread that may need to write to the data,
	though this should only ever need to be the thread owning the
	object created with this template, it's here just in case
	********************/
	std::mutex tsLock;
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
	data = std::make_shared<std::vector<T>>(*(new std::vector<T>));
	clusters = std::make_shared<std::vector<T>>(*(new std::vector<T>));
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

	data.reset(new std::vector<T>);
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
		// TODO: After every Nth input_data, add another min T value for the assigned cluster
		//		 Need to keep up with this padded cluster assignment everywhere else
		//		 that the data is traversed
		for (unsigned int i = 0; i < input_size; i++)
		{
			data->push_back(input_data[i]);

			if (i % input_stride == input_stride - 1)
				data->push_back(std::numeric_limits<T>::min());
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
	log << "Data points: ";

	unsigned int count = 0;
	for(auto it=data->begin(); it!=data->end(); it++)
	{
		if( !(count % stride) )
			log << std::endl;

		if (count == 0)
			log << *it << "\t";

		if ((count % (stride + 1)) != 0)
			log << *it << "\t";

		count++;
	}
#endif

	return (unsigned int)data->size();
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

	std::vector<T>::iterator it = data->begin();

	unsigned int counter = 0;
	while (it != data->end())
	{
		unsigned int i = counter % stride;

		// TODO: Still has a bug in the filter here
		if (i != stride - 1)
		{
			if (*it > ub[i])
				ub[i] = *it;

			if (*it < lb[i])
				lb[i] = *it;
		}
		
		it++;
		counter++;
	}

	// We should now have a lower and upper bound for every dimension in
	// the data, based on traversing all the data

	// For every cluster...
	for (unsigned int idx_c = 0; idx_c < number_of_clusters; idx_c++)
	{
		// ... and every dimension of every cluster...
		for (unsigned int idx_s = 0; idx_s < stride; idx_s++)
		{
			// ...put a random value into the clusters vector that is between the
			// lower bound and upper bound of this particular dimension
			// Using fmod from cmath as modulo is not defined for float
			T val = (T)(std::fmod(rand(), (ub[idx_s] - lb[idx_s]))) + lb[idx_s];
			clusters->push_back(val);
		}
	}

	// TODO: Test for minimum safe distance...

#ifdef _DEBUG
	log << std::endl << std::endl;
	log << "Cluster starting positions:" << std::endl;

	unsigned int it_c_cnt = 0;
	unsigned int cluster_index = 0;
	log << "Cluster " << cluster_index << ":" << std::endl;
	log << "     ";

	std::vector<T>::iterator it_c = clusters->begin();
	while (it_c != clusters->end())
	{
		log << *it_c;

		it_c_cnt++;

		if (!(it_c_cnt % stride) && it_c_cnt < clusters->size())
		{
			cluster_index++;

			log << std::endl;
			log << "Cluster " << cluster_index << ":" << std::endl;
			log << "     ";
		}
		else if (it_c_cnt < clusters->size())
			log << ", ";

		it_c++;
	}

	log << std::endl;
#endif

}

#endif // _TS_CLUSTERS_H