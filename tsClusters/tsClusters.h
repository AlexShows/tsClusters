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
TODO List: 
	-Migrate to a better logger (need timestamps and such)
	-Test data filling functions, bounds check, error handling
	-Add cluster centers data (and starting positions?)
	-Add cluster starting position assignments (minimum safe distance checks, etc.)
	-Add logging support
	-Add cluster assignment
	-Add centroid computation
	-Add logical processor counting
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
private:
	/* We need a shared_ptr here so we can traverse the data across
	multiple threads */
	std::shared_ptr<std::vector<T>> data;
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
};

/*
Default constructor
*/
template <typename T> tsClusters<T>::tsClusters()
{
	// By default we create this shared pointer, but we don't know the stride yet
	// until the data is filled
	data = std::make_shared<std::vector<T>>(*(new std::vector<T>));
	stride = 0;
	number_of_clusters = 0;

#ifdef _DEBUG
	log.open("debug.log", std::fstream::out);
	log << "Initialized tsClusters debug log.";
#endif
}

template <typename T> tsClusters<T>::~tsClusters()
{
#ifdef _DEBUG
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
TODO: Exception handling
*/
template <typename T> unsigned int tsClusters<T>::fill_data_array(T* input_data, unsigned int input_size, unsigned int input_stride)
{
	if(!input_data || !input_size || !input_stride)
		return 0;

	std::lock_guard<std::mutex> lock(tsLock);

	for(unsigned int i=0; i<input_size; i++)
		data->push_back(input_data[i]);
	
	number_of_clusters = input_stride;

#ifdef _DEBUG
	unsigned int count = 0;
	for(auto it=data->begin(); it!=data->end(); it++)
	{
		if(!(count%input_stride))
			log << std::endl;

		log << *it << "\t";
		count++;
	}
#endif

	return data->size();
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

#endif // _TS_CLUSTERS_H