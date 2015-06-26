// kate: replace-tabs off; indent-width 4; indent-mode normal
// vim: ts=4:sw=4:noexpandtab
/*

Copyright (c) 2010--2012,
François Pomerleau and Stephane Magnenat, ASL, ETHZ, Switzerland
You can contact the authors at <f dot pomerleau at gmail dot com> and
<stephane at magnenat dot net>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ETH-ASL BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "DataPointsFiltersImpl.h"
#include "PointMatcherPrivate.h"
#include "MatchersImpl.h"
#include "Functions.h"

#include <algorithm>
#include <boost/format.hpp>

// Eigenvalues
#include "Eigen/QR"

using namespace std;
using namespace PointMatcherSupport;

// IdentityDataPointsFilter
template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::IdentityDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::IdentityDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
}

template struct DataPointsFiltersImpl<float>::IdentityDataPointsFilter;
template struct DataPointsFiltersImpl<double>::IdentityDataPointsFilter;


// RemoveNaNDataPointsFilter
template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::RemoveNaNDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::RemoveNaNDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
	const int nbPointsIn = cloud.features.cols();

	int j = 0;
	for (int i = 0; i < nbPointsIn; ++i)
	{
		const BOOST_AUTO(colArray, cloud.features.col(i).array());
		const BOOST_AUTO(hasNaN, !(colArray == colArray).all());
		if (!hasNaN)
		{
			cloud.setColFrom(j, cloud, i);
			j++;
		}
	}

	cloud.conservativeResize(j);
}

template struct DataPointsFiltersImpl<float>::RemoveNaNDataPointsFilter;
template struct DataPointsFiltersImpl<double>::RemoveNaNDataPointsFilter;


// MaxDistDataPointsFilter
// Constructor
template<typename T>
DataPointsFiltersImpl<T>::MaxDistDataPointsFilter::MaxDistDataPointsFilter(const Parameters& params):
	DataPointsFilter("MaxDistDataPointsFilter", MaxDistDataPointsFilter::availableParameters(), params),
	dim(Parametrizable::get<unsigned>("dim")),
	maxDist(Parametrizable::get<T>("maxDist"))
{
}

// Compute
template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::MaxDistDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::MaxDistDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
	if (dim >= cloud.features.rows() - 1)
	{
		throw InvalidParameter(
			(boost::format("MaxDistDataPointsFilter: Error, filtering on dimension number %1%, larger than authorized axis id %2%") % dim % (cloud.features.rows() - 2)).str());
	}

	const int nbPointsIn = cloud.features.cols();
	const int nbRows = cloud.features.rows();

	int j = 0;
	if(dim == -1) // Euclidean distance
	{
		for (int i = 0; i < nbPointsIn; i++)
		{
			const T absMaxDist = anyabs(maxDist);
			if (cloud.features.col(i).head(nbRows-1).norm() < absMaxDist)
			{
				cloud.setColFrom(j, cloud, i);
				j++;
			}
		}
	}
	else // Single-axis distance
	{
		for (int i = 0; i < nbPointsIn; i++)
		{
			if ((cloud.features(dim, i)) < maxDist)
			{
				cloud.setColFrom(j, cloud, i);
				j++;
			}
		}
	}

	cloud.conservativeResize(j);
}

template struct DataPointsFiltersImpl<float>::MaxDistDataPointsFilter;
template struct DataPointsFiltersImpl<double>::MaxDistDataPointsFilter;


// MinDistDataPointsFilter
// Constructor
template<typename T>
DataPointsFiltersImpl<T>::MinDistDataPointsFilter::MinDistDataPointsFilter(const Parameters& params):
	DataPointsFilter("MinDistDataPointsFilter", MinDistDataPointsFilter::availableParameters(), params),
	dim(Parametrizable::get<unsigned>("dim")),
	minDist(Parametrizable::get<T>("minDist"))
{
}

// Compute
template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::MinDistDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::MinDistDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
	if (dim >= cloud.features.rows() - 1)
		throw InvalidParameter((boost::format("MinDistDataPointsFilter: Error, filtering on dimension number %1%, larger than feature dimensionality %2%") % dim % (cloud.features.rows() - 2)).str());

	const int nbPointsIn = cloud.features.cols();
	const int nbRows = cloud.features.rows();

	int j = 0;
	if(dim == -1) // Euclidean distance
	{
		const T absMinDist = anyabs(minDist);
		for (int i = 0; i < nbPointsIn; i++)
		{
			if (cloud.features.col(i).head(nbRows-1).norm() > absMinDist)
			{
				cloud.setColFrom(j, cloud, i);
				j++;
			}
		}
	}
	else // Single axis distance
	{
		for (int i = 0; i < nbPointsIn; i++)
		{
			if ((cloud.features(dim, i)) > minDist)
			{
				cloud.setColFrom(j, cloud, i);
				j++;
			}
		}
	}

	cloud.conservativeResize(j);

}

template struct DataPointsFiltersImpl<float>::MinDistDataPointsFilter;
template struct DataPointsFiltersImpl<double>::MinDistDataPointsFilter;


// BoundingBoxDataPointsFilter
// Constructor
template<typename T>
DataPointsFiltersImpl<T>::BoundingBoxDataPointsFilter::BoundingBoxDataPointsFilter(const Parameters& params):
	DataPointsFilter("BoundingBoxDataPointsFilter", BoundingBoxDataPointsFilter::availableParameters(), params),
	xMin(Parametrizable::get<T>("xMin")),
	xMax(Parametrizable::get<T>("xMax")),
	yMin(Parametrizable::get<T>("yMin")),
	yMax(Parametrizable::get<T>("yMax")),
	zMin(Parametrizable::get<T>("zMin")),
	zMax(Parametrizable::get<T>("zMax")),
	removeInside(Parametrizable::get<bool>("removeInside"))
{
}

// Compute
template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::BoundingBoxDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::BoundingBoxDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
	const int nbPointsIn = cloud.features.cols();
	const int nbRows = cloud.features.rows();

	int j = 0;
	for (int i = 0; i < nbPointsIn; i++)
	{
		bool keepPt = false;
		const Vector point = cloud.features.col(i);

		// FIXME: improve performance by using Eigen array operations
		const bool x_in = (point(0) > xMin && point(0) < xMax);
		const bool y_in = (point(1) > yMin && point(1) < yMax);
		const bool z_in = (point(2) > zMin && point(2) < zMax) || nbRows == 3;
		const bool in_box = x_in && y_in && z_in;

		if(removeInside)
			keepPt = !in_box;
		else
			keepPt = in_box;

		if(keepPt)
		{
			cloud.setColFrom(j, cloud, i);
			j++;
		}
	}

	cloud.conservativeResize(j);
}

template struct DataPointsFiltersImpl<float>::BoundingBoxDataPointsFilter;
template struct DataPointsFiltersImpl<double>::BoundingBoxDataPointsFilter;


// MaxQuantileOnAxisDataPointsFilter
// Constructor
template<typename T>
DataPointsFiltersImpl<T>::MaxQuantileOnAxisDataPointsFilter::MaxQuantileOnAxisDataPointsFilter(const Parameters& params):
	DataPointsFilter("MaxQuantileOnAxisDataPointsFilter", MaxQuantileOnAxisDataPointsFilter::availableParameters(), params),
	dim(Parametrizable::get<unsigned>("dim")),
	ratio(Parametrizable::get<T>("ratio"))
{
}

// Compute
template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::MaxQuantileOnAxisDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::MaxQuantileOnAxisDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
	if (int(dim) >= cloud.features.rows())
		throw InvalidParameter((boost::format("MaxQuantileOnAxisDataPointsFilter: Error, filtering on dimension number %1%, larger than feature dimensionality %2%") % dim % cloud.features.rows()).str());

	const int nbPointsIn = cloud.features.cols();
	const int nbPointsOut = nbPointsIn * ratio;

	// build array
	vector<T> values;
	values.reserve(cloud.features.cols());
	for (int x = 0; x < cloud.features.cols(); ++x)
		values.push_back(cloud.features(dim, x));

	// get quartiles value
	nth_element(values.begin(), values.begin() + (values.size() * ratio), values.end());
	const T limit = values[nbPointsOut];

	// copy towards beginning the elements we keep
	int j = 0;
	for (int i = 0; i < nbPointsIn; i++)
	{
		if (cloud.features(dim, i) < limit)
		{
			assert(j <= i);
			cloud.setColFrom(j, cloud, i);
			j++;
		}
	}
	assert(j <= nbPointsOut);

	cloud.conservativeResize(j);

}

template struct DataPointsFiltersImpl<float>::MaxQuantileOnAxisDataPointsFilter;
template struct DataPointsFiltersImpl<double>::MaxQuantileOnAxisDataPointsFilter;


// MaxDensityDataPointsFilter
// Constructor
template<typename T>
DataPointsFiltersImpl<T>::MaxDensityDataPointsFilter::MaxDensityDataPointsFilter(const Parameters& params):
	DataPointsFilter("MaxDensityDataPointsFilter", MaxDensityDataPointsFilter::availableParameters(), params),
	maxDensity(Parametrizable::get<T>("maxDensity"))
{
}

// Compute
template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::MaxDensityDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::MaxDensityDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
	typedef typename DataPoints::View View;
	typedef typename DataPoints::ConstView ConstView;

	// Force densities to be computed
	if (!cloud.descriptorExists("densities"))
	{
		throw InvalidField("MaxDensityDataPointsFilter: Error, no densities found in descriptors.");
	}

	const int nbPointsIn = cloud.features.cols();
	View densities = cloud.getDescriptorViewByName("densities");
	const T lastDensity = densities.maxCoeff();
	const int nbSaturatedPts = (densities.cwise() == lastDensity).count();

	// fill cloud values
	int j = 0;
	for (int i = 0; i < nbPointsIn; i++)
	{
		const T density(densities(0,i));
		if (density > maxDensity)
		{
			const float r = (float)std::rand()/(float)RAND_MAX;
			float acceptRatio = maxDensity/density;

			// Handle saturation value of density
			if (density == lastDensity)
			{
				acceptRatio = acceptRatio * (1-nbSaturatedPts/nbPointsIn);
			}

			if (r < acceptRatio)
			{
				cloud.setColFrom(j, cloud, i);
				j++;
			}
		}
		else
		{
			cloud.setColFrom(j, cloud, i);
			j++;
		}
	}

	cloud.conservativeResize(j);
}

template struct DataPointsFiltersImpl<float>::MaxDensityDataPointsFilter;
template struct DataPointsFiltersImpl<double>::MaxDensityDataPointsFilter;


// SurfaceNormalDataPointsFilter
// Constructor
template<typename T>
DataPointsFiltersImpl<T>::SurfaceNormalDataPointsFilter::SurfaceNormalDataPointsFilter(const Parameters& params):
	DataPointsFilter("SurfaceNormalDataPointsFilter", SurfaceNormalDataPointsFilter::availableParameters(), params),
	knn(Parametrizable::get<int>("knn")),
	epsilon(Parametrizable::get<T>("epsilon")),
	keepNormals(Parametrizable::get<bool>("keepNormals")),
	keepDensities(Parametrizable::get<bool>("keepDensities")),
	keepEigenValues(Parametrizable::get<bool>("keepEigenValues")),
	keepEigenVectors(Parametrizable::get<bool>("keepEigenVectors")),
	keepMatchedIds(Parametrizable::get<bool>("keepMatchedIds"))
{
}

// Compute
template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::SurfaceNormalDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::SurfaceNormalDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
	typedef typename DataPoints::View View;
	typedef typename DataPoints::Label Label;
	typedef typename DataPoints::Labels Labels;
	typedef typename MatchersImpl<T>::KDTreeMatcher KDTreeMatcher;
	typedef typename PointMatcher<T>::Matches Matches;

	const int pointsCount(cloud.features.cols());
	const int featDim(cloud.features.rows());
	const int descDim(cloud.descriptors.rows());

	// Validate descriptors and labels
	int insertDim(0);
	for(unsigned int i = 0; i < cloud.descriptorLabels.size(); i++)
		insertDim += cloud.descriptorLabels[i].span;
	if (insertDim != descDim)
		throw InvalidField("SurfaceNormalDataPointsFilter: Error, descriptor labels do not match descriptor data");

	// Reserve memory for new descriptors
	const int dimNormals(featDim-1);
	const int dimDensities(1);
	const int dimEigValues(featDim-1);
	const int dimEigVectors((featDim-1)*(featDim-1));
	//const int dimMatchedIds(knn);

	boost::optional<View> normals;
	boost::optional<View> densities;
	boost::optional<View> eigenValues;
	boost::optional<View> eigenVectors;
	boost::optional<View> matchedValues;

	Labels cloudLabels;
	if (keepNormals)
		cloudLabels.push_back(Label("normals", dimNormals));
	if (keepDensities)
		cloudLabels.push_back(Label("densities", dimDensities));
	if (keepEigenValues)
		cloudLabels.push_back(Label("eigValues", dimEigValues));
	if (keepEigenVectors)
		cloudLabels.push_back(Label("eigVectors", dimEigVectors));
	cloud.allocateDescriptors(cloudLabels);

	if (keepNormals)
		normals = cloud.getDescriptorViewByName("normals");
	if (keepDensities)
		densities = cloud.getDescriptorViewByName("densities");
	if (keepEigenValues)
		eigenValues = cloud.getDescriptorViewByName("eigValues");
	if (keepEigenVectors)
		eigenVectors = cloud.getDescriptorViewByName("eigVectors");
	// TODO: implement keepMatchedIds
//	if (keepMatchedIds)
//	{
//		cloud.allocateDescriptor("normals", dimMatchedIds);
//		matchedValues = cloud.getDescriptorViewByName("normals");
//	}

	// Build kd-tree
	Parametrizable::Parameters param;
	boost::assign::insert(param) ( "knn", toParam(knn) );
	boost::assign::insert(param) ( "epsilon", toParam(epsilon) );
	KDTreeMatcher matcher(param);
	matcher.init(cloud);

	Matches matches(typename Matches::Dists(knn, pointsCount), typename Matches::Ids(knn, pointsCount));
	matches = matcher.findClosests(cloud);

	// Search for surrounding points and compute descriptors
	int degenerateCount(0);
	for (int i = 0; i < pointsCount; ++i)
	{
		// Mean of nearest neighbors (NN)
		Matrix d(featDim-1, knn);
		for(int j = 0; j < int(knn); j++)
		{
			const int refIndex(matches.ids(j,i));
			d.col(j) = cloud.features.block(0, refIndex, featDim-1, 1);
		}

		const Vector mean = d.rowwise().sum() / T(knn);
		const Matrix NN = d.colwise() - mean;

		const Matrix C(NN * NN.transpose());
		Vector eigenVa = Vector::Identity(featDim-1, 1);
		Matrix eigenVe = Matrix::Identity(featDim-1, featDim-1);
		// Ensure that the matrix is suited for eigenvalues calculation
		if(keepNormals || keepEigenValues || keepEigenVectors)
		{
			if(C.fullPivHouseholderQr().rank()+1 >= featDim-1)
			{
				const Eigen::EigenSolver<Matrix> solver(C);
				eigenVa = solver.eigenvalues().real();
				eigenVe = solver.eigenvectors().real();
			}
			else
			{
				//std::cout << "WARNING: Matrix C needed for eigen decomposition is degenerated. Expected cause: no noise in data" << std::endl;
				++degenerateCount;
			}
		}

		if(keepNormals)
			normals->col(i) = computeNormal(eigenVa, eigenVe);
		if(keepDensities)
			(*densities)(0, i) = computeDensity(NN);
		if(keepEigenValues)
			eigenValues->col(i) = eigenVa;
		if(keepEigenVectors)
			eigenVectors->col(i) = serializeEigVec(eigenVe);
	}
	if (degenerateCount)
	{
		LOG_WARNING_STREAM("WARNING: Matrix C needed for eigen decomposition was degenerated in " << degenerateCount << " points over " << pointsCount << " (" << float(degenerateCount)*100.f/float(pointsCount) << " %)");
	}

}

template<typename T>
typename PointMatcher<T>::Vector DataPointsFiltersImpl<T>::SurfaceNormalDataPointsFilter::computeNormal(const Vector eigenVa, const Matrix eigenVe)
{
	// Keep the smallest eigenvector as surface normal
	int smallestId(0);
	T smallestValue(numeric_limits<T>::max());
	for(int j = 0; j < eigenVe.cols(); j++)
	{
		if (eigenVa(j) < smallestValue)
		{
			smallestId = j;
			smallestValue = eigenVa(j);
		}
	}

	return eigenVe.col(smallestId);
}

template<typename T>
typename PointMatcher<T>::Vector DataPointsFiltersImpl<T>::SurfaceNormalDataPointsFilter::serializeEigVec(const Matrix eigenVe)
{
	// serialize row major
	const int eigenVeDim = eigenVe.cols();
	Vector output(eigenVeDim*eigenVeDim);
	for(int k=0; k < eigenVe.cols(); k++)
	{
		output.segment(k*eigenVeDim, eigenVeDim) =
			eigenVe.row(k).transpose();
	}

	return output;
}

template<typename T>
T DataPointsFiltersImpl<T>::SurfaceNormalDataPointsFilter::computeDensity(const Matrix NN)
{
	//volume in meter
	T volume = (4./3.)*M_PI*std::pow(NN.colwise().norm().maxCoeff(), 3);

	//volume in decimeter
	//T volume = (4./3.)*M_PI*std::pow(NN.colwise().norm().maxCoeff()*10.0, 3);
	//const T minVolume = 4.18e-9; // minimum of volume of one millimeter radius
	//const T minVolume = 0.42; // minimum of volume of one centimeter radius (in dm^3)

	//if(volume < minVolume)
	//	volume = minVolume;

	return T(NN.cols())/(volume);
}

template struct DataPointsFiltersImpl<float>::SurfaceNormalDataPointsFilter;
template struct DataPointsFiltersImpl<double>::SurfaceNormalDataPointsFilter;


// SamplingSurfaceNormalDataPointsFilter

// Constructor
template<typename T>
DataPointsFiltersImpl<T>::SamplingSurfaceNormalDataPointsFilter::SamplingSurfaceNormalDataPointsFilter(const Parameters& params):
	DataPointsFilter("SamplingSurfaceNormalDataPointsFilter", SamplingSurfaceNormalDataPointsFilter::availableParameters(), params),
	ratio(Parametrizable::get<T>("ratio")),
	knn(Parametrizable::get<int>("knn")),
	samplingMethod(Parametrizable::get<int>("samplingMethod")),
	maxBoxDim(Parametrizable::get<T>("maxBoxDim")),
	averageExistingDescriptors(Parametrizable::get<bool>("averageExistingDescriptors")),
	keepNormals(Parametrizable::get<bool>("keepNormals")),
	keepDensities(Parametrizable::get<bool>("keepDensities")),
	keepEigenValues(Parametrizable::get<bool>("keepEigenValues")),
	keepEigenVectors(Parametrizable::get<bool>("keepEigenVectors"))
{
}

// Compute
template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::SamplingSurfaceNormalDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::SamplingSurfaceNormalDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
	typedef Matrix Features;
	typedef typename DataPoints::View View;
	typedef typename DataPoints::Label Label;
	typedef typename DataPoints::Labels Labels;

	const int pointsCount(cloud.features.cols());
	const int featDim(cloud.features.rows());
	const int descDim(cloud.descriptors.rows());

	int insertDim(0);
	if (averageExistingDescriptors)
	{
		// TODO: this should be in the form of an assert
		// Validate descriptors and labels
		for(unsigned int i = 0; i < cloud.descriptorLabels.size(); i++)
			insertDim += cloud.descriptorLabels[i].span;
		if (insertDim != descDim)
			throw InvalidField("SamplingSurfaceNormalDataPointsFilter: Error, descriptor labels do not match descriptor data");
	}

	// Compute space requirement for new descriptors
	const int dimNormals(featDim-1);
	const int dimDensities(1);
	const int dimEigValues(featDim-1);
	const int dimEigVectors((featDim-1)*(featDim-1));

	// Allocate space for new descriptors
	Labels cloudLabels;
	if (keepNormals)
		cloudLabels.push_back(Label("normals", dimNormals));
	if (keepDensities)
		cloudLabels.push_back(Label("densities", dimDensities));
	if (keepEigenValues)
		cloudLabels.push_back(Label("eigValues", dimEigValues));
	if (keepEigenVectors)
		cloudLabels.push_back(Label("eigVectors", dimEigVectors));
	cloud.allocateDescriptors(cloudLabels);

	// we keep build data on stack for reentrant behaviour
	View cloudExistingDescriptors(cloud.descriptors.block(0,0,cloud.descriptors.rows(),cloud.descriptors.cols()));
	BuildData buildData(cloud.features, cloud.descriptors);

	// get views
	if (keepNormals)
		buildData.normals = cloud.getDescriptorViewByName("normals");
	if (keepDensities)
		buildData.densities = cloud.getDescriptorViewByName("densities");
	if (keepEigenValues)
		buildData.eigenValues = cloud.getDescriptorViewByName("eigValues");
	if (keepEigenVectors)
		buildData.eigenVectors = cloud.getDescriptorViewByName("eigVectors");
	// build the new point cloud
	buildNew(
		buildData,
		0,
		pointsCount,
		cloud.features.rowwise().minCoeff(),
		cloud.features.rowwise().maxCoeff()
	);

	// Bring the data we keep to the front of the arrays then
	// wipe the leftover unused space.
	std::sort(buildData.indicesToKeep.begin(), buildData.indicesToKeep.end());
	int ptsOut = buildData.indicesToKeep.size();
	for (int i = 0; i < ptsOut; i++){
		int k = buildData.indicesToKeep[i];
		assert(i <= k);
		cloud.features.col(i) = cloud.features.col(k);
		if (cloud.descriptors.rows() != 0)
			cloud.descriptors.col(i) = cloud.descriptors.col(k);
		if(keepNormals)
			buildData.normals->col(i) = buildData.normals->col(k);
		if(keepDensities)
			(*buildData.densities)(0,i) = (*buildData.densities)(0,k);
		if(keepEigenValues)
			buildData.eigenValues->col(i) = buildData.eigenValues->col(k);
		if(keepEigenVectors)
			buildData.eigenVectors->col(i) = buildData.eigenVectors->col(k);
	}
	cloud.features.conservativeResize(Eigen::NoChange, ptsOut);
	cloud.descriptors.conservativeResize(Eigen::NoChange, ptsOut);

	// warning if some points were dropped
	if(buildData.unfitPointsCount != 0)
		LOG_INFO_STREAM("  SamplingSurfaceNormalDataPointsFilter - Could not compute normal for " << buildData.unfitPointsCount << " pts.");
}

template<typename T>
size_t argMax(const typename PointMatcher<T>::Vector& v)
{
	T maxVal(0);
	size_t maxIdx(0);
	for (int i = 0; i < v.size(); ++i)
	{
		if (v[i] > maxVal)
		{
			maxVal = v[i];
			maxIdx = i;
		}
	}
	return maxIdx;
}

template<typename T>
void DataPointsFiltersImpl<T>::SamplingSurfaceNormalDataPointsFilter::buildNew(BuildData& data, const int first, const int last, const Vector minValues, const Vector maxValues) const
{
	const int count(last - first);
	if (count <= int(knn))
	{
		// compute for this range
		fuseRange(data, first, last);
		// TODO: make another filter that creates constant-density clouds,
		// typically by stopping recursion after the median of the bounding cuboid
		// is below a threshold, or that the number of points falls under a threshold
		return;
	}

	// find the largest dimension of the box
	const int cutDim = argMax<T>(maxValues - minValues);

	// compute number of elements
	const int rightCount(count/2);
	const int leftCount(count - rightCount);
	assert(last - rightCount == first + leftCount);

	// sort, hack std::nth_element
	std::nth_element(
		data.indices.begin() + first,
		data.indices.begin() + first + leftCount,
		data.indices.begin() + last,
		CompareDim(cutDim, data)
	);

	// get value
	const int cutIndex(data.indices[first+leftCount]);
	const T cutVal(data.features(cutDim, cutIndex));

	// update bounds for left
	Vector leftMaxValues(maxValues);
	leftMaxValues[cutDim] = cutVal;
	// update bounds for right
	Vector rightMinValues(minValues);
	rightMinValues[cutDim] = cutVal;

	// recurse
	buildNew(data, first, first + leftCount, minValues, leftMaxValues);
	buildNew(data, first + leftCount, last, rightMinValues, maxValues);
}

template<typename T>
void DataPointsFiltersImpl<T>::SamplingSurfaceNormalDataPointsFilter::fuseRange(BuildData& data, const int first, const int last) const
{
	const int colCount(last-first);
	const int featDim(data.features.rows());

	// build nearest neighbors list
	Matrix d(featDim-1, colCount);
	for (int i = 0; i < colCount; ++i)
		d.col(i) = data.features.block(0,data.indices[first+i],featDim-1, 1);
	const Vector box = d.rowwise().maxCoeff() - d.rowwise().minCoeff();
	const T boxDim(box.maxCoeff());
	// drop box if it is too large
	if (boxDim > maxBoxDim)
	{
		data.unfitPointsCount += colCount;
		return;
	}
	const Vector mean = d.rowwise().sum() / T(colCount);
	const Matrix NN = (d.colwise() - mean);

	// compute covariance
	const Matrix C(NN * NN.transpose());
	Vector eigenVa = Vector::Identity(featDim-1, 1);
	Matrix eigenVe = Matrix::Identity(featDim-1, featDim-1);
	// Ensure that the matrix is suited for eigenvalues calculation
	if(keepNormals || keepEigenValues || keepEigenVectors)
	{
		if(C.fullPivHouseholderQr().rank()+1 >= featDim-1)
		{
			const Eigen::EigenSolver<Matrix> solver(C);
			eigenVa = solver.eigenvalues().real();
			eigenVe = solver.eigenvectors().real();
		}
		else
		{
			data.unfitPointsCount += colCount;
			return;
		}
	}

	Vector normal;
	if(keepNormals)
		normal = SurfaceNormalDataPointsFilter::computeNormal(eigenVa, eigenVe);

	T densitie = 0;
	if(keepDensities)
		densitie = SurfaceNormalDataPointsFilter::computeDensity(NN);

	//if(keepEigenValues) nothing to do

	Vector serialEigVector;
	if(keepEigenVectors)
		serialEigVector = SurfaceNormalDataPointsFilter::serializeEigVec(eigenVe);

	// some safety check
	if(data.descriptors.rows() != 0)
		assert(data.descriptors.cols() != 0);

	// Filter points randomly
	if(samplingMethod == 0)
	{
		for(int i=0; i<colCount; i++)
		{
			const float r = (float)std::rand()/(float)RAND_MAX;
			if(r < ratio)
			{
				// Keep points with their descriptors
				int k = data.indices[first+i];
				// Mark the indices which will be part of the final data
				data.indicesToKeep.push_back(k);

				// Build new descriptors
				if(keepNormals)
					data.normals->col(k) = normal;
				if(keepDensities)
					(*data.densities)(0,k) = densitie;
				if(keepEigenValues)
					data.eigenValues->col(k) = eigenVa;
				if(keepEigenVectors)
					data.eigenVectors->col(k) = serialEigVector;

			}
		}
	}
	else
	{

		int k = data.indices[first];
		// Mark the indices which will be part of the final data
		data.indicesToKeep.push_back(k);
		data.features.col(k).topRows(featDim-1) = mean;
		data.features(featDim-1, k) = 1;

		if(data.descriptors.rows() != 0)
		{
			// average the existing descriptors
			if (averageExistingDescriptors)
			{
				Vector mergedDesc(Vector::Zero(data.descriptors.rows()));
				for (int i = 0; i < colCount; ++i)
					mergedDesc += data.descriptors.col(data.indices[first+i]);
				mergedDesc /= T(colCount);
				data.descriptors.col(k) = mergedDesc;
			}
			// else just keep the first one
		}

		// Build new descriptors
		if(keepNormals)
			data.normals->col(k) = normal;
		if(keepDensities)
			(*data.densities)(0,k) = densitie;
		if(keepEigenValues)
			data.eigenValues->col(k) = eigenVa;
		if(keepEigenVectors)
			data.eigenVectors->col(k) = serialEigVector;

	}

}

template struct DataPointsFiltersImpl<float>::SamplingSurfaceNormalDataPointsFilter;
template struct DataPointsFiltersImpl<double>::SamplingSurfaceNormalDataPointsFilter;

// OrientNormalsDataPointsFilter
// Constructor
template<typename T>
DataPointsFiltersImpl<T>::OrientNormalsDataPointsFilter::OrientNormalsDataPointsFilter(const Parameters& params):
	DataPointsFilter("OrientNormalsDataPointsFilter", OrientNormalsDataPointsFilter::availableParameters(), params),
	towardCenter(Parametrizable::get<bool>("towardCenter"))
{
}

// OrientNormalsDataPointsFilter
// Compute
template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::OrientNormalsDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::OrientNormalsDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
	if (!cloud.descriptorExists("normals"))
		throw InvalidField("OrientNormalsDataPointsFilter: Error, cannot find normals in descriptors.");
	if (!cloud.descriptorExists("observationDirections"))
		throw InvalidField("OrientNormalsDataPointsFilter: Error, cannot find observation directions in descriptors.");

	BOOST_AUTO(normals, cloud.getDescriptorViewByName("normals"));
	const BOOST_AUTO(observationDirections, cloud.getDescriptorViewByName("observationDirections"));
	assert(normals.rows() == observationDirections.rows());
	for (int i = 0; i < cloud.features.cols(); i++)
	{
		// Check normal orientation
		const Vector vecP = observationDirections.col(i);
		const Vector vecN = normals.col(i);
		const double scalar = vecP.dot(vecN);

		// Swap normal
		if(towardCenter)
		{
			if (scalar < 0)
				normals.col(i) = -vecN;
		}
		else
		{
			if (scalar > 0)
				normals.col(i) = -vecN;
		}
	}

}

template struct DataPointsFiltersImpl<float>::OrientNormalsDataPointsFilter;
template struct DataPointsFiltersImpl<double>::OrientNormalsDataPointsFilter;


// RandomSamplingDataPointsFilter
// Constructor
template<typename T>
DataPointsFiltersImpl<T>::RandomSamplingDataPointsFilter::RandomSamplingDataPointsFilter(const Parameters& params):
	DataPointsFilter("RandomSamplingDataPointsFilter", RandomSamplingDataPointsFilter::availableParameters(), params),
	prob(Parametrizable::get<double>("prob"))
{
}

// Constructor
template<typename T>
DataPointsFiltersImpl<T>::RandomSamplingDataPointsFilter::RandomSamplingDataPointsFilter(const std::string& className, const ParametersDoc paramsDoc, const Parameters& params):
	DataPointsFilter(className, paramsDoc, params),
	prob(Parametrizable::get<double>("prob"))
{
}

// Compute
template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::RandomSamplingDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::RandomSamplingDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
	const int nbPointsIn = cloud.features.cols();

	int j = 0;
	for (int i = 0; i < nbPointsIn; i++)
	{
		const float r = (float)std::rand()/(float)RAND_MAX;
		if (r < prob)
		{
			cloud.setColFrom(j, cloud, i);
			j++;
		}
	}

	cloud.conservativeResize(j);
}

template struct DataPointsFiltersImpl<float>::RandomSamplingDataPointsFilter;
template struct DataPointsFiltersImpl<double>::RandomSamplingDataPointsFilter;


// MaxPointCountDataPointsFilter
// Constructor
template<typename T>
DataPointsFiltersImpl<T>::MaxPointCountDataPointsFilter::MaxPointCountDataPointsFilter(const Parameters& params):
	RandomSamplingDataPointsFilter("MaxPointCountDataPointsFilter", MaxPointCountDataPointsFilter::availableParameters(), params),
	maxCount(Parametrizable::get<unsigned>("maxCount"))
{
}

// Compute
template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::MaxPointCountDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::MaxPointCountDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
	if (unsigned(cloud.features.cols()) <= maxCount)
		return;
	RandomSamplingDataPointsFilter::inPlaceFilter(cloud);
}

template struct DataPointsFiltersImpl<float>::MaxPointCountDataPointsFilter;
template struct DataPointsFiltersImpl<double>::MaxPointCountDataPointsFilter;


// FixStepSamplingDataPointsFilter
// Constructor
template<typename T>
DataPointsFiltersImpl<T>::FixStepSamplingDataPointsFilter::FixStepSamplingDataPointsFilter(const Parameters& params):
	DataPointsFilter("FixStepSamplingDataPointsFilter", FixStepSamplingDataPointsFilter::availableParameters(), params),
	startStep(Parametrizable::get<unsigned>("startStep")),
	endStep(Parametrizable::get<unsigned>("endStep")),
	stepMult(Parametrizable::get<double>("stepMult")),
	step(startStep)
{
	LOG_INFO_STREAM("Using FixStepSamplingDataPointsFilter with startStep=" << startStep << ", endStep=" << endStep << ", stepMult=" << stepMult);
}


template<typename T>
void DataPointsFiltersImpl<T>::FixStepSamplingDataPointsFilter::init()
{
	step = startStep;
}

// Compute
template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::FixStepSamplingDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::FixStepSamplingDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
	const int iStep(step);
	const int nbPointsIn = cloud.features.cols();
	const int phase(rand() % iStep);

	int j = 0;
	for (int i = phase; i < nbPointsIn; i += iStep)
	{
		cloud.setColFrom(j, cloud, i);
		j++;
	}

	cloud.conservativeResize(j);

	const double deltaStep(startStep * stepMult - startStep);
	step *= stepMult;
	if (deltaStep < 0 && step < endStep)
		step = endStep;
	if (deltaStep > 0 && step > endStep)
		step = endStep;

}

template struct DataPointsFiltersImpl<float>::FixStepSamplingDataPointsFilter;
template struct DataPointsFiltersImpl<double>::FixStepSamplingDataPointsFilter;


// ShadowDataPointsFilter
// Constructor
template<typename T>
DataPointsFiltersImpl<T>::ShadowDataPointsFilter::ShadowDataPointsFilter(const Parameters& params):
	DataPointsFilter("ShadowDataPointsFilter", ShadowDataPointsFilter::availableParameters(), params),
	eps(sin(Parametrizable::get<T>("eps")))
{
	//waring: maxAngle is change to sin(maxAngle)!
}


// Compute
template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::ShadowDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::ShadowDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
	// Check if normals are present
	if (!cloud.descriptorExists("normals"))
	{
		throw InvalidField("ShadowDataPointsFilter, Error: cannot find normals in descriptors");
	}

	const int dim = cloud.features.rows();

	const BOOST_AUTO(normals, cloud.getDescriptorViewByName("normals"));
	int j = 0;

	for(int i=0; i < cloud.features.cols(); i++)
	{
		const Vector normal = normals.col(i).normalized();
		const Vector point = cloud.features.block(0, i, dim-1, 1).normalized();

		const T value = anyabs(normal.dot(point));

		if(value > eps) // test to keep the points
		{
			cloud.features.col(j) = cloud.features.col(i);
			cloud.descriptors.col(j) = cloud.descriptors.col(i);
			j++;
		}
	}

	cloud.conservativeResize(j);
}

template struct DataPointsFiltersImpl<float>::ShadowDataPointsFilter;
template struct DataPointsFiltersImpl<double>::ShadowDataPointsFilter;


// SimpleSensorNoiseDataPointsFilter
// Constructor
template<typename T>
DataPointsFiltersImpl<T>::SimpleSensorNoiseDataPointsFilter::SimpleSensorNoiseDataPointsFilter(const Parameters& params):
	DataPointsFilter("SimpleSensorNoiseDataPointsFilter", SimpleSensorNoiseDataPointsFilter::availableParameters(), params),
	sensorType(Parametrizable::get<unsigned>("sensorType")),
	gain(Parametrizable::get<T>("gain"))
{
  std::vector<string> sensorNames = boost::assign::list_of ("Sick LMS-1xx")("Hokuyo URG-04LX")("Hokuyo UTM-30LX")("Kinect / Xtion")("Sick Tim3xx");
	if (sensorType >= sensorNames.size())
	{
		throw InvalidParameter(
			(boost::format("SimpleSensorNoiseDataPointsFilter: Error, sensorType id %1% does not exist.") % sensorType).str());
	}

	LOG_INFO_STREAM("SimpleSensorNoiseDataPointsFilter - using sensor noise model: " << sensorNames[sensorType]);
}


// SimpleSensorNoiseDataPointsFilter
// Compute
template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::SimpleSensorNoiseDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::SimpleSensorNoiseDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
	cloud.allocateDescriptor("simpleSensorNoise", 1);
	BOOST_AUTO(noise, cloud.getDescriptorViewByName("simpleSensorNoise"));

	switch(sensorType)
	{
	case 0: // Sick LMS-1xx
	{
		noise = computeLaserNoise(0.012, 0.0068, 0.0008, cloud.features);
		break;
	}
	case 1: // Hokuyo URG-04LX
	{
		noise = computeLaserNoise(0.028, 0.0013, 0.0001, cloud.features);
		break;
	}
	case 2: // Hokuyo UTM-30LX
	{
		noise = computeLaserNoise(0.018, 0.0006, 0.0015, cloud.features);
		break;
	}
	case 3: // Kinect / Xtion
	{
		const int dim = cloud.features.rows();
		const Matrix squaredValues(cloud.features.topRows(dim-1).colwise().norm().array().square());
		noise = squaredValues*(0.5*0.00285);
		break;
	}
  case 4: // Sick Tim3xx
  {
    noise = computeLaserNoise(0.004, 0.0053, -0.0092, cloud.features);
    break;
  }
	default:
		throw InvalidParameter(
			(boost::format("SimpleSensorNoiseDataPointsFilter: Error, cannot compute noise for sensorType id %1% .") % sensorType).str());
	}

}

template<typename T>
typename PointMatcher<T>::Matrix DataPointsFiltersImpl<T>::SimpleSensorNoiseDataPointsFilter::computeLaserNoise(const T minRadius, const T beamAngle, const T beamConst, const Matrix features)
{
	typedef typename Eigen::Array<T, 2, Eigen::Dynamic> Array2rows;

	const int nbPoints = features.cols();
	const int dim = features.rows();

	Array2rows evalNoise = Array2rows::Constant(2, nbPoints, minRadius);
	evalNoise.row(0) =  beamAngle * features.topRows(dim-1).colwise().norm();
	evalNoise.row(0) += beamConst;

	return evalNoise.colwise().maxCoeff();

}


template struct DataPointsFiltersImpl<float>::SimpleSensorNoiseDataPointsFilter;
template struct DataPointsFiltersImpl<double>::SimpleSensorNoiseDataPointsFilter;


// ObservationDirectionDataPointsFilter
// Constructor
template<typename T>
DataPointsFiltersImpl<T>::ObservationDirectionDataPointsFilter::ObservationDirectionDataPointsFilter(const Parameters& params):
	DataPointsFilter("ObservationDirectionDataPointsFilter", ObservationDirectionDataPointsFilter::availableParameters(), params),
	centerX(Parametrizable::get<T>("x")),
	centerY(Parametrizable::get<T>("y")),
	centerZ(Parametrizable::get<T>("z"))
{
}

// Compute
template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::ObservationDirectionDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::ObservationDirectionDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
	const int dim(cloud.features.rows() - 1);
	if (dim != 2 && dim != 3)
	{
		throw InvalidField(
			(boost::format("ObservationDirectionDataPointsFilter: Error, works only in 2 or 3 dimensions, cloud has %1% dimensions.") % dim).str()
		);
	}

	Vector center(dim);
	center[0] = centerX;
	center[1] = centerY;
	if (dim == 3)
		center[2] = centerZ;

	cloud.allocateDescriptor("observationDirections", dim);
	BOOST_AUTO(observationDirections, cloud.getDescriptorViewByName("observationDirections"));

	for (int i = 0; i < cloud.features.cols(); i++)
	{
		// Check normal orientation
		const Vector p(cloud.features.block(0, i, dim, 1));
		observationDirections.col(i) = center - p;
	}

}

template struct DataPointsFiltersImpl<float>::ObservationDirectionDataPointsFilter;
template struct DataPointsFiltersImpl<double>::ObservationDirectionDataPointsFilter;

// VoxelGridDataPointsFilter
template <typename T>
DataPointsFiltersImpl<T>::VoxelGridDataPointsFilter::VoxelGridDataPointsFilter() :
	vSizeX(1),
	vSizeY(1),
	vSizeZ(1),
	useCentroid(true),
	averageExistingDescriptors(true) {}

template <typename T>
DataPointsFiltersImpl<T>::VoxelGridDataPointsFilter::VoxelGridDataPointsFilter(const Parameters& params) :
DataPointsFilter("VoxelGridDataPointsFilter", VoxelGridDataPointsFilter::availableParameters(), params),
		vSizeX(Parametrizable::get<T>("vSizeX")),
		vSizeY(Parametrizable::get<T>("vSizeY")),
		vSizeZ(Parametrizable::get<T>("vSizeZ")),
		useCentroid(Parametrizable::get<bool>("useCentroid")),
		averageExistingDescriptors(Parametrizable::get<bool>("averageExistingDescriptors"))
{

}

template <typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::VoxelGridDataPointsFilter::filter(const DataPoints& input)
{
    DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

template <typename T>
void DataPointsFiltersImpl<T>::VoxelGridDataPointsFilter::inPlaceFilter(DataPoints& cloud)
{
    const unsigned int numPoints(cloud.features.cols());
	const int featDim(cloud.features.rows());
	const int descDim(cloud.descriptors.rows());

	assert (featDim == 3 || featDim == 4);

	int insertDim(0);
	if (averageExistingDescriptors)
	{
		// TODO: this should be in the form of an assert
		// Validate descriptors and labels
		for(unsigned int i = 0; i < cloud.descriptorLabels.size(); i++)
			insertDim += cloud.descriptorLabels[i].span;
		if (insertDim != descDim)
			throw InvalidField("VoxelGridDataPointsFilter: Error, descriptor labels do not match descriptor data");
	}

	// TODO: Check that the voxel size is not too small, given the size of the data

	// Calculate number of divisions along each axis
	Vector minValues = cloud.features.rowwise().minCoeff();
	Vector maxValues = cloud.features.rowwise().maxCoeff();

    T minBoundX = minValues.x() / vSizeX;
    T maxBoundX = maxValues.x() / vSizeX;
    T minBoundY = minValues.y() / vSizeY;
    T maxBoundY = maxValues.y() / vSizeY;
    T minBoundZ = 0;
    T maxBoundZ = 0;

    if (featDim == 4)
    {
        minBoundZ = minValues.z() / vSizeZ;
        maxBoundZ = maxValues.z() / vSizeZ;
    }

    // number of divisions is total size / voxel size voxels of equal length + 1
    // with remaining space
    unsigned int numDivX = 1 + maxBoundX - minBoundX;
    unsigned int numDivY = 1 + maxBoundY - minBoundY;;
    unsigned int numDivZ = 0;

    // If a 3D point cloud
    if (featDim == 4 )
        numDivZ = 1 + maxBoundZ - minBoundZ;

    unsigned int numVox = numDivX * numDivY;
    if ( featDim == 4)
        numVox *= numDivZ;

    // Assume point cloud is randomly ordered
    // compute a linear index of the following type
    // i, j, k are the component indices
    // nx, ny number of divisions in x and y components
    // idx = i + j * nx + k * nx * ny
    std::vector<unsigned int> indices(numPoints);

    // vector to hold the first point in a voxel
    // this point will be ovewritten in the input cloud with
    // the output value

    std::vector<Voxel>* voxels;

    // try allocating vector. If too big return error
    try {
    	voxels = new std::vector<Voxel>(numVox);
    } catch (std::bad_alloc&) {
    	throw InvalidParameter((boost::format("VoxelGridDataPointsFilter: Memory allocation error with %1% voxels.  Try increasing the voxel dimensions.") % numVox).str());
    }


    for (unsigned int p = 0; p < numPoints; p++ )
    {
        unsigned int i = floor(cloud.features(0,p)/vSizeX - minBoundX);
        unsigned int j = floor(cloud.features(1,p)/vSizeY- minBoundY);
        unsigned int k = 0;
        unsigned int idx;
        if ( featDim == 4 )
        {
            k = floor(cloud.features(2,p)/vSizeZ - minBoundZ);
            idx = i + j * numDivX + k * numDivX * numDivY;
        }
        else
        {
            idx = i + j * numDivX;
        }

        unsigned int pointsInVox = (*voxels)[idx].numPoints + 1;

        if (pointsInVox == 1)
        {
            (*voxels)[idx].firstPoint = p;
        }

        (*voxels)[idx].numPoints = pointsInVox;

        indices[p] = idx;

    }

    // store which points contain voxel position
    std::vector<unsigned int> pointsToKeep;

    // Store voxel centroid in output
    if (useCentroid)
    {
        // Iterate through the indices and sum values to compute centroid
        for (unsigned int p = 0; p < numPoints ; p++)
        {
            unsigned int idx = indices[p];
            unsigned int firstPoint = (*voxels)[idx].firstPoint;

            // If this is the first point in the voxel, leave as is
            // if not sum up this point for centroid calculation
            if (firstPoint != p)
            {
            	// Sum up features and descriptors (if we are also averaging descriptors)

            	for (int f = 0; f < (featDim - 1); f++ )
            	{
            		cloud.features(f,firstPoint) += cloud.features(f,p);
            	}

            	if (averageExistingDescriptors) {
            		for (int d = 0; d < descDim; d++)
            		{
            			cloud.descriptors(d,firstPoint) += cloud.descriptors(d,p);
            		}
            	}
            }
        }

        // Now iterating through the voxels
        // Normalize sums to get centroid (average)
        // Some voxels may be empty and are discarded
        for(unsigned int idx = 0; idx < numVox; idx++)
        {
            unsigned int numPoints = (*voxels)[idx].numPoints;
            unsigned int firstPoint = (*voxels)[idx].firstPoint;
            if(numPoints > 0)
            {
                for ( int f = 0; f < (featDim - 1); f++ )
                    cloud.features(f,firstPoint) /= numPoints;

                if (averageExistingDescriptors) {
                	for ( int d = 0; d < descDim; d++ )
                		cloud.descriptors(d,firstPoint) /= numPoints;
                }

                pointsToKeep.push_back(firstPoint);
            }
        }
    }
    else
    {
    	// Although we don't sum over the features, we may still need to sum the descriptors
    	if (averageExistingDescriptors)
    	{
    		// Iterate through the indices and sum values to compute centroid
    		for (unsigned int p = 0; p < numPoints ; p++)
    		{
    			unsigned int idx = indices[p];
    			unsigned int firstPoint = (*voxels)[idx].firstPoint;

    			// If this is the first point in the voxel, leave as is
    			// if not sum up this point for centroid calculation
    			if (firstPoint != p)
    			{
    				for (int d = 0; d < descDim; d++ )
    				{
    					cloud.descriptors(d,firstPoint) += cloud.descriptors(d,p);
    				}
    			}
    		}
    	}

        for (unsigned int idx = 0; idx < numVox; idx++)
        {
            unsigned int numPoints = (*voxels)[idx].numPoints;
            unsigned int firstPoint = (*voxels)[idx].firstPoint;

            if (numPoints > 0)
            {
                // get back voxel indices in grid format
                // If we are in the last division, the voxel is smaller in size
                // We adjust the center as from the end of the last voxel to the bounding area
                unsigned int i = 0;
                unsigned int j = 0;
                unsigned int k = 0;
                if (featDim == 4)
                {
                    k = idx / (numDivX * numDivY);
                    if (k == numDivZ)
                        cloud.features(3,firstPoint) = maxValues.z() - (k-1) * vSizeZ/2;
                    else
                        cloud.features(3,firstPoint) = k * vSizeZ + vSizeZ/2;
                }

                j = (idx - k * numDivX * numDivY) / numDivX;
                if (j == numDivY)
                    cloud.features(2,firstPoint) = maxValues.y() - (j-1) * vSizeY/2;
                else
                    cloud.features(2,firstPoint) = j * vSizeY + vSizeY / 2;

                i = idx - k * numDivX * numDivY - j * numDivX;
                if (i == numDivX)
                    cloud.features(1,firstPoint) = maxValues.x() - (i-1) * vSizeX/2;
                else
                    cloud.features(1,firstPoint) = i * vSizeX + vSizeX / 2;

                // Descriptors : normalize if we are averaging or keep as is
                if (averageExistingDescriptors) {
                	for ( int d = 0; d < descDim; d++ )
                		cloud.descriptors(d,firstPoint) /= numPoints;
                }

                pointsToKeep.push_back(firstPoint);
            }
        }

    }

    // deallocate memory for voxels information
    delete voxels;

    // Move the points to be kept to the start
    // Bring the data we keep to the front of the arrays then
	// wipe the leftover unused space.
	std::sort(pointsToKeep.begin(), pointsToKeep.end());
	int numPtsOut = pointsToKeep.size();
	for (int i = 0; i < numPtsOut; i++){
		int k = pointsToKeep[i];
		assert(i <= k);
		cloud.features.col(i) = cloud.features.col(k);
		if (cloud.descriptors.rows() != 0)
			cloud.descriptors.col(i) = cloud.descriptors.col(k);
	}
	cloud.features.conservativeResize(Eigen::NoChange, numPtsOut);

	if (cloud.descriptors.rows() != 0)
		cloud.descriptors.conservativeResize(Eigen::NoChange, numPtsOut);
}

template struct DataPointsFiltersImpl<float>::VoxelGridDataPointsFilter;
template struct DataPointsFiltersImpl<double>::VoxelGridDataPointsFilter;


// CutAtDescriptorThresholdDataPointsFilter
// Constructor
template<typename T>
DataPointsFiltersImpl<T>::CutAtDescriptorThresholdDataPointsFilter::CutAtDescriptorThresholdDataPointsFilter(const Parameters& params):
	DataPointsFilter("CutAtDescriptorThresholdDataPointsFilter", CutAtDescriptorThresholdDataPointsFilter::availableParameters(), params),
	descName(Parametrizable::get<std::string>("descName")),
	useLargerThan(Parametrizable::get<bool>("useLargerThan")),
	threshold(Parametrizable::get<T>("threshold"))
{
}

// Compute
template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::CutAtDescriptorThresholdDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::CutAtDescriptorThresholdDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
	// Check field exists
	if (!cloud.descriptorExists(descName))
	{
		throw InvalidField("CutAtDescriptorThresholdDataPointsFilter: Error, field not found in descriptors.");
	}

	const int nbPointsIn = cloud.features.cols();
	typename DataPoints::View values = cloud.getDescriptorViewByName(descName);

	// fill cloud values
	int j = 0;
	if (useLargerThan)
	{
		for (int i = 0; i < nbPointsIn; i++)
		{
			const T value(values(0,i));
			if (value <= threshold)
			{
				cloud.setColFrom(j, cloud, i);
				j++;
			}
		}
	}
	else
	{
		for (int i = 0; i < nbPointsIn; i++)
		{
			const T value(values(0,i));
			if (value >= threshold)
			{
				cloud.setColFrom(j, cloud, i);
				j++;
			}
		}
	}
	cloud.conservativeResize(j);
}

template struct DataPointsFiltersImpl<float>::CutAtDescriptorThresholdDataPointsFilter;
template struct DataPointsFiltersImpl<double>::CutAtDescriptorThresholdDataPointsFilter;


// InterestPointSamplingDataPointsFilter

// Here we compare first voxel > second voxel
// Empty voxels will always be considered less than a non-empty voxel.
// If both contain points, we compare the dispersion ratio.
template <typename T>
bool greater_than_voxel (
	const typename DataPointsFiltersImpl<T>::InterestPointSamplingDataPointsFilter::Voxel first,
	const typename DataPointsFiltersImpl<T>::InterestPointSamplingDataPointsFilter::Voxel second)
{
	// First empty == second empty
	if(first.pointIndices.empty() && second.pointIndices.empty()) return false;

	// Non-empty voxels > Empty voxels
	if(not first.pointIndices.empty() &&     second.pointIndices.empty()) return true;
	if(    first.pointIndices.empty() && not second.pointIndices.empty()) return false;

	// If they're both not empty, sort by the dispersion ratio
	return first.dispersionRatio > second.dispersionRatio;
}

template <typename T>
DataPointsFiltersImpl<T>::InterestPointSamplingDataPointsFilter::InterestPointSamplingDataPointsFilter() :
	vSizeX(1),
	vSizeY(1),
	vSizeZ(1),
	perc(1),
	keepInteresting(true) {}

template <typename T>
DataPointsFiltersImpl<T>::InterestPointSamplingDataPointsFilter::InterestPointSamplingDataPointsFilter(const Parameters& params) :
DataPointsFilter("InterestPointSamplingDataPointsFilter", InterestPointSamplingDataPointsFilter::availableParameters(), params),
		vSizeX(Parametrizable::get<T>("vSizeX")),
		vSizeY(Parametrizable::get<T>("vSizeY")),
		vSizeZ(Parametrizable::get<T>("vSizeZ")),
		perc(Parametrizable::get<T>("perc")),
		keepInteresting(Parametrizable::get<bool>("keepInteresting"))
{

}

template <typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::InterestPointSamplingDataPointsFilter::filter(const DataPoints& input)
{
    DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

template <typename T>
void DataPointsFiltersImpl<T>::InterestPointSamplingDataPointsFilter::inPlaceFilter(DataPoints& cloud)
{

	// Assert cloud contains a 3D point cloud
	assert (cloud.features.rows() == 4);

    const unsigned int numPoints(cloud.features.cols());

	// TODO: Check that the voxel size is not too small, given the size of the data

	// Calculate number of divisions along each axis
	Vector minValues = cloud.features.rowwise().minCoeff();
	Vector maxValues = cloud.features.rowwise().maxCoeff();

    T minBoundX = minValues.x() / vSizeX;
    T maxBoundX = maxValues.x() / vSizeX;
    T minBoundY = minValues.y() / vSizeY;
    T maxBoundY = maxValues.y() / vSizeY;
    T minBoundZ = minValues.z() / vSizeZ;
    T maxBoundZ = maxValues.z() / vSizeZ;

    // number of divisions is total size / voxel size voxels of equal length + 1
    // with remaining space
    unsigned int numDivX = 1 + maxBoundX - minBoundX;
    unsigned int numDivY = 1 + maxBoundY - minBoundY;;
    unsigned int numDivZ = 1 + maxBoundZ - minBoundZ;


    unsigned int numVox = numDivX * numDivY * numDivZ;

    // Assume point cloud is randomly ordered
    // compute a linear index of the following type
    // i, j, k are the component indices
    // nx, ny number of divisions in x and y components
    // idx = i + j * nx + k * nx * ny
    // std::vector<unsigned int> indices(numPoints);

    // vector to hold the first point in a voxel
    // this point will be ovewritten in the input cloud with
    // the output value

    std::vector<Voxel>* voxels;

    // try allocating vector. If too big return error
    try {
    	voxels = new std::vector<Voxel>(numVox);
    } catch (std::bad_alloc&) {
    	throw InvalidParameter((boost::format("InterestPointSamplingDataPointsFilter: Memory allocation error with %1% voxels.  Try increasing the voxel dimensions.") % numVox).str());
    }


    for (unsigned int p = 0; p < numPoints; p++ )
    {
        unsigned int i = floor(cloud.features(0,p)/vSizeX - minBoundX);
        unsigned int j = floor(cloud.features(1,p)/vSizeY - minBoundY);
        unsigned int k = floor(cloud.features(2,p)/vSizeZ - minBoundZ);
        unsigned int idx = i + j * numDivX + k * numDivX * numDivY; // Vox index

        Vector3 point = cloud.features.col(p).head(3);

		(*voxels)[idx].sumOuterProd += point * point.transpose();
		(*voxels)[idx].sumPoints += point;
		(*voxels)[idx].pointIndices.push_back(p);

        // indices[p] = idx;

    }

    // Now iterate through voxes computing dispersion ratio

    // NOTE: maybe it's better to use adjoint matrix to get the eigenvalues, using something from the code snipped below?
    //    Matrix3 centered = mat.topRows(3).colwise() - mat.topRows(3).rowwise().mean();
	//    Matrix3 cov = centered.adjoint() * centered;
	//    Eigen::SelfAdjointEigenSolver<Matrix3> eig(cov);

    for(unsigned int idx = 0; idx < numVox; idx++)
    {
        unsigned int numPointsInVox = (*voxels)[idx].pointIndices.size();
        if(numPointsInVox > 0)
        {
        	// Compute the covariance matrix of the voxel
        	// TODO: Is this prone to subject cancelation when using floats? Maybe this and the covariance values should always be double?
        	Matrix3 cov = ((*voxels)[idx].sumOuterProd / (T)numPointsInVox) - ((*voxels)[idx].sumPoints/(T)numPointsInVox)*((*voxels)[idx].sumPoints/(T)numPointsInVox).transpose();

        	// Solve to get eigenvalues and eigenvectors
        	Eigen::EigenSolver<Matrix3> eig(cov);

        	// Get min and max eigenvalues
        	T minEigValue = eig.eigenvalues().real().minCoeff(), maxEigValue = eig.eigenvalues().real().maxCoeff();

        	// Compute the dispersion ration
        	(*voxels)[idx].dispersionRatio = maxEigValue/minEigValue;

        }
    }

    // Sort voxes by decreasing dispersionRatio (empty voxels will be at the end)
    std::sort((*voxels).begin(), (*voxels).end(), greater_than_voxel<T>);

    // Gather all points to keep
    std::vector<unsigned int> pointsToKeep;
    if(keepInteresting)
    {
    	for(unsigned int idx = 0; idx < numVox; idx++)
    	{
			std::copy ((*voxels)[idx].pointIndices.begin(),(*voxels)[idx].pointIndices.end(),back_inserter(pointsToKeep));
			if( ((T)pointsToKeep.size() / (T)numPoints) >= perc ) break;
    	}
    }
    else
    {
		for(unsigned int idx = 0; idx < numVox; idx++)
    	{
			std::copy ((*voxels)[(numVox-1)-idx].pointIndices.begin(),(*voxels)[(numVox-1)-idx].pointIndices.end(),back_inserter(pointsToKeep));
			if( ((T)pointsToKeep.size() / (T)numPoints) >= 1 - perc ) break;
    	}
    }

    // deallocate memory for voxels information
    delete voxels;

    // Move the points to be kept to the start
    // Bring the data we keep to the front of the arrays then
	// wipe the leftover unused space.
	std::sort(pointsToKeep.begin(), pointsToKeep.end());
	int numPtsOut = pointsToKeep.size();
	for (int i = 0; i < numPtsOut; i++){
		int k = pointsToKeep[i];
		assert(i <= k);
		cloud.features.col(i) = cloud.features.col(k);
		if (cloud.descriptors.rows() != 0)
			cloud.descriptors.col(i) = cloud.descriptors.col(k);
	}
	cloud.features.conservativeResize(Eigen::NoChange, numPtsOut);

	if (cloud.descriptors.rows() != 0)
		cloud.descriptors.conservativeResize(Eigen::NoChange, numPtsOut);

}

template struct DataPointsFiltersImpl<float>::InterestPointSamplingDataPointsFilter;
template struct DataPointsFiltersImpl<double>::InterestPointSamplingDataPointsFilter;

// SegmentGroundDataPointsFilter
// Constructor
template<typename T>
DataPointsFiltersImpl<T>::SegmentGroundDataPointsFilter::SegmentGroundDataPointsFilter(const Parameters& params):
	DataPointsFilter("SegmentGroundDataPointsFilter", SegmentGroundDataPointsFilter::availableParameters(), params),
		seedRadius    (Parametrizable::get<T>   ("seedRadius")),
		seedBase      (Parametrizable::get<T>   ("seedBase")),
		gpLengthScale (Parametrizable::get<T>   ("gpLengthScale")),
		gpSignalVar   (Parametrizable::get<T>   ("gpSignalVar")),
		gpNoiseVar    (Parametrizable::get<T>   ("gpNoiseVar")),
		zGroundVar    (Parametrizable::get<T>   ("zGroundVar")),
		thModel       (Parametrizable::get<T>   ("thModel")),
		thData        (Parametrizable::get<T>   ("thData")),
		keepGround    (Parametrizable::get<bool>("keepGround")),
		keepObjects   (Parametrizable::get<bool>("keepObjects")),
		keepUnknown   (Parametrizable::get<bool>("keepUnknown"))
{
}

template<typename T>
typename PointMatcher<T>::DataPoints DataPointsFiltersImpl<T>::SegmentGroundDataPointsFilter::filter(
	const DataPoints& input)
{
	DataPoints output(input);
	inPlaceFilter(output);
	return output;
}

// In-place filter
template<typename T>
void DataPointsFiltersImpl<T>::SegmentGroundDataPointsFilter::inPlaceFilter(
	DataPoints& cloud)
{
	const unsigned int numPoints(cloud.features.cols());

	// Get seed points (they'll be the first inliers) and put them in the begining of the cloud
	unsigned int seed_offset = 0;
	for (unsigned int curr_offset = 0; curr_offset < numPoints; curr_offset++ )
	{
		if( cloud.features.col(curr_offset)(2) <= seedBase &&
		    cloud.features.col(curr_offset).head(2).norm() <= seedRadius )
		{
			if(curr_offset != seed_offset)
			{
				// swap curr_offset and seed_offset
				cloud.features.col(curr_offset).swap(cloud.features.col(seed_offset));
				// swap descriptors if any
				if (cloud.descriptors.rows() != 0)
					cloud.descriptors.col(curr_offset).swap(cloud.descriptors.col(seed_offset));
			}
			seed_offset++;
		}
	}

	// Main loop of gp-insac
	unsigned int inliers_size = 0;
	unsigned int new_inliers_size = seed_offset;
	unsigned int outliers_offset, outliers_size, unknowns_offset, unknowns_size;
	while(new_inliers_size > 0)
	{
		// DEBUG
		cout << "new_inliers_size = " << new_inliers_size << endl;

		// 1) Generate the current model from current inliers
		inliers_size += new_inliers_size;
		Matrix K(inliers_size,inliers_size); K << createKernel(cloud.features.block(0,0,2,inliers_size), cloud.features.block(0,0,2,inliers_size));

		// 2) Get the estimated mean and variance for the test points
		unsigned int test_size = cloud.features.cols() - inliers_size;
		// use cholesky to get Linv of (K + gpNoiseVar*I)
		Matrix L( (K + gpNoiseVar*Matrix::Identity(inliers_size,inliers_size)).llt().matrixL() );
		Matrix Linv( L.inverse() );
		// compute the estimated mean values for the test points
		Matrix Kstar(inliers_size,test_size); Kstar << createKernel(cloud.features.block(0,0,2,inliers_size),cloud.features.block(0,inliers_size,2,test_size));
		Vector fstar(test_size); fstar << Kstar.transpose() * (Linv.transpose() * (Linv * cloud.features.row(2).segment(inliers_size,test_size).transpose()));
		// compute the estimated variance for the test points
		Matrix v(Linv * Kstar);
		Matrix Kstar2(test_size,test_size); Kstar2 << createKernel(cloud.features.block(0,inliers_size,2,test_size),cloud.features.block(0,inliers_size,2,test_size));
		Vector Vstar((Kstar2 - v.transpose() * v).diagonal());

		// 3) Classify the test points into inliers, outliers or unknown
		// At this point the all the inliers will be at the begining of the cloud
		// and we consider all test points unknown. Then we sort the cloud putting
		// the new inliers together with the old inliers, and the outliers at the
		// end of the cloud. So, what we're doing here is going from:
		//
        // ,- inliers_offset (always zero)
		// |                   ,- unknowns_offset(=inliers_size)          ,- outliers_offset(=inliers_size+test_size)
 		// v                   v                                          v
		// ---------------------------------------------------------------
		// |      inliers     |         unknowns (=test points)          |
		// ---------------------------------------------------------------
		//                                                                ^- outliers_size(=0)
		//                     ^--------unknowns_size(=test_size)--------^
		//                     ^-new_inliers_size(=0)
		// ^---inliers_size---^
		//
		// to:
		//
		// ,- inliers_offset (always zero)
		// |                       ,- unknowns_offset    ,- outliers_offset
 		// v                       v                     v
		// ---------------------------------------------------------------
		// |      inliers     !   |       unknowns      |    outliers    |
		// ---------------------------------------------------------------
		//                                               ^-outliers_size-^
		//                         ^----unknowns_size---^
		//                     ^--^
		//                       `new_inliers_size
		// ^---inliers_size-- ^

		// reset class pointers and counters
		// inliers_offset is always zero.
		new_inliers_size = 0;
		unknowns_offset = inliers_size;
		unknowns_size = test_size;
		outliers_offset = inliers_size + unknowns_size;
		outliers_size = 0;
		unsigned i = 0;
		while(i < outliers_offset - inliers_size)
		{
			if(sqrt(Vstar(i)) < thModel)
			{
				double mahal_dist = (cloud.features.row(2).segment(inliers_size,test_size)(i) - fstar(i))/sqrt(zGroundVar + Vstar(i));
				if(mahal_dist < thData)
				{
					if(new_inliers_size != i)
					{
						// swap features
						cloud.features.col(inliers_size + new_inliers_size).swap(cloud.features.col(inliers_size + i));
						// swap descriptors if any
						if (cloud.descriptors.rows() != 0)
							cloud.descriptors.col(inliers_size + new_inliers_size).swap(cloud.descriptors.col(inliers_size + i));
					}
					// adjust pointers and counters
					new_inliers_size++;
					unknowns_offset++;
					unknowns_size--;
					i++;
				}
				else
				{
					// point is obstacle; add to outliers
					outliers_offset--;
					if(inliers_size + i != outliers_offset)
					{
						// swap features
						cloud.features.col(inliers_size + i).swap(cloud.features.col(outliers_offset));
						// swap descriptors if any
						if (cloud.descriptors.rows() != 0)
							cloud.descriptors.col(inliers_size + i).swap(cloud.descriptors.col(outliers_offset));
					}
					outliers_size++;
					unknowns_size--;
				}
			}
			else
			{
				// Point is unknown; just advance to next point
				i++;
			}
		}
	}

	// DEBUG: Save point clouds of each class
	if(inliers_size > 0)
	{
		DataPoints ground_cloud(cloud.features.block(0,0,cloud.features.rows(),inliers_size), cloud.featureLabels);
		ground_cloud.save("ground.pcd");
	}

	if(outliers_size > 0)
	{
		DataPoints object_cloud(cloud.features.block(0,outliers_offset,cloud.features.rows(),outliers_size), cloud.featureLabels);
		object_cloud.save("object.pcd");
	}

	if(unknowns_size > 0)
	{
		DataPoints unknown_cloud(cloud.features.block(0,unknowns_offset,cloud.features.rows(),unknowns_size), cloud.featureLabels);
		unknown_cloud.save("unknown.pcd");
	}

	// Crop the point cloud based on the flags keepGround, keepObjects and keepUnknown
	unsigned int curr_offset = 0;
	if(keepGround)
	{
		// All ground points should be already at the begining of the cloud. Just adjust curr_offset
		curr_offset += inliers_size;
	}
	if(keepUnknown)
	{
		if(unknowns_size != 0 && curr_offset != unknowns_offset)
		{
			cloud.features.block(0,curr_offset,cloud.features.rows(),unknowns_size) = cloud.features.block(0,unknowns_offset,cloud.features.rows(),unknowns_size);
			if (cloud.descriptors.rows() != 0)
				cloud.descriptors.block(0,curr_offset,cloud.descriptors.rows(),unknowns_size) = cloud.descriptors.block(0,unknowns_offset,cloud.descriptors.rows(),unknowns_size);
		}
		curr_offset += unknowns_size;
	}
	if(keepObjects)
	{
		if(outliers_size != 0 && curr_offset != outliers_offset)
		{
			cloud.features.block(0,curr_offset,cloud.features.rows(),outliers_size) = cloud.features.block(0,outliers_offset,cloud.features.rows(),outliers_size);
			if (cloud.descriptors.rows() != 0)
				cloud.descriptors.block(0,curr_offset,cloud.descriptors.rows(),outliers_size) = cloud.descriptors.block(0,outliers_offset,cloud.descriptors.rows(),outliers_size);
		}
		curr_offset += outliers_size;
	}

	cloud.features.conservativeResize(Eigen::NoChange, curr_offset);
	if (cloud.descriptors.rows() != 0)
		cloud.descriptors.conservativeResize(Eigen::NoChange, curr_offset);

}

template<typename T>
typename DataPointsFiltersImpl<T>::SegmentGroundDataPointsFilter::Matrix
DataPointsFiltersImpl<T>::SegmentGroundDataPointsFilter::createKernel(
	const Matrix & Xp, const Matrix & Xq)
{
	Matrix K(Xp.cols(), Xq.cols());
	Matrix aux(Xq.rows(),Xq.cols());
	for (unsigned int i = 0; i < Xp.cols(); i++)
	{
		aux = -Xq;
		aux.colwise() += Xp.col(i);
		aux *= 1/gpLengthScale;
		K.row(i) = gpSignalVar*Eigen::exp(-(1/2)*(aux.row(1).array().pow(2) + aux.row(2).array().pow(2)));
	}

	return K;
}


template struct DataPointsFiltersImpl<float>::SegmentGroundDataPointsFilter;
template struct DataPointsFiltersImpl<double>::SegmentGroundDataPointsFilter;
