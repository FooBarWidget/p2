#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>

#include <iostream>
#include <stdio.h>

using namespace std;

template<unsigned int nQuantiles = 1>
class P2 {
private:
	static const unsigned int nMarkers = 2 + nQuantiles * 3;

	unsigned int count;
	int positions[nMarkers];            // n
	double heights[nMarkers];           // q
	double increments[nMarkers];        // dn'
	double desiredPositions[nMarkers];  // n'

	void initializeState() {
		memset(positions, 0, sizeof(positions));
		memset(heights, 0, sizeof(heights));
		memset(increments, 0, sizeof(increments));
		memset(desiredPositions, 0, sizeof(desiredPositions));
		increments[0] = 0.0;
		increments[1] = 1.0;
	}

	// An insertion sort implementation for sorting small arrays.
	// Faster than std::sort() which likely uses quicksort or mergesort and
	// is more suitable for larger arrays.
	template<typename ValueType>
	static void smallsort(ValueType *array, unsigned short size) {
		for (unsigned short i = 0; i < size; i++) {
			int j = i;
			while (j > 0 && array[j] < array[j - 1]) {
				ValueType temp = array[j];
				array[j] = array[j - 1];
				array[j - 1] = temp;
				j--;
			}
		}
	}

	double parabolic(int i, int d) const {
		return heights[i]
			+ d
			/ double(positions[i + 1] - positions[i - 1])
			* ((positions[i] - positions[i - 1] + d)
				* (heights[i + 1] - heights[i])
				/ (positions[i + 1] - positions[i])
				+ (positions[i + 1] - positions[i] - d)
				* (heights[i] - heights[i - 1])
				/ (positions[i] - positions[i - 1])
			);
	}

	double linear(int i, int d) const {
		return heights[i] + d * (heights[i + d] - heights[i]) / (positions[i + d] - positions[i]);
	}

	int sign(double value) const {
		if (value >= 0.0) {
			return 1.0;
		} else {
			return -1.0;
		}
	}

	void prepareAlgorithmInitialization(double value) {
		heights[count] = value;
		count++;
		if (count == nMarkers) {
			initializeAlgorithm();
		}
	}

	void initializeAlgorithm() {
		smallsort(heights, nMarkers);
		for (unsigned i = 0; i < nMarkers; i ++) {
			positions[i] = i + 1;
		}
	}

	void runAlgorithm(double value) {
		unsigned int cellIndex;
		unsigned int i;

		// Algorithm step B.1
		if (value < heights[0]) {
			heights[0] = value;
			cellIndex = 1;
		} else if (value >= heights[nMarkers - 1]) {
			heights[nMarkers - 1] = value;
			cellIndex = nMarkers - 1;
		} else {
			cellIndex = 1;
			for (i = 1; i < nMarkers; i++) {
				if (value < heights[i]) {
					cellIndex = i;
					break;
				}
			}
		}

		// Algorithm step B.2
		for (i = cellIndex; i < nMarkers; i++) {
			positions[i]++;
			desiredPositions[i] = desiredPositions[i] + increments[i];
		}
		for (i = 0; i < cellIndex; i++) {
			desiredPositions[i] = desiredPositions[i] + increments[i];
		}

		// Algorithm step B.3
		for (i = 1; i < nMarkers - 1; i++) {
			double d = desiredPositions[i] - positions[i];
			if ((d >=  1.0 && positions[i + 1] - positions[i] > 1)
			 || (d <= -1.0 && positions[i - 1] - positions[i] < -1.0))
			{
				double newq = parabolic(i, sign(d));
				if (heights[i - 1] < newq && newq < heights[i + 1]) {
					heights[i] = newq;
				} else {
					heights[i] = linear(i, sign(d));
				}
				positions[i] += sign(d);
			}
		}
	}

public:
	P2()
		: count(0)
	{
		initializeState();
	}

	P2(double quantile)
		: count(0)
	{
		assert(nQuantiles == 1);
		initializeState();
		setQuantile(0, quantile);
		finalizeQuantiles();
	}

	void setQuantile(unsigned int index, double quantile) {
		assert(index < nQuantiles);
		increments[2 + index * 3 + 0] = quantile;
		increments[2 + index * 3 + 1] = quantile / 2.0;
		increments[2 + index * 3 + 2] = (1.0 + quantile) / 2.0;
	}

	void finalizeQuantiles() {
		smallsort(increments, nMarkers);
		for (unsigned int i = 0; i < nMarkers; i++) {
			desiredPositions[i] = (nMarkers - 1) * increments[i] + 1;
		}
	}

	void add(double value) {
		if (count < nMarkers) {
			prepareAlgorithmInitialization(value);
		} else {
			runAlgorithm(value);
		}
	}

	double result() {
		return result(increments[2]);
	}

	double result(double quantile) {
		if (count < nMarkers) {
			unsigned int closest = 1;
			smallsort(heights, count);
			for (unsigned int i = 2; i < count; i++) {
				if (fabs(double(i) / count - quantile) < fabs(double(closest) / nMarkers - quantile)) {
					closest = i;
				}
			}
			return heights[closest];
		} else {
			// Figure out which quantile is the one we're looking for by nearest increment.
			unsigned int closest = 1;
			for (unsigned int i = 2; i < nMarkers - 1; i ++) {
				if (fabs(increments[i] - quantile) < fabs(increments[closest] - quantile)) {
					closest = i;
				}
			}
			return heights[closest];
		}
	}
};

/* int
main() {
	P2<> p2;
	p2.setQuantile(0, 0.96);
	p2.finalizeQuantiles();
	p2.add(1);
	p2.add(5);
	p2.add(9);
	p2.add(2);
	p2.add(4);
	printf("%f\n", p2.result());
	return 0;
}
*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    // We will calculate the 10, 50, and 90th percentile in several different ways to compare accuracy
    // p2_10, p2_50 and p2_90 will each be used to calculate a single percentile
	P2<> p2_10(0.1), p2_50(0.5), p2_90(0.9);
    // multi will be used to simultaneously calculate the 10, 50 and 90 percentiles
	P2<3> multi;
	FILE *fi;
	double d;
	char buf[30];

	multi.setQuantile(0, 0.1);
	multi.setQuantile(1, 0.5);
	multi.setQuantile(2, 0.9);
	multi.finalizeQuantiles();

	if( argc != 2 ) {
		cout << "No file specified" << endl;
		exit(1);
	}

	fi = fopen(argv[1], "r");

	if( fi == NULL ) {
		printf( "Failed to open %s\n", argv[1] );
		exit(1);
	}

	int iround = 0;

	// Read data from a file, one value per line. Add it into each of the p2_t trackers
	while(true) {
		fgets(buf, 30, fi);
        if(feof(fi)) break;
		d = strtod(buf, NULL);

		//printf("--------------------------> round %d\n", iround++);
		p2_10.add( d );
		p2_50.add( d );
		p2_90.add( d );
		multi.add( d );
	}
	fclose(fi);

	// Print out the results from each of the trackers at the 10, 50 and 90 percentiles
	printf("%g %g %g\n", p2_10.result( ), p2_50.result( ), p2_90.result( ) );
	printf("%g %g %g\n", multi.result( 0.1 ), multi.result( 0.5 ), multi.result( 0.9 ) );
}
