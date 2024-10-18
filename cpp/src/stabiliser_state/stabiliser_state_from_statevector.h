#ifndef _FAST_STABILISER_STABILISER_STATE_FROM_VECTOR_H
#define _FAST_STABILISER_STABILISER_STATE_FROM_VECTOR_H

#include <complex>

#include "stabiliser_state.h"

namespace fst
{
	/// Convert a state vector of complex amplitudes into a stabiliser state object.
	///
	/// Assuming valid is faster, but will result in undefined behaviour if the state vector is not in fact a
	/// valid stabaliser state
	Stabiliser_State stabiliser_from_statevector(const std::vector<std::complex<float>> &statevector, bool assume_valid = false);

	/// ;)
	Stabiliser_State stab_in_the_dark(const std::vector<std::complex<float>> &statevector);

	/// Test wheter a state vector of complex amplitudes corresponds to a stabiliser state.
	bool is_stabiliser_state(const std::vector<std::complex<float>> &statevector);
}

#endif