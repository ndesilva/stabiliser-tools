#ifndef _FAST_STABILISER_STABILISER_STATE_FROM_VECTOR_H
#define _FAST_STABILISER_STABILISER_STATE_FROM_VECTOR_H

#include <span>
#include <complex>
#include <optional>

#include "stabiliser_state.h"

namespace fst
{
	/// Convert a state vector of complex amplitudes into a stabiliser state object.
	/// 
	/// Assuming valid is faster, but will result in undefined behaviour if the state vector is not in fact a
	/// valid stabaliser state
	std::optional<Stabiliser_State> make_stabalizer_assume_valid( const std::span<const std::complex<float>> statevector );
	std::optional<Stabiliser_State> make_stabalizer_and_validate( const std::span<const std::complex<float>> statevector );
}

#endif