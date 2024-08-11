import benchmarking.generators as gs
import benchmarking.brute_force.stab_state_check as brute

import sys
PATH_TO_LIBRARY = './build/ninja-multi-vcpkg/cpp/src/Release'
PATH_TO_STIM_MOCK = './build/ninja-multi-vcpkg/cpp/benchmarking/stim/Release'
sys.path.extend([PATH_TO_LIBRARY, PATH_TO_STIM_MOCK])
import fast as fst
import stim_mock as sm

configs = [
    {
        "pre_string": "converting S1 to efficient rep",
        "title": "Converting a random stabiliser state vector (S1) to an efficient representation",
        "functions_to_time": [
            # brute.stab_to_xmatr,
            sm.circuit_from_statevector,
            fst.stabiliser_state_from_statevector
        ],
        "function_strings": [
            # "brute force",
            "stim",
            "our method"
        ],
        "generation_types": [
            gs.random_stab_state_with_assump,
            gs.random_stab_state
        ],
        "generation_strings": [
            "random stab state with assump",
            "random stab state without assump"
        ]
    },

    {
        "pre_string": "converting S1 to efficient rep",
        "title": "Converting a stabiliser state vector (S1) to an efficient representation, "
                 "extremal cases",
        "functions_to_time": [
            # brute.stab_to_xmatr,
            sm.circuit_from_statevector,
            fst.stabiliser_state_from_statevector
        ],
        "function_strings": [
            # "brute force",
            "stim",
            "our method"
        ],
        "generation_types": [
            gs.best_case_stab_state_with_assump,
            gs.worst_case_stab_state
        ],
        "generation_strings": [
            "best case stab state with assump",
            "worst case stab state without assump"
        ]
    },

    {
        "pre_string": "non-stab reject",
        "title": "Rejecting a non-stabiliser state",
        "functions_to_time": [
            # brute.stab_to_xmatr,
            sm.circuit_from_statevector,
            fst.is_stabiliser_state
        ],
        "function_strings": [
            # "brute force",
            "stim",
            "our method"
        ],
        "generation_types": [
            gs.worst_case_almost_stab_state
        ],
        "generation_strings": [
            "worst case non-stab state"
        ]
    }
]
