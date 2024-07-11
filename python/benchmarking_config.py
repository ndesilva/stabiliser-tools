import benchmarking.generators as gs

configs = [
    {
        "title": "Converting a random state vector (S1) to a quadratic form triple (S2)",
        "functions_to_time": [],
        "function_strings": [
            "brute force",
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
        "title": "Converting a state vector (S1) to a quadratic form triple (S2), "
                 "extremal cases",
        "functions_to_time": [],
        "function_strings": [
            "brute force",
            "stim",
            "our method"
        ],
        "generation_types": [],
        "generation_strings": [
            "best case stab state with assump",
            "worst case stab state without assump"
        ]
    },

    {
        "title": "Rejecting a non-stabiliser state",
        "functions_to_time": [],
        "function_strings": [
            "brute force",
            "stim",
            "our method"
        ],
        "generation_types": [],
        "generation_strings": [
            "worst case non-stab state"
        ]
    }
]
