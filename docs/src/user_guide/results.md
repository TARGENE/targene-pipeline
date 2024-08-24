# Understanding TarGene's Outputs

```json
{
    "TMLE_GLM_GLM": {
        "n": 485,
        "estimates": [
            {
                "std": 15.622429857504532,
                "n": 485,
                "type": "TMLEstimate",
                "IC": [],
                "estimate": -0.6186142282443571,
                "estimand": {
                    "outcome_extra_covariates": [
                        "Cheese intake",
                        "Number of vehicles in household"
                    ],
                    "type": "ATE",
                    "treatment_values": {
                        "1:15944765:G:A": {
                            "control": 0,
                            "case": 1
                        }
                    },
                    "outcome": "Body mass index (BMI)",
                    "treatment_confounders": {
                        "1:15944765:G:A": [
                            "PC1",
                            "PC2",
                            "PC3",
                            "PC4",
                            "PC5",
                            "PC6"
                        ]
                    }
                }
            },
            {
                "std": 10.820703990065084,
                "n": 485,
                "type": "TMLEstimate",
                "IC": [],
                "estimate": -0.404886039735767,
                "estimand": {
                    "outcome_extra_covariates": [
                        "Cheese intake",
                        "Number of vehicles in household"
                    ],
                    "type": "ATE",
                    "treatment_values": {
                        "1:15944765:G:A": {
                            "control": 1,
                            "case": 2
                        }
                    },
                    "outcome": "Body mass index (BMI)",
                    "treatment_confounders": {
                        "1:15944765:G:A": [
                            "PC1",
                            "PC2",
                            "PC3",
                            "PC4",
                            "PC5",
                            "PC6"
                        ]
                    }
                }
            }
        ],
        "type": "TMLE.JointEstimate",
        "cov": [
            [244.06031465264908, -71.3344199317966],
            [-71.3344199317966, 117.08763484061052]
        ],
        "estimand": {
            "type": "TMLE.JointEstimand",
            "args": [
                {
                    "outcome_extra_covariates": [
                        "Cheese intake",
                        "Number of vehicles in household"
                    ],
                    "type": "ATE",
                    "treatment_values": {
                        "1:15944765:G:A": {
                            "control": 0,
                            "case": 1
                        }
                    },
                    "outcome": "Body mass index (BMI)",
                    "treatment_confounders": {
                        "1:15944765:G:A": [
                            "PC1",
                            "PC2",
                            "PC3",
                            "PC4",
                            "PC5",
                            "PC6"
                        ]
                    }
                },
                {
                    "outcome_extra_covariates": [
                        "Cheese intake",
                        "Number of vehicles in household"
                    ],
                    "type": "ATE",
                    "treatment_values": {
                        "1:15944765:G:A": {
                            "control": 1,
                            "case": 2
                        }
                    },
                    "outcome": "Body mass index (BMI)",
                    "treatment_confounders": {
                        "1:15944765:G:A": [
                            "PC1",
                            "PC2",
                            "PC3",
                            "PC4",
                            "PC5",
                            "PC6"
                        ]
                    }
                }
            ]
        }
    },
    "OSE_GLM_GLM": {
        "n": 485,
        "estimates": [
            {
                "std": 15.695012511050656,
                "n": 485,
                "type": "OSEstimate",
                "IC": [],
                "estimate": -0.5009846779180457,
                "estimand": {
                    "outcome_extra_covariates": [
                        "Cheese intake",
                        "Number of vehicles in household"
                    ],
                    "type": "ATE",
                    "treatment_values": {
                        "1:15944765:G:A": {
                            "control": 0,
                            "case": 1
                        }
                    },
                    "outcome": "Body mass index (BMI)",
                    "treatment_confounders": {
                        "1:15944765:G:A": [
                            "PC1",
                            "PC2",
                            "PC3",
                            "PC4",
                            "PC5",
                            "PC6"
                        ]
                    }
                }
            },
            {
                "std": 10.820968353562348,
                "n": 485,
                "type": "OSEstimate",
                "IC": [],
                "estimate": -0.4036203473738083,
                "estimand": {
                    "outcome_extra_covariates": [
                        "Cheese intake",
                        "Number of vehicles in household"
                    ],
                    "type": "ATE",
                    "treatment_values": {
                        "1:15944765:G:A": {
                            "control": 1,
                            "case": 2
                        }
                    },
                    "outcome": "Body mass index (BMI)",
                    "treatment_confounders": {
                        "1:15944765:G:A": [
                            "PC1",
                            "PC2",
                            "PC3",
                            "PC4",
                            "PC5",
                            "PC6"
                        ]
                    }
                }
            }
        ],
        "type": "TMLE.JointEstimate",
        "cov": [
            [246.33341772203678, -71.41251770729579],
            [-71.41251770729579, 117.09335610879779]
        ],
        "estimand": {
            "type": "TMLE.JointEstimand",
            "args": [
                {
                    "outcome_extra_covariates": [
                        "Cheese intake",
                        "Number of vehicles in household"
                    ],
                    "type": "ATE",
                    "treatment_values": {
                        "1:15944765:G:A": {
                            "control": 0,
                            "case": 1
                        }
                    },
                    "outcome": "Body mass index (BMI)",
                    "treatment_confounders": {
                        "1:15944765:G:A": [
                            "PC1",
                            "PC2",
                            "PC3",
                            "PC4",
                            "PC5",
                            "PC6"
                        ]
                    }
                },
                {
                    "outcome_extra_covariates": [
                        "Cheese intake",
                        "Number of vehicles in household"
                    ],
                    "type": "ATE",
                    "treatment_values": {
                        "1:15944765:G:A": {
                            "control": 1,
                            "case": 2
                        }
                    },
                    "outcome": "Body mass index (BMI)",
                    "treatment_confounders": {
                        "1:15944765:G:A": [
                            "PC1",
                            "PC2",
                            "PC3",
                            "PC4",
                            "PC5",
                            "PC6"
                        ]
                    }
                }
            ]
        }
    }
}
```