# Capstone Project Requirement Details
## Module: Primer Collection
### Description:
This module identifies all possible upstream and downstream primers for amplifying a given target region using [PCR](https://en.wikipedia.org/wiki/Polymerase_chain_reaction)). Given a string of nucleotides and the position of the target region within the string, the upstream primers are collected from the sequence to the left of the target and the downstream primers are collected from the sequence to the right. The primers will represent all possible k-mers (over a user-specified range of k values) in the string. A number of characteristics is calculated for each primer and they are stored in a database.
### Input:
#### Required from user:
1. A string of nucleotide sequence (in 5' to 3' orientation) for each desired amplicon (can include [IUPAC ambiguous base codes](https://droog.gs.washington.edu/parc/images/iupac.html)). The sequence should include the target region and flanking upstream and downstream primer regions. The length of the sequence will generally range from 100-2000 characters.
    - **Interface**:
        - *ideal:* Should be able to upload sequences either from a local [fasta](https://en.wikipedia.org/wiki/FASTA_format) file and/or allow them to be pasted into a text area. When loading from a fasta, the fasta header should be automatically filled in as the unique sequence identifier but it should still be editable. When pasting into a text area, each sequence should be added separately so there should be a add next sequence button that will provide a new text area for each additional sequence. 
        - *basic* Allow only one of these methods
    - **Validation**: No unexpected characters. Should only contain (A, C, T, G, M, R, W, S, Y, K, V, H, D, B, N, X) and X should only be used to identify regions in the primer zone to avoid (see below).
    - **Hints**:
        - It needs to be clear to the user that they need to provide a continuous sequence that includes their desired upstream and downstream primer regions, not just their target amplicon.
1. Target region identification
    - **Interface**:
        - *ideal:* Provide input fields for entering the start and end positions of the target region within each sequence and highlight the selected region in the sequence when complete. and/or allow user to highlight the desired region in the sequence and parse the range automatically.
        - *basic:* Provide input fields for entering the start and end positions of the target region for each sequence.
    - **Validation** The range must be bounded by the size of the provided sequence and should provide enough sequence (greater than 50 nucleotides) on both sides of the target region. Start should be less than end. There must be a target region for each sequence and it should include at least one nucleotide.
    - **Hints**:
        Describe how to use the interface (for example: let user know if we are indexing from 0 or 1, indexing from 1 may be more intuitive for users but on the back end it needs to be 0-indexed). Validator hints should provide explanation for any errors. 
1. Regions in the primer zone to avoid (optional). If a user wants a region of the primer zone to be avoided in primer design, they may replace the characters in that region of the sequence with 'X's.
    - **Interface**
        - *ideal:* Highlight any X's in the sequence and identify them as regions that will be excluded from primer design.
        - *basic:* Notify user that they have excluded regions in a particular sequence.
    - **Validation**: There must be at least 50 continuous nucleotides left in both the upstream and downstream primer regions after sites are excluded.
1. Unique identifier for each sequence
    - **Interface**:
        - *ideal:* This should be a text field that is automatically filled in when loading sequences from a fasta file.
        - *basic:* A text field that is filled in by user.
    - **Validation**: Up to 120 characters not containing any tabs or commas.
    - **Hints**:
        - Validator hint that warns if there are duplicate names or invalid names (with reason).
1. Size range for the primers.
    - **Interface:**
        - *ideal:* A way to indicate the range of primer sizes (range/slider or min/max fields) only allow range between 8-40.
        - *basic:* Fields for min and max that validate if outside of range.
    - **Default:** Min=18, Max=22
    - **Validation**: Must be between 8 and 40 and min must be smaller than max.
    - **Hints**: If max is higher than default warn that longer primers will have difficulty binding at PCR annealing temperature. If min is lower than default warn that shorter primers may not have adequate binding specificity.
1. Paths to background sequences (in fasta format) that should be used in calculating primer specificity. (optional)
    - **Interface**:
        - *ideal:* Query on-line database for sequence data and retrieval (via ftp) and accept paths to local fasta files.
        - *basic:* Take paths to local fasta files.
    - **Validation**: Check that file paths are valid
    - **Hints:** User should provide sequence data for any background DNA that may be in the PCR. This allows the program to determine if a primer is specific to the target region or if it will bind in many places in the DNA resulting in poor amplification.
1. PCR Concentration Parameters for calculating melting temperature (in Millimolar)
    - **Default:**
        1. Na=50
        1. K=0
        1. Tris=0
        1. Mg=0 
        1. dNTPs=0
    - **Interface**: In a delineated section dedicated to PCR concentration parameters, provide fields for each value with defaults prefilled.
    - **Validation**: They must be positive values
    - **Hints**: The units should be clearly labeled (mM) next to the input field.

#### Pseudoschema for primer collection module input json
```
{
    primer_collection:{
        params:{
            pcr_salts:{
                na: {
                    type: number,
                    default: 50},
                k: {
                    type: number,
                    default: 0},
                tris: {
                    type: number,
                    default: 0}
                mg: {
                    type: number,
                    default: 0}
                dntps: {
                    type: number,
                    default: 0}
            },
            background_seq: {
                type: array,
                required: false
                description: paths to background sequences
            }
        }
    },
    sequences: {
        description: sequence objects
        required: at least one
            seq_id: {
                seq: {
                    type: string,
                    required: true},
                target_start: {
                    type: number,
                    required: true},
                target_end: {
                    type: number,
                    required: true},
                primer_len_range: {
                    min: {
                        type: number,
                        default: 18},
                    max: {
                        type: number,
                        default: 22}
                }
            }
        }
    }
```

#### Pseudoschema for primer collection module output json
```
{ ...
    sequences: {
        description: sequence objects
        required: at least one
        seq_id: {
        ...
            primers: {
                forward: {
                    type: primer objects
                    description: primers for upstream region,
                    primer_id: {
                        seq_id: {
                            type: string,
                            description: Name of the target sequence from which the primer originates
                        }
                        seq: {
                            type: string,
                            description: primer sequence
                        },
                        flank: {
                            type: boolean,
                            description: 0 for upstream primer 1 for downstream primer
                        },
                        start: {
                            type: number,
                            description: index of starting position of primer in sequence
                        },
                        tm: {
                            type: number,
                            description: melting temperature of primer
                        },
                        gc: {
                            type: number,
                            description: percent GC content
                        },
                        length: {
                            type: number,
                            description: length of primer sequence
                        },
                        degenerate: {
                            type: number,
                            description: number of degenerate/ambiguous bases in sequence
                        },
                        homopolymers: {
                            description: runs of the same nuclotide longer than 3,
                            percent: {
                                type: number,
                                description: percent of primer that is part of homopolymer
                            },
                            run: {
                                type: number,
                                description: length of the longest run of homopolymer
                            }
                        },
                        dimerization: {
                            description: how likely the primer is to bind to itself,
                            percent: {
                                type: number,
                                description: percent of nucletides that self-hybridize
                            },
                            run: {
                                type: number,
                                description: length of longest run of complementary bases 
                            },
                            median: {
                                type: number,
                                description: median length of complementary bases
                            },
                            hairpin: {
                                type: boolean,
                                description: if 3' end hybridizes to self forming secondary structure that blocks amplification
                            }
                        },
                        specificity: {
                            description: how many times primer appears in background DNA,
                            mm0: {
                                type: number,
                                description: with no mismatches
                            },
                            mm1: {
                                type: number,
                                description: with 1 mismatch
                            },
                            mm2: {
                                type: number,
                                description: with 2 mismatches
                            },
                            mm3: {
                                type: number,
                                description: with 3 mismatches
                            }
                        }
                    }
                },
                reverse: {
                    type: primer objects
                    description: primers for downstream region
                        ...
                }
            }
        }
    }
}
```
### Output:
Summary stats per sequence per flank (ideally graphical) for primer data. 
Include the total number of primer options for each. At this point user can see if any of the targets have very poor primer options.
To mitigate this, the user should be able to:
1. Remove a target completely
    - This will require the sequence objects for the removed targets to be removed from the json
    - Primer collection does not need to be rerun
1. Add a new target/modify an existing target
    - Primer collection will need to be rerun for these targets and data will be merged with previous json
1. Change any of the parameters
    - Primer collection will need to be rerun for all targets, previous json can be passed to primer collection to reduce the amount of repeated work.

## Module: Primer Scoring
This module scores groups of primer for each target and each flank based on user provided weights for different features.
### Input:
#### Required from user:
1. Optimal melting temperature (°C) for the reaction.
    - **Interface:** Input field 
        - **Validation:** Number between 40-80°C
        - **Default:** 55°C
1. Optimal range for GC content (%)
    - **Interface:** Input field or slider/range for min and max GC content.
        - **Validation:** Numbers between 0 and 100, min <= max
        - **Defaults:** 40%-60%
1. Relative importance of features
    - Tm: minimize the distance of primers from optimal melting temperature
    - GC: minimize distance of primer GC% out of optimal range
    - Homopolymers: minimize the amount of homopolymer in primer
    - dimerization: minimize the amount of self dimerization in primer
    - specificity: minimize the number of times primers appear in background
    - degenerate: minimize the number of degenerate/ambiguous bases in the primer
    - **Interface:** A slider or similar indicating the relative importance of each of these features when determining the final ranking of primers.
    - **Validation:** Each value can range from 0 to 100. 
    - **Default:** 1 for each feature.
    - **Hints:** A value of zero means that the feature will not be considered at all during scoring.
#### Pseudoschema for primer scoring module input json
```
{...
    primer_scores: {
        params: {
            tm_opt: {
                type: number,
                default: 55
            },
            gc_opt: {
                min: {
                    type: number,
                    default: 40
                },
                max: {
                    type: number,
                    default: 60
                }
            },
            weights: {
                tm: {
                    type: number,
                    default: 1
                },
                gc: {
                    type: number,
                    default: 1
                },
                homopolymer: {
                    type: number,
                    default: 1
                },
                dimer: {
                    type: number,
                    default: 1
                },
                specificity: {
                    type: number,
                    default: 1
                },
                degenerate: {
                    type: number,
                    default: 1
                }
            }
        }
    }
}
```
#### Pseudoschema for primer scoring module output json
```
{ ...
    sequences: {
        description: sequence objects
        required: at least one
        seq_id: {
        ...
            primers: {
                forward: {
                    type: primer objects
                    description: primers for upstream region,
                    primer_id: {
                        ...
                        score: {
                            type: number,
                            description: relative score within primers from same target and same flank (lower is better)
                        },
                        rank: {
                            type: number,
                            description: rank relative to primers from same target and same flank
                        }
                    }
                },
                reverse: {
                    ...
                }
            }
        }
    }
}
```
## Module: Set Optimization
### Input:
#### Required from user:
1. Max number of iterations for optimization.
    - **Interface:** input field
    - **Default:** 100
1. Optimum amplicon size range. This will help pick primer pairs that will produce an amplicon in this size range. 
    - **Interface** min/max input field. If none is provided, amplicon size will not be considered. User should be able to set a different range for each target sequence. (optional)
1. Max distance of target from a primer. If the amplicons are going to be sequenced, the target should be with the sequencing platforms read length. This can be set in several ways. (optional)
    1. Distance from forward primer
    2. Distance from reverse primer
    3. Distance from at least one primer
    4. Distance from both primers
    - **Interface** A radial input selecting one of these options and a input field for the distance.
    - **Default:** None
1. Background Primers. Allow user to add existing primers that can included in the optimization. (optional)
    - **Interface** Provide a text area to paste sequences in fasta format or do file upload.
1. Relative importance of features
    - Similar melting temperatures: minimizes the variation in melting temperatures for primers in the set
    - Primer scores: minimizes the individual primer scores for primers in the set to favor high ranking primers
    - Cross dimerization: minimizes the number of cross interactions between all primers
    - Amplicon size: minimizes distance outside of range of optimal amplicon size
    - Target distance from primer: minimizes distance of primer from target
    - **Interface:** A slider or similar indicating the relative importance of each of these features when optimizing sets of primers.
    - **Validation:** Each value can range from 0 to 100. 
    - **Default:** 1 for each feature, 0 for amplicon size and target distance from primer if no values are selected,
    - **Hints:** A value of zero means that the feature will not be considered at all during scoring.
1. Include primers/exclude primers. User should be able to manually remove and add certain primers.
    - **Interface"** 
        - *ideal* By default, the top 10% of primers from each target and each flank will be moved forward to be considered in the primer set optimization. However, the user should be able to modify this value for each primer group. There should be some indication of how many primers will be included from each set and it should update anytime the percent changes. Additionally, the user should be able to open up a table for each primer group that shows each of the primers in order of rank (best to worst) and summarizes their features. Within this table, the primers that are in the currently selected top X% should have marked check boxes indicating that they are to be included. The user should then be able to remove any of the checked boxes and check any additional primers they want to be included.
        - *basic* Just allow the user to modify what % of top primers move forward.

#### Pseudoschema for primer set optimization module input json
```
{...
    set_optimization: {
        params: {
            iter: {
                type: number,
                default: 100
            },
            amp_size: {
                min: {
                    type: number,
                    default: null
                },
                max: {
                    type: number,
                    default: null
                }
            },
            target_distance: {
                forward: {
                    type: number,
                    default: null
                },
                reverse: {
                    type: number,
                    default: null
                },
                any: {
                    type: number,
                    default: null
                },
                both: {
                    type: number,
                    default: null
                }
            },
            background: {
                type: background primer objects
                primer_id: {
                    seq: {
                        type: string
                    }
                }
            },
            weights: {
                tm: {
                    type: number,
                    default: 1
                },
                scores: {
                    type: number,
                    default: 1
                },
                cross_dimerization: {
                    type: number,
                    default: 1
                },
                size: {
                    type: number,
                    default: 1
                },
                target_dist: {
                    type: number,
                    default: 1
                }
            },
            include: {
                seq_id: {
                    forward: {
                        primer_ids: [
                            type: array
                            description: array of primer ids that should be included
                        ]},
                    reverse: {
                        ...
                    }
                },
                ...
            }
        }
    }
}
```
#### Pseudoschema for primer set optimization module output json
```
{ ...
    set_optimization: {
        result: [
            {seq_id: {forward: primer_id, reverse: primer_id}}
            ...
        ]

        }
}        
```
### Output:
*basic* provide a summary table of features for each selected primer. Provide the primers in fasta format that can be copied or saved locally.
*ideal* provide some additional graphics
    1. For each target sequence, illustrate the location of the primers and the target sequences and provide the expected size and sequence of the amplicon, allow the amplicon sequences to be saved to a file.
    2. Provide average, min, and max melting temperature for the primers and offer suggested annealing temperature for PCR.
    3. Show alignments of top X number of strong hybridizations.
    4. Color code extreme values in summary table
Have two separate options to go back at this point:
    1. rerun from scratch
        - go back one step and default to including top X% of primer options
    2. rerun with modifications
        - here is where the user can keep optimized primers for some targets but further optimize primers for other targets. In this case, the value for the included primers for the non-modified targets should include just the primer from the previous run.



