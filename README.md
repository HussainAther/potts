# No Pseudolikelihood for Me: Training Potts Models with Contrastive Divergence for Protein Design

The current 250 million known protein sequences are the result of the generative process evolution. Similar sequences that have likely evolved from the same ancestral sequence, also known as homologs, are clustered in families stored in the Pfam database. Potts models trained on a Pfam family can predict a common protein 3D structure to the sequences in that family, predict the beneficial or detrimental effect of mutations, and recently be sampled to generate novel sequences for protein design.

A protein is a string of amino acids, also called residues, that in cells can adopt a 3-dimensional fold. These folds mean that two residues far apart in the sequence can be nearby in 3D space. One assumption underlying all protein science is that proteins structures prefer to adopt low energy states. Interactions with neighbouring atoms should be locally favourable. Thus, for mutations to not disrupt protein folds, residues nearby in 3D space tend to coevolve, i.e. a residue’s evolutionary trajectory is dependent on nearby residues. Potts models are attractive for protein modelling because of its inductive bias. In particular the pairwise term can capture coevolutionary information in sequences.


## Project Deliverable Goals

On Pfam family and HHblits result of 6EQD, 4EB0, 6QZ4 do:
1. Train a Potts model with bmDCA, pseudolikelihood, and contrastive divergence. Show that negative log-likelihood is indeed being minimized.
2. Evaluate model: goal is to have high log-likelihood for real data, low log-likelihood for detrimental mutations/shuffled sequences/nonsensical sequences/fragmented sequences. Plot results in three bar graphs, one for each training method.
3. Draw 10000 samples from models trained with the three training methods. Report the log-likelihood for each sample.
4. Visualize sequences with a sequence logo, compute average sequence identity to training sequences.

Nice-to-haves
1. Using an independent homology predictor, compute the fraction of sampled sequences that are homologous to the wildtype from each training method
2. Model the structure of some generated sequences from each training method and compute the Rosetta energy
3. Introduce a temperature parameter in the Potts model to control sequence diversity
4. Try another training strategy: No MCMC for me
5. Try another discrete sampling strategy: Oops I took a gradient

## Potts Model and Protein Covariation: https://ronlevygroup.cst.temple.edu/courses/2017_spring/chem5412/lectures/Lecture4_SBII_2017.pdf

+ HMMs can struggle with sequence identities lower than 20%, which can still have similar structure (e.g. fibronectin)
+ History: coevolution implies structural contacts
+ Can use mutual information or statistical coupling, but these model local correlations
+ Potts model can model indirect correlations, beyond pairwise (A → B → C could mean A and C are also correlated)
+ Start from maximum entropy distribution for p(sequence), then add constraints (first order and second order amino acid frequencies in data), solve using Lagrange multipliers
+ Rearrangement gives a Boltzmann energy functional form, where energies comes from single and double pairs from data
+ Applications
+ + Contact prediction, structure prediction from Direct Coupling Analysis from contact maps
+ + Free energy landscapes
+ + Fitness landscapes

## Inverse Statistical Physics of Protein Sequences: A Key Issues Review: https://arxiv.org/abs/1703.01222

+ Use of sequences to infer Boltzmann distributions
+ + Homology detection
+ + Traditionally use profile models, which are Potts models with only single-site terms
+ + HMMs are also profile models but with insertion/deletion states
+ Protein contact map prediction
+ + Pairwise correlations motivated from hypothesis: residues in contact coevolve
+ + Just using correlations does not work because of chaining effects i → j → k
+ + Potts models work well at picking out direct couplings from indirect correlations
Mutation effect, fitness landscapes
Jij represent epistatic contributions
But function could depend non-linearly on difference in Potts energies
Protein design
Sampling from Potts models may be enough to design new sequences!
Requires model to be fully generative: samples able to produce statistics in MSA, even beyond pairwise
A Potts model can be fitted for each protein family (MSA of homologs)
For a protein of length 263, you have about 14.5 million parameters!
Regularization is therefore needed since MSAs typically have 100 to 10,000 seqs
Also add pseudocounts with a Dirichlet prior
Inference
