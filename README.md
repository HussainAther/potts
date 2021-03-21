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
+ + ![](https://github.com/HussainAther/potts/blob/main/img/structuralcontacts.png)
+ Can use mutual information or statistical coupling, but these model local correlations
+ Potts model can model indirect correlations, beyond pairwise (A → B → C could mean A and C are also correlated)
+ Start from maximum entropy distribution for p(sequence), then add constraints (first order and second order amino acid frequencies in data), solve using Lagrange multipliers
+ Rearrangement gives a Boltzmann energy functional form, where energies comes from single and double pairs from data
+ + ![](https://github.com/HussainAther/potts/blob/main/img/boltzmannenergy.png)
+ Applications
  + Contact prediction, structure prediction from Direct Coupling Analysis from contact maps
  + Free energy landscapes
  + Fitness landscapes

## Inverse Statistical Physics of Protein Sequences: A Key Issues Review: https://arxiv.org/abs/1703.01222

+ Use of sequences to infer Boltzmann distributions
  + Homology detection
  + Traditionally use profile models, which are Potts models with only single-site terms
  + HMMs are also profile models but with insertion/deletion states
+ Protein contact map prediction
  + Pairwise correlations motivated from hypothesis: residues in contact coevolve
  + Just using correlations does not work because of chaining effects i → j → k
  + Potts models work well at picking out direct couplings from indirect correlations
+ Mutation effect, fitness landscapes
  + Jij represent epistatic contributions
  + But function could depend non-linearly on difference in Potts energies
+ Protein design
  + Sampling from Potts models may be enough to design new sequences!
  + Requires model to be fully generative: samples able to produce statistics in MSA, even beyond pairwise
+ A Potts model can be fitted for each protein family (MSA of homologs)
  + For a protein of length 263, you have about 14.5 million parameters!
  + Regularization is therefore needed since MSAs typically have 100 to 10,000 seqs
  + Also add pseudocounts with a Dirichlet prior
+ Inference
  + Boltzmann machine learning
  + + ![](https://github.com/HussainAther/potts/blob/main/img/boltzmannml.png)
  + Guaranteed to converge for small enough ε 
  + Runtime is still slow for sequence lengths > 200 
  + Approximation methods like Gaussian Approximation, Mean field approximation and Psudolikelihood Maximization does not accurately reproduce empirical frequencies!
  + But they are still able to predict contact maps if substantial coevolution signal exists. 
+ Applications
  + Contact map prediction
  + Potts models can be extended for two families to identify protein-protein interactions!
  + Predicting viral escape mutations
+ Design applications
  + Mutation effect
  + How far away a designed sequence is from an MSA of functional sequences
  + Seems that different inference techniques does not affect sequence-function prediction performance
+ Outlook
  + Amino acids are still represented as abstract symbols--perhaps should take advantage of their physicochemical properties or representations of them.
 
## Protein 3D Structure Computed from Evolutionary Sequence Variation: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0028766

+ EVfold
+ Follows the maximum entropy principle
  + Simultaneously want a uniform distribution (maximum entropy) and one that explains the data well
  + The Potts model probability distribution form is actually a solution to this
+ Mutual information between residue pairs don’t work well
  + Can’t resolve transitive correlations (A → B, B → C but A → C is false)
+ Input is the MSA from a Pfam! More coverage the more accurate
+ + ![](https://github.com/HussainAther/potts/blob/main/img/evmutationfitting.png)

## EVMutation: https://marks.hms.harvard.edu/evmutation/

+ Important note about DMS: measured is organism fitness, which may have a nonlinear relationship with protein function!
+ Fit Potts models on sequence family MSAs, infer effect of mutation on 21 different protein DMS experiments
  + Uses L2 regularization: justification is the number of parameters is larger than the number of data examples
+ Training data
  + For each sequence, get an MSA using jackhmmer, searching all of UniRef100
  + Ensure the number of sequences is >=10L
+ Delete positions of MSA with more than 30% gaps, remove sequences that align to less than half of the target sequence

## An evolution-based model for designing chorismate mutase enzymes: https://science.sciencemag.org/content/369/6502/440.abstract

+ Basically uses a Potts model trained on MSA of 1259 natural homologs of CM enzyme (part of Phe and Tyr synthesis pathways)
  + Code: https://github.com/ranganathanlab/bmDCA 
+ Low sequence identities of designs (~25%) but many have appreciable activity, some even exceeds activity of wildtype
+ Probabilities under the model correspond well with function of CM enzyme, below a threshold very few sequences function
  + But just sampling from the PSSM derived from MSA (only first order frequencies are modelled here), though high probability under model, performs poorly.
  + ![](https://github.com/HussainAther/potts/blob/main/img/probabiltiescm.png)
+ Three-residue correlations exist in the Potts model, even though it was only trained on first and second order correlations!
+ + ![](https://github.com/HussainAther/potts/blob/main/img/threeresidue.png)
+ Estimate the size of sequence space that can form functional enzymes
  + e^{entropy}
+ Extensions
  + New experimentally characterized sequences can be directly added to the MSA, ensure diversity because they already have low seq identity
+ Limitations
  + Not denovo--heavily relies on the existence of function already in nature

## Energy-based Out-of-distribution Detection: https://arxiv.org/abs/2010.03759

+ For classification problems, instead of relying on softmax confidence score, use an energy score
  + Don’t need to train an EBM generative model for this!
+ Ultimate goal: reduce false positives, and may provide insight into where to explore
+ Introduced an energy-based classifier by including an energy-based regularization term in the loss
+ + ![](https://github.com/HussainAther/potts/blob/main/img/Lenergy.png)
  + Explicitly penalizes high E for in-distribution and low E for out-of-distribution
  + Improves false positives at TPR95% on all datasets tested without sacrificing classification accuracy


