Does a structure predictor notice a mutation?

CAIR Collective, UCLA — Structural Biology and Molecular Biophysics


QUESTION

AlphaFold2 returns a near wild-type structure when a single residue is
substituted, including for substitutions known to destabilize the protein. The
explanation usually offered is that the model reads a multiple sequence
alignment of thousands of homologs, so one changed position is outweighed by the
alignment around it.

That explanation is testable. ESMFold folds from a single sequence through a
protein language model and uses no alignment at all. If alignment depth is the
reason for the insensitivity, a predictor without an alignment should respond
where an alignment-based one does not.

We ask whether it does, and whether any response it produces is physically
meaningful or numerical noise.


BACKGROUND

The insensitivity itself is established. Published analyses report that
AlphaFold2 cannot reliably estimate the structural or functional consequence of
a point mutation, that RMSD and confidence scores are discordant for
destabilizing variants, and that changes in confidence show no correlation with
measured fluorescence across a large set of GFP mutants. A dissenting preprint
argues that mutation-induced change can be recovered, but only by manipulating
the alignment, which concedes that the alignment is the controlling variable.

What is missing from that literature is the control. These studies compare a
predicted wild type to a predicted mutant without establishing what the
predictor does when nothing is changed at all, so the effect sizes they report
have no baseline to be measured against. None of them tests the alignment
explanation directly by removing the alignment.


HYPOTHESES

H1. Alignment dominance explains the insensitivity. A single-sequence predictor
deforms measurably more than an alignment-based one for the same substitution.

H0. The insensitivity is intrinsic to structure prediction from sequence.
Neither predictor exceeds its own run-to-run variation.

H2. Any motion produced is physically meaningful. Displacement concentrates
near the substituted residue and scales with measured fitness cost.

H2-null. The motion is numerical noise, delocalized across the chain and
unrelated to fitness.

A null result is a result, provided the reference frames below are measured.


REFERENCE FRAMES

A displacement between a predicted wild type and a predicted mutant means
nothing in isolation. Four distances are needed, and each answers a different
question.

Replicate folds of the same sequence give the run-to-run noise floor: what the
predictor does when nothing has changed.

A predicted wild type against its experimental structure gives the accuracy
floor: how far the model sits from reality before any mutation is introduced.

A wild type coordinate set carrying the substituted side chain with no
repacking gives the trivial mutant: the structure a model would return if it
registered the substitution but did not respond to it.

Experimental wild type and mutant structures, where matched pairs exist, give
the ground truth deformation: how far a real protein actually moves. This frame
matters because the premise that a model should move is itself an assumption.
Most point substitutions perturb side chain packing without displacing a
backbone, and if real structures also barely move, then insensitivity is in part
correct behavior rather than failure. Matched experimental pairs carry a
survivorship bias that must be stated: variants severe enough to prevent folding
are not represented, because they were never crystallized.


APPROACH

Panel. Small, well-characterized proteins with published deep mutational
scanning data and experimental reference structures: ubiquitin, the PSD95 PDZ3
domain, avGFP, and TEM-1 beta-lactamase. Short chains, dense variant coverage,
and folds that are recognizable on sight.

Sampling. Substitutions are selected stratified on measured fitness and on
burial, rather than at random. Random sampling supports reporting an average and
little else; stratification is what makes the dependence on structural context
answerable, and it cannot be retrofitted after selection.

Prediction. Wild type and mutant sequences are folded with ESMFold and with
ColabFold, with replicates of every wild type.

Measurement. Structures are superposed and compared along four axes.

Per-residue displacement locates where a model moved, which global RMSD cannot.
Global RMSD is reported for comparability with prior work, not as the finding.

Contact maps record which residue pairs lie within a distance cutoff. Packing
rearrangement can leave RMSD nearly unchanged while altering contacts, so a
substitution that opens a cavity should break contacts around the site. A model
that swaps a residue without repacking leaves the contact map intact, which is
the behavior under examination.

Local solvent accessibility measures burial and its response. Accessibility at
the substituted position is dominated by the amino acid swap itself rather than
by any structural response, so change is measured over the surrounding shell and
against the trivial mutant instead.

Confidence change is the model's own estimate of local accuracy, and is the one
signal prior work reports as responsive. It is a statement about certainty, not
about energy, and a decrease does not evidence destabilization. Confidence is
compared within a predictor and across predictors only as a difference, since
the two scores come from differently trained heads.


EVALUATION

Results are reported per protein and never pooled across proteins, paired by
substitution when predictors are compared, with bootstrap confidence intervals.
Every claimed effect is reported against the noise floor and the accuracy floor.


VERIFICATION

Residue numbering in deep mutational scanning tables, in UniProt, and in
deposited structures disagree routinely, usually by a constant offset arising
from a cleaved signal peptide or an expression tag. The mismatch raises no
error. It yields a complete set of plausible structures answering a question
that was never asked. Every mutant sequence is therefore verified in code to
carry the intended residue at the intended position, with the wild-type residue
matching its source table, before anything is folded.


EXTENSION

If alignment depth governs the insensitivity, then progressively subsampling the
alignment should restore sensitivity in proportion to depth. That experiment
converts a comparison between two tools into a statement about why either
behaves as it does, and the work described here is what makes it possible.


REPOSITORY

notebooks/   one lab notebook per member
data/        raw and derived data; contents are not tracked
src/         shared code: variant handling and structural geometry

Findings, row counts, numbering offsets and the reasons behind decisions belong
in data/README.md, recorded alongside the work that produced them.
