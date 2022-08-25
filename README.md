# RDK_vMMM

### Random dot kinematogram data analysis using a von Mises mixture model

This code belongs to the paper "Psychophysics and computational modeling of feature-continuous motion perception" by [Töpfer, Barbieri](https://twitter.com/DecisionGuy) et al. (2022), accepted for publication in the *Journal of Vision*. It consists of MATLAB scripts processing behavioral data, also available in this repository, and generating results and figures, as they appear in the paper.

- Paper: TBA
- Data: https://github.com/JoramSoch/RDK_vMMM/tree/master/data
- Toolbox: https://github.com/JoramSoch/RDK_vMMM/tree/master/tools/vMMM


## Requirements

This code was written by [Felix Töpfer](https://www.researchgate.net/profile/Felix-Toepfer) (main folder) and [Joram Soch](https://orcid.org/0000-0002-8879-5666) (sub-folder `tools/vMMM/`).

The data analysis was developed and run in [MATLAB R2019b](https://de.mathworks.com/help/matlab/release-notes.html). No further toolboxes are required.


## Instructions

To re-analyze the behavioral data, proceed as follows:
1. Clone the repository to some folder on your computer.
2. Move to that folder and run `BEHAVIOR_ANALYSIS.m`.
3. Modify [line 21](https://github.com/JoramSoch/RDK_vMMM/blob/master/BEHAVIOR_ANALYSIS.m#L21) (see [line 22](https://github.com/JoramSoch/RDK_vMMM/blob/master/BEHAVIOR_ANALYSIS.m#L22)) and re-run the script.

This should generate Figures 3/4, 5A/B, 6/7 and 9-16 from the paper, depending on whether you have chosen the response method to be `meth = 'trackball'` (trackball, see step 2) or `meth = 'bar'` (rotating bar, see step 3). Note that the panels of Figure 3/4 appear in separate windows, as each panel consists of several subplots. These windows are labeled as "Scatterplot of trial-wise responses (`stimulus-type`: `coherence-level`)" where `stimulus_type` is either "TM" (transparent motion), "BM" (Brownian motion) or "WM" (white noise motion) and `coherence-level` varies between 0, 12.5, 25, 50 and 100%.
