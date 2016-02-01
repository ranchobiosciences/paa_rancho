----------------------------------------------------------------------------
PAA (Protein Array Analyzer) package with modifications by Rancho BioSciences
----------------------------------------------------------------------------
       John Obenauer, Ph.D., and Ivan Grishagin, Ph.D., February 2016

The Bioconductor package PAA, created by Michael Turewicz at Ruhr University Bochum,
provides data analysis functions for protein microarrays.  In our analysis of Thermo Fisher
ProtoArrays at Rancho BioSciences, we needed functions for robust linear model (RLM)
normalization and M-statistic calculations for determining differentially expressed proteins.
These two methods are recommended by Mark Gerstein's lab in their paper introducing the
ProtoArray platform:

    Sboner A, Karpikov A, Chen G, Smith M, Dawn M, Freeman-Cook L, Schweitzer B, Gerstein MB.
    Robust-linear-model normalization to reduce technical variability in functional protein
    microarrays.  Journal of Proteome Research 8:5451-5464, 2009.

Gerstein's data is publicly available and is provided as test data for our version of PAA.
They used a mixture of six positive control serum samples (PCS) containing antibodies against
the antigens Ro/SS-A, La/SS-B, Smith antigen, RNP complex, DNA topoisomerase I (Scl-70), and
Jo-1.  They made dilutions of this positive control mixture using one normal serum sample with
no known pathology.  In the sample information file we prepared for this data set, we used the
names PCS000, PCS025, PCS050, PCS075, and PCS100 for the samples containing 0%, 25%, 50%, 75%,
or 100% of the positive control sample mixture.  These dilutions were used in the paper cited
above to show the improved correlations obtained with RLM normalization compared to raw data
or data normalized other ways.

The files we provide are described as follows:

PAA_Rancho_v21.R -- This is our version of the PAA package, containing all the functions needed
for analyzing Thermo Fisher ProtoArrays.  The functions should also work on other protein arrays
but have not been tested on any others.

example.R -- This is an example R script showing how to call the PAA functions for a typical
analysis and find differentially expressed proteins.

README -- This file, explaining what the other files contain.

Users can obtain the Rancho version of PAA from our GitHub page:

https://github.com/ranchobiosciences/paa_rancho

We are also submitting these new functions for inclusion in future versions of the PAA
package in Bioconductor in order to make the functions easier to obtain.  Rancho's modified 
PAA code is available under the open source BSD license in order to be compatible with the 
parent PAA package.

Questions about the Rancho-modified PAA package should be directed to John Obenauer
(john@ranchobiosciences.com) or Ivan Grishagin (ivan@ranchobiosciences.com).

