---
title: "Functional redundancy of Amazonian dung beetles confers community-level resistance to primary forest disturbance"

author: "Cássio Alencar Nunes, Jos Barlow, Filipe França, Erika Berenguer, Ricardo R. C. Solar, Julio Louzada, Rafael Leitão, Laís Maia, Victor H. F. Oliveira, Rodrigo Fagundes Braga, Fernando Z. Vaz-de-Mello, Emma J. Sayer"

date: "April 2021"
output: html_document
---

<style>
body {
text-align: justify}
</style>

**e-mail**: cassioalencarnunes@gmail.com


#### Code repository for the paper: Nunes et al. Functional redundancy of Amazonian dung beetles confers community-level resistance to primary forest disturbance. Biotropica (2021).


In this repository we included codes and data that we used to run the Rarity Index calculation and the simulations of species loss. As we explain in the manuscript all these analyses were adapted from Leitão et al. (2016).

1. You may find the codes to run the Rarity Index calculation in the "rarity_index.html" file.

2. The R script "indices_functions.R" provides the codes of the functions used to calculate Functional Richness (FRic), Functional Specialization (FSpe) and Functional Originality (FOri) with different levels of species extinctions. The functions calculate functional indices in scenarios losing rarest species first, most common species first and losing species randomly at the local scale (see details of Methods below or in the manuscript). Functions were mainly designed by Leitão et al. (2016).

3. The R script "script_simulations.R" provides the codes to run simulations of species loss at regional and local scales for FRic, FSpe and FOri indices.


## Abstract of the manuscript

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Tropical forest biodiversity is being threatened by human activities, and species losses during forest disturbance can compromise important ecosystem functions and services. We assessed how species losses due to tropical forest disturbance affect community functional structure, using Amazonian dung beetles as a model group. We collected empirical data from 106 forest transects and used simulated extinction scenarios to determine how species loss influences community structure at regional and local scales. Although functional and taxonomic community metrics were largely unaffected by primary forest disturbance, they differed markedly between primary and secondary forests. However, our extinction scenarios demonstrated scale-dependence of species losses, whereby functional structure only eroded with species extinction at the local scale. Hence, we extend the spatial insurance hypothesis by demonstrating that landscape-scale functional redundancy offsets the impact of local species losses and confers community-level resistance to primary forest disturbance.


## Part of the Methods section as it is in the manuscript


### 2.2.2. Rarity index

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; We used data from all 161 transects to define species rarity by combining estimates of local abundance, geographical range, and habitat breadth of each dung beetle species. The local abundance of a dung beetle species $_i$ (LA$_i$) was calculated as the mean number of individuals in all transects where the species occurred. The geographical range (GR$_i$) was estimated as the area (ha) inside the smallest polygon joining the outermost sites in which the species occurred using QGIS software (QGIS Development Team 2017). If the species occurred in < 3 transects, the GR$_i$ was considered as the sum of the area within a 1-km radius of the central point of each transect. To estimate the habitat breadth (HB) of a dung beetle species $_i$, we used the “tolerance” metric from Outlying Mean Index analysis (“ade4” package; Dray & Dufour 2007). The Outlying Mean Index is a measure of the species’ niche breadth relative to the niche space of the region, and the tolerance metric describes the spatial variance of the niche across measured environmental conditions or resources (Dolédec et al. 2000). We used tree species richness, aboveground biomass, understorey density, canopy openness, soil texture and elevation (Appendix A1 – section A) to estimate the habitat breadth of each dung beetle species. We calculated the rarity index (RI) for each species using LA, GR and HB, following Leitão et al. (2016). We first log-transformed each metric and then standardized the data by dividing each value by the maximum value across all species, to give values between 0 and 1 for each metric. We also accounted for the degree of dependence among the three metrics by weighting each by its correlation with the other two. 
The rarity index for a species $_i$ (RI$_i$) was calculated by the following formula:

$$RI_i = \frac{\left[ \left( {LA_i \times w_{la}}\right) + \left( {GR_i \times w_{gr}}\right) + \left( {HB_i \times w_{hb}}\right) \right]}{\left( w_{la} + w_{gr} + w_{hb} \right)}$$

where $w_{la}$, $w_{gr}$ and $w_{hb}$ are the weighting parameters for local abundance, geographical range and habitat breadth, respectively. The weighing parameter for each metric $_x$  was calculated by the following formula:

$$w_x = \frac{1}{2} + \left[ \left( \frac{1 - |r_{x1}|}{2} \right) + \left( \frac{1 - |r_{x2}|}{2} \right) \right] $$

where $r_{x1}$ and $r_{x2}$ are the Pearson’s correlation coefficients between the given metric x and each of the other two metrics. Values of RI$_i$ range between 0-1, whereby the rarest species have values close to 0 and the most common species have values close to 1; for clarity of discussion, we thus refer to species with RI < 0.5 as ‘rare’ and to those with RI > 0.5 as ‘common’. For six species that were only collected in 2016, we calculated GR from the 37 transects sampled in that year. Since we were interested in testing hypotheses related to species rarity, we excluded from the analyses seven species that are known not to be attracted to mammal dung and are therefore undersampled by baited pitfall traps: five species of the genus *Anomiopus*, as well as *Bdelyrus paraensis* and *Dendropaemon* aff. *refulgens* (18 individuals in total).


### 2.4. Simulating species loss

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; We simulated scenarios of species loss to assess the consequences of possible extinctions on the functional structure of dung beetle communities. We ran simulations for extinction scenarios at a local (transects) and a regional scale (pool of species), using dung beetle communities from the undisturbed primary forest (12 transects), and assessing the outcomes of the scenarios from the change in the three functional indices.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For the regional simulations, we assessed three scenarios: 1) “rarest first” in which we sequentially removed species from the pool, from the rarest species (lowest RI values) to the most common species (highest RI values), and recalculated the three indices after each removal; 2) “common first” in which we sequentially removed species from the most common to the rarest; 3) “null scenario” in which we randomly removed species from the pool by shuffling the order of species removal 1000 times. We then evaluated the level of functional erosion (the decline in the values of the three functional indices) for each scenario, comparing the outcome of the first two scenarios against the null scenario.

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; The local simulations also included the same three scenarios, but we calculated FSpe and FOri indices for nine levels of species loss (from 10% to 90%), and we calculated FRic for seven levels (10% to 70%), to ensure there were sufficient species per transect to calculate the convex hull volume of the functional space in all scenarios (four species minimum; Mouillot, Graham, et al. 2013). We compared the values of functional indices resulting from the three scenarios using Friedman paired tests.


### References

DOLÉDEC, S., D. CHESSEL, and C. GIMARET-CARPENTIER. 2000. Niche separation in community analysis: a new method. Ecology 81: 2914–2927.

DRAY, S., and A.-B. DUFOUR. 2007. The ade4 Package: Implementing the Duality Diagram for Ecologists. J. Stat. Softw. 22. Available at: http://www.jstatsoft.org/v22/i04/.

LEITÃO, R. P., J. ZUANON, S. VILLÉGER, S. E. WILLIAMS, C. BARALOTO, C. FORTUNEL, F. P. MENDONÇA, and D. MOUILLOT. 2016. Rare species contribute disproportionately to the functional structure of species assemblages. Proc. R. Soc. B Biol. Sci. 283: 20160084.

MOUILLOT, D., N. A. J. GRAHAM, S. VILLÉGER, N. W. H. MASON, and D. R. BELLWOOD. 2013. A functional approach reveals community responses to disturbances. Trends Ecol. Evol. 28: 167–177.

QGIS DEVELOPMENT TEAM. 2017. QGIS Geographic Information System. Open Source Geospatial Foundation Project. Available at: http://qgis.org.

