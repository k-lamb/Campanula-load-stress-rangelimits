### All code required to run analyses in 'Genetic drift constrains range expansion into novel habitats'
##### Code banks are specific to the analyses:
######  - genetic load
    - calculate genetic load from heterosis data (./load/all_het_dat_prep.R)
    - determine if cohort effects influence estimation of load/heterosis (./load/cohort_equality_crosses_per_pop.R)
    - bootstrap estimates of genetic load to determine if reproductive isolation is more or less likely to be present (./load/all_het_bootstrap_dat_prep.R)
    - final analysis (./load/het_analysis_final.R)
######  - local adaptation (using relativized fitness data)
    - data preparation (./local_adaptation/*/*_dat_prep.R)
    - running relative models (./local_adaptation/*/*_rel_abs_LA.R)
######  - mapping (generic mapping and spatial analysis)
    - plot maps in Fig. 1a,b with SDM results (./maps/elev_map.R)
    - quantify environmental gradients Fig. 1c (./maps/gradient_quant.R)
    - plot population suitability Fig. 1d (./maps/pop_suitability.R)
    - plot final spatial data used in SDM (./maps/ensemble_datplot.R)
    - map SDMs trimmed for connectivity (./maps/plotWestApp_together.R)
######  - SDM
    - generate pseudo-absences and perform spatial thinning and down-sampling (./SDM/dat_prep_fix.R)
    - test Random Forest and Maximum Entropy modeling strategies (./SDM/models/maxent_model.R & ./SDM/models/rf_model.R)
######  - load x adaptation (and stress)
    - bring together all above data for unified analysis (./adaptation_load/*_het_la_dat_prep.R) 
    - integrate above analyses into a unified modeling framework (./adaptation_load/elev_lat_ensemble.R)
