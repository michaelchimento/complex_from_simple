# Data

Please find below a short description of each Rda used for analysis.

df_solves_cleaned.Rda - this contains all solving activity at the puzzleboxes for all 3 stages of the experiment.

| variable name       | description     |
| :------------- | ----------: |
| ABSDATE | integer day of experiment   |
| UNIT   | puzzlebox ID |
| SITE  | patch within subpopulation |
| ACTION | component (e.g. dial_r) performed at puzzlebox. Two consecutive components make up a complex solve. |
| ID  | PIT code ID |
| RING  | ring number |
| SPECIES  | species code |
| SEX  | M or F when known |
| YRBORN  | year born when known |
| AGE  | age of bird when recorded, when known |
| DEMO  | binary variable of whether bird was a trained demonstrator |
| PATCH  | Subpopulation |
| time_stamp  | full posix timestamp |
| experiment  | experimental stage |
| solution  | solution recorded |
| solution category  | whether solution is dial, slide or complex |
| TTS  | time-to-solve (seconds) |
| ind_solve_index  | running tally of solutions by bird |
| ind_solve_index_byexp  | running tally of solutions by bird within experiment |
| ind_solve_index_bytype  | running tally of solutions by bird within solution |
| ind_solve_index_bytype_byexp  | running tally of solutions by bird within solution and experiment |
| ind_solve_index_bycat  | running tally of solutions by bird within solution category |
| ind_solve_index_bycat_byexp  | running tally of solutions by bird within solution category and experiment |
| total_solutions_byexp  | total nr of solutions produced by bird within experiment |
| total_solutions  | total nr of solutions produced by bird over all 3 stages |

df_visitors_all.Rda - record of all individual great tits who visited the puzzlebox

df_visitors_nonsolvers.Rda - record of all individual great tits who visited the puzzlebox but never produced any total_solutions

df_sim_{ind, site, subpop}.Rda - simulated data for figure 5

NBDA_input_data.Rda - lists of input data needed for NBDA analysis by SW.
