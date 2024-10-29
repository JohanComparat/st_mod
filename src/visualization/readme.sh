
cd $GIT_STMOD/src/visualization
conda activate stmod
#
#
# DECOMPOSITION ALL = CEN + SAT
#
#
python show_cen_sat_decomp.py
python show_cen_sat_decomp_MS.py
#
#
# DECOMPOSITION CEN = SUM components : hot gas + AGN + XRB
#
#
python show_profile_w_agn_xrb.py
python show_profile_w_agn_xrb_MS.py
#
# Blue Cloud - Red Sequence
python show_profile_w_agn_xrb_BC_RS.py
python show_profile_w_agn_xrb_MS_BC_RS.py
#
# plots the scaling relation (intermediate result to check)
#
python plot_LX_MS.py
python plot_clusters_scaling_relations_M500c_LX.py

