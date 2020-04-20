import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import flopy
import pyemu

t_d = "template"
nam_file = "xsec.nam"
def plot_domain():
	m = flopy.modflow.Modflow.load(nam_file,model_ws=t_d,check=False)
	fig,axes = plt.subplots(2,1,figsize=(8,8))
	

def run_prior_sweep():
	pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
	pst.pestpp_options["ies_num_reals"] = 1000
	pst.control_data.noptmax = -1
	pr_pst = "pest_prior.pst"
	pst.write(os.path.join(t_d,pr_pst))
	pyemu.os_utils.start_workers(t_d,pr_pst,"pestpp-ies",num_workers=20,master_dir="master_prior")

if __name__ == "__main__":
	#plot_domain()
	run_prior_sweep()
