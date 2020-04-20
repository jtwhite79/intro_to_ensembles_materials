import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle as rect
import flopy
import pyemu

t_d = "template"
nam_file = "xsec.nam"
def plot_domain():
	pyemu.os_utils.run("mfnwt {0}".format(nam_file),cwd=t_d)
	m = flopy.modflow.Modflow.load(nam_file,model_ws=t_d,check=False)
	x = np.arange(m.ncol) + 1
	harr = np.loadtxt(os.path.join(t_d,nam_file.replace(".nam",".hds")))
	fig,axes = plt.subplots(2,1,figsize=(8,4))
	axes[0].plot(x,harr[0,:],marker='.',color='b')
	axes[1].plot(x,harr[1,:],marker='.',color='b')

	axes[0].set_xticks(x)
	axes[1].set_xticks(x)
	axes[0].set_ylim(0,harr.max()*1.2)
	axes[1].set_ylim(0, harr.max()*1.2)

	for i in range(m.ncol):
		t = m.dis.top.array[0,i]
		b = m.dis.botm.array[0,0,i]
		c = "0.5"
		if i == 0:
			c = "g"
		if i == m.ncol -1:
			c = "m"
		x = i + 0.5
		r1 = rect((x, b), 1.0, t - b, facecolor=c, edgecolor="k", alpha=0.5)
		r2 = rect((x, b), 1.0, t - b, facecolor=c, edgecolor="k", alpha=0.5)
		axes[0].add_patch(r1)
		axes[1].add_patch(r2)
		axes[0].text(i+1,(t-b)/2.0,str(i+1),ha="center",va="center")
		axes[1].text(i + 1, (t - b) / 2.0, str(i + 1), ha="center", va="center")

	ox = [4,6]
	oy = harr[0,[3,5]]
	axes[0].scatter(ox,oy,marker="^",color="b",s=55)

	axes[0].annotate("water level obs",
					 xy=(ox[0], oy[0]), xycoords='data',
					 xytext=(5, 3), textcoords='data',
					 size=8, va="center", ha="center",
					 bbox=dict(boxstyle="round4", fc="w"),
					 arrowprops=dict(arrowstyle="-|>",
									 connectionstyle="arc3,rad=0.5",
									 fc="w"))
	axes[0].annotate("water level obs",
					 xy=(ox[1], oy[1]), xycoords='data',
					 xytext=(5, 3), textcoords='data',
					 size=8, va="center", ha="center",
					 bbox=dict(boxstyle="round4", fc="w"),
					 arrowprops=dict(arrowstyle="-|>",
									 connectionstyle="arc3,rad=-0.5",
									 fc="w"))

	ox = 8
	oy = harr[0, 7]
	axes[0].scatter(ox, oy, marker="^", color="r", s=55,zorder=10)
	axes[0].annotate("water level\nforecast",
					 xy=(ox, oy), xycoords='data',
					 xytext=(7, 4), textcoords='data',
					 size=8, va="center", ha="center",
					 bbox=dict(boxstyle="round4", fc="w"),
					 arrowprops=dict(arrowstyle="-|>",
									 connectionstyle="arc3,rad=-0.5",
									 fc="w"))

	ox = 8
	oy = harr[1, 7]
	axes[1].scatter(ox, oy, marker="^", color="r", s=55,zorder=10)
	axes[1].annotate("water level\nforecast",
					 xy=(ox, oy), xycoords='data',
					 xytext=(7, 4.7), textcoords='data',
					 size=8, va="center", ha="center",
					 bbox=dict(boxstyle="round4", fc="w"),
					 arrowprops=dict(arrowstyle="-|>",
									 connectionstyle="arc3,rad=-0.5",
									 fc="w"))


	axes[0].text(0.5,-0.5,"specified\nwater level",ha="left",va="center")
	axes[0].text(10.5, -0.5, "specified\ninflow", ha="right", va="center")
	axes[1].text(0.5, -0.5, "specified\nwater level", ha="left", va="center")
	axes[1].text(10.5, -0.5, "specified\ninflow", ha="right", va="center")

	axes[0].set_ylabel("water level")
	axes[1].set_ylabel("water level")
	axes[0].set_xticks([])
	axes[1].set_xticks([])
	axes[0].set_xlabel("active model cell")
	axes[1].set_xlabel("active model cell")

	axes[0].set_title("A) stress period 1: history",loc="left")
	axes[1].set_title("B) stress period 2: future",loc="left")

	axes[0].annotate("0.5 $\\frac{L^3}{T}$",
                  xy=(10.0, 1.0), xycoords='data',
                  xytext=(8.75, 1.95), textcoords='data',
                  size=8, va="center", ha="center",
                  bbox=dict(boxstyle="round4", fc="w"),
                  arrowprops=dict(arrowstyle="-|>",
                                  connectionstyle="arc3,rad=-0.5",
                                  fc="w"))
	axes[1].annotate("1.0 $\\frac{L^3}{T}$",
					 xy=(10.0, 1.0), xycoords='data',
					 xytext=(8.75, 1.95), textcoords='data',
					 size=8, va="center", ha="center",
					 bbox=dict(boxstyle="round4", fc="w"),
					 arrowprops=dict(arrowstyle="-|>",
									 connectionstyle="arc3,rad=-0.5",
									 fc="w"))
	axes[0].set_xlim(0.5,(m.ncol)+0.5)
	axes[1].set_xlim(0.5, (m.ncol) + 0.5)
	plt.tight_layout()
	plt.show()

	

def run_prior_sweep():
	pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
	pst.pestpp_options["ies_num_reals"] = 100000
	pst.control_data.noptmax = -1
	pr_pst = "pest_prior.pst"
	pst.write(os.path.join(t_d,pr_pst))
	pyemu.os_utils.start_workers(t_d,"pestpp-ies",pr_pst,num_workers=50,master_dir="master_prior")

def run_prior_sweep():
	pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
	pst.control_data.noptmax = -1
	pst.write(os.path.join(t_d,"pest.pst"))
	pyemu.os_utils.start_workers(t_d,"pestpp-glm","pest.pst",num_workers=1,master_dir="master_jcb")

if __name__ == "__main__":
	plot_domain()
	#run_prior_sweep()
