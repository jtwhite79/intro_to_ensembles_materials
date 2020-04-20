import os
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
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
	pst.pestpp_options["ies_num_reals"] = 10000
	pst.control_data.noptmax = -1
	pr_pst = "pest_prior.pst"
	pst.write(os.path.join(t_d,pr_pst))
	pyemu.os_utils.start_workers(t_d,"pestpp-ies",pr_pst,num_workers=50,master_dir="master_prior")

def run_jco():
	pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
	pst.control_data.noptmax = -1
	pst.write(os.path.join(t_d,"pest.pst"))
	pyemu.os_utils.start_workers(t_d,"pestpp-glm","pest.pst",num_workers=1,master_dir="master_jcb")


def jco_vs_en(num_reals=1000):
	pst = pyemu.Pst(os.path.join(t_d,"pest.pst"))
	jco = pyemu.Jco.from_binary(os.path.join("master_jcb","pest.jcb")).to_dataframe()
	# pen = pyemu.ParameterEnsemble.from_csv(pst=pst,filename=os.path.join("master_prior","pest_prior.0.par.csv"))
	# oen = pyemu.ParameterEnsemble.from_csv(pst=pst, filename=os.path.join("master_prior", "pest_prior.0.obs.csv"))
	pen = pd.read_csv(os.path.join("master_prior", "pest_prior.0.par.csv"),index_col=0)
	oen = pd.read_csv(os.path.join("master_prior", "pest_prior.0.obs.csv"),index_col=0)
	obs = pst.nnz_obs_names
	obs.extend([o for o in pst.obs_names if "_08" in o])
	print(obs)
	print(pen.shape)
	with PdfPages("jco_vs_en.pdf") as pdf:
		for par in pen.columns:
			print(par)
			fig,axes = plt.subplots(1,len(obs),figsize=(10,2))
			for ax,o in zip(axes,obs):
				slope = jco.loc[o,par]
				pv = pen.loc[:,par].copy()
				ub = pst.parameter_data.loc[par,"parubnd"]
				lb = pst.parameter_data.loc[par, "parlbnd"]

				pv[pv>=ub] = np.NaN
				pv[pv <= lb] = np.NaN
				pv = pv.dropna()
				pv = pv.apply(np.log10)
				ov = oen.loc[pv.index,o]
				pv = pv.iloc[:num_reals]
				ov = ov.iloc[:num_reals]
				ax.scatter(pv,ov,marker='.',s=1,color="0.5")
				df = pd.DataFrame({"pv": pv, "ov": ov})
				cc = df.corr().iloc[1, 0]
				cc *= ov.std()
				cc /= pv.std()
				xlim,ylim = ax.get_xlim(),ax.get_ylim()
				m,b = np.polyfit(pv,ov,1)
				x = np.linspace(xlim[0],xlim[1],100)
				ax.plot(x,b + (x*slope),"r",lw=2.0,label="jacobian")
				ax.plot(x,b + (x*cc),"k",lw=2.0,dashes=(1,1),label="ensemble")

				print(b, m,cc)
				print(cc)
				ax.set_title("par:{0}, obs:{1}\nsen:{2:4.2E} vs CC:{3:4.2E}".format(par,o,slope,cc))
				ax.set_xlabel("parameter value")
				ax.set_ylabel("observation value")
				ax.legend(loc="upper left")

			plt.tight_layout()
			pdf.savefig()
			plt.close(fig)



def example_jco_vs_en():
	np.random.seed(12356)
	fig,ax = plt.subplots(1,1,figsize=(10,5))
	x = np.linspace(0,10,100)
	slope1 = -0.2
	slope2 = 0.0
	y = 3 + x*slope1
	x_en = x.copy()  # + np.random.normal(0.0,1.0,size=x.size)
	y_en = y + np.random.normal(0.0, 2.1, size=x.size) * 0.5
	ylim = [y_en.min() * 0.5,y_en.max() * 1.3]

	ax.plot(x,y,"r",label="jacobian")
	ax.set_xlabel("parameter value")
	ax.set_ylabel("observation value")
	ax.legend(loc="upper right")
	ax.annotate("finite-diff gradient={0}\n(jacobian entry)".format(slope1),
				xy=(x[0], y[0]), xycoords='data',
				xytext=(3, 4.5), textcoords='data',
				size=8, va="center", ha="center",
				bbox=dict(boxstyle="round4", fc="w"),
				arrowprops=dict(arrowstyle="-|>",
								connectionstyle="arc3,rad=0.5",
								fc="w"))
	ax.set_ylim(ylim)
	plt.savefig("example1.pdf")


	ax.scatter(x_en,y_en,marker=".",color="0.5",s=3)
	ax.annotate("ensemble pairs\n(1 model run each)",
				xy=(x_en[-30], y_en[-30]), xycoords='data',
				xytext=(8, 5.5), textcoords='data',
				size=8, va="center", ha="center",
				bbox=dict(boxstyle="round4", fc="w"),
				arrowprops=dict(arrowstyle="-|>",
								connectionstyle="arc3,rad=0.5",
								fc="w"))
	ax.set_ylim(ylim)
	plt.savefig("example2.pdf")
	m, b = np.polyfit(x_en, y_en, 1)
	y = b + (x * m)
	ax.plot(x, y, "k", label="ensemble", dashes=(1, 1))
	ax.legend(loc="upper right")
	ax.annotate("empirical gradient={0:3.2f}\n(from ensemble)".format(m),
				xy=(x[-10], y[-10]), xycoords='data',
				xytext=(8, 0.5), textcoords='data',
				size=8, va="center", ha="center",
				bbox=dict(boxstyle="round4", fc="w"),
				arrowprops=dict(arrowstyle="-|>",
								connectionstyle="arc3,rad=0.5",
								fc="w"))
	ax.set_ylim(ylim)
	plt.savefig("example3.pdf")
	plt.show()



def set_truth():
	pst = pyemu.Pst(os.path.join(t_d, "pest.pst"))
	jco = pyemu.Jco.from_binary(os.path.join("master_jcb", "pest.jcb")).to_dataframe()
	pen = pyemu.ParameterEnsemble.from_csv(pst=pst,filename=os.path.join("master_prior","pest_prior.0.par.csv"))
	oen = pyemu.ObservationEnsemble.from_csv(pst=pst, filename=os.path.join("master_prior", "pest_prior.0.obs.csv"))
	phi = oen.phi_vector
	phi.sort_values(inplace=True)
	idx = phi.index[-100]
	pst.observation_data.loc[:,"obsval"] = oen.loc[idx,pst.obs_names]
	pst.write(os.path.join(t_d,"pest.pst"))


def run_ies(num_reals,noptmax):
	pst = pyemu.Pst(os.path.join(t_d, "pest.pst"))
	pst.pestpp_options["ies_num_reals"] = num_reals
	pst.control_data.noptmax = noptmax
	pst_name = "pest_ies.pst"
	pst.write(os.path.join(t_d,pst_name))
	m_d = "master_ies_{0}_{1}".format(num_reals,noptmax)
	pyemu.os_utils.start_workers(t_d, "pestpp-ies", pst_name, num_workers=min(num_reals,20),
								 master_dir=m_d)
	return m_d



def plot_ies_result(m_d):


	pst = pyemu.Pst(os.path.join(m_d,"pest_ies.pst"))
	obs = pst.nnz_obs_names
	obs.extend([o for o in pst.obs_names if "_08" in o])
	oe = pd.read_csv(os.path.join(m_d,"pest_ies.{0}.obs.csv".format(pst.control_data.noptmax)),index_col=0)
	fig,axes = plt.subplots(1,len(obs),figsize=(10,4))
	for o,ax in zip(obs,axes):
		ax.set_title(o,loc="left")
		oe.loc[:,o].hist(ax=ax)
		v = pst.observation_data.loc[o,"obsval"]
		ax.plot([v,v],ax.get_ylim(),"r")
	plt.show()


if __name__ == "__main__":
	#plot_domain()
	#run_prior_sweep()
	#run_jco()
	#jco_vs_en()
	#example_jco_vs_en()
	#set_truth()
	m_d = run_ies(4,7)
	plot_ies_result(m_d)

