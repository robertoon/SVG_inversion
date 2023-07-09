import os
import shutil
import numpy as np
import pandas as pd
import pyemu


#<-- change this to where the t2 bin is on your system...
path_2_t2 = os.path.join("..","tephra2/tephra2_2020")

def setup_interface(org_dir,rg_file,out_file,wind_file_2_run):
	"""setup a pest interface for the tephra2 files in 
	`org_dir`.  The model+pest files are setup in the 
	`org_dir`+"_pest" dir.  The wind files are also included.
	
	Args:
		org_dir : original tephra2 input files
		rg_file : the file in `org_dir` that contains the ranges for ESPs to estimate
		out_file : the text output file from tephra2 that has the observed tephra thickness
		wind_file_2_run: the wind db file to use in the tephra2 sim.  This is only a temporary fix
			since we arent estimating the wind profile yet...

	Note:
		`wind_db` is assumed to exist in `org_dir`.  Wind parameters are setup for the files found in `wind_db`
		but these parameters are setup to write wind files at the upper directory level (`wind_db` is removed from 
		the new dir)

	"""

	rg_file = os.path.join(org_dir,rg_file)
	out_file = os.path.join(org_dir,out_file)
	assert os.path.exists(rg_file)
	assert os.path.exists(out_file)

	new_dir = org_dir + "_pest"

	if os.path.exists(new_dir):
		shutil.rmtree(new_dir)
	shutil.copytree(org_dir,new_dir)

	w_dir = os.path.join(new_dir,"wind_db")
	wind_files = os.listdir(w_dir)
	assert wind_file_2_run in wind_files

	cmd_str = "./tephra2_2020 tephra2.in grid.txt {0} >{1}".format(wind_file_2_run,os.path.split(out_file)[-1]) 
	with open(os.path.join(new_dir,"forward_run.py"),'w') as f:
		f.write("import os\n")
		f.write("os.system('{0}')\n".format(cmd_str))
	shutil.copy2(path_2_t2,os.path.join(new_dir,os.path.split(path_2_t2)[-1]))
	#pyemu.os_utils.run(cmd_str,cwd=new_dir)

	f_org = open(rg_file)

	#first the T2 config
	in_file = os.path.join(new_dir,'tephra2.in')
	tpl_file = os.path.join(new_dir,'tephra2.in.tpl')

	f_in = open(in_file,'w')
	f_tpl = open(tpl_file,'w')
	f_tpl.write("ptf ~\n")
	pdata = {}
	for line in f_org:
		if len(line.strip()) == 0 or line.strip()[0] == "#":
			f_tpl.write(line)
			f_in.write(line)
		else:
			raw = line.strip().split()
			if len(raw) == 2:
				f_tpl.write(line)
				f_in.write(line)
			else:
				pname = raw[0].strip().lower()
				lb = float(raw[1])
				ub = float(raw[2])
				line = raw[0] + " ~  " + pname + "   ~\n"
				f_tpl.write(line)
				if "median" in pname:
					mn = (ub+lb)/2.0
				else:
					mn = 10**((np.log10(ub)+np.log10(lb))/2.0)
				pdata[pname] = [mn,lb,ub]
				f_in.write("{0}  {1}\n".format(raw[0],(ub+lb)/2.0))
	f_org.close()
	f_tpl.close()
	f_in.close()

	tpl_files = [tpl_file]
	in_files = [in_file]

	# now the wind
	
	for w_file in os.listdir(w_dir):
		shutil.copy2(os.path.join(w_dir,w_file),os.path.join(new_dir,w_file))
		w_file = os.path.join(new_dir,w_file)
		tpl_file = w_file+".tpl"
		f_in = open(w_file,'r')
		f_tpl = open(tpl_file,'w')
		f_tpl.write("ptf ~\n")
		for i,line in enumerate(f_in):
			raw = line.strip().split()
			h = float(raw[0])
			d = float(raw[1])
			v = float(raw[2])
			dname = os.path.split(w_file)[1] + "_height:{0:4.1f}_dir".format(h)
			pdata[dname] = [d,d-100,d+100]
			vname = os.path.split(w_file)[1] + "_height:{0:4.1f}_vel".format(h)
			pdata[vname] = [v,v*0.5,v*1.5]
			f_tpl.write("{0}  ~  {1}   ~  ~   {2}  ~\n".format(raw[0],dname,vname))
			

		f_tpl.close()
		f_in.close()
		tpl_files.append(tpl_file)
		in_files.append(tpl_file.replace(".tpl",""))

	shutil.rmtree(w_dir)

	# now the obs
	odata = {}
	out_file = os.path.join(new_dir,os.path.split(out_file)[1])
	ins_file = out_file+".ins"
	f_out = open(out_file,'r')
	f_ins = open(ins_file,'w')
	f_ins.write("pif ~\nl1\n")
	f_out.readline() #header
	for line in f_out:
		raw = line.strip().split()
		oname = "mass_east:{0}_north:{1}".format(raw[0],raw[1])
		f_ins.write("l1 w w w !{0}!\n".format(oname))
		odata[oname] = float(raw[-1])
	f_out.close()
	f_ins.close()

	pst = pyemu.Pst.from_io_files(tpl_files,in_files,[ins_file],[out_file],pst_path=".")
	pst.model_command = "python forward_run.py"
	par = pst.parameter_data
	par.loc[:,"parval1"] = [pdata[p][0] for p in pst.par_names]
	par.loc[:,"parlbnd"] = [pdata[p][1] for p in pst.par_names]
	par.loc[:,"parubnd"] = [pdata[p][2] for p in pst.par_names]
	par.loc[:,"pargp"] = par.parnme
	par.loc[par.parnme.apply(lambda x: "_vel" in x or "_dir" in x),"partrans"] = "fixed"
	par.loc[par.parnme.apply(lambda x: wind_file_2_run in x),"partrans"] = "none"
	par.loc[par.parnme.apply(lambda x: "_vel" in x or "_dir" in x),"pargp"] = par.loc[par.parnme.apply(lambda x: "_vel" in x or "_dir" in x),"parnme"].apply(lambda x: x.split("_")[0])
	
	par.loc[par.parnme.str.contains(wind_file_2_run),"partrans"] = "none"
	par.loc["median_grainsize","partrans"] = "none"
	par.loc["median_grainsize","parval1"] = 1


	obs = pst.observation_data
	obs.loc[:,"standard_deviation"] = obs.obsval * 0.05
	obs.loc[obs.standard_deviation<0.1,"standard_deviation"] = 0.1
	obs.loc[:,"lower_bound"] = 0.0
	obs.loc[:,"weight"] = 1.0/(obs.obsval * 0.1)
	obs.loc[obs.weight>10.0,"weight"] = 10.0
	obs.loc[obs.weight<0.1,"weight"] = 0.1
	
	
	print(obs.obsval)
	pst.control_data.noptmax = 0
	pst.write(os.path.join(new_dir,"pest.pst"),version=2)

	pyemu.os_utils.run("pestpp-glm pest.pst", cwd=new_dir) ######

########
def run_pestpp(master_dir,tool="pestpp-ies",num_reals=100,noptmax=5):
    """run pestpp-glm in parallel locally"""
    pst = pyemu.Pst(os.path.join("U2_pest", "pest.pst"))
    pst.control_data.noptmax = noptmax
    pst.pestpp_options["ies_num_reals"] = num_reals
    pst.pestpp_options["ies_bad_phi_sigma"] = 1.25
    pst.pestpp_options["panther_agent_freeze_on_fail"] = True
    pst.pestpp_options["par_sigma_range"] = 6
    pst.pestpp_options["ies_multimodal_alpha"] = 0.05
    pst.write(os.path.join("U2_pest", "pest_run.pst"),version=2)
    master_dir += "_{0}reals".format(num_reals)
    pyemu.os_utils.start_workers("U2_pest", tool, "pest_run.pst", num_workers=10,
                                 master_dir=master_dir,worker_root=".")
    return master_dir
#########


# def plot_results(master_dir):

# 	pst = pyemu.Pst(os.path.join(master_dir,"pest_run.pst"))
# 	opr = pd.read_csv(os.path.join(master_dir,"pest_run.0.obs.csv"),index_col=0)
# 	opt = pd.read_csv(os.path.join(master_dir,"pest_run.{0}.obs.csv".format(pst.control_data.noptmax)),index_col=0)


# 	pr = pd.read_csv(os.path.join(master_dir,"pest_run.0.par.csv"),index_col=0)
# 	pt = pd.read_csv(os.path.join(master_dir,"pest_run.{0}.par.csv".format(pst.control_data.noptmax)),index_col=0)
# 	import matplotlib.pyplot as plt
# 	from matplotlib.backends.backend_pdf import PdfPages
# 	with PdfPages(os.path.join(master_dir,"results.pdf")) as pdf:

# 		fig,ax = plt.subplots(1,1,figsize=(5,5))
# 		ax.hist(np.log10(pyemu.ObservationEnsemble(pst,opr).phi_vector),bins=20,facecolor="0.5",density=True,alpha=0.5)
# 		ax.hist(np.log10(pyemu.ObservationEnsemble(pst,opt).phi_vector),bins=20,facecolor="b",density=True,alpha=0.5)
# 		ax.set_yticks([])
# 		ax.set_title("phi",loc="left")
# 		plt.tight_layout()
# 		pdf.savefig()
# 		plt.close(fig)


# 		fig,ax = plt.subplots(1,1,figsize=(5,5))
# 		obs = pst.observation_data.loc[pst.nnz_obs_names,:]
# 		for oname,oval in zip(obs.obsnme,obs.obsval):
# 			ax.scatter([oval]*opr.shape[0],opr.loc[:,oname],marker=".",c="0.5",alpha=0.5)
# 			ax.scatter([oval]*opt.shape[0],opt.loc[:,oname],marker=".",c="b",alpha=0.5)
# 		# mn = min(ax.get_xlim()[0],ax.get_ylim()[0])
# 		# mx = max(ax.get_xlim()[1],ax.get_ylim()[1])
# 		mn = opt.min().min()
# 		mx = opt.max().max()
# 		ax.plot([mn,mx],[mn,mx],"k--",lw=2.5)
# 		ax.set_xlim(mn,mx)
# 		ax.set_ylim(mn,mx)
# 		plt.tight_layout()
# 		ax.grid()
# 		pdf.savefig()
# 		plt.close(fig)

# 		for oname,oval,weight in zip(obs.obsnme,obs.obsval,obs.weight):
# 			fig,ax = plt.subplots(1,1,figsize=(5,5))
# 			#mn = min(opr.loc[:,oname].values.min(),opt.loc[:,oname].values.min())
# 			#mx = max(opr.loc[:,oname].max(),opt.loc[:,oname].max())
# 			#bins = np.linspace(mn,mx,100)
# 			bins = 20
			
# 			ax.hist(opr.loc[:,oname].values,bins=bins,density=True,facecolor="0.5",alpha=0.5)
# 			ax.hist(opt.loc[:,oname].values,bins=bins,density=True,facecolor="b",alpha=0.5)
# 			ylim = ax.get_ylim()
# 			ax.plot([oval,oval],ylim,"r-",lw=2.0)
# 			#ax.set_xlim(0,oval*3)
# 			ax.set_title("{0}, weight:{1:3.1f}".format(oname,weight),loc="left")
# 			plt.tight_layout()
# 			pdf.savefig()
# 			plt.close(fig)


# 		for pname in pst.adj_par_names:

# 			partrans = pst.parameter_data.loc[pname,"partrans"]
# 			fig,ax = plt.subplots(1,1,figsize=(5,5))
# 			vals = pr.loc[:,pname].values.copy()
# 			if partrans == "log":
# 				vals = np.log10(vals)
# 				ax.set_title("$log_{10}$ " + pname,loc="left")
				
# 			else:
# 				ax.set_title(pname,loc="left")
# 			ax.hist(vals,bins=20,density=True,facecolor="0.5",alpha=0.5)
# 			vals = pt.loc[:,pname].values.copy()
# 			if partrans == "log":
# 				vals = np.log10(vals)
				
# 			ax.hist(vals,bins=20,density=True,facecolor="b",alpha=0.5)
# 			ax.set_yticks([])
			
			
# 			plt.tight_layout()
# 			pdf.savefig()
# 			plt.close(fig)


if __name__ == "__main__":
	

	org_dir = "U2"
	rg_file = "inversion_ranges_u2.txt"
	out_file = "svg_u2_data.txt"

	# setup_interface(org_dir,rg_file,out_file,"wind0")
	# m_d = run_pestpp("master_dir_wind0",num_reals=1000,noptmax=10)
	#plot_results(m_d)

	# setup_interface(org_dir,rg_file,out_file,"wind1")
	# m_d = run_pestpp("master_dir_wind1",num_reals=1000,noptmax=10)
	# plot_results(m_d)

	# setup_interface(org_dir,rg_file,out_file,"wind2")
	# m_d = run_pestpp("master_dir_wind2",num_reals=1000,noptmax=10)
	#plot_results(m_d)

	setup_interface(org_dir,rg_file,out_file,"wind3")
	m_d = run_pestpp("master_dir_wind3",num_reals=1000,noptmax=10)
 	#plot_results(m_d)

	# setup_interface(org_dir,rg_file,out_file,"wind4")
	# m_d = run_pestpp("master_dir_wind4",num_reals=1000,noptmax=10)
	#plot_results(m_d)

	# setup_interface(org_dir,rg_file,out_file,"wind5")
	# m_d = run_pestpp("master_dir_wind5",num_reals=1000,noptmax=10)
	# plot_results(m_d)
	
	# setup_interface(org_dir,rg_file,out_file,"wind6")
	# m_d = run_pestpp("master_dir_wind6",num_reals=1000,noptmax=10)
	# plot_results(m_d)

	# setup_interface(org_dir,rg_file,out_file,"wind7")
	# m_d = run_pestpp("master_dir_wind7",num_reals=1000,noptmax=10)
	# plot_results(m_d)

	# setup_interface(org_dir,rg_file,out_file,"wind8")
	# m_d = run_pestpp("master_dir_wind8",num_reals=1000,noptmax=10)
	# plot_results(m_d)

	# setup_interface(org_dir,rg_file,out_file,"wind9")
	# m_d = run_pestpp("master_dir_wind9",num_reals=1000,noptmax=10)
	# plot_results(m_d)

	# setup_interface(org_dir,rg_file,out_file,"wind10")
	# m_d = run_pestpp("master_dir_wind10",num_reals=1000,noptmax=10)
	# plot_results(m_d)

	# setup_interface(org_dir,rg_file,out_file,"wind11")
	# m_d = run_pestpp("master_dir_wind11",num_reals=1000,noptmax=10)
	# plot_results(m_d)

	# setup_interface(org_dir,rg_file,out_file,"wind12")
	# m_d = run_pestpp("master_dir_wind12",num_reals=1000,noptmax=10)
	# plot_results(m_d)

	# setup_interface(org_dir,rg_file,out_file,"wind13")
	# m_d = run_pestpp("master_dir_wind13",num_reals=1000,noptmax=10)
	# plot_results(m_d)

	# setup_interface(org_dir,rg_file,out_file,"wind14")
	# m_d = run_pestpp("master_dir_wind14",num_reals=1000,noptmax=10)
	# plot_results(m_d)

	# setup_interface(org_dir,rg_file,out_file,"wind15")
	# m_d = run_pestpp("master_dir_wind15",num_reals=1000,noptmax=10)
	# plot_results(m_d)