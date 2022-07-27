import numpy as np
import pandas as pd
import os
import math
class convert:
	def __init__(self, rawfilename):
		self.rawfilename = rawfilename
		self.trktohist()

	def isfloat(self, val):
		try:
			float(val)
			return True
		except ValueError:
			return False

	def rightname(self, val):
		retval = True
		if self.isfloat(val)==False:
			if val[0].isnumeric()==False:
				retval = False
		return retval
		
	def eformat(self, f, prec, exp_digits):
		s = "%.*e"%(prec, f)
		mantissa, exp = s.split('e')
		# add 1 to digits as 1 is taken by sign +/-
		return "%se%+0*d"%(mantissa, exp_digits+1, int(exp))

	def trktohist(self):
		filename = self.rawfilename[:-4]
		file1=open(self.rawfilename,'r')
		lines = file1.readlines()
		Neutrino=False
		for line in lines:
			if line[1:9]=='Neutrino':
				Neutrino=True

		if Neutrino==False:
			#df = pd.read_table(self.rawfilename,sep='\s+',skiprows=13,names=['Model #', 'shells', 'AGE', 'log(L/Lsun)', 'log(R/Rsun)', 'log(g)', 'log(Teff)', 'Mconv. core', 'Mconv. env.', 'Rconv. env.','M He core', 'Xenv', 'Zenv'])
		
			namearray = ['Model #', 'shells', 'AGE', 'log(L/Lsun)', 'log(R/Rsun)', 'log(g)', 'log(Teff)', 'Mconv. core', 'Mconv. env.', 'Rconv. env.','M He core', 'Xenv', 'Zenv',
			'Luminosity: ppI', 'Luminosity: ppII', 'Luminosity: ppIII', 'Luminosity: CNO', 'Luminosity: triple-alpha', 'Luminosity: He-C', 'gravity', 'neutrinos (old)', '% Grav. eng.', 'Itot', 
			'Central: log(T)', 'Central: log(RHO)', 'Central: log(P)', 'Central: BETA', 'Central: ETA', 'Central: X', 'Central: Z', 'Central: H', 'Central: H shell midpoint', 'Central: H shell mass', 'Central: T and rho at base of c.z.',
			 'Central Abundances: He3', 'Central Abundances: C12', 'Central Abundances: C13', 'Central Abundances: N14', 'Central Abundances: N15', 'Central Abundances: O16', 'Central Abundances: O17', 'Central Abundances: O18',
			 'Surface Abundances: He3', 'Surface Abundances: C12', 'Surface Abundances: C13', 'Surface Abundances: N14', 'Surface Abundances: N15', 'Surface Abundances: O16', 'Surface Abundances: O17','Surface Abundances: O18']

			#clean bad value
			#cols = df.columns.values.tolist()
			#df[cols] = df[cols].apply(pd.to_numeric, errors='coerce', axis=1)
			#df.replace('NaN',0)
			arraysize = (len(lines[13:])//5)*len(namearray)
			xss = np.char.split(lines[13:])
			data = np.nan_to_num(pd.to_numeric(np.array([x for xs in xss for x in xs]),errors='coerce'))
			dataarray = data[:arraysize].reshape((len(lines[13:])//5, len(namearray)))

			dp = pd.DataFrame(dataarray,columns=namearray)

		else:
			#df = pd.read_table(self.rawfilename,sep='\s+',skiprows=14,names=['Model #', 'shells', 'AGE', 'log(L/Lsun)', 'log(R/Rsun)', 'log(g)', 'log(Teff)', 'Mconv. core', 'Mconv. env.', 'Rconv. env.','M He core', 'Xenv', 'Zenv'])
		
			namearray = ['Model #', 'shells', 'AGE', 'log(L/Lsun)', 'log(R/Rsun)', 'log(g)', 'log(Teff)', 'Mconv. core', 'Mconv. env.', 'Rconv. env.','M He core', 'Xenv', 'Zenv',
			'Luminosity: ppI', 'Luminosity: ppII', 'Luminosity: ppIII', 'Luminosity: CNO', 'Luminosity: triple-alpha', 'Luminosity: He-C', 'gravity', 'neutrinos (old)', '% Grav. eng.', 'Itot', 
			'Central: log(T)', 'Central: log(RHO)', 'Central: log(P)', 'Central: BETA', 'Central: ETA', 'Central: X', 'Central: Z', 'Central: H', 'Central: H shell midpoint', 'Central: H shell mass', 'Central: T and rho at base of c.z.',
			 'Central Abundances: He3', 'Central Abundances: C12', 'Central Abundances: C13', 'Central Abundances: N14', 'Central Abundances: N15', 'Central Abundances: O16', 'Central Abundances: O17', 'Central Abundances: O18',
			 'Surface Abundances: He3', 'Surface Abundances: C12', 'Surface Abundances: C13', 'Surface Abundances: N14', 'Surface Abundances: N15', 'Surface Abundances: O16', 'Surface Abundances: O17','Surface Abundances: O18',
			 'Neutrinos: pp', 'Neutrinos: pep', 'Neutrinos: hep', 'Neutrinos: Be7', 'Neutrinos: B8', 'Neutrinos: N13', 'Neutrinos: O15', 'Neutrinos: F17','Neutrinos: Cl37 flux', 'Neutrinos: Ga71 flux']
			
			#clean bad value
			arraysize = (len(lines[14:])//6)*len(namearray)
			xss = np.char.split(lines[14:])
			data = np.nan_to_num(pd.to_numeric(np.array([x for xs in xss for x in xs]),errors='coerce'))
			dataarray = data[:arraysize].reshape((len(lines[14:])//6, len(namearray)))

			dp = pd.DataFrame(dataarray,columns=namearray)

		#idxloc = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
		idxloc = np.arange(1,23).tolist()
		hdrloc= ['star_age','star_mass','star_mdot','he_core_mass','c_core_mass','log_L','log_LH','log_LHe','log_Teff','log_R','log_g','surface_h1','surface_he3','surface_he4','surface_c12','surface_o16','log_center_T','log_center_Rho','center_gamma','center_h1','center_he4','center_c12']
		lendp = len(dp)
		listzip = [[]*22]*lendp
		for i in range(lendp):
			listzip[i] = [(dp['AGE'].values[i])*1e9,int(self.rawfilename[1:5])/1000,0,dp['M He core'].values[i],0, dp['log(L/Lsun)'].values[i],np.log10((dp['Luminosity: ppI'].values[i] + dp['Luminosity: ppII'].values[i] + dp['Luminosity: ppIII'].values[i] + dp['Luminosity: CNO'].values[i])/(3.839e33) + 1e-10), np.log10((dp['Luminosity: triple-alpha'].values[i])/(3.839e33) + 1e-10),dp['log(Teff)'].values[i],dp['log(R/Rsun)'].values[i],dp['log(g)'].values[i],dp['Xenv'].values[i],dp['Surface Abundances: He3'].values[i],1 - dp['Xenv'].values[i] - dp['Zenv'].values[i],dp['Surface Abundances: C12'].values[i], dp['Surface Abundances: O16'].values[i], dp['Central: log(T)'].values[i],dp['Central: log(RHO)'].values[i],0,dp['Central: X'].values[i],1 - dp['Central: X'].values[i] - dp['Central: Z'].values[i],dp['Central Abundances: C12'].values[i]]
		#listzip = list(zip(dp['AGE'],lendp*[dp['Mconv. env.'].values[0]],lendp*[0],dp['M He core'],lendp*[0], dp['log(L/Lsun)'],np.log10(dp['Luminosity: ppI'] + dp['Luminosity: ppII'] + dp['Luminosity: ppIII'] + dp['Luminosity: CNO'] + 1e-10), np.log10(dp['Luminosity: triple-alpha'] + 1e-10),dp['log(Teff)'],dp['log(R/Rsun)'],dp['log(g)'],dp['Xenv'],dp['Surface Abundances: He3'],1 - dp['Xenv'] - dp['Zenv'],dp['Surface Abundances: C12'], dp['Surface Abundances: O16'], dp['Central: log(T)'],dp['Central: log(RHO)'],lendp*[0],dp['Central: X'],1 - dp['Central: X'] - dp['Central: Z'],dp['Central Abundances: C12']))
		#listzip1 = [['']*len(listzip[0])]*len(listzip)
		for i in range(len(listzip)):
			for j in range(len(listzip[i])):
				try:
					listzip[i][j] = str(self.eformat(listzip[i][j], 16, 3))
				except:
					listzip[i][j] = listzip[i-1][j]
		dfloc=pd.DataFrame(listzip, columns = hdrloc,dtype=str)
		data3 = dfloc.to_string(col_space=40, index=False)
		with open("{}.data".format(filename), 'w') as fp:
			fp.write('                                       1                                        2                                        3                                        4                                        5                                        6                                        7                                        8\n')
			fp.write('                          version_number                                 compiler                                    build                         MESA_SDK_version                             math_backend                                     date                                burn_min1                                burn_min2\n')
			fp.write('                                 "15140"                               "gfortran"                                 "10.2.0"                   "x86_64-linux-20.12.1"                                 "CRMATH"                               "20210328"                  5.0000000000000000E+001                  1.0000000000000000E+003\n')
			fp.write('\n')
			fp.write('                                       1                                        2                                        3                                        4                                        5                                        6                                        7                                        8                                        9                                       10                                       11                                       12                                       13                                       14                                       15                                       16                                       17                                       18                                       19                                       20                                       21                                       22 \n')
			fp.write(data3)
