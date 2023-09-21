#This code is used to convert the standard output
#from iso_D to the format that can be used in rml analysis
#We assume the output has columns from my_history_columns_basic.list

import pandas as pd
def __init__(self, iso_path, wrt_path):
    #define columns width
    colspecs = [(0,5)]
    for i in range(23):
        colspecs.append((5 + i*32, 37 + i*32))
    df = pd.read_fwf(iso_path, colspecs=colspecs,skiprows=10)
    df1 = df[['# EEP','isochrone_age_yr','initial_mass','log_L', 'log_R']].copy()
    df2 = df1[(df1['# EEP'] != '# num') & (df1['# EEP'] != '#   1') & (df1['# EEP'] != '# EEP')].astype(float).copy()
    df2.to_csv(wrt_path,index=False,header=['EEP', 'Age', 'MMs','LogLLs','LogRRs'])