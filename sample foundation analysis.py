import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats


# In[4]:


def independence_test(drugUSE, drugNMU, df):
    drug = [n-ca[drugUSE].sum(), int(ca[drugUSE].sum()-ca[drugNMU].sum()), int(ca[drugNMU].sum())]
    #print(drug)

    gender = [int(ca.loc[ca['DEM_GENDER'] == 1, 'DEM_GENDER'].sum()), int(ca.loc[ca['DEM_GENDER'] == 2, 'DEM_GENDER'].sum()/2)]
    #print(gender)

    genderDrugObserved = [[ca.loc[(ca['DEM_GENDER'] == 1) & (ca[drugUSE] == 0), drugUSE].count(), 
                               ca.loc[(ca['DEM_GENDER'] == 1) & (ca[drugNMU] == 0), drugNMU].count(), 
                               ca.loc[(ca['DEM_GENDER'] == 1) & (ca[drugNMU] == 1), drugNMU].count()],
                              [ca.loc[(ca['DEM_GENDER'] == 2) & (ca[drugUSE] == 0), drugUSE].count(), 
                               ca.loc[(ca['DEM_GENDER'] == 2) & (ca[drugNMU] == 0), drugNMU].count(), 
                               ca.loc[(ca['DEM_GENDER'] == 2) & (ca[drugNMU] == 1), drugNMU].count()]]
    #print(genderDrugObserved)

    genderDrugExpected = []
    for i in range(len(gender)):
        new=[]
        for j in range(len(drug)):
            new.append(gender[i]*drug[j]/n)
        genderDrugExpected.append(new)
    #print(genderDrugExpected)

    chisq = 0
    for i in range(len(gender)):
        for j in range(len(drug)):
            chisq += (genderDrugObserved[i][j] - genderDrugExpected[i][j]) * (genderDrugObserved[i][j] - genderDrugExpected[i][j]) / genderDrugExpected[i][j]
    #print(chisq)

    #print('p-value: ', 1 - stats.chi2.cdf(chisq, df))
    
    return (1 - stats.chi2.cdf(chisq, df))




