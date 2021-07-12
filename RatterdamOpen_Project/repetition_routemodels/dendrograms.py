# -*- coding: utf-8 -*-
"""
Created on Sun Jul 11 19:13:46 2021

@author: whockei1
"""

rownames = []

for i in range(df.shape[0]):
    n = f"{df.iloc[i]['dirM1']}{df.iloc[i]['dirC']}{df.iloc[i]['dirP1']}"
    rownames.append(n)
    

savepath = 'E:\\Ratterdam\\R_data_repetition\\210624-onward_turnLM\\models_21-07-10\\'
    
for unitname in df['unit'].unique():
    print(unitname)
    cdf = df.loc[df['unit']==unitname]
    for field in cdf['field'].unique():
        fdf = cdf.loc[cdf['field']==field]
        Z = linkage(np.reshape(list(fdf["rate"]),(fdf["rate"].shape[0],1)),'single')
        d = dendrogram(Z, leaf_rotation=90,leaf_font_size=10,labels=list(fdf['turnCode']))
        plt.title(f"{unitname} Field {field} Dendrogram method=single")
        
        newname = unitname.split("\\")[0]+unitname.split("\\")[1]+f"_{field}_dendrogram"
        fig=plt.gcf()
        fig.set_size_inches(15,12)
        plt.savefig(savepath+newname+".png",dpi=300)
        plt.close()
        