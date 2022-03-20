

#######################################
## script for the generation of figures
#######################################


#######################################
## import the libraries
#######################################

import pandas as pd
import numpy as np
import plotly.express as px
import os


#######################################
## Load the data
#######################################

## load the epidemiological data
epidemiological_data = pd.read_csv(
    snakemake.input['subgroup_epidemiological_data_all'],
    sep = ","
)

#######################################
## Analysis of the data
#######################################

#print(epidemiological_data.columns.values)
#exit()


histology_counting = epidemiological_data.loc[:,['cluster', 'classification2']] 
histology_counting = histology_counting.replace(np.nan, 'unknown', regex = True) # replace the NAN value by unknown
histology_counting = histology_counting.groupby(by = ['cluster','classification2']).size().reset_index(name='counts')

print(histology_counting.loc[:, 'classification2'].values)


fig = px.sunburst(
    histology_counting,
    path = ['cluster', 'classification2'],
    values = 'counts',
    color = 'classification2'
    #color_discrete_map = {
    #    'glioblastoma':'darkred',
    #    'oligodendroglioma':'gold',
    #    'astrocytoma':'darkblue',
    #    'unknown':'grey'
    #}
)

colors=['#EF553B',
         '#636EFA',
         '#00CC96',
         '#AB63FA',
         '#FFA15A',
         '#19D3F3',
         '#FF6692',
         '#B6E880',
         '#FF97FF',
         '#FECB52']

D = histology_counting['classification2'].unique()

for i, m in enumerate(D):
    fig.add_annotation(dict(font=dict(color=colors[i],size=14),
                                        x=0.8,
                                        y=1-(i/10),
                                        showarrow=False,
                                        text=D[i],
                                        textangle=0,
                                        xanchor='left',
                                        xref="paper",
                                        yref="paper"))

fig.update_layout(
    paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)'
)



fig.update_layout(
    uniformtext = dict(minsize=18)
)


fig.show()

fig.write_image("data/test.svg")


exit()


