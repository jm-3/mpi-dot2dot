#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import json
import matplotlib.pyplot as plt
import matplotlib.ticker as tckr
from os import listdir
import re
from matplotlib.axis import Axis  
import statistics
import os
import numpy as np
import seaborn as sns
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("dir", help="Directory with benchmark results (form: dataset)", nargs="+")
args = parser.parse_args()

parameter_values=[]
fig, (axes1, axes2) = plt.subplots(nrows=1, ncols=2)
df = pd.DataFrame(columns=['Alg', 'procs', 'mem'])

for dataset in args.dir:
        dataset_name = os.path.basename(dataset)
        axes1.set_title(dataset_name)
        axes2.set_title(dataset_name)
        algnames = listdir(dataset)
        for algname in algnames:
                if (re.match("secuencial.*", algname)):
                        print("Skiping " + algname) 
                        continue
                confs = listdir(dataset+"/"+algname)
                parameter_values = [re.search("threads-(\d+)",b).group(1) for b in confs]
                times_mean = []
                mem_mean = []
                for conf in confs:
                        tests = listdir(dataset+"/"+algname+"/"+conf)
                        times_values = []
                        mem_values = []
                        for test in tests:
                                try:
                                        with open(dataset + "/" + algname + "/" + conf + "/" + test + "/sacct.txt") as f:
                                                lines = f.read().splitlines()
                                                last_line = lines[-1]
                                                fields = last_line.split()
                                                times_values.append(int(fields[4]))
                                                mem_values.append(int(fields[7][:-1]))
                                                result = df.append(pd.Series([algname , int(conf[8:]), int(fields[7][:-1])], index=df.columns),ignore_index=True)
                                                df = result
                                except:
                                        print("Falta test: " + dataset + "/" + algname + "/" + conf + "/" + test)
                        if(len(times_values) > 1):
                                times_mean.append(statistics.mean(times_values))
                                mem_mean.append(statistics.mean(mem_values))
                        else:
                                times_mean.append(0)
                                mem_mean.append(0)


                parameter_values = list(map(int,parameter_values))
                times_mean = [x for _,x in sorted(zip(parameter_values,times_mean))]
                mem_mean = [x/1024 for _,x in sorted(zip(parameter_values,mem_mean))]
                parameter_values.sort()

                time_base = times_mean[0]
                speedup = [time_base / number for number in times_mean]

                print("parameter " + str(parameter_values))
                print("medias " + str(times_mean))
                print("speedup " + str(speedup))
                print("memoria " + str(mem_mean))


                axes1.plot(parameter_values,
                        speedup,
                        marker="o",
                        label=algname
                )

                #ind = np.arange(len(mem_mean))  # the x locations for the groups
                #width = 0.35  # the width of the bars
                #axes2.bar(ind - width/2,
                #       mem_mean,
                #       width,
                #       label=algname)   

""" add caso ideal """
axes1.plot([1, 384],
         [1, 384],
         "r-",
        linestyle="dashed",
         marker="",
         label="Ideal")



axes1.set_xlabel("nº hebras/procesos")
axes1.set_ylabel("speedup")
#axes.set_xlim([0,28])
#axes.set_ylim([0,28])

#axes2.set_xlabel("nº hebras/procesos")
#axes2.set_ylabel("RAM (GB)")

sns.barplot(y="mem", x="procs", hue="Alg", data=df)


#axes2.set_xticks(ind)
#axes2.set_xticklabels(('1', '2', '4', '8', '12', '24', '48', '96', '192', '384'))
#axes2.set_yscale('log')

axes1.legend()
axes2.legend()

print(df)

plt.tight_layout()
plt.show()
