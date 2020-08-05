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

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("dir", help="Directory with benchmark results (form: dataset)", nargs="+")
args = parser.parse_args()

parameter_values=[]
fig, axes = plt.subplots()

for dataset in args.dir:
        dataset_name = os.path.basename(dataset)
        axes.set_title(dataset_name)
        algnames = listdir(dataset)
        for algname in algnames:
                if (re.match("secuencial.*", algname)):
                        print("Skiping " + algname) 
                        continue
                confs = listdir(dataset+"/"+algname)
                parameter_values = [re.search("threads-(\d+)",b).group(1) for b in confs]
                times_mean = []
                for conf in confs:
                        tests = listdir(dataset+"/"+algname+"/"+conf)
                        times_values = []
                        for test in tests:
                                try:
                                        with open(dataset + "/" + algname + "/" + conf + "/" + test + "/sacct.txt") as f:
                                                lines = f.read().splitlines()
                                                last_line = lines[-1]
                                                fields = last_line.split()
                                                times_values.append(int(fields[4]))
                                except:
                                        print("Falta test: " + dataset + "/" + algname + "/" + conf + "/" + test)
                        if(len(times_values) > 1):
                                times_mean.append(statistics.mean(times_values))
                        else:
                                times_mean.append(0)


                parameter_values = list(map(int,parameter_values))
                times_mean = [x for _,x in sorted(zip(parameter_values,times_mean))]
                parameter_values.sort()

                time_base = times_mean[0]
                speedup = [time_base / number for number in times_mean]

                print("parameter " + str(parameter_values))
                print("medias " + str(times_mean))
                print("speedup " + str(speedup))


                axes.plot(parameter_values,
                        speedup,
                        marker="o",
                        label=algname
                )        

""" add caso ideal """
axes.plot([1, 384],
         [1, 384],
         "r-",
        linestyle="dashed",
         marker="",
         label="Ideal")



axes.set_xlabel("nº hebras/procesos")
axes.set_ylabel("speedup")
#axes.set_xlim([0,28])
#axes.set_ylim([0,28])

axes.legend()

plt.tight_layout()
plt.show()
