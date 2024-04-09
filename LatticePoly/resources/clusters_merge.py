import os
import numpy as np
import json
full_forks=np.zeros((2000,50))
repl_dir = "/Xnfs/physbiochrom/ddasaro/year3/October23/cohesion_S/cohesion_mode_0_loops_0.0006_slow_Jsister_100"
#repl_dir = "./"
for k,folder in enumerate(os.listdir(repl_dir)):
    if folder.startswith('.')==False and  folder.endswith('gz')==False and folder.endswith('cool')==False and folder.endswith('.res')==False :
        for file_name in os.listdir(repl_dir+"/"+folder):
            if file_name.endswith('.json') and file_name.endswith('cool')==False and file_name.endswith('.res')==False :
                file_path = os.path.join(repl_dir+"/"+folder, file_name)
                if  os.path.getsize(file_path) != 0:
                    print(file_path)
                    cluster = json.load(open(file_path))
                    for i in range(2000):
                        newclu=list(filter(lambda a: a != -1, cluster[str(i)]))
                        for e in newclu:
                            full_forks[i,e]+=1
                            #print(e)
                  
np.savetxt("/Xnfs/physbiochrom/ddasaro/year3/October23/cohesion_S/cohesion_mode_0_loops_0.0006_slow_Jsister_100/clusters.res",full_forks)
