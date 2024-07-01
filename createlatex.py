import sys
import numpy as np


loaded_data = np.load("./potential_loop/data.npz")

nsites_array = loaded_data["nsites_array"]
nsources_array = loaded_data["nsources_array"]
nmesh_array = loaded_data["nmesh_array"]
noise_coeff_array = loaded_data["noise_coeff_array"]
finit_array = loaded_data["finit_array"]
noise_level_array = loaded_data["noise_level_array"]
normgpinit_array = loaded_data["normgpinit_array"]
ffinal_array = loaded_data["ffinal_array"]
normgpfinal_array  = loaded_data["normgpfinal_array"]
erroropt_array  = loaded_data["erroropt_array"]
errorinit_array  = loaded_data["errorinit_array"]
flagsol_array = loaded_data["flagsol_array"]
iter_array = loaded_data["iter_array"]
numevalf_array = loaded_data["numevalf_array"]
CPU_time_array = loaded_data["CPU_time_array"]


#def createLatex(nsites_list, nsources_list,nmesh_list, noise_coeff_list, error_opt, error_init, product_dict):

with open("./potential_loop/results.tex", "w") as f:
    f.write(r"""
\documentclass{article}
\usepackage{graphicx}
\usepackage{array}
\usepackage{subcaption}
%\usepackage{fullpage}
\usepackage{verbatim}
\usepackage[margin=0in]{geometry}
\begin{document}
\begin{table}[htb]
\centering
\begin{tabular}{|m{2cm}|m{2cm}|m{5cm}|m{5cm}|m{5cm}|}
\hline
data & results & ground truth & initialization & optimized diagram \\
\hline
""")

    for k in range(len(nsites_array)): 

        nsites, nsources, nmesh, noise_coeff, finit, noise_level, normgpinit, ffinal, normgpfinal, erroropt, errorinit, flagsol, iter, numevalf, CPU_time = nsites_array[k], nsources_array[k], nmesh_array[k], noise_coeff_array[k], finit_array[k], noise_level_array[k], normgpinit_array[k], ffinal_array[k], normgpfinal_array[k], erroropt_array[k], errorinit_array[k], flagsol_array[k], iter_array[k], numevalf_array[k], CPU_time_array[k]

        k = k + 1

        #[finit, normgpinit, ffinal, normgpfinal, noise_level, error_opt, error_init, error_flag] = product_dict[f'{nsites}_{nsources}_{nmesh}_{noise_coeff}']

        input_folder = './'+str(nsites)+'_'+str(nsources)+'_'+str(nmesh)+'_'+str(int(1000.0*noise_coeff))+'/'+str(1)+'/'
        f.write(r""" nsites =  """ + str(nsites)+ r""" \newline nsources =  """ + str(nsources)+r""" \newline nmesh = """ + str(nmesh) +r""" \newline noiselevel = """ + str(round(noise_level,5))+r""" \newline flagsol = """ + str(flagsol)+r""" \newline iter = """ + str(iter)+r""" \newline numevalf = """ + str(numevalf)+r""" \newline CPU time (s) = """ + str(round(CPU_time, 2))+r"""  
        & erroropt =  """ + str(round(erroropt,5))+ r""" \newline errorinit = """ + str(round(errorinit,5))+ r""" \newline f initial = """ +\
            str(round(finit,5)) +r"""\newline g initial = """ + str(round(normgpinit,5))+r"""\newline f final = """ + str(round(ffinal,5)) +r"""\newline g final = """ + str(round(normgpfinal,5))  +r"""&
        \begin{minipage}{.4\textwidth}\includegraphics[width=0.4\textwidth]{"""+ input_folder  + r"""ground.png} \end{minipage} & 
        \begin{minipage}{.4\textwidth}\includegraphics[width=0.4\textwidth]{"""+ input_folder  + r"""init.png} \end{minipage}& 
        \begin{minipage}{.4\textwidth}\includegraphics[width=0.4\textwidth]{"""+ input_folder  + r"""optimized.png} \end{minipage}\\
        \hline
            """)
        if k % 4 == 0:
            f.write(r"""
            \end{tabular}
            \end{table}
            \newpage
            \begin{table}[htb]
            \centering
            \begin{tabular}{|m{2cm}|m{2cm}|m{5cm}|m{5cm}|m{5cm}|}
            \hline
            data & results & ground & init & optimized \\
            \hline
                """)

    f.write(r"""
    \end{tabular}
    \end{table}
    \end{document}""")
    
    

print('Latex file results.tex was created')