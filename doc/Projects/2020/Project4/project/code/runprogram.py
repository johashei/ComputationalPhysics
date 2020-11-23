import subprocess,os
"""
All the runs needed to generate the data used in the final report
"""
# Make a directory to store datafiles
dir = "datafiles"
if not os.path.isdir("./"+dir):
    os.mkdir(dir)
'''
print("Compiling")
subprocess.run(['g++-10','-fopenmp','-O3','-c','main.cpp','not_para.cpp','Ising.cpp','-std=c++11'])
subprocess.run(['g++-10','-fopenmp','-O3','-o','main_omp','main.o','Ising.o'])
subprocess.run(['g++-10','-O3','-o','main_sng','not_para.o','Ising.o'])


print("Working on 4d) and 4e) : L=20, MCS = 1e6")
T_list = ['1.0',    '1.0',      '2.4',      '2.4']
S_list = ['aligned','random',   'aligned',  'random']
F_list = ['d1a.txt','d1r.txt',  'd24a.txt', 'd24r.txt']
for T,S,F in zip(T_list,S_list,F_list):
    cmd = ['./main_sng','20','1e6',T,S,F]
    print('\t',' '.join(cmd))
    subprocess.run(cmd)
    subprocess.run(['mv',F,dir])
'''

print("Working on 4f) MCS = 3e6")
L_list = [  '40',       '40',       '60',       '60',       '80',       '80',       '100',       '100']
F_list = [  'f40ZC.txt','f40ZX.txt','f60ZC.txt','f60ZX.txt','f80ZC.txt','f80ZX.txt','f100ZC.txt','f100ZX.txt']
Ts_list = [ '2.278',    '2.317',    '2.278',    '2.294',    '2.277',    '2.286',    '2.270',     '2.286']
Tf_list = [ '2.290',    '2.323',    '2.290',    '2.306',    '2.283',    '2.298',    '2.282',     '2.298']
dT_list = [ '0.004',    '0.002',    '0.004',    '0.004',    '0.002',    '0.004',    '0.004',     '0.004']
for L,F,Ts,Tf,dT in zip(L_list,F_list,Ts_list,Tf_list,dT_list):
    cmd = ['./main_omp',L,'3e6','1e5',Ts,Tf,dT,'default',F]
    print('\t',' '.join(cmd))
    subprocess.run(cmd)
    subprocess.run(['mv',F,dir])
