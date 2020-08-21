import numpy
import subprocess
import re
import os


root_dir = os.getcwd()
print('Current directory is' + root_dir)

#Pl = ['0.0', '0.002', '0.004', '0.006','0.008', '0.010', '0.012', '0.014','0.016', '0.018', \
#      '0.020', '0.022', '0.024', '0.026','0.028', '0.030', '0.032', '0.034','0.036', '0.038', \
#      '0.040', '0.042', '0.044', '0.046','0.048', '0.050', '0.052', '0.054','0.056', '0.058', \
#      '0.060', '0.062', '0.064', '0.066','0.068', '0.070', '0.072', '0.074','0.076', '0.078', \
#      '0.080', '0.082', '0.084', '0.086','0.088', '0.090', '0.092', '0.094','0.096', '0.098']

#Pl = ['0.1', '0.102', '0.104', '0.106','0.108', '0.110', '0.112', '0.114','0.116', '0.118', \
#      '0.120', '0.122', '0.124', '0.126','0.128', '0.130', '0.132', '0.134','0.136', '0.138', \
#      '0.140', '0.142', '0.144', '0.146','0.148', '0.150', '0.152', '0.154','0.156', '0.158', \
#      '0.160', '0.162', '0.164', '0.166','0.168', '0.170', '0.172', '0.174','0.176', '0.178', \
#      '0.180', '0.182', '0.184', '0.186','0.188', '0.190', '0.192', '0.194','0.196', '0.198']

b_star = ['0.0','0.5','1.0','1.5','2.0','2.5'] #,'6.0','7.0','8.0']#,'0.4','0.5','0.6','0.7'] #['1.0','2.0','4.0','8.0', '16.0','32.0','64.0','128.0']

k_star = ['0.0'] #['1.0','1.5','2.0','3.0','4.0','5.0']

b_star_count = 20
for i in b_star:
  b_star_count = b_star_count + 1
  k_star_count = 0
  for j in k_star:
    k_star_count = k_star_count + 1
    
    if b_star_count < 10:
      extra_str = '0'
    else:
      extra_str = ''

    # Make directory
    dirstr = extra_str + str(b_star_count) + '_' + str(k_star_count)

    command = 'mkdir ' + dirstr
    mkd = subprocess.Popen(command, shell = True)
    mkd.wait()

    # Copy code
    command = 'cp ./code/* ./' + dirstr
    mkd = subprocess.Popen(command, shell = True)
    mkd.wait()

    # Enter directory
    print('Changing to directory ' + dirstr)
    os.chdir(dirstr)

    # Set up strings to replace
    f = open('qg_jacob20_repl.py')
    text = f.read()
    f.close()

    SF = re.compile('REPL1')
    SF2 = re.compile('REPL2')
    ############################
    ## Basic constants
    new = SF.sub('b_star = ' + i,text)
    new = SF2.sub('k_star = ' + j, new)
    f = open('qg_jacob20_repl.py', 'w')
    f.write(new)
    f.close()


#     # Set up strings to replace
#    f = open('qg_jacob20_restart_repl.py')
#    text = f.read()
#    f.close()
#
#    SF = re.compile('REPL1')
#    SF2 = re.compile('REPL2')
#    ############################
#    ## Basic constants
#    new = SF.sub('b_star = ' + i,text)
#    new = SF2.sub('k_star = ' + j, new)
#    f = open('qg_jacob20_restart_repl.py', 'w')
#    f.write(new)
#    f.close()
 
    ###### Edit start_dedalus
    f = open('start_dedalus.sh')
    text = f.read()
    f.close()

    SF = re.compile('REPL1')
    new = SF.sub('#PBS -d /project/oce/bachman/MITgcm/dedalus/two_layer_beta/' + dirstr + '/',text)   ##
    f = open('start_dedalus.sh', 'w')
    f.write(new)
    f.close()

    f = open('start_dedalus.sh')
    text = f.read()
    f.close()
 
    SF = re.compile('REPL2')
    new = SF.sub('mpiexec -v -n 16 -display-map --hostfile $PBS_NODEFILE /home/bachman/.linuxbrew/bin/python3 /project/oce/bachman/MITgcm/dedalus/two_layer_beta/' + dirstr + '/qg_jacob20_repl.py',text)

    f = open('start_dedalus.sh', 'w')
    f.write(new)
    f.close()
 
#    ###### Edit restart_dedalus
#    f = open('restart_dedalus.sh')
#    text = f.read()
#    f.close()
#
#    SF = re.compile('REPL1')
#    new = SF.sub('#PBS -d /project/oce/bachman/MITgcm/dedalus/QG/20_LAYERS/' + dirstr + '/',text)   ##
#    f = open('restart_dedalus.sh', 'w')
#    f.write(new)
#    f.close()
#
#    f = open('restart_dedalus.sh')
#    text = f.read()
#    f.close()
#
#    SF = re.compile('REPL2')
#    new = SF.sub('mpiexec -v -n 16 -display-map --hostfile $PBS_NODEFILE /home/bachman/.linuxbrew/bin/python3 /project/oce/bachman/MITgcm/dedalus/QG/20_LAYERS/' + dirstr + '/qg_jacob20_restart_repl.py',text)
#
#    f = open('restart_dedalus.sh', 'w')
#    f.write(new)
#    f.close()


    # Submit job
    command = 'qsub start_dedalus.sh'
    mkd = subprocess.Popen(command, shell = True)
    mkd.wait()

    os.chdir(root_dir)
