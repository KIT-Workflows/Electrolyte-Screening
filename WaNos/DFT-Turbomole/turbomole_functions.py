import os,glob,sys,yaml
import subprocess
import numpy as np

#function definitions

def utf8_enc(var):
    if sys.version_info >= (3, 0): return var.encode('utf-8')
    else: return var

def utf8_dec(var):
    if sys.version_info >= (3, 0): return var.decode('utf-8')
    else: return var

def humo_lumo_number_orbt(file_name, string_to_search):
    #Searching for the given string in file along with its line numbers
    line_number = 0
    list_of_results = []
    # Opening the file in read only mode
    with open(file_name, 'r') as read_obj:
        # Reading all lines in the file one by one by iterating the file
        for line in read_obj:
            # checking each line, if the line contains the string
            if string_to_search in line:
                # If it contains the string, then add the line number & line as a tuple in the list
                list_of_results.append((line_number, line.rstrip()))
            line_number += 1

    # Return list of tuples containing line numbers and lines where string is found

    with open(file_name, 'r') as file:
        content = file.readlines()

    line_number = 1 + list_of_results[0][0]

    homo_line = int(content[line_number].split()[2])
    lumo_line = int(content[line_number+1].split()[2])
    
    return homo_line, lumo_line

def make_define_str(settings,coord):
    
    if settings['use old mos']:
        with open('old_results/rendered_wano.yml') as infile:
            old_settings = yaml.full_load(infile)
        same_basis = old_settings['Basis set']['Basis set type'] == settings['basis set']

    if not (settings['use old mos'] and same_basis):
        define_string='\n%s\na %s\n'%(settings['title'],coord)
        #implement symmetry
        if settings['int coord']: define_string+='ired\n*\n'
        else: define_string+='*\nno\n'

        define_string+='b all %s\n*\n'%(settings['basis set']) # add options for basis sets (different for different atoms)

        if not settings['use old mos']:
            define_string+='eht\n\n'
            define_string+='%i\n'%(settings['charge'])
            #implement symmetry
            if settings['multiplicity'] < 3: define_string+='\n\n\n'
            else: define_string+='n\nu %i\n*\n\n'%(settings['multiplicity']-1)
        else: 
            define_string+='use old_results/control\n\n'

        define_string+='scf\niter\n%i\n\n'%(settings['scf iter'])
        if settings['use ri']: define_string+='ri\non\nm %i\n\n'%(settings['ricore'])
        if settings['functional'] != 'None': define_string+='dft\non\nfunc %s\ngrid %s\n\n'%(settings['functional'],settings['grid size'])
        if settings['disp'] != 'off': define_string+='dsp\n%s\n\n'%(settings['disp'])

        if settings['tddft']:
            define_string+='ex\n'
            if settings['multiplicity'] > 1: define_string+='urpa\n*\n'
            else: define_string+='rpa%s\n*\n'%(settings['exc state type'][0].lower())
            define_string+='a %i\n*\n'%(settings['num exc states'])
            define_string+='*\n\n'
            #implement symmetry
        define_string+='*\n'
    
    else:
        output_files=['alpha','auxbasis','basis','beta','control','hessapprox','mos']
        for filename in output_files:
            if os.path.isfile('old_results/%s'%(filename)): os.system('cp old_results/%s .'%(filename))
        os.system('cp coord_0 coord')
        define_string='\n\n\n\n\n'
        if settings['use ri']: define_string+='ri\non\nm %i\n\n'%(settings['ricore'])
        else: define_string+='ri\noff\n\n'
        if settings['functional'] != 'None': define_string+='dft\non\nfunc %s\ngrid %s\n\n'%(settings['functional'],settings['grid size'])
        else: define_string+='dft\noff\n\n'
        define_string+='dsp\n%s\n\n'%(settings['disp'])
        
        if not settings['tddft']:
            if old_settings['Type of calculation']['Excited states calculation']:
                os.system('sed -i \'s/#$max/$max/g\' control')
                for dg in 'soes','scfinstab','rpacor','denconv': os.system('kdg %s'%(dg))

        elif not old_settings['Type of calculation']['Excited states calculation']:
            define_string+='ex\n'
            if settings['multiplicity'] > 1: define_string+='urpa\n*\n'
            else: define_string+='rpa%s\n*\n'%(settings['exc state type'][0].lower())
            define_string+='a %i\n*\n'%(settings['num exc states'])
            define_string+='*\n\n'
            #implement symmetry
        else: 
            define_string+='ex\n'
            if old_settings['Type of calculation']['TDDFT options']['Type of excited states'] != settings['exc state type']:
                if settings['multiplicity'] > 1: define_string+='urpa\n'
                else: define_string+='rpa%s\n'%(settings['exc state type'][0].lower())
            define_string+='*\n'
            define_string+='a %i\n*\n'%(settings['num exc states'])
            define_string+='*\n\n'

        define_string+='*\n'

    return define_string

def inputprep(program,input_string):
    outfilename='%s.out'%(program)
    process=subprocess.Popen([program],stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate(input=utf8_enc(input_string))

    with open(outfilename,'w') as outfile: outfile.write(utf8_dec(out))
    if not 'normally' in utf8_dec(err).split():
       print('An error occured when running %s:'%(program))
       with open(outfilename) as infile:
           lines=infile.readlines()
       for line in lines: print(line)

def single_point_calc(settings,tmp=False):
    if settings['use ri']: 
        scf_program='ridft'
    else: 
        scf_program='dscf'

    if tmp: 
        suffix = '_tmp'
    elif settings['opt']: 
        suffix = '_0'
    else: 
        suffix = ''

    output = scf_program+suffix+'.out'

    num_iter,done = 0,False
    while not done:
        run_turbomole(scf_program,output)
        os.system('eiger > eiger.out')
        done,err = check_scf(output)
        if not done:
            if err == 'not converged':
                num_iter += settings['scf iter']
                if num_iter > settings['max scf iter']:
                    print('SCF not converged in maximum number of iterations (%i)'%(settings['max scf iter']))
                    exit(0)
            elif err == 'negative HLG':
                print('Attention: negative HOMO-LUMO gap found - please check manually')
                exit(0)

    if settings['tddft']:
        run_turbomole('escf','escf%s.out'%(suffix))
        done,err = check_escf(output)
        if not done:
            print('Problem with escf calculation found - please check manually')
            exit(0)

def aoforce():
    run_turbomole('aoforce','aoforce.out')

def plt_homo_lumo_orbt():
    run_turbomole('nohup riper -proper','riper.out')

def hyper_pol():
    run_turbomole('escf','escf.out')

def jobex(settings):
    options = ''
    if settings['use ri']: options += ' -ri'
    if settings['tddft']: options += ' -ex %i'%(settings['opt exc state'])
    options += ' -c %i'%(settings['opt cyc'])
    
    num_cycles,done = 0,False

    while not done:
        run_turbomole('jobex%s'%(options))
        os.system('eiger > eiger.out')
        done,err = check_opt()

        if not done:
            if err == 'opt not converged':
                num_cycles += settings['opt cyc']
                if num_cycles > settings['max opt cyc']:
                     print('Structure optimisation not converged in maximum number of cycles (%i)'%(settings['max opt cyc']))
                     exit(0)
            elif err.startswith('scf problem'):
                num_cycles += int(err.split()[-1])
                single_point_calc(settings,tmp=True)
            else:
                print('An error occurred during the structure optimisation - please check manually')
                exit(0)
            #implement other errors

        elif not settings['tddft']:
           scf_done,err = check_scf('job.last')
           if not scf_done:
               print('The Structure optimisations converged but the last SCF run showed an error.')
               exit(0)

def run_turbomole(command,outfile = None):
    if outfile == None: outfile=command.split()[0]+'.out'
    
    proc_cmd = ['nohup']+command.split()
    if not command.startswith('jobex'):
        if outfile == None: outfile = command+'.out'
        proc_cmd += ['>','outfile']

    with open(outfile,'w') as tm_out:
        tm_process = subprocess.Popen(proc_cmd,stdout=tm_out,stderr=subprocess.PIPE)
        out, err = tm_process.communicate()

    err = utf8_dec(err)

    # if not 'normally' in err.split('\n')[-2].split():
    #     print('Error while running %s:'%(command.split()[0]))
    #     print(err)
    #     exit(0)

def check_scf(output_file):
    conv = False
    with open(output_file,'r') as infile:
        for line in infile.readlines():
            if 'convergence criteria cannot be satisfied' in line:
                conv = False
                break
            if 'convergence criteria satisfied' in line:
                conv = True
                break

    if conv == False:
        return False, 'not converged'
    else:
        with open('eiger.out','r') as infile:
            for line in infile.readlines():
                if 'Gap' in line:
                    hlg=float(line.split()[-2])
                    break
    if hlg < 0: return False, 'negative HLG'
    else: return True,None

def check_escf(output_file):
    conv = False
    with open(output_file,'r') as infile:
        for line in infile.readlines():
            if 'all done' in line:
                conv = True
                break
    return conv, None

def check_opt():
   
    conv = os.path.isfile('GEO_OPT_CONVERGED')

    if not conv:
        if os.path.isfile('GEO_OPT_RUNNING'): err = 'jobex did not end properly'
        else:
            with open('GEO_OPT_FAILED','r') as infile:
                for line in infile.readlines():
                    if 'OPTIMIZATION DID NOT CONVERGE' in line:
                        err = 'opt not converged'
                        break
                    if 'your energy calculation did not converge' in line:
                        err = 'scf problem during step nr. %s'%(glob.glob('job.[123456789]')[0].split('.')[-1])
                        break
    else: err = '' # check for problems with last run here, e.g. neg. HLG
    
    return conv, err 

# Get the hyper polarizability  
def get_hyper_polarizability():
    with open("escf.out", "r") as infile:
        escf_data = infile.readlines()

    # Electronic dipole hyperpolarizability
    begin_hyper_pol = int([i for i, line in enumerate(escf_data) if "Electronic dipole hyperpolarizability" in line][0]) + 4
    # read beta_x
    beta = np.zeros((3,3,3))
    lines = escf_data[begin_hyper_pol:begin_hyper_pol+9]
    line_split = [line.split() for line in lines]
    for i in range(3):
        beta[i,0,0] = float(line_split[0][i+(i+1)])
        beta[i,1,0] = float(line_split[1][i+(i+1)])
        beta[i,2,0] = float(line_split[2][i+(i+1)])
        beta[i,0,1] = float(line_split[3][i+(i+1)])
        beta[i,1,1] = float(line_split[4][i+(i+1)])
        beta[i,2,1] = float(line_split[5][i+(i+1)])
        beta[i,0,2] = float(line_split[6][i+(i+1)])
        beta[i,1,2] = float(line_split[7][i+(i+1)])
        beta[i,2,2] = float(line_split[8][i+(i+1)])
    return beta

# get the dipole_moment
def get_dipole_moment():
    with open("control", "r") as infile:
        control_data = infile.readlines()
    begin_dipole = int(
        [i for i, line in enumerate(control_data) if "$dipole" in line][0]) + 1
    line_split = control_data[begin_dipole].split()

    dipole = np.zeros((3))
    for i in range(3):
        dipole[i] = float(line_split[i+(i+1)])
    return dipole



# def man_occ(occ_vec, settings):
#     if settings['multiplicity'] == 1:
#         instring='\n\n\ny\neht\n0\n\n'
#         define_process=subprocess.Popen(['define'],stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         define_out=define_process.communicate(input=instring)[0].decode('utf-8').split('\n')
#         for line in define_out:
#             if 'NUMBER OF ELECTRONS IN YOUR MOLECULE IS' in line:
#                 n_occ=int(line.split()[-1])/2
#                 break
#         instring+='c 1-%i\n'%(n_occ-len(occ_vec))
#         for i in len(occ_vec):
#             if occ_vec[i] == 1: instring+='c %i\n'%(n_occ-len(occ_vec)+1+i)
#         instring+='*\n\n*\n'   

#     else: raise Error("manual occupation not implemented for open-shell systems")

