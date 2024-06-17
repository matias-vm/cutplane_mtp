#!/usr/bin/python

###############################################################################
#
# This code was written and is being maintained by Matias Villagra,
# PhD Student in Operations Research @ Columbia, supervised by 
# Daniel Bienstock.
#
# Please report any bugs or issues (for sure there will be) to
#                         mjv2153@columbia.edu
#
# Oct 2023
###############################################################################

import sys
import os
import numpy as np
import time
import reader
from myutils import breakexit
from versioner import *
from log import danoLogger
from cutplane_mtp import gocutplane

def read_config(log, filename):

    log.joint("reading config file " + filename + "\n")

    try:
        f = open(filename, "r")
        lines = f.readlines()
        f.close()
    except:
        log.stateandquit("cannot open file", filename)
        sys.exit("failure")

    casefilename                 = 'NONE'
    lpfilename                   = ""
    lpfilename_cuts              = ""
    linenum                      = 0

    i2                           = 0

    jabrcuts                     = 0
    i2cuts                       = 0
    losscuts                     = 0
    limitcuts                    = 0
    objective_cuts               = 0    
    loud_cuts                    = 0



    most_violated_fraction_loss  = 1
    
    jabr_inequalities            = 0
    i2_inequalities              = 0
    limit_inequalities           = 0
    loss_inequalities            = 0
    qloss_inequalities           = 0



    threshold                    = 1e-5
    threshold_i2                 = 1e-1
    tolerance                    = 1e-5
    threshold_objcuts            = 1e-5
    threshold_dotprod            = 5e-1
    most_violated_fraction_jabr  = 1
    dropjabr                     = 0
    droploss                     = 0
    dropi2                       = 0
    cut_age_limit                = 20
    solver_method                = 2
    primal_bound                 = 'NONE'
    crossover                    = 0
    timelimit                    = 'NONE'
    max_rounds                   = 100
    cut_analysis                 = 0

    most_violated_fraction_i2    = 1
    mincut                       = 0
    mincut_active                = 0
    mincut_reactive              = 0
    mincut_switch                = 0

    dographics                   = 0
    
    obbt                         = 0

    hybrid                       = 0

    jabr_validity                = 0
    i2_validity                  = 0
    loss_validity                = 0
    qloss_validity               = 0
    limit_validity               = 0

    getsol                       = 0

    fixflows                     = 0
    fixcs                        = 0
    fix_tolerance                = 1e-5
    ampl_sol                     = 0
    matpower_sol                 = 0
    writeACsol                   = 0
    writesol                     = 0
    writesol_wtol                = 0
    newsol                       = 0

    FeasibilityTol               = 1e-5


    linear_objective             = 0

    writecuts                    = 0
    addcuts                      = 0

    droplimit                    = 0
    most_violated_fraction_limit = 1
    threshold_limit              = 1e-5


    writelps                     = 0

    ftol                         = 1e-3
    ftol_iterates                = 5
    max_time                     = 200

    losstest                     = 0

    parallel_check               = 0
    
    T                            = 2
    nperturb                     = 0.01
    uniform                      = 0
    uniform_drift                = 0.02
    uniform2                     = 0
    uniform3                     = 0
    uniform4                     = 0
    uniform5                     = 0
    uniform6                     = 0        
    arpae2                       = 0
    arpae                        = 0                
    pglib_reverse                = 0
    writelastLP                  = 0

    barconvtol                   = 1e-6
    feastol                      = 1e-6
    opttol                       = 1e-6
    
    while linenum < len(lines):
        thisline = lines[linenum].split()
        if len(thisline) > 0:

            if thisline[0] == 'casefilename':
                casefilename = thisline[1]

            elif thisline[0] == 'lpfilename':
                lpfilename = thisline[1]

            elif thisline[0] == 'lpfilename_cuts':
                lpfilename_cuts = thisline[1]
            
            elif thisline[0] == 'jabr_inequalities':
                jabr_inequalities = 1

            elif thisline[0] == 'cut_analysis':
                cut_analysis = 1
                
            elif thisline[0] == 'dropjabr':
                dropjabr     = 1

            elif thisline[0] == 'droploss':
                droploss = 1

            elif thisline[0] == 'dropi2':
                dropi2 = 1                

            elif thisline[0] == 'limit_inequalities':
                limit_inequalities = 1

            elif thisline[0] == 'solver_method':
                solver_method = int(thisline[1])

            elif thisline[0] == 'objective_cuts':
                objective_cuts = 1

            elif thisline[0] == 'linear_objective':
                linear_objective = 1

            elif thisline[0] == 'i2cuts':
                i2cuts           = 1

            elif thisline[0] == 'most_violated_fraction_i2':
                most_violated_fraction_i2 = float(thisline[1])                
                
            elif thisline[0] == 'cut_age_limit':
                cut_age_limit = int(thisline[1])

            elif thisline[0] == 'threshold':
                threshold = float(thisline[1])

            elif thisline[0] == 'threshold_i2':
                threshold_i2 = float(thisline[1])
     
            elif thisline[0] == 'threshold_objcuts':
                threshold_objcuts = float(thisline[1])

            elif thisline[0] == 'threshold_dotprod':
                threshold_dotprod = float(thisline[1])
                parallel_check    = 1
                
            elif thisline[0] == 'tolerance':
                tolerance = float(thisline[1])

            elif thisline[0] == 'mincut':
                mincut          = 1
                mincut_active   = 1

            elif thisline[0] == 'mincut_reactive':
                mincut          = 1
                mincut_active   = 0
                mincut_reactive = 1

            elif thisline[0] == 'mincut_switch':
                mincut_switch   = 1

            elif thisline[0] == 'most_violated_fraction_jabr':
                most_violated_fraction_jabr = float(thisline[1])

            elif thisline[0] == 'primal_bound':
                primal_bound = float(thisline[1])

            elif thisline[0] == 'crossover':
                crossover = int(thisline[1])

            elif thisline[0] == 'max_time':
                max_time  = float(thisline[1])

            elif thisline[0] == 'max_rounds':
                max_rounds = int(thisline[1])

            elif thisline[0] == 'dographics':
                dographics = 1

            elif thisline[0] == 'obbt':
                obbt = 1

            elif thisline[0] == 'hybrid':
                hybrid = 1

            elif thisline[0] == 'i2_inequalities':
                i2_inequalities = 1

            elif thisline[0] == 'jabr_validity':
                jabr_validity = 1

            elif thisline[0] == 'i2_validity':
                i2_validity = 1

            elif thisline[0] == 'loss_validity':
                loss_validity = 1

            elif thisline[0] == 'qloss_validity':
                qloss_validity = 1

            elif thisline[0] == 'fixflows':
                fixflows = 1

            elif thisline[0] == 'ampl_sol':
                ampl_sol = 1

            elif thisline[0] == 'fixcs':
                fixcs    = 1

            elif thisline[0] == 'fix_tolerance':
                fix_tolerance = float(thisline[1])

            elif thisline[0] == 'getsol':
                getsol     = 1

            elif thisline[0] == 'loss_inequalities':
                loss_inequalities  = 1

            elif thisline[0] == 'qloss_inequalities':
                qloss_inequalities  = 1

            elif thisline[0] == 'most_violated_fraction_loss':
                most_violated_fraction_loss = float(thisline[1])

            elif thisline[0] == 'FeasibilityTol':
                FeasibilityTol = float(thisline[1])

            elif thisline[0] == 'writecuts':
                writecuts  = 1

            elif thisline[0] == 'addcuts':
                addcuts     = 1
                
            elif thisline[0] == 'fromscratch':
                addcuts     = 0

            elif thisline[0] == 'droplimit':
                droplimit  = 1

            elif thisline[0] == 'most_violated_fraction_limit':
                most_violated_fraction_limit    = 1

            elif thisline[0] == 'threshold_limit':
                threshold_limit     = float(threshold_limit)

            elif thisline[0] == 'limit_validity':
                limit_validity  = 1

            elif thisline[0] == 'limitcuts':
                limitcuts    = 1

            elif thisline[0] == 'writelps':
                writelps     = 1

            elif thisline[0] == 'ftol':
                ftol          = float(thisline[1])

            elif thisline[0] == 'ftol_iterates':
                ftol_iterates = int(thisline[1])

            elif thisline[0] == 'matpower_sol':
                matpower_sol  = 1

            elif thisline[0] == 'ampl_sol':
                ampl_sol    = 1

            elif thisline[0] == 'writeACsol':
                writeACsol    = 1

            elif thisline[0] == 'jabrcuts':
                jabrcuts      = 1

            elif thisline[0] == 'losscuts':
                losscuts      = 1

            elif thisline[0] == 'loud_cuts':
                loud_cuts     = 1

            elif thisline[0] == 'i2':
                i2            = 1

            elif thisline[0] == 'writesol':
                writesol      = 1
                writesol_wtol = 0

            elif thisline[0] == 'writesol_wtol':
                if len(thisline) == 1:
                    writesol_wtol = 1e-5
                else:
                    writesol_wtol = float(thisline[1])
                writesol = 0

            elif thisline[0] == 'newsol':
                newsol      = 1                    

            elif thisline[0] == 'losstest':
                losstest    = 1

            elif thisline[0] == 'T':
                T           = int(thisline[1])

            elif thisline[0] == 'nperturb':
                nperturb    = float(thisline[1])
                uniform     = 0
                uniform2    = 0
                uniform3    = 0
                uniform4    = 0
                uniform5    = 0
                uniform6    = 0
                arpae       = 0                
                arpae2      = 0
                pglib_reverse = 0

            elif thisline[0] == 'uniform':
                nperturb    = 0
                uniform     = 1
                uniform2    = 0
                uniform3    = 0
                uniform4    = 0
                uniform5    = 0
                uniform6    = 0
                arpae       = 0                                
                arpae2      = 0                
                pglib_reverse = 0
                
            elif thisline[0] == 'uniform_drift':
                uniform_drift = float(thisline[1])
                
            elif thisline[0] == 'uniform2':
                nperturb    = 0
                uniform     = 0
                uniform2    = 1
                uniform3    = 0
                uniform4    = 0
                uniform5    = 0                
                uniform6    = 0
                arpae       = 0                                
                arpae2      = 0
                pglib_reverse = 0

            elif thisline[0] == 'uniform3':
                nperturb    = 0
                uniform     = 0
                uniform2    = 0
                uniform3    = 1
                uniform4    = 0
                uniform5    = 0
                uniform6    = 0
                arpae       = 0                                
                arpae2      = 0
                pglib_reverse = 0                
                
            elif thisline[0] == 'uniform4':
                nperturb    = 0
                uniform     = 0
                uniform2    = 0
                uniform3    = 0
                uniform4    = 1
                uniform5    = 0
                uniform6    = 0                
                arpae       = 0                                
                arpae2      = 0
                pglib_reverse = 0                
                
            elif thisline[0] == 'pglib_reverse':
                nperturb    = 0
                uniform     = 0
                uniform2    = 0
                uniform3    = 0
                uniform4    = 0
                uniform5    = 0
                uniform6    = 0
                arpae       = 0                
                arpae2      = 0
                pglib_reverse = 1
                
            elif thisline[0] == 'arpae':
                nperturb    = 0
                uniform     = 0
                uniform2    = 0
                uniform3    = 0
                uniform4    = 0
                uniform5    = 0
                uniform6    = 0
                arpae       = 1                
                arpae2      = 0
                pglib_reverse = 0                

            elif thisline[0] == 'arpae2':
                nperturb    = 0
                uniform     = 0
                uniform2    = 0
                uniform3    = 0
                uniform4    = 0
                uniform5    = 0
                uniform6    = 0
                arpae       = 0                
                arpae2      = 1
                pglib_reverse = 0

            elif thisline[0] == 'uniform5':
                nperturb    = 0
                uniform     = 0
                uniform2    = 0
                uniform3    = 0
                uniform4    = 0
                uniform5    = 1
                uniform6    = 0                                                
                arpae2      = 0
                pglib_reverse = 0

            elif thisline[0] == 'uniform6':
                nperturb    = 0
                uniform     = 0
                uniform2    = 0
                uniform3    = 0
                uniform4    = 0
                uniform5    = 0
                uniform6    = 1                                                
                arpae2      = 0
                pglib_reverse = 0                                                

            elif thisline[0] == 'barconvtol':
                barconvtol    = float(thisline[1])

            elif thisline[0] == 'feastol':
                feastol       = float(thisline[1])

            elif thisline[0] == 'opttol':
                opttol        = float(thisline[1])

            elif thisline[0] == 'writelastLP':
                writelastLP   = 1
                
            elif thisline[0] == 'END':
                break
                
            else:
                sys.exit("main_mtp: illegal input " + thisline[0] + " bye")


        linenum += 1

    all_data                    = {}
    all_data['casefilename']    = casefilename
    casename = all_data['casename']        = casefilename.split('data/')[1].split('.m')[0] 
    all_data['lpfilename']      = all_data['casename'] + '.lp'
    all_data['lpfilename_cuts'] = all_data['casename'] + '_cuts.lp'    

    if len(lpfilename):
        all_data['lpfilename']  = lpfilename
    
    if len(lpfilename_cuts):
        all_data['lpfilename_cuts'] = lpfilename_cuts
    

    all_data['jabr_inequalities']           = jabr_inequalities
    all_data['limit_inequalities']          = limit_inequalities
    all_data['i2_inequalities']             = i2_inequalities


    all_data['max_rounds']      = max_rounds
    all_data['solver_method']   = solver_method
    all_data['cut_analysis']    = cut_analysis
    
    #cuts
    all_data['initial_threshold']      = threshold
    all_data['threshold']              = threshold
    all_data['tolerance']              = tolerance
    
    #parallel cuts
    all_data['threshold_dotprod'] = threshold_dotprod
    
    all_data['tight_cuts'] = {}
    all_data['tight_cuts_fraction'] = {}
    
    all_data['cut_age_limit'] = cut_age_limit

    all_data['linear_objective']   = linear_objective
    all_data['hybrid']             = hybrid

    if linear_objective or hybrid:
        all_data['objective_cuts']     = objective_cuts
        all_data['obj_cuts']           = {}
        all_data['num_objective_cuts'] = 0
        all_data['threshold_objcuts']  = threshold_objcuts

    all_data['jabrcuts'] = jabrcuts
    if jabrcuts:
        all_data['most_violated_fraction_jabr'] = most_violated_fraction_jabr
        all_data['jabr_cuts']                   = {}
        all_data['num_jabr_cuts_rnd']           = {}
        all_data['num_jabr_cuts_added']         = 0
        all_data['num_jabr_cuts_dropped']       = 0
        all_data['dropped_jabrs']               = []    
        all_data['jabr_cuts_info']              = {}
        all_data['max_error_jabr']              = 0
        all_data['total_jabr_dropped']          = 0

    all_data['ID_jabr_cuts']         = 0 # Before it was define iff jabrcuts =1
    all_data['num_jabr_cuts']        = 0 # Before it was define iff jabrcuts =1
    all_data['dropjabr']             = dropjabr
    all_data['jabr_validity']        = jabr_validity

    all_data['limitcuts'] = limitcuts
    if limitcuts:
        all_data['most_violated_fraction_limit'] = most_violated_fraction_limit
        all_data['limit_cuts']                   = {}
        all_data['num_limit_cuts_rnd']           = {}
        all_data['num_limit_cuts_added']         = 0
        all_data['num_limit_cuts_dropped']       = 0
        all_data['limit_cuts_info']              = {}
        all_data['threshold_limit']              = threshold_limit
        all_data['dropped_limit']                = []
        all_data['max_error_limit']              = 0
        all_data['total_limit_dropped']          = 0

    all_data['ID_limit_cuts']   = 0 # Before it was defnd iff limitcuts = 1
    all_data['num_limit_cuts']  = 0 # Before it was defnd iff limitcuts = 1
    all_data['droplimit']       = droplimit
    all_data['limit_validity']  = limit_validity    


    all_data['i2cuts']                              = i2cuts
    all_data['i2']                                  = i2

    if (all_data['i2cuts'] == 1) or (all_data['i2_inequalities']):
        all_data['i2'] = 1
        
    if i2cuts:
        all_data['most_violated_fraction_i2'] = most_violated_fraction_i2
        all_data['i2_cuts']                   = {}
        all_data['num_i2_cuts_rnd']           = {}
        all_data['num_i2_cuts_added']         = 0
        all_data['num_i2_cuts_dropped']       = 0
        all_data['i2_cuts_info']              = {}
        all_data['threshold_i2']              = threshold_i2
        all_data['dropped_i2']                = []
        all_data['max_error_i2']              = 0
        all_data['total_i2_dropped']          = 0

    all_data['ID_i2_cuts']        = 0 # Before it was defnd iff i2cuts = 1
    all_data['num_i2_cuts']       = 0 # Before it was defnd iff i2cuts = 1
    all_data['dropi2']            = dropi2
    all_data['i2_validity']       = i2_validity

    ###########
    
    all_data['primal_bound']                  = primal_bound
    all_data['crossover']                     = crossover
    all_data['max_time']                      = max_time

    all_data['mincut']                        = mincut
    all_data['mincut_active']                 = mincut_active
    all_data['mincut_reactive']               = mincut_reactive
    all_data['mincut_switch']                 = mincut_switch

    all_data['dographics']                    = dographics
    
    if all_data['dographics']:
        all_data['coordsfilename'] = None

    all_data['obbt']                          = obbt
    if all_data['obbt']:
        all_data['obbt_itcount'] = 0
        all_data['mincut_jabrs'] = []
        all_data['mincut_i2']    = []


    all_data['FeasibilityTol']                = FeasibilityTol
    all_data['fixflows']                      = fixflows
    all_data['fixcs']                         = fixcs
    all_data['fix_tolerance']                 = fix_tolerance
    all_data['matpower_sol']                  = matpower_sol
    all_data['ampl_sol']                      = ampl_sol
    all_data['writeACsol']                    = writeACsol


    all_data['loss_inequalities']  = loss_inequalities
    all_data['qloss_inequalities'] = qloss_inequalities
    all_data['losscuts']           = losscuts

    if loss_inequalities or losscuts:
        all_data['loss_cuts']                   = {}
        all_data['num_loss_cuts']               = 0
        all_data['num_loss_cuts_added']         = 0
        all_data['num_loss_cuts_dropped']       = 0
        all_data['dropped_loss']                = []
        all_data['most_violated_fraction_loss'] = most_violated_fraction_loss

        all_data['max_error_loss']              = 0

    all_data['droploss']                      = droploss
    all_data['loss_validity']                 = loss_validity


    all_data['qloss_validity']                 = qloss_validity


    all_data['writecuts']                     = writecuts
    all_data['addcuts']                       = addcuts
    all_data['writelps']                      = writelps

    all_data['ftol']                          = ftol
    all_data['ftol_iterates']                 = ftol_iterates

    all_data['loud_cuts']                     = loud_cuts
    all_data['writesol']                      = writesol
    all_data['writesol_wtol']                 = writesol_wtol
    all_data['newsol']                        = newsol

    all_data['losstest']                      = losstest
    all_data['parallel_check']                = parallel_check
    all_data['T']                             = T
    all_data['nperturb']                      = nperturb    
    all_data['uniform']                       = uniform
    all_data['uniform_drift']                 = uniform_drift    
    all_data['uniform2']                      = uniform2
    all_data['uniform3']                      = uniform3
    all_data['uniform4']                      = uniform4
    all_data['uniform5']                      = uniform5
    all_data['uniform6']                      = uniform6
    all_data['arpae']                         = arpae    
    all_data['arpae2']                        = arpae2
    all_data['pglib_reverse']                   = pglib_reverse

    all_data['barconvtol']                    = barconvtol
    all_data['feastol']                       = feastol
    all_data['opttol']                        = opttol
    all_data['writelastLP']                   = writelastLP    


    all_data['duals']                         = {}
    all_data['dual_diff']                     = {}

    casetype = ''

    if all_data['nperturb']:
        casetype = 'gaussian'
        log.joint('case ' + all_data['casename'] + ' multi-period ' + str(T)
                  + ' ' + casetype + '\n')
    elif all_data['uniform']:
        casetype = 'uniform_'+str(uniform_drift)
        log.joint('case ' + all_data['casename'] + ' multi-period ' + str(T)
                  + ' ' + casetype + '\n')
    elif all_data['uniform2']:
        casetype = 'uniform2'
        log.joint('case ' + all_data['casename'] + ' multi-period ' + str(T)
                  + ' ' + casetype + '\n')
    elif all_data['uniform3']:
        casetype = 'uniform3'
        log.joint('case ' + all_data['casename'] + ' multi-period ' + str(T)
                  + ' ' + casetype + '\n')
    elif all_data['uniform4']:
        casetype = 'uniform4'
        log.joint('case ' + all_data['casename'] + ' multi-period ' + str(T)
                  + ' ' + casetype + '\n')
    elif all_data['uniform5']:
        casetype = 'uniform5'
        log.joint('case ' + all_data['casename'] + ' multi-period ' + str(T)
                  + ' ' + casetype + '\n')
    elif all_data['uniform6']:
        casetype = 'uniform6'
        log.joint('case ' + all_data['casename'] + ' multi-period ' + str(T)
                  + ' ' + casetype + '\n')                
    elif all_data['pglib_reverse']:
        casetype = 'pglib_reverse'
        log.joint('case ' + all_data['casename'] + ' multi-period ' + str(T)
                  + ' ' + casetype + '\n')
    elif all_data['arpae']:
        casetype = 'arpae'
        log.joint('case ' + all_data['casename'] + ' multi-period ' + str(T)
                  + ' ' + casetype + '\n')                        
    elif all_data['arpae2']:
        casetype = 'arpae2'
        log.joint('case ' + all_data['casename'] + ' multi-period ' + str(T)
                  + ' ' + casetype + '\n')                
        
    all_data['casetype'] = casetype

    if all_data['nperturb']:
        loadsfilename = '../../ampl_aopcf/data/mtploads/' + casename + '_mtploads_' + str(T) + '_n1.txt'
    elif all_data['uniform']:
        loadsfilename = '../../ampl_acopf/data/mtploads/' + casename + '_mtploads_' + str(T) + '_u' + str(all_data['uniform_drift']) + '.txt'
    elif all_data['uniform2']:
        loadsfilename = '../../ampl_acopf/data/mtploads/' + casename + '_mtploads_' + str(T) + '_u2.txt'
    elif all_data['uniform3']:
        loadsfilename = '../../ampl_acopf/data/mtploads/' + casename + '_mtploads_' + str(T) + '_u3.txt'
    elif all_data['uniform4']:
        loadsfilename = '../../ampl_acopf/data/mtploads/' + casename + '_mtploads_' + str(T) + '_u4.txt'
    elif all_data['uniform5']:
        loadsfilename = '../../ampl_acopf/data/mtploads/' + casename + '_mtploads_' + str(T) + '_u5_' + str(all_data['uniform_drift']) + '.txt'
    elif all_data['uniform6']:
        loadsfilename = '../../ampl_acopf/data/mtploads/' + casename + '_mtploads_' + str(T) + '_u6_' + str(all_data['uniform_drift']) + '.txt'           
    elif all_data['pglib_reverse']:
        loadsfilename = '../data/mtploads/' + casename + '_mtploads_' + str(T) + '_' + str(all_data['uniform_drift']) + '_pglib.txt'
    elif all_data['arpae']:
        loadsfilename = '../../ampl_acopf/data/mtploads/' + casename + '_mtploads_' + str(T) + '_arpae.txt'                        
    elif all_data['arpae2']:
        loadsfilename = '../../ampl_acopf/data/mtploads/' + casename + '_mtploads_' + str(T) + '_arpae2.txt'                

    all_data['loadsfilename'] = loadsfilename
        
    return all_data
        

if __name__ == '__main__':
    if len(sys.argv) > 4:
        print ('Usage: main_mtp.py file.conf [sols] [logfile]\n')
        exit(0)

    T0 = time.time()
    
    if len(sys.argv) == 4:
        sols      = sys.argv[2] + '/'
        mylogfile = sys.argv[3]
    else:
        sols      = ""
        mylogfile = "main.log"

    
    log = danoLogger(mylogfile)
        
    log.joint('\n  *********************************************************\n')
    log.joint(' ***********************************************************\n')
    log.joint(' ****                                                   ****\n')
    log.joint(' ****    Initializing AC-OPF Cutting-plane algorithm    ****\n')
    log.joint(' ****                                                   ****\n')
    log.joint(' ***********************************************************\n')
    log.joint('  *********************************************************\n\n')

    stateversion(log)
    
    all_data              = read_config(log,sys.argv[1])
    all_data['T0']        = T0
    all_data['sols']      = sols
    all_data['mylogfile'] = mylogfile
    all_data['datetime']  = mylogfile.strip("CPexp_").strip(".log")

    readcode = reader.readcase(log,all_data,all_data['casefilename'])

    gocutplane(log,all_data)
    
    log.closelog()
    
