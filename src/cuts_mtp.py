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


from myutils import breakexit
from log import danoLogger
import time
import math
from json import dumps
from gurobipy import *
import numpy as np
from numpy import linalg as LA


def compute_normal(*args):
    n = len(args)
    v = np.zeros(n,dtype = 'float')
    for i in range(n):
        v[i] = args[i]

    norm = LA.norm(v)

    if norm == 0:
        log.joint(' BUG: check vector (cut)\n')
        brekaexit('check cut')
    return v / norm


def compute_coeffs_cutnorm(cft,sft,cff,ctt):

    v = np.zeros(4,dtype = 'float')
    
    cutnorm   = math.sqrt( (2 * cft)**2 + (2 * sft)**2 + (cff - ctt)**2 )
    coeff_cft = 4 * cft
    coeff_sft = 4 * sft
    coeff_cff = cff - ctt - cutnorm
    coeff_ctt = - (cff - ctt) - cutnorm

    return v,cutnorm

# def new_cut(coeff_cft,coeff_sft,coeff_cff,coef__ctt):

#     #cutexp += coeff_cft * cvar[branch] + coeff_sft * svar[branch] + coeff_cff * cvar[buses[c\
# ount_of_f]] + coeff_ctt * cvar[buses[count_of_t]]

#     #if coeff_cft
#     #cutnorm = math.sqrt( (2 * cft)**2 + (2 * sft)**2 + (cff - ctt)**2 )
#     coeff_cft = 4 * cft
#     coeff_sft = 4 * sft
#     coeff_cff = cff - ctt - cutnorm
#     coeff_ctt = - (cff - ctt) - cutnorm

#     return v,cutnorm

def loss_cuts(log,all_data):

    log.joint('\n')
    log.joint(' **** Loss-cuts ****\n')

    themodel     = all_data['themodel']
    IDtoCountmap = all_data['IDtoCountmap']
    branches     = all_data['branches']
    Pvar_f       = all_data['Pvar_f']
    Pvar_t       = all_data['Pvar_t']
    Pfvalues     = all_data['Pfvalues']
    Ptvalues     = all_data['Ptvalues']

    # cvar         = all_data['cvar']
    # svar         = all_data['svar']
    # cvalues      = all_data['cvalues']
    # svalues      = all_data['svalues']

    #add with care, once a loss is added, is it possible to be violated again? (feasibility tolerance..!)
    #which is better, loss ineqs with fundamental vars or flow vars?
    
    rnd                  = all_data['round']
    loss_cuts            = all_data['loss_cuts']
    threshold            = all_data['threshold']
    violated             = {}
    violated_count       = 0
    
    FeasibilityTol       = all_data['FeasibilityTol']

    #loss-analysis
    #loss_cuts_info = all_data['loss_cuts_info']
    #loss_cuts_info_updated = all_data['loss_cuts_info_updated']
    
    all_data['NO_loss_violated'] = 0
    
    log.joint(' checking for violations of active loss inequalities ... \n')
    for branch in branches.values():

        if branch.r < 0: #!!!
            continue 

        f = branch.f
        t = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]
        violation = - ( Pfvalues[branch] + Ptvalues[branch] )
        if violation > threshold:
            violated_count  += 1
            violated[branch] = violation

    if violated_count == 0:
        all_data['NO_loss_violated'] = 1
        log.joint(' all violations below threshold\n' )
        log.joint(' no more cuts to add for current threshold\n' )
        return None
    
    log.joint('  number violated loss inequalities ' + str(violated_count) + '\n')

    log.joint(' sorting most violated loss inequalities ... \n')

    num_selected        =  math.ceil(violated_count * all_data['most_violated_fraction_loss'] )
    most_violated       = dict(sorted(violated.items(), key = lambda x: x[1], reverse = True)[:num_selected])
    most_violated_count = 0

    log.joint(' adding loss inequalities ... \n')

    for branch in most_violated.keys():

        branchcount = branch.count
    
        if (most_violated_count == 0):
            most_violated_branch = branchcount
            max_error = most_violated[branch]

        most_violated_count += 1
        violation            = most_violated[branch]
        f                    = branch.f
        t                    = branch.t
        Pf                   = Pfvalues[branch]
        Pt                   = Ptvalues[branch]

        if all_data['loud_cuts']:
            log.joint('  --> new loss inequality\n')
            log.joint('  branch ' + str(branchcount) + ' f ' + str(f) + ' t ' 
                      + str(t) + ' violation ' + str(violation) + '\n')
            log.joint('  values ' + ' Pf ' + str(Pf) + ' Pt ' + str(Pt) + '\n' )

        #sanity check
        if all_data['loss_validity']:
            
            sol_Pf     = all_data['sol_Pfvalues'][branch]
            sol_Pt     = all_data['sol_Ptvalues'][branch]
            violation = - (sol_Pf + sol_Pt)
         
            if violation > FeasibilityTol:
                log.joint('  WARNING, the loss inequality associated to branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' is violated by the AC solution!\n')
                log.joint('  violation ' + str(violation) + '\n')
                log.joint('  values (AC solution) ' + ' Pf ' + str(sol_Pf) + ' Pt ' + str(sol_Pt) + '\n')
                breakexit('check!')
            else:
                log.joint('  AC solution satisfies loss inequality at branch ' + str(branch.count) + ' with slack ' + str(violation) + '\n')
                 
        loss_cuts[branch] = (rnd,violation,threshold)        
        lossexp = LinExpr()
        constrname = "loss_ineq_"+str(branchcount)+"_"+str(f)+"_"+str(t)
        lossexp += Pvar_f[branch] + Pvar_t[branch]
        themodel.addConstr(lossexp >= 0, name = constrname)
        
    log.joint('  number loss-ineq added ' + str(most_violated_count) + '\n')
    log.joint('  max error (abs) ' + str(max_error) + ' at branch ' + str(most_violated_branch) + '\n' )

    all_data['most_violated_branch_loss'] = most_violated_branch
    all_data['max_error_loss']            = max_error
    all_data['num_loss_cuts']            += most_violated_count
    all_data['num_loss_cuts_added']       = most_violated_count

    if all_data['droploss'] and (all_data['round'] >= all_data['cut_age_limit']):
        t0_drop = time.time()
        drop_loss(log,all_data)
        t1_drop = time.time()

        log.joint('  time spent on drop loss ' + str(t1_drop - t0_drop) + '\n')


def drop_loss(log,all_data):
    
    #log.joint('\n')
    #log.joint(' **** drop Loss cuts ****\n')
    log.joint(' dropping old and slack loss inequalities ...\n')

    themodel      = all_data['themodel']
    branches      = all_data['branches']
    buses         = all_data['buses']
    IDtoCountmap  = all_data['IDtoCountmap']
    cut_age_limit = all_data['cut_age_limit']
    threshold     = all_data['threshold']
    rnd           = all_data['round']
    loss_cuts     = all_data['loss_cuts']
    dropped_loss  = all_data['dropped_loss']
    
    drop_loss     = []

    for branch in loss_cuts.keys():
        cut         = loss_cuts[branch]
        cut_rnd     = cut[0]
        cut_age     = rnd - cut_rnd

        if (cut_age == 0):
            log.joint(' branch ' + str(branch.count) + 'skipped \n')
            continue #cut new born, skip                                                                                    
        branchcount = branch.count 
        f           = branch.f
        t           = branch.t
        constrname  = "loss_ineq_"+str(branchcount)+"_"+str(f)+"_"+str(t)
        constr      = themodel.getConstrByName(constrname)
        slack       = constr.getAttr("slack")

        if ( slack < - threshold ) and (cut_age > cut_age_limit):
            drop_loss.append(branch)
            themodel.remove(themodel.getConstrByName(constrname))
            if all_data['loud_cuts']:
                log.joint('  --> removed loss-cut\n')
                log.joint('  the cut ' + str(key)
                          + ' was removed from the model\n')
            
            
    num_drop_loss = len(drop_loss)
    
    if num_drop_loss:
        all_data['num_loss_cuts_dropped'] = num_drop_loss
        all_data['num_loss_cuts']        -= num_drop_loss

        dropped_loss.extend(drop_loss)
        for branch in drop_loss:
            dropped_loss.remove(branch)
            loss_cuts.pop(branch)

        all_data['dropped_loss'] = dropped_loss
        all_data['loss_cuts']    = loss_cuts
        log.joint('\n the loss cuts in drop_loss list were removed from dict loss_cuts\n' )

    else:
        all_data['num_loss_cuts_dropped'] = 0
        log.joint('\n no loss cuts were dropped this round\n')        


def i2_cuts(log,all_data):
        
    log.joint('\n')
    log.joint(' **** i2-cuts ****\n')

    themodel       = all_data['themodel']
    IDtoCountmap   = all_data['IDtoCountmap']
    buses          = all_data['buses']
    branches       = all_data['branches']
    FeasibilityTol = all_data['FeasibilityTol']
    T              = all_data['T']
    
    cvar      = all_data['cvar']
    Pvar_f    = all_data['Pvar_f']
    Qvar_f    = all_data['Qvar_f']
    i2var_f   = all_data['i2var_f']

    #### CHANGE Do this once!
    alphadic  = all_data['alphadic']
    goodi2s   = {}
    for branch in branches.values():
        if alphadic[branch] < all_data['alpha_threshold']:
            goodi2s[branch.count] = branch

    #####
    
    Pfvalues  = all_data['Pfvalues']
    Qfvalues  = all_data['Qfvalues']
    cvalues   = all_data['cvalues']
    i2fvalues = all_data['i2fvalues']

    rnd                      = all_data['round']
    i2_cuts                  = all_data['i2_cuts']
    i2_cuts_info             = all_data['i2_cuts_info']
    num_cuts                 = all_data['ID_i2_cuts']
    threshold                = all_data['threshold_i2']
    violated                 = {}
    violated_count           = 0
    num_i2_cuts_added        = 0

    
    
    # New code
    most_violated = {}
    
    for k in range(T):
        violated[k] = {}

    ##
        
    all_data['NO_i2_cuts_violated'] = 0
    
    log.joint(' checking for violations of i2 inequalities ... \n')
    for branch in goodi2s.values():
        
        f = branch.f
        t = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]

        for k in range(T):
            violation = Pfvalues[k][branch]*Pfvalues[k][branch] + Qfvalues[k][branch]*Qfvalues[k][branch] - cvalues[k][buses[count_of_f]] * i2fvalues[k][branch]
            #violation_f = Pfvalues[branch]*Pfvalues[branch] + Qfvalues[branch]*Qfvalues[branch] - cvalues[buses[count_of_f]] * i2fvalues[branch]
            #violation_t = Ptvalues[branch]*Ptvalues[branch] + Qtvalues[branch]*Qtvalues[branch] - cvalues[buses[count_of_t]] * i2tvalues[branch]
            #violation = max(violation_f,violation_t)
            if violation > threshold:
                if all_data['loud_cuts']:
                    log.joint('  violation ' + str(violation) + ' at branch '
                              + str(branch.count) + ' time-period ' + str(k)  + ' f '
                              + str(f) + ' t ' + str(t) + ' i2 value '
                              + str(i2fvalues[k][branch]) + '\n')
                violated_count  += 1
                violated[k][branch] = violation

    if violated_count == 0:
        all_data['NO_i2_cuts_violated'] = 1
        log.joint(' all i2 violations below threshold\n' )
        log.joint(' no more i2 cuts to add for current threshold\n' )
        return None

    log.joint('  number violated i2 ineqs ' + str(violated_count) + '\n')

    log.joint(' sorting most violated i2-envelope cuts ...\n')

    num_selected        =  math.ceil(violated_count * all_data['most_violated_fraction_i2'] )

    for k in range(T):
        most_violated[k] = dict(sorted(violated[k].items(),
                                       key = lambda x: x[1],
                                       reverse = True)[:num_selected])

    log.joint(' computing i2-envelope cuts ... \n')

    most_violated_count  = 0
    most_violated_branch = 'none'
    max_error            = -1

    for k in range(T):
        numkcuts = 0
        for branch in most_violated[k].keys():
            if (most_violated_count == 0):
                if most_violated[k][branch] > max_error:
                    most_violated_branch = 'branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' k ' + str(k)
                    max_error = most_violated[k][branch]

                
            f                    = branch.f
            t                    = branch.t
            count_of_f           = IDtoCountmap[f]
            count_of_t           = IDtoCountmap[t]
            Pft                  = Pfvalues[k][branch]
            Qft                  = Qfvalues[k][branch]
            cff                  = cvalues[k][buses[count_of_f]]
            i2ft                 = i2fvalues[k][branch]

            cutnorm    = math.sqrt( (2 * Pft)**2 + (2 * Qft)**2 + (cff - i2ft)**2 )
            coeff_Pft  = 4 * Pft
            coeff_Qft  = 4 * Qft
            coeff_cff  = cff - i2ft - cutnorm
            coeff_i2ft = - (cff - i2ft) - cutnorm

            if all_data['parallel_check']:
                if parallel_check_i2(log,all_data,branch,coeff_Pft,coeff_Qft,
                                     coeff_cff,coeff_i2ft,k):
                    continue

            most_violated_count += 1
            violation            = most_violated[k][branch]
            cutid                = num_cuts + most_violated_count

            if all_data['loud_cuts']:
                log.joint('  --> new i2-cut\n')
                log.joint('  branch ' + str(branch.count) + ' time-period ' + str(k)
                          + ' f ' + str(f) + ' t ' + str(t) + ' violation '
                          + str(violation) + ' cut id ' + str(cutid) + '\n' )
                log.joint('  values ' + ' Pft ' + str(Pft) + ' Qft ' + str(Qft) 
                          + ' cff ' + str(cff) + ' i2ft ' + str(i2ft) + '\n' )
                log.joint('  LHS coeff ' + ' Pft ' + str(coeff_Pft) + ' Qft ' 
                          + str(coeff_Qft) + ' cff ' + str(coeff_cff) + ' i2ft ' 
                          + str(coeff_i2ft) + '\n' )
                log.joint('  cutnorm ' + str(cutnorm) + '\n')

            # Sanity check
            if all_data['i2_validity']:
            
                sol_Pf        = all_data['sol_Pfvalues'][k][branch]
                sol_Qf        = all_data['sol_Qfvalues'][k][branch]
                sol_c         = all_data['sol_cvalues'][k][branch]
                sol_s         = all_data['sol_svalues'][k][branch]
                sol_cbusf     = all_data['sol_cvalues'][k][buses[count_of_f]]
                sol_cbust     = all_data['sol_cvalues'][k][buses[count_of_t]]
                sol_i2f       = computei2value(log,all_data,branch,sol_c,sol_s,sol_cbusf,sol_cbust)
                violation     = coeff_Pft * sol_Pf + coeff_Qft * sol_Qf + coeff_cff * sol_cbusf + coeff_i2ft * sol_i2f
                relviolation  = violation / ( ( coeff_Pft**2 + coeff_Qft**2 + coeff_cff**2 + coeff_i2ft**2 )**0.5 ) 

                if relviolation > FeasibilityTol:
                    log.joint('  WARNING, the loss inequality associated to branch '
                              + str(branch.count) + ' time-period ' + str(k) +' f '
                              + str(branch.f) + ' t ' + str(branch.t)
                              + ' is violated by the AC solution!\n')
                    log.joint('  violation ' + str(violation) + '\n')
                    log.joint('  relative violation ' + str(relviolation) + '\n')
                    log.joint('  values (AC solution) ' + ' Pft ' + str(sol_Pf) + ' Qft ' + str(sol_Qf) + ' cff ' + str(sol_cbusf) + ' i2ft ' + str(sol_i2f) + '\n' )
                    breakexit('check!')
                else:
                    log.joint('  AC solution satisfies loss inequality at branch '
                              + str(branch.count) + ' with slack ' + str(violation) + '\n')
        
        
            i2_cuts[(cutid,branch.count)] = (rnd,threshold,k)
            i2_cuts_info[k][branch][cutid]   = (rnd,violation,coeff_Pft,coeff_Qft,coeff_cff,coeff_i2ft,threshold,cutid)

        
            cutexp     = LinExpr()
            constrname = "i2_cut_"+str(cutid)+"_"+str(branch.count)+"r_"+str(rnd)+"k_"+str(k)+"_"+str(f)+"_"+str(t)
            cutexp    += coeff_Pft * Pvar_f[k][branch] + coeff_Qft * Qvar_f[k][branch] + coeff_cff * cvar[k][buses[count_of_f]] + coeff_i2ft * i2var_f[k][branch]

            themodel.addConstr(cutexp <= 0, name = constrname)
            numkcuts += 1

        log.joint('  time-period = ' + str(k) + ' : i2-envelope cuts'
                  + ' added '  + str(numkcuts) + '\n')
    log.joint('  number i2-envelope cuts added ' + str(most_violated_count)
              + '\n')
    log.joint('  max error (abs) ' + str(max_error) + ' at '
              + str(most_violated_branch) + '\n' )


    all_data['most_violated_branch_i2'] = most_violated_branch
    all_data['max_error_i2']            = max_error
    all_data['ID_i2_cuts']             += most_violated_count
    all_data['num_i2_cuts_added']       = most_violated_count
    all_data['num_i2_cuts']            += most_violated_count
    all_data['num_i2_cuts_rnd'][rnd]    = most_violated_count


    if all_data['dropi2']:
        if all_data['addcuts']:
            t0_drop = time.time()
            drop_i2(log,all_data)
            t1_drop = time.time()
            log.joint('  time spent on drop i2 ' + str(t1_drop - t0_drop) + '\n')
        elif all_data['round'] >= all_data['cut_age_limit']: ######
            t0_drop = time.time()
            drop_i2(log,all_data)
            t1_drop = time.time()
            log.joint('  time spent on drop i2 ' + str(t1_drop - t0_drop) + '\n')


def drop_i2(log,all_data):
    
    #log.joint('\n')
    #log.joint(' **** drop i2-envelope cuts ****\n')
    log.joint(' dropping old and slack i2-envelope cuts ...\n')

    themodel             = all_data['themodel']
    branches             = all_data['branches']
    buses                = all_data['buses']
    IDtoCountmap         = all_data['IDtoCountmap']
    current_rnd          = all_data['round']
    cut_age_limit        = all_data['cut_age_limit']
    i2_cuts              = all_data['i2_cuts']
    i2_cuts_info         = all_data['i2_cuts_info']
    dropped_i2           = all_data['dropped_i2']
    drop_i2              = []
    T                    = all_data['T']
    num_i2_cuts_dropped  = all_data['num_i2_cuts_dropped'] 

    for key in i2_cuts.keys():  #key = (cutid,branch.count), value = (rnd,threshold)
        cut     = i2_cuts[key] 
        cut_rnd = cut[0]
        cut_age = current_rnd - cut_rnd

        if cut_age <= cut_age_limit: ###### fix
            continue #baby cuts, skip!
        
        cutid         = key[0]
        branchid      = key[1]
        branch        = branches[branchid]
        f             = branch.f
        t             = branch.t
        count_of_f    = IDtoCountmap[f]
        count_of_t    = IDtoCountmap[t]
        cut_threshold = cut[1]
        k             = cut[2]
        
        constrname    = "i2_cut_"+str(cutid)+"_"+str(branchid)+"r_"+str(cut_rnd)+"k_"+str(k)+"_"+str(f)+"_"+str(t)
        constr        = themodel.getConstrByName(constrname)
        slack         = - constr.getAttr("slack")

        if ( slack < - cut_threshold ):
            drop_i2.append(key)
            themodel.remove(themodel.getConstrByName(constrname))
            #if all_data['loud_cuts']:
            #    log.joint('  --> removed i2-cut\n')
            #    log.joint('  cut ' + str(key) + ' branch ' 
            #              + str(branchid) + ' f ' + str(f) + ' t ' + str(t) 
            #              + ' was removed from the model\n')
    
    num_drop_i2       = len(drop_i2)
    
    if num_drop_i2:
        all_data['num_i2_cuts_dropped'] = num_drop_i2
        all_data['num_i2_cuts']        -= num_drop_i2
        all_data['total_i2_dropped']   += num_drop_i2

        dropped_i2.extend(drop_i2)

        
        ##### drop i2-envelope cuts from i2_cuts_info
        for key in drop_i2:
            cutid    = key[0]
            branchid = key[1]
            k        = i2_cuts[key][2] 
            cuts_branch = i2_cuts_info[k][branches[branchid]]
            if all_data['loud_cuts']:
                log.joint('  --> i2-cut removed\n')
                log.joint('  cutid ' + str(cutid) + ' branchid '
                          + str(branchid) + ' k ' + str(k) + '\n')
                log.joint(' coeff: Pft ' + str(cuts_branch[cutid][2])
                          + ' Qft ' + str(cuts_branch[cutid][3]) + ' cff '
                          + str(cuts_branch[cutid][4]) + ' i2ft '
                          + str(cuts_branch[cutid][5]) + ' rnd '
                          + str(cuts_branch[cutid][0]) + '\n')            
            cuts_branch.pop(cutid)
        log.joint('  cuts in drop_i2 list were removed from i2_cuts_info\n')

        ##### drop i2-envelope cuts from i2_cuts        
        for key in drop_i2:
            i2_cuts.pop(key)
        log.joint('  the cuts in drop_i2 list were removed from dict i2_cuts\n')
        
    else:
        all_data['num_i2_cuts_dropped'] = 0
        log.joint('  no i2-envelope cuts were dropped this round\n')
    

def limit_cuts(log,all_data):
        
    log.joint('\n')
    log.joint(' **** Limit-cuts ****\n')  #should we add f->t and t->f? if r small, |Pf| = |Pt| for sure, but Qf and Qt might differ 

    themodel         = all_data['themodel']
    IDtoCountmap     = all_data['IDtoCountmap']
    buses            = all_data['buses']
    branches         = all_data['branches']
    Pvar_f           = all_data['Pvar_f']
    Qvar_f           = all_data['Qvar_f']
    Pvar_t           = all_data['Pvar_t']
    Qvar_t           = all_data['Qvar_t']
    FeasibilityTol   = all_data['FeasibilityTol']
    T                = all_data['T']
    
    Pfvalues = all_data['Pfvalues']
    Qfvalues = all_data['Qfvalues']
    Ptvalues = all_data['Ptvalues']
    Qtvalues = all_data['Qtvalues']
    
    rnd                   = all_data['round']
    limit_cuts            = all_data['limit_cuts']    
    limit_cuts_info       = all_data['limit_cuts_info']
    num_cuts              = all_data['ID_limit_cuts']
    num_limit_cuts_added  = 0
    threshold             = all_data['threshold']
    violated              = {}
    violated_count        = 0

    # New code
    most_violated = {}

    for k in range(T):
        violated[k] = {}

    ##
    
    all_data['NO_limit_cuts_violated'] = 0
    
    log.joint(' checking for violations of limit inequalities ... \n')
    for branch in branches.values():
        for k in range(T):
            violation_f = Pfvalues[k][branch]*Pfvalues[k][branch] + Qfvalues[k][branch]*Qfvalues[k][branch] - branch.limit**2
            violation_t = Ptvalues[k][branch]*Ptvalues[k][branch] + Qtvalues[k][branch]*Qtvalues[k][branch] - branch.limit**2
            if violation_f > threshold:
                violated_count += 1
                violated[k][branch] = (violation_f,'f')
            elif violation_t > threshold:
                violated_count += 1
                violated[k][branch] = (violation_t,'t')    

    if violated_count == 0:
        all_data['NO_limit_cuts_violated'] = 1
        log.joint(' all limit violations below threshold\n' )
        log.joint(' no more limit-cuts to add for current threshold\n' )
        return None
    
    log.joint('  number violated limits ' + str(violated_count) + '\n')
    
    log.joint(' sorting most violated limit-envelope cuts ...\n')

    num_selected =  math.ceil(violated_count * all_data['most_violated_fraction_limit'] )

    for k in range(T):
        most_violated[k] = dict(sorted(violated[k].items(),
                                       key = lambda x: x[1][0],
                                       reverse = True)[:num_selected])

    log.joint(' computing limit-envelope cuts ... \n')

    most_violated_count  = 0
    most_violated_branch = 'none'
    max_error            = -1

    for k in range(T):
        numkcuts = 0
        for branch in most_violated[k].keys():
            if (most_violated_count == 0):
                if most_violated[k][branch][0] > max_error:
                    most_violated_branch = 'branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' k ' + str(k)
                    max_error = most_violated[k][branch][0]
  

            from_or_to = most_violated[k][branch][1]
    
            if from_or_to == 'f':
                Pval                  = Pfvalues[k][branch]
                Qval                  = Qfvalues[k][branch]
            elif from_or_to == 't':
                Pval                 = Ptvalues[k][branch]
                Qval                 = Qtvalues[k][branch]
            else:
                log.joint(' we have a bug\n')
                breakexit('look for bug')

            f          = branch.f
            t          = branch.t
            count_of_f = IDtoCountmap[f]
            count_of_t = IDtoCountmap[t]
            violation  = most_violated[k][branch][0]
            u          = branch.limit
            u2         = u**2

            if violation + u2 < 1e-05:
                log.joint('  limit cuts: check branch\n') 
                breakexit('check')  

            t0         = 1 / (u * math.sqrt(violation + u2))
            coeff_P    = t0 * Pval
            coeff_Q    = t0 * Qval
            z          = 1

            if all_data['parallel_check']:
                if parallel_check_limit(log,all_data,branch,coeff_P,coeff_Q,
                                        from_or_to,k):
                    continue

            # We add the cut
            most_violated_count += 1
            cutid                = num_cuts + most_violated_count

            if all_data['loud_cuts']:
                log.joint('  --> new cut\n')
                log.joint('  branch ' + str(branch.count) + ' time-period '
                          + str(k) + ' f ' + str(f) + ' t '  + str(t)
                          + ' violation ' + str(violation) + ' cut id '
                          + str(cutid) + '\n')
                if from_or_to == 'f':
                    log.joint('  values ' + ' Pft ' + str(Pval) + ' Qft ' 
                              + str(Qval) + '\n')
                    log.joint('  LHS coeff ' + ' Pft ' + str(coeff_P) + ' Qft '
                              + str(coeff_Q) + ' RHS ' + str(z) + '\n')
                elif from_or_to == 't':
                    log.joint('  values ' + ' Ptf ' + str(Pval) + ' Qtf ' 
                              + str(Qval) + '\n')
                    log.joint('  LHS coeff ' + ' Ptf ' + str(coeff_P) + ' Qtf '
                              + str(coeff_Q) + ' RHS ' + str(z) + '\n')
        
            #sanity check
            if all_data['limit_validity']:
            
                if from_or_to == 'f':
                    sol_Pval           = all_data['sol_Pfvalues'][k][branch]
                    sol_Qval           = all_data['sol_Qfvalues'][k][branch]
                elif from_or_to == 't':
                    sol_Pval           = all_data['sol_Ptvalues'][k][branch]
                    sol_Qval           = all_data['sol_Qtvalues'][k][branch]

                slack    = coeff_P * sol_Pval + coeff_Q * sol_Qval - z

                if slack > FeasibilityTol:
                    log.joint('  this cut is not valid!\n')
                    log.joint('  branch ' + str(branch.count) + ' time-period '
                          + str(k) + ' f ' + str(f) + ' t '  + str(t)
                          + ' cut id ' + str(cutid) + '\n')
                    log.joint('  violation ' + str(slack) + '\n')
                    log.joint('  values (a primal bound)' + ' P '
                              + str(sol_Pval) + ' Q ' + str(sol_Qval)
                              + ' branch limit ' + str(u) + '\n')
                    breakexit('check!')
                else:
                    log.joint('  valid cut at branch ' + str(branch.count)
                              + ' with slack ' + str(slack) + '\n')
        
            limit_cuts[(cutid,branch.count)]  = (rnd,threshold,from_or_to,k)
            limit_cuts_info[k][branch][cutid] = (rnd,violation,coeff_P,coeff_Q,threshold,cutid,from_or_to)
        
            cutexp = LinExpr()

            if from_or_to == 'f':
                constrname = "limit_cut_" + str(cutid) + "_" + str(branch.count) + "r_" + str(rnd) + "k_" + str(k) + "_" + str(f) + "_" + str(t)
                cutexp += coeff_P * Pvar_f[k][branch] + coeff_Q * Qvar_f[k][branch]
            elif from_or_to == 't':
                constrname = "limit_cut_" + str(cutid) + "_" + str(branch.count) + "r_" + str(rnd) + "k_" + str(k) + "_" + str(t) + "_" + str(f)
                cutexp += coeff_P * Pvar_t[k][branch] + coeff_Q * Qvar_t[k][branch]
            else:
                log.joint('  we have a bug\n')
                breakexit('look for bug')
            
            themodel.addConstr(cutexp <= z, name = constrname)
            numkcuts += 1
        log.joint('  time-period = ' + str(k) + ' : limit-envelope cuts'
                  + ' added ' + str(numkcuts) + '\n')
    log.joint('  number limit-envelope cuts added '
              + str(most_violated_count) + '\n')
    log.joint('  max error (abs) ' + str(max_error) + ' at '
              + str(most_violated_branch) + '\n' )

    all_data['most_violated_branch_limit']   = most_violated_branch
    all_data['max_error_limit']              = max_error
    all_data['ID_limit_cuts']               += most_violated_count
    all_data['num_limit_cuts_added']         = most_violated_count
    all_data['num_limit_cuts']              += most_violated_count
    all_data['num_limit_cuts_rnd'][rnd]      = most_violated_count


    if all_data['droplimit']:
        if all_data['addcuts']:
            t0_drop = time.time()
            drop_limit(log,all_data)
            t1_drop = time.time()
            log.joint('  time spent on drop limit ' + str(t1_drop - t0_drop) + '\n')
        elif all_data['round'] >= all_data['cut_age_limit']: ######
            t0_drop = time.time()
            drop_limit(log,all_data)
            t1_drop = time.time()
            log.joint('  time spent on drop limit ' + str(t1_drop - t0_drop) + '\n')


def drop_limit(log,all_data):
    
    #log.joint('\n')
    #log.joint(' **** drop Limit-envelope cuts ****\n')
    log.joint(' dropping old and slack limit-envelope cuts ...\n')

    themodel                = all_data['themodel']
    branches                = all_data['branches']
    buses                   = all_data['buses']
    IDtoCountmap            = all_data['IDtoCountmap']
    current_rnd             = all_data['round']
    cut_age_limit           = all_data['cut_age_limit']
    limit_cuts              = all_data['limit_cuts']
    limit_cuts_info         = all_data['limit_cuts_info']
    drop_limit              = []
    dropped_limit           = all_data['dropped_limit']
    num_limit_cuts_dropped  = all_data['num_limit_cuts_dropped']
    T                       = all_data['T']

    for key in limit_cuts.keys():
        cut     = limit_cuts[key]   #        limit_cuts[(cutid,branch.count)] = (rnd,threshold,from_or_to)
        cut_rnd = cut[0]
        cut_age = current_rnd - cut_rnd

        if cut_age <= cut_age_limit: ###fix this ws
            continue #baby cuts, skip!
        
        cutid         = key[0]
        branchid      = key[1]
        branch        = branches[branchid]
        f             = branch.f
        t             = branch.t
        count_of_f    = IDtoCountmap[f]
        count_of_t    = IDtoCountmap[t]
        cut_threshold = cut[1]
        from_or_to    = cut[2]
        k             = cut[3]

        if from_or_to == 'f':
            constrname = "limit_cut_" + str(cutid) + "_" + str(branch.count) + "r_" + str(cut_rnd) + "k_" + str(k) + "_" + str(f) + "_" + str(t)
        elif from_or_to == 't':
            constrname = "limit_cut_" + str(cutid) + "_" + str(branch.count) + "r_" + str(cut_rnd) + "k_" + str(k) + "_" + str(t) + "_" + str(f)
        else:
            log.joint(' we have a bug\n')
            breakexit('look for bug')

        constr        = themodel.getConstrByName(constrname)
        slack         = - constr.getAttr("slack")

        if ( slack < - cut_threshold ):
            drop_limit.append(key)
            themodel.remove(themodel.getConstrByName(constrname))
            #if all_data['loud_cuts']:
            #    log.joint(' --> removed limit-cut\n')
            #    log.joint(' the cut ' + str(key)
            #              + ' was removed from the model\n')
    
    num_drop_limit = len(drop_limit)
    if num_drop_limit:
        all_data['num_limit_cuts_dropped'] = num_drop_limit
        all_data['num_limit_cuts']        -= num_drop_limit
        all_data['total_limit_dropped']   += num_drop_limit

        dropped_limit.extend(drop_limit)

        ##### drop limit-envelope cuts from limit_cuts_info
        for key in drop_limit:
            cutid      = key[0]
            branchid   = key[1]
            from_or_to = limit_cuts[key][2]
            k          = limit_cuts[key][3]
            cuts_branch = limit_cuts_info[k][branches[branchid]]
            if all_data['loud_cuts']:
                log.joint('  --> limit-cut removed\n')
                log.joint('  cutid ' + str(cutid) + ' branchid '
                          + str(branchid) + ' k ' + str(k) + '\n')
                if from_or_to == 'f':
                    log.joint(' coeff: Pft ' + str(cuts_branch[cutid][2])
                              + ' Qft ' + str(cuts_branch[cutid][3])
                              + ' rnd ' + str(cuts_branch[cutid][0]) + '\n')
                elif from_or_to == 't':
                    log.joint(' coeff: Ptf ' + str(cuts_branch[cutid][2])
                              + ' Qtf ' + str(cuts_branch[cutid][3])
                              + ' rnd ' + str(cuts_branch[cutid][0]) + '\n') 
            cuts_branch.pop(cutid)
        log.joint('  cuts in drop_limit list were removed from limit_cuts_info\n')

        ##### drop limit-envelope cuts from limit_cuts        
        for key in drop_limit:
            limit_cuts.pop(key)

        log.joint('  the cuts in drop_limit list were removed from dict limit_cuts\n' )        
    else:
        all_data['num_limit_cuts_dropped'] = 0
        log.joint('  no limit-envelope cuts were dropped this round\n')


def jabr_cuts(log,all_data):
        
    log.joint('\n')
    log.joint(' **** Jabr-cuts ****\n')

    themodel       = all_data['themodel']
    IDtoCountmap   = all_data['IDtoCountmap']
    buses          = all_data['buses']
    branches       = all_data['branches']
    cvar           = all_data['cvar']
    svar           = all_data['svar']
    FeasibilityTol = all_data['FeasibilityTol']
    cvalues        = all_data['cvalues']
    svalues        = all_data['svalues']
    T              = all_data['T']
    
    rnd                 = all_data['round']
    jabr_cuts           = all_data['jabr_cuts']
    jabr_cuts_info      = all_data['jabr_cuts_info']
    num_cuts            = all_data['ID_jabr_cuts']
    num_jabr_cuts_added = 0
    threshold           = all_data['threshold']
    violated            = {}
    violated_count      = 0

    ##new code
    most_violated       = {}
    
    for k in range(T):
        violated[k] = {}
    ##
        
    all_data['NO_jabrs_violated'] = 0
    
    t0_violation = time.time()

    log.joint(' checking for violations of Jabr inequalities ... \n')

    for branch in branches.values():
        f = branch.f
        t = branch.t
        count_of_f = IDtoCountmap[f]
        count_of_t = IDtoCountmap[t]

        for k in range(T):
            violation = cvalues[k][branch]*cvalues[k][branch] + svalues[k][branch]*svalues[k][branch] - cvalues[k][buses[count_of_f]]*cvalues[k][buses[count_of_t]]
            if violation > threshold:
                violated_count += 1
                violated[k][branch] = violation

    if violated_count == 0:
        all_data['NO_jabrs_violated'] = 1
        log.joint(' all violations below threshold\n' )
        log.joint(' no more cuts to add for current threshold\n' )
        return None

    t1_violation = time.time()

    log.joint('  number violated Jabrs ' + str(violated_count) + '\n')  
    log.joint('  time spent on violation Jabrs '
              + str(t1_violation - t0_violation) + '\n')


    log.joint(' sorting most violated Jabr-envelope cuts ... \n')

    t0_mostviol = time.time()

    num_selected = math.ceil( violated_count *
                             all_data['most_violated_fraction_jabr'] )

    for k in range(T):
        most_violated[k] = dict(sorted(violated[k].items(),
                                       key = lambda x: x[1],
                                       reverse = True)[:num_selected])
        
    t1_mostviol = time.time()

    log.joint('  time spent sorting most violated Jabrs '
              + str(t1_mostviol - t0_mostviol) + '\n')

    log.joint(' computing Jabr-envelope cuts ... \n')
    t0_compute  = time.time()

    most_violated_count  =  0
    most_violated_branch = 'none'
    max_error            = -1
    
    for k in range(T):
        numkcuts = 0
        for branch in most_violated[k].keys():
            if (most_violated_count == 0):
                if most_violated[k][branch] > max_error:
                    most_violated_branch = 'branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' k ' + str(k)
                    max_error = most_violated[k][branch]

            f          = branch.f
            t          = branch.t
            count_of_f = IDtoCountmap[f]
            count_of_t = IDtoCountmap[t]
            cft        = cvalues[k][branch]
            sft        = svalues[k][branch]
            cff        = cvalues[k][buses[count_of_f]]
            ctt        = cvalues[k][buses[count_of_t]]

            cutnorm = math.sqrt( (2 * cft)**2 + (2 * sft)**2 + (cff - ctt)**2 )
            coeff_cft = 4 * cft
            coeff_sft = 4 * sft
            coeff_cff = cff - ctt - cutnorm
            coeff_ctt = - (cff - ctt) - cutnorm

            if all_data['parallel_check']:
                if parallel_check(log,all_data,branch,coeff_cft,coeff_sft,
                                  coeff_cff,coeff_ctt,k):
                    continue

            #we add the cut
            most_violated_count += 1
            violation            = most_violated[k][branch]
            cutid                = num_cuts + most_violated_count

            if all_data['loud_cuts']:
                log.joint('  --> new cut\n')
                log.joint('  branch ' + str(branch.count) + ' time-period '
                          + str(k) + ' f ' + str(f) + ' t ' + str(t)
                          + ' violation ' + str(violation) + ' cut id '
                          + str(cutid) + '\n' )
                log.joint('  values ' + ' cft ' + str(cft) + ' sft ' + str(sft)
                          + ' cff ' + str(cff) + ' ctt ' + str(ctt) + '\n' )
                log.joint('  LHS coeff ' + ' cft ' + str(coeff_cft) + ' sft ' 
                          + str(coeff_sft) + ' cff ' + str(coeff_cff) + ' ctt '
                          + str(coeff_ctt) + '\n' )
                log.joint('  cutnorm ' + str(cutnorm) + '\n')

            # Sanity check
            if all_data['jabr_validity']:
                sol_c     = all_data['sol_cvalues'][k][branch]
                sol_s     = all_data['sol_svalues'][k][branch]
                sol_cbusf = all_data['sol_cvalues'][k][buses[count_of_f]]
                sol_cbust = all_data['sol_cvalues'][k][buses[count_of_t]]
                slack     = coeff_cft * sol_c + coeff_sft * sol_s + coeff_cff * sol_cbusf + coeff_ctt * sol_cbust

                if slack > FeasibilityTol:
                    log.joint('  this cut is not valid!\n')
                    log.joint('  violation ' + str(slack) + '\n')
                    log.joint('  values (a primal bound)' + ' cft '
                              + str(sol_c)  + ' sft ' + str(sol_s) + ' cff '
                              + str(sol_busf)  + ' ctt ' + str(sol_cbust)
                              + '\n' )
                    breakexit('check!')
                else:
                    log.joint('  valid cut at branch ' + str(branch.count)
                              + ' with slack ' + str(slack) + '\n')
        
            jabr_cuts[(cutid,branch.count)]  = (rnd,threshold,k)
            jabr_cuts_info[k][branch][cutid] = (rnd,violation,coeff_cft,
                                                coeff_sft,coeff_cff,coeff_ctt,
                                                threshold,cutid)
        
            cutexp = LinExpr()
            constrname = "jabr_cut_"+str(cutid)+"_"+str(branch.count)+"r_"+str(rnd)+"k_"+str(k)+"_"+str(f)+"_"+str(t)
            cutexp += coeff_cft * cvar[k][branch] + coeff_sft * svar[k][branch] + coeff_cff * cvar[k][buses[count_of_f]] + coeff_ctt * cvar[k][buses[count_of_t]]

            themodel.addConstr(cutexp <= 0, name = constrname)
            numkcuts += 1
        log.joint('  time-period = ' + str(k) + ' : Jabr-envelope cuts'
                  + ' added ' + str(numkcuts) + '\n')

    log.joint('  number Jabr-envelope cuts added ' + str(most_violated_count)
              + '\n')
    log.joint('  max error (abs) ' + str(max_error) + ' at '
              + str(most_violated_branch) + '\n' )

    t1_compute = time.time()

    log.joint('  time spent on computing Jabrs '
              + str(t1_compute - t0_compute) + '\n')

    all_data['most_violated_branch_jabr'] = most_violated_branch
    all_data['max_error_jabr']            = max_error
    all_data['ID_jabr_cuts']             += most_violated_count
    all_data['num_jabr_cuts_added']       = most_violated_count
    all_data['num_jabr_cuts']            += most_violated_count
    all_data['num_jabr_cuts_rnd'][rnd]    = most_violated_count

    
    if all_data['dropjabr']:
        if all_data['addcuts']:
            t0_drop = time.time()
            drop_jabr(log,all_data)
            t1_drop = time.time()
            log.joint('  time spent on drop Jabrs ' + str(t1_drop - t0_drop)
                      + '\n')
        elif all_data['round'] >= all_data['cut_age_limit']: ######
            t0_drop = time.time()
            drop_jabr(log,all_data)
            t1_drop = time.time()
            log.joint('  time spent on drop Jabrs ' + str(t1_drop - t0_drop)
                      + '\n')

def drop_jabr(log,all_data):
    

    log.joint(' dropping old and slack Jabr-envelope cuts ...\n')

    themodel               = all_data['themodel']
    branches               = all_data['branches']
    buses                  = all_data['buses']
    IDtoCountmap           = all_data['IDtoCountmap']
    current_rnd            = all_data['round']
    cut_age_limit          = all_data['cut_age_limit']
    jabr_cuts              = all_data['jabr_cuts']
    jabr_cuts_info         = all_data['jabr_cuts_info']
    drop_jabrs             = []
    num_jabr_cuts_dropped  = all_data['num_jabr_cuts_dropped'] 
    dropped_jabrs          = all_data['dropped_jabrs']
    T                      = all_data['T']
    
    for key in jabr_cuts.keys():
        cut     = jabr_cuts[key] 
        cut_rnd = cut[0]
        cut_age = current_rnd - cut_rnd 

        if cut_age <= cut_age_limit: #fix this... we are not cleaning up precomputed cuts...
            continue #baby cuts, skip! 

        cutid         = key[0]
        branchid      = key[1]
        branch        = branches[branchid]
        f             = branch.f
        t             = branch.t
        count_of_f    = IDtoCountmap[f]
        count_of_t    = IDtoCountmap[t]
        cut_threshold = cut[1]
        k             = cut[2]

        constrname    = "jabr_cut_"+str(cutid)+"_"+str(branchid)+"r_"+str(cut_rnd)+"k_"+str(k)+"_"+str(f)+"_"+str(t)
        constr        = themodel.getConstrByName(constrname)
        slack         = - constr.getAttr("slack")
        
        if ( slack < - cut_threshold ):
            drop_jabrs.append(key)
            themodel.remove(themodel.getConstrByName(constrname))
            #if all_data['loud_cuts']:
            #    log.joint('  --> removed Jabr-cut\n')
            #    log.joint('  the cut ' + str(key)
            #              + ' was removed from the model\n')
    
    num_drop_jabrs = len(drop_jabrs)

    if num_drop_jabrs:
        all_data['num_jabr_cuts_dropped'] = num_drop_jabrs
        all_data['num_jabr_cuts']        -= num_drop_jabrs
        all_data['total_jabr_dropped']   += num_drop_jabrs

        dropped_jabrs.extend(drop_jabrs)

        ##### drop Jabr-envelope cuts from jabr_cuts_info                
        for key in drop_jabrs:
            cutid       = key[0]
            branchid    = key[1]
            k           = jabr_cuts[key][2]
            cuts_branch = jabr_cuts_info[k][branches[branchid]]
            if all_data['loud_cuts']:
                log.joint('  --> Jabr-cut removed\n')
                log.joint('  cutid ' + str(cutid) + ' branchid '
                          + str(branchid) + ' k ' + str(k) + '\n')
                log.joint(' coeff: cft ' + str(cuts_branch[cutid][2])
                          + ' sft ' + str(cuts_branch[cutid][3]) + ' cff '
                          + str(cuts_branch[cutid][4]) + ' ctt '
                          + str(cuts_branch[cutid][5]) + ' rnd '
                          + str(cuts_branch[cutid][0]) + '\n')
            cuts_branch.pop(cutid)
    
        log.joint('  cuts in drop_jabrs list were removed from jabr_cuts_info\n')

        ##### drop Jabr-envelope cuts from jabr_cuts        
        for key in drop_jabrs:
            jabr_cuts.pop(key)

        log.joint('  the cuts in drop_jabrs list were removed from dict jabr_cuts\n')
        
    else:
        all_data['num_jabr_cuts_dropped'] = 0
        log.joint('  no Jabr-envelope cuts were dropped this round\n')


def add_cuts(log,all_data):

    if '_b' in all_data['casename']:
        original_casename = all_data['casename'][:len(all_data['casename']) - 2]
    elif ('_n_5_5' in all_data['casename'] or '_n_0_5' in all_data['casename'] 
          or '_n_1_1' in all_data['casename'] or '_pline' in all_data['casename']):
        original_casename = all_data['casename'][:len(all_data['casename']) - 6]
    #elif '_line' in all_data['casename']:
    #    original_casename = all_data['casename'][:len(all_data['casename']) - 5]
    #elif '_line5' in all_data['casename']:
    #    original_casename = all_data['casename'][:len(all_data['casename']) - 6]
    else:
        original_casename = all_data['casename']

    filename = 'newcuts/cuts_' + original_casename + '.txt' ###cuts | newcuts
    log.joint(" opening file with cuts " + filename + "\n")

    try:
        thefile = open(filename, "r") 
        lines = thefile.readlines()
    except:
        log.stateandquit(" cannot open file", filename)
        sys.exit("failure")

    numlines  = len(lines)
    theround  = lines[0].split()[3]
    firstline = lines[1].split()
    jabr      = 1
    i2        = 0
    
    log.joint(' loading cuts from round ' + str(theround) + '\n')

    if firstline[0] == '#Jabr-envelope':
        numjabrcuts = int(firstline[3])
        all_data['addcuts_numjabrcuts'] = numjabrcuts
        log.joint(' number of Jabr-envelope cuts in file = ' + str(numjabrcuts) + '\n')
    elif firstline[0] == '#i2-envelope':
        jabr      = 0
        i2        = 1
        numi2cuts = int(firstline[3])
        all_data['addcuts_numi2cuts'] = numi2cuts
        log.joint(' no Jabr-envelope cuts to add\n')
        log.joint(' number of i2-envelope cuts in file = ' + str(numi2cuts) + '\n')
    elif firstline[0] == '#limit-envelope':
        jabr         = 0
        i2           = 0
        numlimitcuts = int(firstline[3])
        all_data['addcuts_numlimitcuts'] = numlimitcuts
        log.joint(' no Jabr nor i2-envelope cuts to add\n')
        log.joint(' number of limit-envelope cuts in file = ' + str(numlimitcuts) + '\n')
    else:
        log.joint(' no cuts added\n')
        return None

    linenum = 2

    themodel       = all_data['themodel']
    buses          = all_data['buses']
    branches       = all_data['branches']
    IDtoCountmap   = all_data['IDtoCountmap']
    cvar           = all_data['cvar']
    svar           = all_data['svar']
    Pvar_f         = all_data['Pvar_f']
    Qvar_f         = all_data['Qvar_f']
    Pvar_t         = all_data['Pvar_t']
    Qvar_t         = all_data['Qvar_t']
    T              = all_data['T']
    FeasibilityTol = all_data['FeasibilityTol']
    alphadic       = all_data['alphadic']

    if all_data['i2']:
        i2var_f    = all_data['i2var_f']
        numi2added = 0
    else:
        numi2cuts  = 0

        
    while linenum < numlines: 
        thisline = lines[linenum].split()
        if thisline[0] == '#i2-envelope' and jabr:
            if all_data['i2']:
                numi2cuts = int(thisline[3])
                all_data['addcuts_numi2cuts'] = numi2cuts
                log.joint(' number of i2-envelope cuts in file = '
                          + str(numi2cuts) + '\n')
            linenum += 1
            jabr = 0
            i2   = 1
            continue
            

        elif thisline[0] == '#limit-envelope' and i2: #i2 before
            numlimitcuts = int(thisline[3])
            all_data['addcuts_numlimitcuts'] = numlimitcuts
            log.joint(' number of limit-envelope cuts in file = ' + str(numlimitcuts) 
                      + '\n')
            linenum += 1
            i2       = 0
            continue

        elif all_data['i2'] == 0 and i2:
            linenum += 1
            continue
            

        elif jabr:
            branchid  = int(thisline[1])
            f         = int(thisline[3])
            t         = int(thisline[5])
            cutid     = int(thisline[7])
            rnd       = int(thisline[9])
            coeff_cft = float(thisline[15])
            coeff_sft = float(thisline[17])
            coeff_cff = float(thisline[19])
            coeff_ctt = float(thisline[21])


            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' 
                          + str(branchid) + ' f ' + str(f) + ' t ' + str(t) 
                          + ' was turned OFF\n')
                linenum += 1
                continue

            branch     = branches[branchid]
            count_of_f = IDtoCountmap[f] 
            count_of_t = IDtoCountmap[t]

            if ( (branchid != branch.count) or f != (branch.f) 
                 or (t != branch.t) ):

                log.joint(' branchid ' + str(branchid) + ' branch.count ' 
                          + str(branch.count) + ' f ' + str(f) + ' branch.f ' 
                          + str(branch.f) + ' t ' + str(t) + ' branch.t ' 
                          + str(branch.t) + '\n')
                breakexit('bug')

            if all_data['loud_cuts']:
                log.joint(' --> new Jabr-envelope cut for every time period\n')
                log.joint(' branch ' + str(branchid) + ' f ' + str(f) + ' t ' 
                          + str(t) + ' cutid ' + str(cutid) + '\n' )
                log.joint(' LHS coeff ' + ' cft ' + str(coeff_cft) + ' sft ' 
                          + str(coeff_sft) + ' cff ' + str(coeff_cff) 
                          + ' ctt ' + str(coeff_ctt) + '\n' )

            
            if all_data['jabr_validity']:
                for k in range(T):
                    sol_c       = all_data['sol_cvalues'][k][branch]
                    sol_s       = all_data['sol_svalues'][k][branch]
                    sol_cbusf   = all_data['sol_cvalues'][k][buses[count_of_f]]
                    sol_cbust   = all_data['sol_cvalues'][k][buses[count_of_t]]
                    sol_viol    = coeff_cft * sol_c + coeff_sft * sol_s + coeff_cff * sol_cbusf + coeff_ctt * sol_cbust
                    sol_relviol = sol_viol / ( ( coeff_cft**2 + coeff_sft**2 + coeff_cff**2 + coeff_ctt**2 )**0.5 ) 

                    if sol_relviol > FeasibilityTol:
                        log.joint(' WARNING, the Jabr-envelope cut associated '
                                  + 'to branch ' + str(branch.count) + ' k '
                                  + str(k) + ' f ' + str(branch.f)
                                  + ' t ' + str(branch.t)
                                  + ' is violated by the AC solution!\n')
                        log.joint(' violation ' + str(sol_viol) + '\n')
                        log.joint(' relative violation ' + str(sol_relviol)
                                  + '\n')
                        log.joint(' values (AC solution)' + ' cft '
                                  + str(sol_c) + ' sft ' + str(sol_s)
                                  + ' cff ' + str(sol_busf) + ' ctt '
                                  + str(sol_cbust) + '\n' )
                        breakexit('check!')
                    else:
                        log.joint(' AC solution satisfies Jabr inequality '
                                  + 'at branch ' + str(branch.count)
                                  + ' with slack ' + str(sol_relviol) + '\n')
                    
            for k in range(T):
                cutexp = LinExpr()
                constrname = "jabr_cut_"+str(cutid)+"_"+str(branchid)+"r_"+str(rnd)+"k_"+str(k)+"_"+str(f)+"_"+str(t)
                cutexp += coeff_cft * cvar[k][branch] + coeff_sft * svar[k][branch] + coeff_cff * cvar[k][buses[count_of_f]] + coeff_ctt * cvar[k][buses[count_of_t]]

                themodel.addConstr(cutexp <= 0, name = constrname)
                
            linenum += 1

        elif (jabr == 0) and i2 and all_data['i2']:

            branchid   = int(thisline[1])
            f          = int(thisline[3])
            t          = int(thisline[5])
            cutid      = int(thisline[7])
            coeff_Pft  = float(thisline[15])
            coeff_Qft  = float(thisline[17])
            coeff_cff  = float(thisline[19])
            coeff_i2ft = float(thisline[21])
            alpha      = all_data['alphadic'][branches[branchid]]

            #print(thisline)
            #log.joint(' cutid ' + str(cutid) + '\n')
            #log.joint(' alpha ' + str(alpha) + '\n')
            #log.joint(' branch ' + str(branch) + ' id ' + str(branch.count)
            #  + ' alpha ' + str(alpha) + '\n')
            
            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' 
                          + str(branchid) + ' f ' + str(f) + ' t ' + str(t) 
                          + 'was turned OFF\n')
                linenum += 1
                continue

            if alpha >= all_data['alpha_threshold']:
                log.joint(' stability: branch ' + str(branchid) + ' f '
                          + str(f) + ' t ' + str(t) + ' was skipped\n')
                linenum += 1
                continue
            
            branch     = branches[branchid]
            count_of_f = IDtoCountmap[f] 
            count_of_t = IDtoCountmap[t]

            if ( (branchid != branch.count) or f != (branch.f) 
                 or (t != branch.t) ):
                breakexit('there might be bug')

            if all_data['loud_cuts']:
                log.joint(' --> new i2-envelope cut for every time period\n')
                log.joint(' branch ' + str(branchid) + ' f ' + str(f) + ' t ' 
                          + str(t) + ' cutid ' + str(cutid) + '\n' )
                log.joint(' LHS coeff ' + ' Pft ' + str(coeff_Pft) + ' Qft ' 
                          + str(coeff_Qft) + ' cff ' + str(coeff_cff) 
                          + ' i2ft ' + str(coeff_i2ft) + '\n' )
            
            if all_data['i2_validity']:
                sol_Pf           = all_data['sol_Pfvalues'][branch]
                sol_Qf           = all_data['sol_Qfvalues'][branch]
                sol_c            = all_data['sol_cvalues'][branch]
                sol_s            = all_data['sol_svalues'][branch]
                sol_cbusf        = all_data['sol_cvalues'][buses[count_of_f]]
                sol_cbust        = all_data['sol_cvalues'][buses[count_of_t]]
                sol_i2f          = computei2value(log,all_data,branch,sol_c,sol_s,sol_cbusf,sol_cbust)
                sol_viol    = coeff_Pft * sol_Pf + coeff_Qft * sol_Qf + coeff_cff * sol_cbusf + coeff_i2ft * sol_i2f
                sol_relviol = sol_viol / ( ( coeff_Pft**2 + coeff_Qft**2 + coeff_cff**2 + coeff_i2ft**2 )**0.5 ) 

                if sol_relviol > FeasibilityTol:
                    log.joint(' WARNING, the i2-envelope cut associated to branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' is violated by the AC solution!\n')
                    log.joint(' violation ' + str(sol_viol) + '\n')
                    log.joint(' relative violation ' + str(sol_relviol) + '\n')
                    log.joint(' values (AC solution) ' + ' Pft ' + str(sol_Pf) + ' Qft ' + str(sol_Qf) + ' cff ' + str(sol_cbusf) + ' i2ft ' + str(sol_i2f) + '\n' )
                    breakexit('check!')
                else:
                    log.joint(' AC solution satisfies i2 inequality at branch ' + str(branch.count) + ' with slack ' + str(sol_relviol) + '\n')

            for k in range(T):
                #log.joint(' are we here? ' + str(alpha) + '\n')
                #log.joint(' branch ' + str(branch) + ' id '
                #          + str(branch.count) + ' k ' + str(k) + '\n')
                #log.joint(' i2 ' + str(i2var_f[k][branch]) + '\n')
                #breakexit('c')
                
                cutexp     = LinExpr()
                constrname = "i2_cut_"+str(cutid)+"_"+str(branch.count)+"r_"+str(rnd)+"k_"+str(k)+"_"+str(f)+"_"+str(t)
                cutexp    += coeff_Pft * Pvar_f[k][branch] + coeff_Qft * Qvar_f[k][branch] + coeff_cff * cvar[k][buses[count_of_f]] + coeff_i2ft * i2var_f[k][branch]

                themodel.addConstr(cutexp <= 0, name = constrname)
                numi2added += 1
                
            linenum += 1


        elif (jabr == 0) and (i2 == 0):
            branchid   = int(thisline[1])
            f          = int(thisline[3])
            t          = int(thisline[5])
            cutid      = int(thisline[7])
            if thisline[14] == 'Pft':
                from_or_to = 'f'
            elif thisline[14] == 'Ptf':
                from_or_to = 't'
            else:
                log.joint(' look for a bug\n')
                breakexit('bug')

            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' 
                          + str(branchid) + ' f ' + str(f) + ' t ' + str(t) 
                          + 'was turned OFF\n')
                linenum += 1
                continue

            coeff_P    = float(thisline[15])
            coeff_Q    = float(thisline[17])

            branch     = branches[branchid]
            count_of_f = IDtoCountmap[f] 
            count_of_t = IDtoCountmap[t]

            if ( (branchid != branch.count) or f != (branch.f) 
                 or (t != branch.t) ):
                breakexit('there might be bug')

            if all_data['loud_cuts']:
                log.joint(' --> new limit-envelope cut for every time period\n')
                log.joint(' branch ' + str(branchid) + ' f ' + str(f) + ' t ' 
                          + str(t) + ' cutid ' + str(cutid) + '\n' )
                if from_or_to == 'f':
                    log.joint(' LHS coeff ' + ' Pft ' + str(coeff_P) 
                              + ' Qft ' + str(coeff_Q) + '\n')
                elif from_or_to == 't':
                    log.joint(' LHS coeff ' + ' Ptf ' + str(coeff_P) 
                              + ' Qtf ' + str(coeff_Q) + '\n')
            
            
            if all_data['limit_validity']:
                if from_or_to == 'f':
                    sol_P           = all_data['sol_Pfvalues'][branch]
                    sol_Q           = all_data['sol_Qfvalues'][branch]
                elif from_or_to == 't':
                    sol_P           = all_data['sol_Ptvalues'][branch]
                    sol_Q           = all_data['sol_Qtvalues'][branch]
                sol_viol    = coeff_P * sol_P + coeff_Q * sol_Q - 1

                if sol_viol > FeasibilityTol:
                    log.joint(' WARNING, the limit-envelope cut associated to branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' is violated by the AC solution!\n')
                    log.joint(' violation ' + str(sol_viol) + '\n')
                    if from_or_to == 'f':
                        log.joint(' values (AC solution) ' + ' Pft ' + str(sol_P) + ' Qft ' + str(sol_Q) + '\n')
                    elif from_or_to == 't':
                        log.joint(' values (AC solution) ' + ' Ptf ' + str(sol_P) + ' Qtf ' + str(sol_Q) + '\n')
                    breakexit('check!')
                else:
                    log.joint(' AC solution satisfies limit inequality at branch ' + str(branch.count) + ' with slack ' + str(sol_viol) + '\n')


            for k in range(T):
                cutexp     = LinExpr()
                if from_or_to == 'f':
                    constrname = "limit_cut_"+str(cutid)+"_"+str(branch.count)+"r_"+str(rnd)+"k_"+str(k)+"_"+str(f)+"_"+str(t)
                    cutexp    += coeff_P * Pvar_f[k][branch] + coeff_Q * Qvar_f[k][branch]
                elif from_or_to == 't':
                    constrname = "limit_cut_"+str(cutid)+"_"+str(branch.count)+"r_"+str(rnd)+"_"+str(f)+"_"+str(t)+"_"+str(k)
                    cutexp    += coeff_P * Pvar_t[k][branch] + coeff_Q * Qvar_t[k][branch]

                themodel.addConstr(cutexp <= 1, name = constrname)
            linenum += 1

    all_data['num_jabr_cuts']  = numjabrcuts
    all_data['num_i2_cuts']    = numi2cuts
    all_data['num_limit_cuts'] = numlimitcuts
    all_data['ID_jabr_cuts']   = numjabrcuts
    all_data['ID_i2_cuts']     = numi2cuts
    all_data['ID_limit_cuts']  = numlimitcuts


def write_cuts(log,all_data):

    #filename = 'cuts_' + all_data['casename'] + '.txt' #change to newcuts when running bash script
    filename = 'cuts/cuts_' + all_data['casename'] + '.txt' #change to newcuts when running bash script
    log.joint(" opening file with cuts " + filename + "\n")

    try:
        thefile = open(filename, "w+") #a+ if we want to append cuts of different rounds
        lines = thefile.readlines()
    except:
        log.stateandquit(" cannot open file", filename)
        sys.exit("failure")
    
    jabr_cuts = all_data['jabr_cuts_info']
    rnd       = all_data['round'] 

    log.joint(' writing down to ' + filename + ' jabr-envelope cuts in round ' + str(rnd) + '\n')

    thefile.write('current round = ' + str(all_data['round']) + '\n')
    thefile.write('#Jabr-envelope cuts = ' + str(all_data['num_jabr_cuts']) + '\n')
    for branch in jabr_cuts.keys():
        branchcount = branch.count
        f           = branch.f
        t           = branch.t
        branch_cuts = jabr_cuts[branch]
        
        for cutid in branch_cuts.keys():
            cut = branch_cuts[cutid]
            if len(cut) == 0:
                continue
            cut_info = 'branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' cutid ' + str(cutid) + ' round ' + str(cut[0]) + ' violation ' + str(cut[1]) + ' threshold ' + str(cut[6]) + ' cft ' + str(cut[2]) + ' sft ' + str(cut[3]) + ' cff ' + str(cut[4]) + ' ctt ' + str(cut[5]) + '\n'
            thefile.write(cut_info)
        
    log.joint(' done writing down jabr-envelope cuts\n')

    if all_data['i2cuts']:
        log.joint(' writing down to ' + filename + ' i2-envelope cuts in round ' + str(rnd) + '\n')
        thefile.write('#i2-envelope cuts = ' + str(all_data['num_i2_cuts']) + '\n')
        
        i2_cuts = all_data['i2_cuts_info']
        for branch in i2_cuts.keys(): #BUG found
            branchcount = branch.count
            f           = branch.f
            t           = branch.t
            branch_cuts = i2_cuts[branch]
        
            for cutid in branch_cuts.keys():
                cut = branch_cuts[cutid]
                if len(cut) == 0:
                    continue
                cut_info = 'branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' cutid ' + str(cutid) + ' round ' + str(cut[0]) + ' violation ' + str(cut[1]) + ' threshold ' + str(cut[6]) + ' Pft ' + str(cut[2]) + ' Qft ' + str(cut[3]) + ' cff ' + str(cut[4]) + ' i2ft ' + str(cut[5]) + '\n'
                thefile.write(cut_info)        

        log.joint(' done writing down i2-envelope cuts\n')


    if all_data['limitcuts']:
        log.joint(' writing down to ' + filename + ' limit-envelope cuts in round ' + str(rnd) + '\n')
        thefile.write('#limit-envelope cuts = ' + str(all_data['num_limit_cuts']) + '\n')
        
        limit_cuts = all_data['limit_cuts_info']
        for branch in limit_cuts.keys():
            branchcount = branch.count
            f           = branch.f    #(rnd,violation,coeff_P,coeff_Q,threshold,cutid,from_or_to)
            t           = branch.t
            branch_cuts = limit_cuts[branch]
        
            for cutid in branch_cuts.keys():
                cut        = branch_cuts[cutid]

                if len(cut) == 0:
                    continue

                from_or_to = cut[6]
                if from_or_to == 'f':
                    cut_info = 'branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' cutid ' + str(cutid) + ' round ' + str(cut[0]) + ' violation ' + str(cut[1]) + ' threshold ' + str(cut[4]) + ' Pft ' + str(cut[2]) + ' Qft ' + str(cut[3]) + '\n'
                elif from_or_to == 't':
                    cut_info = 'branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' cutid ' + str(cutid) + ' round ' + str(cut[0]) + ' violation ' + str(cut[1]) + ' threshold ' + str(cut[4]) + ' Ptf ' + str(cut[2]) + ' Qtf ' + str(cut[3]) + '\n'
                else:
                    log.joint(' we have a bug\n')
                    breakexit('look for the bug')

                thefile.write(cut_info)        

        log.joint(' done writing down limit-envelope cuts\n')


    thefile.close()
    #breakexit('check file')

def parallel_check(log,all_data,branch,coeff_cft,coeff_sft,coeff_cff,coeff_ctt,k): 

    threshold         = all_data['threshold']
    threshold_dotprod = all_data['threshold_dotprod']    
    jabr_cuts_info    = all_data['jabr_cuts_info']
    cuts_branch       = jabr_cuts_info[k][branch] 
  
    if all_data['loud_cuts']:
        log.joint('\n -- parallel check wrt previous Jabr-envelope cuts at branch ' + str(branch.count) +'\n')

    if len(cuts_branch) == 0:
        if all_data['loud_cuts']:
            log.joint(' first Jabr-envelope cut, we add it\n')
        return 0 
    else:
        v = compute_normal(coeff_cft,coeff_sft,coeff_cff,coeff_ctt)
        if all_data['loud_cuts']:
            log.joint(' LHS coeffs potential cut ' + ' cft ' + str(coeff_cft) 
                      + ' sft ' + str(coeff_sft) + ' cff ' + str(coeff_cff) 
                      + ' ctt ' + str(coeff_ctt) + '\n' )

        for cut in cuts_branch.values():
            cutid         = cut[7]
            cut_coeff_cft = cut[2]
            cut_coeff_sft = cut[3]
            cut_coeff_cff = cut[4]
            cut_coeff_ctt = cut[5]
            w = compute_normal(cut_coeff_cft,cut_coeff_sft,cut_coeff_cff,
                               cut_coeff_ctt)
            dotprod = np.dot(v,w)
            angle = np.arccos(dotprod)
            angle_deg = angle * 180 / np.pi
            if all_data['loud_cuts']:
                log.joint(' LHS coeffs of cutid ' + str(cutid) + ' cft ' 
                          + str(cut_coeff_cft) + ' sft ' + str(cut_coeff_sft) 
                          + ' cff ' + str(cut_coeff_cff) + ' ctt ' 
                          + str(cut_coeff_ctt) + '\n' )
                log.joint(' angle (rad) ' + str(angle) + ' angle (deg) ' 
                          + str(angle_deg) +  ' dot-product ' + str(dotprod) 
                          + '\n')
            
            if dotprod > 1 - threshold_dotprod:
                if all_data['loud_cuts']:
                    log.joint(' parallel cut, should not be added\n')
                return 1
            else:
                if all_data['loud_cuts']:
                    log.joint(' cut should be added\n')
                return 0
        
     
def parallel_check_i2(log,all_data,branch,coeff_Pft,coeff_Qft,coeff_cff,coeff_i2ft,k):

    threshold         = all_data['threshold']
    threshold_dotprod = all_data['threshold_dotprod']    
    i2_cuts_info      = all_data['i2_cuts_info']      
    cuts_branch       = i2_cuts_info[k][branch]          

    if all_data['loud_cuts']:
        log.joint('\n -- parallel check wrt previous i2-envelope cuts at branch ' + str(branch.count) +'\n')

    if len(cuts_branch) == 0:
        if all_data['loud_cuts']:
            log.joint(' first i2-envelope cut, we add it\n')
        return 0 
    else:
        v = compute_normal(coeff_Pft,coeff_Qft,coeff_cff,coeff_i2ft)
        if all_data['loud_cuts']:
            log.joint(' LHS coeffs potential cut ' + ' Pft ' + str(coeff_Pft) 
                      + ' Qft ' + str(coeff_Qft) + ' cff ' + str(coeff_cff) 
                      + ' i2ft ' + str(coeff_i2ft) + '\n' )

        for cut in cuts_branch.values():
            cutid      = cut[7]
            cut_coeff_Pft  = cut[2]
            cut_coeff_Qft  = cut[3]
            cut_coeff_cff  = cut[4]
            cut_coeff_i2ft = cut[5]
            w = compute_normal(cut_coeff_Pft,cut_coeff_Qft,cut_coeff_cff,
                               coeff_i2ft)
            dotprod = np.dot(v,w)
            angle = np.arccos(dotprod)
            angle_deg = angle * 180 / np.pi
            if all_data['loud_cuts']:
                log.joint(' LHS coeffs of cutid ' + str(cutid) + ' cft ' 
                          + str(cut_coeff_Pft) + ' sft ' + str(cut_coeff_Qft) 
                          + ' cff ' + str(cut_coeff_cff) + ' ctt ' 
                          + str(cut_coeff_i2ft) + '\n' )
                log.joint(' angle (rad) ' + str(angle) + ' angle (deg) ' 
                          + str(angle_deg) +  ' dot-product ' + str(dotprod) 
                          + '\n')
            if dotprod > 1 - threshold_dotprod:
                if all_data['loud_cuts']:
                    log.joint(' parallel cut, should not be added\n')
                return 1
            else:
                if all_data['loud_cuts']:
                    log.joint(' cut should be added\n')
                return 0


def parallel_check_limit(log,all_data,branch,coeff_P,coeff_Q,from_or_to,k):

    threshold         = all_data['threshold']
    threshold_dotprod = all_data['threshold_dotprod']    
    limit_cuts_info   = all_data['limit_cuts_info']
    cuts_branch       = limit_cuts_info[k][branch]          

    if all_data['loud_cuts']:
        log.joint('\n -- parallel check wrt previous limit-envelope cuts at branch ' + str(branch.count) +'\n')

    if len(cuts_branch) == 0:
        if all_data['loud_cuts']:
            log.joint(' first limit-envelope cut, we add it\n')
        return 0 
    else:        
        v = compute_normal(coeff_P,coeff_Q)
        if all_data['loud_cuts']:
            if from_or_to == 'f':
                log.joint(' LHS coeffs potential cut ' + ' Pft ' 
                          + str(coeff_P) + ' Qft ' + str(coeff_Q) + '\n')
            elif from_or_to == 't':
                log.joint(' LHS coeffs potential cut ' + ' Ptf ' 
                          + str(coeff_P) + ' Qtf ' + str(coeff_Q) + '\n' )

        for cut in cuts_branch.values(): 
            cutid          = cut[5]
            cut_coeff_P    = cut[2]
            cut_coeff_Q    = cut[3]
            cut_from_or_to = cut[6]

            if from_or_to != cut_from_or_to: #we check whether potential cut has the same sense (from/to) as the incumbent cut
                continue

            w         = compute_normal(cut_coeff_P,cut_coeff_Q)
            dotprod   = np.dot(v,w)
            angle     = np.arccos(dotprod)
            angle_deg = angle * 180 / np.pi
            
            if all_data['loud_cuts']:
                if cut_from_or_to == 'f':
                    log.joint(' LHS coeffs of cutid ' + str(cutid) + ' Pft ' 
                              + str(cut_coeff_P) + ' Qft ' + str(cut_coeff_Q) 
                              + '\n')
                elif cut_from_or_to == 't':
                    log.joint(' LHS coeffs of cutid ' + str(cutid) + ' Ptf ' 
                              + str(cut_coeff_P) + ' Qtf ' + str(cut_coeff_Q) 
                              + '\n')
                log.joint(' angle (rad) ' + str(angle) + ' angle (deg) ' 
                          + str(angle_deg) +  ' dot-product ' + str(dotprod) 
                          + '\n')

            if dotprod > 1 - threshold_dotprod:
                if all_data['loud_cuts']:
                    log.joint(' parallel cut, should not be added\n')
                return 1
            else:
                if all_data['loud_cuts']:
                    log.joint(' cut should be added\n')
                return 0
            
        
def objective_cuts(log,all_data):

    log.joint('\n')
    log.joint(' **** objective-cuts ****\n')
    log.joint('\n')

    themodel           = all_data['themodel']
    gens               = all_data['gens']
    GenPvar            = all_data['GenPvar']
    GenPvalues         = all_data['GenPvalues']
    rnd                = all_data['round']
    num_objective_cuts = all_data['num_objective_cuts']
    dicGenPvalues      = all_data['dicGenPvalues']
    threshold          = all_data['threshold_objcuts']
    GenTvar            = all_data['GenTvar']
    
        
    if rnd == 1:
        for gen in gens.values():
            if gen.costdegree == 2 and gen.costvector[0] != 0:
                gen_values = dicGenPvalues[gen]
                current_sol = GenPvalues[gen]
                index = 0
                a = gen.costvector[0]
                b = gen.costvector[1]
                x0 = gen_values[index]
            
                const = 0
                slope = b
                trueRHS = a * current_sol * current_sol + b * current_sol
                approxRHS = const + slope * current_sol
                violation = trueRHS - approxRHS

                #print("genid = {0} gen_values = {1}, index = {2}, violation = {3}".format(gen.count,gen_values,index,violation))
                if violation > threshold:
                    linexpr_gen = LinExpr()
                    current_const = - a * current_sol * current_sol
                    current_slope = 2 * a * current_sol + b
                    linexpr_gen += current_const + current_slope * GenPvar[gen]
                    
                    if all_data['loud_cuts']:
                        log.joint(' ----\n')
                        log.joint(' violated objective-cut at gen id ' 
                                  + str(gen.count) + ' violation ' 
                                  + str(violation) + ' threshold ' 
                                  + str(threshold) +  '\n')
                        log.joint(' LHS coeff violated cut: constant ' 
                                  + str(const) + ' slope ' + str(slope) 
                                  + ' RHS ' + str(approxRHS) + '\n')
                        log.joint(' GenPvalues: old Pgen ' + str(x0) 
                                  + ' new Pgen ' + str(current_sol) + '\n')
                        log.joint(' LHS coeff new cut: constant ' 
                                  + str(current_const) + ' slope ' 
                                  + str(current_slope) + ' RHS ' 
                                  + str(trueRHS) + '\n')

                    themodel.addConstr(linexpr_gen <= GenTvar[gen], name = 'qcost_gen_cut_' + str(rnd) + '_' + str(gen.count) + '_' + str(gen.nodeID) )
                    num_objective_cuts += 1
                    
                    if all_data['loud_cuts']:
                        log.joint(' the objective-cut associated to gen id ' 
                                  + str(gen.count) + ' was added \n' )
            
    else:
        for gen in gens.values():
            if gen.costdegree == 2 and gen.costvector[0] != 0:
                gen_values  = dicGenPvalues[gen]
                current_sol = GenPvalues[gen]
                index_sol   = gen_values.index(current_sol)
                a           = gen.costvector[0]
                b           = gen.costvector[1]

                if index_sol == 0 or index_sol == rnd - 1:
                    if index_sol == 0:
                        index = index_sol + 1
                    else:
                        index = index_sol - 1
                    x0 = gen_values[index]
                    const = - a * x0 * x0
                    slope = 2 * a * x0 + b
                    trueRHS = a * current_sol * current_sol + b * current_sol
                    approxRHS = const + slope * current_sol
                    violation = trueRHS - approxRHS
                    #print(" genid = {0} gen_values = {1}, index = {2}, violation = {3}".format(gen.count,gen_values,index,violation))
                    if violation > threshold:
                        linexpr_gen = LinExpr()
                        current_const = - a * current_sol * current_sol
                        current_slope = 2 * a * current_sol + b                    
                        linexpr_gen += current_const + current_slope * GenPvar[gen]
                        if all_data['loud_cuts']:    
                            log.joint(' ----\n')
                            log.joint(' violated objective-cut at gen id ' 
                                      + str(gen.count) + ' violation ' 
                                      + str(violation) + ' threshold ' 
                                      + str(threshold) +  '\n')
                            log.joint(' LHS coeff violated cut: constant ' 
                                      + str(const) + ' slope ' + str(slope) 
                                      + ' RHS ' + str(approxRHS) + '\n')
                            log.joint(' GenPvalues: old Pgen ' + str(x0) 
                                      + ' new Pgen ' + str(current_sol) + '\n')
                            log.joint(' LHS coeff new cut: constant ' 
                                      + str(current_const) + ' slope ' 
                                      + str(current_slope) + ' RHS ' 
                                      + str(trueRHS) + '\n')

                        themodel.addConstr(linexpr_gen <= GenTvar[gen], name = 'qcost_gen_cut_' + str(rnd) + '_' + str(gen.count) + '_' + str(gen.nodeID) )
                        num_objective_cuts += 1
                        
                        if all_data['loud_cuts']:
                            log.joint(' the objective-cut associated to gen id ' + str(gen.count) + ' was added \n' )
                    
                else:
                    index_left = index_sol - 1
                    index_right = index_sol + 1                    
                    x0_left = gen_values[index_left]
                    x0_right = gen_values[index_right]
                    const_left = - a * x0_left * x0_left
                    const_right = - a * x0_right * x0_right
                    slope_left = 2 * a * x0_left + b
                    slope_right = 2 * a * x0_right + b
                    trueRHS = a * current_sol * current_sol + b * current_sol
                    approxRHS_left = const_left + slope_left * current_sol
                    approxRHS_right = const_right + slope_right * current_sol
                    violation_left = trueRHS - approxRHS_left
                    violation_right = trueRHS - approxRHS_right
                    #print(" genid = {0} gen_values = {1}, index_left = {2}, index_right = {3}, violation_left = {4}, violation_right = {5}".format(gen.count,gen_values,index_left,index_right,violation_left,violation_right))
                    if violation_left <= violation_right:
                        if violation_left > threshold:
                            linexpr_gen = LinExpr()
                            current_const = - a * current_sol * current_sol
                            current_slope = 2 * a * current_sol + b                    
                            linexpr_gen += current_const + current_slope * GenPvar[gen]
                            if all_data['loud_cuts']:
                                log.joint(' -----\n')
                                log.joint(' violated objective-cut at gen id '
                                          + str(gen.count) + ' violation ' 
                                          + str(violation_left) + ' threshold '
                                          + str(threshold) +  '\n')
                                log.joint(' LHS coeff violated cut: constant ' 
                                          + str(const_left) + ' slope ' 
                                          + str(slope_left) + ' RHS ' 
                                          + str(approxRHS_left) + '\n')
                                log.joint(' GenPvalues: old Pgen ' 
                                          + str(x0_left) + ' new Pgen ' 
                                          + str(current_sol) + '\n')
                                log.joint(' LHS coeff new cut: constant ' 
                                          + str(current_const) + ' slope ' 
                                          + str(current_slope) + ' RHS ' 
                                          + str(trueRHS) + '\n')

                            themodel.addConstr(linexpr_gen <= GenTvar[gen], name = 'qcost_gen_cut_' + str(rnd) + '_' + str(gen.count) + '_' + str(gen.nodeID) )
                            num_objective_cuts += 1
                            if all_data['loud_cuts']:
                                log.joint(' the objective-cut associated to gen id ' + str(gen.count) + ' was added \n' )
                    
                    else:
                        if violation_right > threshold:
                            linexpr_gen = LinExpr()
                            current_const = - a * current_sol * current_sol
                            current_slope = 2 * a * current_sol + b                    
                            linexpr_gen += current_const + current_slope * GenPvar[gen]
                            if all_data['loud_cuts']:
                                log.joint(' ----\n')
                                log.joint(' violated objective-cut at gen id ' 
                                          + str(gen.count) + ' violation ' 
                                          + str(violation_right) 
                                          + ' threshold ' + str(threshold) 
                                          + '\n')
                                log.joint(' LHS coeff violated cut: constant ' 
                                          + str(const_right) + ' slope ' 
                                          + str(slope_right) + ' RHS ' 
                                          + str(approxRHS_right) + '\n')
                                log.joint(' GenPvalues: old Pgen ' 
                                          + str(x0_right) + ' new Pgen ' 
                                          + str(current_sol) + '\n')
                                log.joint(' LHS coeff new cut: constant ' 
                                          + str(current_const) + ' slope ' 
                                          + str(current_slope) + ' RHS ' 
                                          + str(trueRHS) + '\n')
                            
                            themodel.addConstr(linexpr_gen <= GenTvar[gen], name = 'qcost_gen_cut_' + str(rnd) + '_' + str(gen.count) + '_' + str(gen.nodeID) )
                            num_objective_cuts += 1
                            
                            if all_data['loud_cuts']:
                                log.joint(' the objective-cut associated to gen id ' + str(gen.count) + ' was added \n' )

            
    all_data['num_objective_cuts'] = num_objective_cuts

    
                     
def cut_analysis(log,all_data):

    log.joint('\n')
    log.joint(' **** Jabr-envelope cuts analysis ****\n')
    log.joint('\n')
    
    themodel = all_data['themodel']
    branches = all_data['branches']
    current_rnd =all_data['round']    

    tight_cuts = all_data['tight_cuts']
    tight_cuts_fraction = all_data['tight_cuts_fraction']
    count_tight = 0
    count_rnd = 1
    
    num_jabr_cuts = all_data['num_jabr_cuts_rnd']
    jabr_cuts_info = all_data['jabr_cuts_info'] #havent fixed this...
    jabr_cuts_info_updated = all_data['jabr_cuts_info_updated'] #havent fixed this... *_updated deleted
    
    
    most_violated_branches = sorted(jabr_cuts_info, key=lambda k: len(jabr_cuts_info[k]), reverse=True)[:10]
    most_violated_branches_updated = sorted(jabr_cuts_info_updated, key=lambda k: len(jabr_cuts_info_updated[k]), reverse=True)[:10]

    
    log.joint(' ---- Top 10 most violated branches\n')
    for branch in most_violated_branches:
        branchid = branch.count
        f = branch.f
        t = branch.t
        cuts_branch = jabr_cuts_info[branch] #dictionary (rnd,violation,coeff_cft,coeff_sft,coeff_cff,coeff_ctt,threshold)
        num_cuts = len(cuts_branch)
        log.joint(' -- Jabr-envelope cuts at branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ', number cuts = '+ str(num_cuts) + '\n')
        for cutid in cuts_branch.keys():            
            log.joint(' cutid ' + str(cutid) + ' added in rnd ' + str(cuts_branch[cutid][0]) + ' violation ' + str(cuts_branch[cutid][1]) + ' threshold ' + str(cuts_branch[cutid][6] ) + '\n')
            log.joint(' LHS coeff: ' + ' cft ' + str(cuts_branch[cutid][2]) + ' sft ' + str(cuts_branch[cutid][3]) + ' cff ' + str(cuts_branch[cutid][4]) + ' ctt ' + str(cuts_branch[cutid][5]) + '\n' )
            log.joint('\n')

    log.joint(' ---- Top 10 most violated branches (updated)\n')
    for branch in most_violated_branches_updated:
        branchid = branch.count
        f = branch.f
        t = branch.t
        cuts_branch = jabr_cuts_info_updated[branch] #dictionary (rnd,violation,coeff_cft,coeff_sft,coeff_cff,coeff_ctt,threshold)
        num_cuts = len(cuts_branch)
        log.joint(' -- Jabr-envelope cuts at branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ', number cuts = '+ str(num_cuts) + '\n')

    log.joint('\n')
            
    log.joint(' ---- ALL Jabr-envelope cuts\n')
    for branch in branches.values():
        branchid = branch.count
        f = branch.f
        t = branch.t
        cuts_branch = jabr_cuts_info_updated[branch] #dictionary (rnd,violation,coeff_cft,coeff_sft,coeff_cff,coeff_ctt,threshold)
        num_cuts = len(cuts_branch)

        log.joint(' -- Jabr-envelope cuts at branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ', number cuts = '+ str(num_cuts) + '\n')

        for cutid in cuts_branch.keys():            
            log.joint(' cutid ' + str(cutid) + ' added in rnd ' + str(cuts_branch[cutid][0]) + ' violation ' + str(cuts_branch[cutid][1]) + ' threshold ' + str(cuts_branch[cutid][6] ) + '\n')
            log.joint(' LHS coeff: ' + ' cft ' + str(cuts_branch[cutid][2]) + ' sft ' + str(cuts_branch[cutid][3]) + ' cff ' + str(cuts_branch[cutid][4]) + ' ctt ' + str(cuts_branch[cutid][5]) + '\n' )
            log.joint('\n')



    #tightness
    for key in all_data['jabr_cuts'].keys():
        cut = all_data['jabr_cuts'][key] 
        cut_rnd = cut[0]
        cut_age = current_rnd - cut_rnd
        if cut_age == 0:
            break #cut new born, skip
        
        cutid = key[0]
        branchid = key[1]
        branch = branches[branchid]
        f = branch.f
        t = branch.t
        cut_threshold = cut[6]
        
        constrname = "jabr_cut_"+str(cutid)+"_"+str(branchid)+"r_"+str(cut_rnd)+"_"+str(f)+"_"+str(t)
        constr = themodel.getConstrByName(constrname)
        slack = - constr.getAttr("slack")
        
        if ( slack >= - cut_threshold ):
            if count_rnd != cut_rnd:
                tight_cuts[count_rnd] = count_tight
                tight_cuts_fraction[count_rnd] = round(count_tight/num_jabr_cuts[count_rnd],3)
                count_tight = 1
                count_rnd = cut_rnd
            else:
                count_tight += 1
            log.joint(' cutid ' + str(cutid) + ' tight at branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' slack ' +str(slack) + '\n' )

    all_data['tight_cuts'] = tight_cuts
    all_data['tight_cuts_fraction'] = tight_cuts_fraction

    log.joint('\n')
    log.joint(' number Jabr-envelope cuts per round ' + dumps(num_jabr_cuts) + ' \n' )
    log.joint(' number tight Jabr-envelope cuts ' + dumps(tight_cuts) + ' \n' )
    log.joint(' % of tight Jabr-envelope cuts ' + dumps(tight_cuts_fraction) + ' \n' )


def computei2value(log,all_data,branch,sol_c,sol_s,sol_cbusf,sol_cbust):
    
    ratio  = branch.ratio
    y      = branch.y
    g      = y.real
    b      = y.imag
    bshunt = branch.bc
    angle  = branch.angle_rad
                                                                                
    #first f
    i2f  = 0
    i2f += (g*g + b*b)/(ratio*ratio) * ( (sol_cbusf/(ratio*ratio)) + sol_cbust - (2/ratio) * ( sol_c * math.cos(angle) + sol_s * math.sin(angle) ) )
    i2f += b*bshunt/(ratio**3) * ( (sol_cbusf/ratio) - (sol_c * math.cos(angle) + sol_s * math.sin(angle) )) 
    i2f += g*bshunt/(ratio**3) * ( sol_s * math.cos(angle) - sol_c * math.sin(angle) )
    i2f += (bshunt*bshunt*sol_cbusf/(4*(ratio**4)) )

    #now t
    #     #expr_t += (g*g + b*b) * ( cvar[buses[count_of_t]] + cvar[buses[count_of_f]]/(ratio*ratio) - (2/ratio) * ( cvar[branch] * math.cos(angle) + svar[branch] * math.sin(angle) ) )
    #     #expr_t += b*bshunt * ( cvar[buses[count_of_t]] - (1/ratio) * (cvar[branch] * math.cos(angle) + svar[branch] * math.sin(angle) ))
    #     #expr_t += (g*bshunt/ratio) * ( svar[branch] * math.cos(angle) - cvar[branch] * math.sin(angle) )
    #     #expr_t += ( bshunt*bshunt*cvar[buses[count_of_t]]/4 )
    # else:
    #     i2f += (g*g + b*b) * ( mp_cbusf + mp_cbust - 2 * mp_c ) + bshunt * ( mp_cbusf * ( (bshunt/4) + b ) - b * mp_c + g * mp_s )
    #     #expr_t += (g*g + b*b) * ( cvar[buses[count_of_f]] + cvar[buses[count_of_t]] - 2 * cvar[branch] ) + bshunt * ( cvar[buses[count_of_t]] * ( (bshunt/4) + b ) - b * cvar[branch] - g * svar[bra 
        
    return i2f



def add_cuts_ws(log,all_data):

    log.joint('\n')
    log.joint(' **** loading precomputed cuts ****\n')
    
    if '_b' in all_data['casename']:
        original_casename = all_data['casename'][:len(all_data['casename']) - 2]
    elif ('_n_5_5' in all_data['casename'] or '_n_0_5' in all_data['casename'] 
          or '_n_1_1' in all_data['casename'] or '_pline' in all_data['casename']):
        original_casename = all_data['casename'][:len(all_data['casename']) - 6]        
    #elif '_line5' in all_data['casename']:
    #    original_casename = all_data['casename'][:len(all_data['casename']) - 6]
    else:
        original_casename = all_data['casename']

    filename = 'newcuts/cuts_' + original_casename + '.txt' ###cuts | newcuts
    log.joint(" opening file with cuts " + filename + "\n")


    try:
        thefile = open(filename, "r") 
        lines = thefile.readlines()
    except:
        log.stateandquit(" cannot open file", filename)
        sys.exit("failure")

    numlines  = len(lines)
    theround  = lines[0].split()[3]
    firstline = lines[1].split()
    jabr      = 1
    
    log.joint(' loading cuts from round ' + str(theround) + '\n')

    if firstline[0] == '#Jabr-envelope':
        numjabrcuts = int(firstline[3])
        log.joint(' number of Jabr-envelope cuts in file = ' + str(numjabrcuts) + '\n')
    elif firstline[0] == '#i2-envelope':
        jabr      = 0
        i2        = 1
        numi2cuts = int(firstline[3])
        log.joint(' no Jabr-envelope cuts to add\n')
        log.joint(' number of i2-envelope cuts in file = ' + str(numi2cuts) + '\n')
    elif firstline[0] == '#limit-envelope':
        jabr         = 0
        i2           = 0
        numlimitcuts = int(firstline[3])
        all_data['addcuts_numlimitcuts'] = numlimitcuts
        log.joint(' no Jabr nor i2-envelope cuts to add\n')
        log.joint(' number of limit-envelope cuts in file = ' + str(numlimitcuts) + '\n')
    else:
        log.joint(' no cuts added\n')
        return None

    linenum = 2

    themodel       = all_data['themodel']
    buses          = all_data['buses']
    branches       = all_data['branches']
    IDtoCountmap   = all_data['IDtoCountmap']
    cvar           = all_data['cvar']
    svar           = all_data['svar']
    Pvar_f         = all_data['Pvar_f']
    Qvar_f         = all_data['Qvar_f']
    Pvar_t         = all_data['Pvar_t']
    Qvar_t         = all_data['Qvar_t']    
    FeasibilityTol = all_data['FeasibilityTol']
    T              = all_data['T']
    alphadic       = all_data['alphadic']

    if all_data['i2']:
        i2var_f    = all_data['i2var_f']

    jabr_cuts       = all_data['jabr_cuts']
    jabr_cuts_info  = all_data['jabr_cuts_info']

    if all_data['i2cuts']:
        i2_cuts         = all_data['i2_cuts']     
        i2_cuts_info    = all_data['i2_cuts_info']
        
    limit_cuts      = all_data['limit_cuts']
    limit_cuts_info = all_data['limit_cuts_info']

    cutid_jabr  = 0
    cutid_i2    = 0
    cutid_limit = 0

    numjabradded  = 0
    numi2added    = 0
    numlimitadded = 0
    
    while linenum < numlines: 
        thisline = lines[linenum].split()
        if thisline[0] == '#i2-envelope' and jabr:
            numi2cuts = int(thisline[3])
            all_data['addcuts_numi2cuts'] = numi2cuts
            log.joint(' number of i2-envelope cuts in file = '
                      + str(numi2cuts) + '\n')
            linenum += 1
            jabr     = 0
            i2       = 1
            continue

        elif thisline[0] == '#limit-envelope' and i2:
            numlimitcuts = int(thisline[3])
            log.joint(' number of limit-envelope cuts in file = '
                      + str(numlimitcuts) + '\n')
            linenum += 1
            i2       = 0
            continue

        elif jabr:
            branchid  = int(thisline[1])
            f         = int(thisline[3])
            t         = int(thisline[5])
            #cutid     = int(thisline[7])
            #rnd       = int(thisline[9])
            rnd       = 0
            violation = float(thisline[11])
            threshold = float(thisline[13])
            coeff_cft = float(thisline[15])
            coeff_sft = float(thisline[17])
            coeff_cff = float(thisline[19])
            coeff_ctt = float(thisline[21])


            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' + str(branchid) 
                          + ' f ' + str(f) + ' t ' + str(t) 
                          + ' was turned OFF\n')
                linenum += 1
                continue

            branch     = branches[branchid]
            count_of_f = IDtoCountmap[f] 
            count_of_t = IDtoCountmap[t]

            if (branchid != branch.count) or f != (branch.f) or (t != branch.t):
                log.joint(' branchid ' + str(branchid) + ' branch.count ' 
                          + str(branch.count) + ' f ' + str(f) + ' branch.f ' 
                          + str(branch.f) + ' t ' + str(t) + ' branch.t ' 
                          + str(branch.t) + '\n')
                breakexit('bug')


            if all_data['loud_cuts']:
                log.joint(' --> new Jabr-envelope cut for every time period\n')
                log.joint(' branch ' + str(branchid) + ' f ' + str(f) + ' t ' 
                          + str(t) + ' rnd ' + str(rnd)
                          + ' cutid ' + str(cutid_jabr) + '\n' )
                log.joint(' LHS coeff ' + ' cft ' + str(coeff_cft) + ' sft ' 
                          + str(coeff_sft) + ' cff ' + str(coeff_cff) 
                          + ' ctt ' + str(coeff_ctt) + '\n' )

            if all_data['jabr_validity']:
                for k in range(T):
                    sol_c       = all_data['sol_cvalues'][k][branch]
                    sol_s       = all_data['sol_svalues'][k][branch]
                    sol_cbusf   = all_data['sol_cvalues'][k][buses[count_of_f]]
                    sol_cbust   = all_data['sol_cvalues'][k][buses[count_of_t]]
                    sol_viol    = coeff_cft * sol_c + coeff_sft * sol_s + coeff_cff * sol_cbusf + coeff_ctt * sol_cbust
                    sol_relviol = sol_viol / ( ( coeff_cft**2 + coeff_sft**2 + coeff_cff**2 + coeff_ctt**2 )**0.5 ) 

                    if sol_relviol > FeasibilityTol:
                        log.joint(' WARNING, the Jabr-envelope cut associated '
                                  + 'to branch ' + str(branch.count) + ' k '
                                  + str(k) + ' f ' + str(branch.f)
                                  + ' t ' + str(branch.t)
                                  + ' is violated by the AC solution!\n')
                        log.joint(' violation ' + str(sol_viol) + '\n')
                        log.joint(' relative violation ' + str(sol_relviol)
                                  + '\n')
                        log.joint(' values (AC solution)' + ' cft '
                                  + str(sol_c) + ' sft ' + str(sol_s)
                                  + ' cff ' + str(sol_busf) + ' ctt '
                                  + str(sol_cbust) + '\n' )
                        breakexit('check!')
                    else:
                        log.joint(' AC solution satisfies Jabr inequality '
                                  + 'at branch ' + str(branch.count)
                                  + ' with slack ' + str(sol_relviol) + '\n')
                    
            
            for k in range(T):
                # Updating dictionaries
                jabr_cuts[(cutid_jabr,branchid)]   = (rnd,threshold,k)
                jabr_cuts_info[k][branch][cutid_jabr] = (rnd,violation,
                                                         coeff_cft,coeff_sft,
                                                         coeff_cff,coeff_ctt,
                                                         threshold,cutid_jabr)
                
                cutexp = LinExpr()
                constrname = "jabr_cut_"+str(cutid_jabr)+"_"+str(branchid)+"r_"+str(rnd)+"k_"+str(k)+"_"+str(f)+"_"+str(t)
                cutexp += (coeff_cft * cvar[k][branch]
                           + coeff_sft * svar[k][branch]
                           + coeff_cff * cvar[k][buses[count_of_f]]
                           + coeff_ctt * cvar[k][buses[count_of_t]])
                
                themodel.addConstr(cutexp <= 0, name = constrname)
                numjabradded += 1
                cutid_jabr += 1                
                
            linenum += 1

        elif (jabr == 0) and i2 and all_data['i2']:
            
            branchid   = int(thisline[1])
            f          = int(thisline[3])
            t          = int(thisline[5])
            #cutid      = int(thisline[7])
            #rnd        = int(thisline[9])
            rnd        = 0
            violation  = float(thisline[11])
            threshold  = float(thisline[13])
            coeff_Pft  = float(thisline[15])
            coeff_Qft  = float(thisline[17])
            coeff_cff  = float(thisline[19])
            coeff_i2ft = float(thisline[21])
            alpha      = all_data['alphadic'][branches[branchid]]
            
            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + 'was turned OFF\n')
                linenum += 1
                continue

            if alpha >= all_data['alpha_threshold']:
                if all_data['loud_cuts']:
                    log.joint(' bad coeffs: branch ' + str(branchid) + ' f '
                              + str(f) + ' t ' + str(t) + ' was skipped\n')
                linenum += 1
                continue
            
            branch     = branches[branchid]
            count_of_f = IDtoCountmap[f] 
            count_of_t = IDtoCountmap[t]

            if (branchid != branch.count) or f != (branch.f) or (t != branch.t):
                breakexit('there might be bug')
            

            if all_data['loud_cuts']:
                log.joint(' --> new i2-envelope cut for every time period\n')
                log.joint(' branch ' + str(branchid) + ' f ' + str(f) + ' t ' 
                          + str(t) + ' cutid ' + str(cutid_i2) + '\n' )
                log.joint(' LHS coeff ' + ' Pft ' + str(coeff_Pft) + ' Qft ' 
                          + str(coeff_Qft) + ' cff ' + str(coeff_cff) 
                          + ' i2ft ' + str(coeff_i2ft) + '\n' )

            if all_data['i2_validity']:
                for k in range(T):
                    
                    sol_Pf    = all_data['sol_Pfvalues'][k][branch]
                    sol_Qf    = all_data['sol_Qfvalues'][k][branch]
                    sol_c     = all_data['sol_cvalues'][k][branch]
                    sol_s     = all_data['sol_svalues'][k][branch]
                    sol_cbusf = all_data['sol_cvalues'][k][buses[count_of_f]]
                    sol_cbust = all_data['sol_cvalues'][k][buses[count_of_t]]
                    sol_i2f   = computei2value(log,all_data,branch,sol_c,
                                               sol_s,sol_cbusf,sol_cbust)   
                    sol_viol    = coeff_Pft * sol_Pf + coeff_Qft * sol_Qf + coeff_cff * sol_cbusf + coeff_i2ft * sol_i2f
                    sol_relviol = sol_viol / ( ( coeff_Pft**2 + coeff_Qft**2 + coeff_cff**2 + coeff_i2ft**2 )**0.5 )                     

                    if sol_relviol > FeasibilityTol:
                        log.joint(' WARNING, the i2-envelope cut associated '
                                  + 'to branch ' + str(branch.count) + ' k '
                                  + str(k) + ' f ' + str(branch.f)
                                  + ' t ' + str(branch.t)
                                  + ' is violated by the AC solution!\n')
                        log.joint(' violation ' + str(sol_viol) + '\n')
                        log.joint(' relative violation ' + str(sol_relviol)
                                  + '\n')
                        log.joint(' values (AC solution) ' + ' Pft '
                                  + str(sol_Pf) + ' Qft ' + str(sol_Qf)
                                  + ' cff ' + str(sol_cbusf) + ' i2ft '
                                  + str(sol_i2f) + '\n' )                 
                        breakexit('check!')
                    else:
                        log.joint(' AC solution satisfies the i2-cut '
                                  + 'at branch ' + str(branch.count)
                                  + ' with slack ' + str(sol_relviol) + '\n')

            for k in range(T):
                # Updating dictionaries            
                i2_cuts[(cutid_i2,branchid)]   = (rnd,threshold,k)
                i2_cuts_info[k][branch][cutid_i2] = (rnd,violation,coeff_Pft,
                                                     coeff_Qft,coeff_cff,
                                                     coeff_i2ft,threshold,cutid_i2)
                    
                cutexp     = LinExpr()
                constrname = "i2_cut_"+str(cutid_i2)+"_"+str(branchid)+"r_"+str(rnd)+"k_"+str(k)+"_"+str(f)+"_"+str(t)
                cutexp    += (coeff_Pft * Pvar_f[k][branch]
                              + coeff_Qft * Qvar_f[k][branch]
                              + coeff_cff * cvar[k][buses[count_of_f]]
                              + coeff_i2ft * i2var_f[k][branch])

                themodel.addConstr(cutexp <= 0, name = constrname)
                numi2added += 1
                cutid_i2   += 1 # updating i2-cut counter
                
            linenum += 1


        elif (jabr == 0) and (i2 == 0):
            branchid   = int(thisline[1])
            f          = int(thisline[3])
            t          = int(thisline[5])
            cutid      = int(thisline[7])
            if thisline[14] == 'Pft':
                from_or_to = 'f'
            elif thisline[14] == 'Ptf':
                from_or_to = 't'
            else:
                log.joint(' look for a bug\n')
                breakexit('bug')

            if branchid not in branches.keys(): # perturbed lines
                log.joint(' we do not add this cut since branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + 'was turned OFF\n')
                linenum += 1
                continue

            #rnd        = int(thisline[9])
            rnd        = 0
            violation  = float(thisline[11])
            threshold  = float(thisline[13])
            coeff_P    = float(thisline[15])
            coeff_Q    = float(thisline[17])

            branch     = branches[branchid]
            count_of_f = IDtoCountmap[f] 
            count_of_t = IDtoCountmap[t]

            if (branchid != branch.count) or f != (branch.f) or (t != branch.t):
                breakexit('there might be bug')

            if all_data['loud_cuts']:
                log.joint(' --> new limit-envelope cut for every time period\n')
                log.joint(' branch ' + str(branchid) + ' f ' + str(f) + ' t ' 
                          + str(t) + ' cutid ' + str(cutid_limit) + '\n' )
                if from_or_to == 'f':
                    log.joint(' LHS coeff ' + ' Pft ' + str(coeff_P) 
                              + ' Qft ' + str(coeff_Q) + '\n')
                elif from_or_to == 't':
                    log.joint(' LHS coeff ' + ' Ptf ' + str(coeff_P) 
                              + ' Qtf ' + str(coeff_Q) + '\n')
            
            if all_data['i2_validity']:
                for k in range(T):
                    
                    sol_Pf    = all_data['sol_Pfvalues'][k][branch]
                    sol_Qf    = all_data['sol_Qfvalues'][k][branch]
                    sol_c     = all_data['sol_cvalues'][k][branch]
                    sol_s     = all_data['sol_svalues'][k][branch]
                    sol_cbusf = all_data['sol_cvalues'][k][buses[count_of_f]]
                    sol_cbust = all_data['sol_cvalues'][k][buses[count_of_t]]
                    sol_i2f   = computei2value(log,all_data,branch,sol_c,
                                               sol_s,sol_cbusf,sol_cbust)   
                    sol_viol    = coeff_Pft * sol_Pf + coeff_Qft * sol_Qf + coeff_cff * sol_cbusf + coeff_i2ft * sol_i2f
                    sol_relviol = sol_viol / ( ( coeff_Pft**2 + coeff_Qft**2 + coeff_cff**2 + coeff_i2ft**2 )**0.5 )                     

                    if sol_relviol > FeasibilityTol:
                        log.joint(' WARNING, the i2-envelope cut associated '
                                  + 'to branch ' + str(branch.count) + ' k '
                                  + str(k) + ' f ' + str(branch.f)
                                  + ' t ' + str(branch.t)
                                  + ' is violated by the AC solution!\n')
                        log.joint(' violation ' + str(sol_viol) + '\n')
                        log.joint(' relative violation ' + str(sol_relviol)
                                  + '\n')
                        log.joint(' values (AC solution) ' + ' Pft '
                                  + str(sol_Pf) + ' Qft ' + str(sol_Qf)
                                  + ' cff ' + str(sol_cbusf) + ' i2ft '
                                  + str(sol_i2f) + '\n' )                 
                        breakexit('check!')
                    else:
                        log.joint(' AC solution satisfies the i2-cut '
                                  + 'at branch ' + str(branch.count)
                                  + ' with slack ' + str(sol_relviol) + '\n')
                    
            if all_data['limit_validity']:
                for k in range(T):
                    if from_or_to == 'f':
                        sol_P   = all_data['sol_Pfvalues'][k][branch]
                        sol_Q   = all_data['sol_Qfvalues'][k][branch]
                    elif from_or_to == 't':
                        sol_P   = all_data['sol_Ptvalues'][k][branch]
                        sol_Q   = all_data['sol_Qtvalues'][k][branch]
                    sol_viol    = coeff_P * sol_P + coeff_Q * sol_Q - 1

                if sol_viol > FeasibilityTol:
                    log.joint(' WARNING, the limit-envelope cut associated '
                              + 'to branch ' + str(branchid) + ' k '
                              + str(k) + ' f ' + str(f) + ' t ' + str(t)
                              + ' is violated by the AC solution!\n')
                    log.joint(' violation ' + str(sol_viol) + '\n')
                    
                    if from_or_to == 'f':
                        log.joint(' values (AC solution) ' + ' Pft ' 
                                  + str(sol_P) + ' Qft ' + str(sol_Q) + '\n')
                    elif from_or_to == 't':
                        log.joint(' values (AC solution) ' + ' Ptf ' 
                                  + str(sol_P) + ' Qtf ' + str(sol_Q) + '\n')
                    breakexit('check!')
                else:
                    log.joint(' AC solution satisfies limit inequality '
                              + 'at branch ' + str(branchid)
                              + ' with slack ' + str(sol_viol) 
                              + '\n')

            for k in range(T):            
                # Updating dictionaries
                limit_cuts[(cutid_limit,branchid)]   = (rnd,threshold,
                                                        from_or_to,k)
                limit_cuts_info[k][branch][cutid_limit] = (rnd,violation,
                                                           coeff_P,coeff_Q,
                                                           threshold,cutid_limit,
                                                           from_or_to)                

                cutexp     = LinExpr()
                if from_or_to == 'f':
                    constrname = "limit_cut_"+str(cutid_limit)+"_"+str(branchid)+"r_"+str(rnd)+"k_"+str(k)+"_"+str(f)+"_"+str(t)
                    cutexp    += (coeff_P * Pvar_f[k][branch]
                                  + coeff_Q * Qvar_f[k][branch])
                elif from_or_to == 't':
                    constrname = "limit_cut_"+str(cutid_limit)+"_"+str(branchid)+"r_"+str(rnd)+"k_"+str(k)+"_"+str(t)+"_"+str(f)
                    cutexp    += (coeff_P * Pvar_t[k][branch]
                                  + coeff_Q * Qvar_t[k][branch])
                themodel.addConstr(cutexp <= 1, name = constrname)
                numlimitadded += 1
                cutid_limit   += 1 # updating limit-cut counter                
            
            linenum += 1
    

    all_data['num_jabr_cuts']  = numjabradded
    all_data['num_i2_cuts']    = numi2added
    all_data['num_limit_cuts'] = numlimitadded
    all_data['ID_jabr_cuts']   = numjabradded
    all_data['ID_i2_cuts']     = numi2added
    all_data['ID_limit_cuts']  = numlimitadded
    
    all_data['addcuts_numjabrcuts']  = numjabradded
    all_data['addcuts_numi2cuts']    = numi2added
    all_data['addcuts_numlimitcuts'] = numlimitadded
        
    all_data['jabr_cuts']       = jabr_cuts
    all_data['jabr_cuts_info']  = jabr_cuts_info
    all_data['i2_cuts']         = i2_cuts
    all_data['i2_cuts_info']    = i2_cuts_info
    all_data['limit_cuts']      = limit_cuts
    all_data['limit_cuts_info'] = limit_cuts_info

    log.joint(' -- number of cuts added and propagated\n')
    log.joint(' Jabr-envelope cuts= ' + str(numjabradded)
              + '\n')
    log.joint(' i2-envelope cuts = ' + str(numi2added)
              + '\n')
    log.joint(' limit-envelope cuts = ' + str(numlimitadded)
              + '\n')    

def add_i2(log,all_data):
        
    log.joint('\n')
    log.joint(' **** new i2 vars ****\n')

    themodel       = all_data['themodel']
    buses          = all_data['buses']
    branches       = all_data['branches']    
    i2var_f        = all_data['i2var_f']
    newi2          = all_data['newi2']

    for branch in newi2:

      branchcount = branch.count
      f           = branch.f
      t           = branch.t
      count_of_f  = IDtoCountmap[f]
      count_of_t  = IDtoCountmap[t]
      bus_f       = buses[count_of_f]

      upperbound_f = branch.limit**2 / (bus_f.Vmin * bus_f.Vmin)

      i2var_f[branch] = themodel.addVar(obj = 0.0, lb = 0, ub = upperbound_f ,
                                        name = "i2_" + str(branch.count) + "_"
                                        + str(f) + "_" + str(t))

    
    all_data['i2var_f']  = i2var_f
    all_data['themodel'] = themodel

    themodel.update()
    log.joint(' model updated with new i2 variables\n')


