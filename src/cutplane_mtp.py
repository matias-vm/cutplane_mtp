###############################################################################
##                                                                           ##
## This code was written and is being maintained by Matias Villagra,         ##
## PhD Student in Operations Research @ Columbia, supervised by              ##
## Daniel Bienstock.                                                         ##
##                                                                           ##
## Please report any bugs or issues (for sure there will be) to              ##
##                         mjv2153@columbia.edu                              ##
##                                                                           ##
## Oct 2023                                                                  ##
###############################################################################

import sys
import math
from log import danoLogger
from gurobipy import *
import numpy as np
import bisect
from myutils import breakexit
import reader
import time
import math
from cuts_mtp import *
from flowdecomp import *
from tools import *
import random
import os
import platform
from random import randint

def gocutplane(log, all_data):

  formulation_start = time.time()
  themodel          = Model("Cutplane")
  buses             = all_data['buses']
  numbuses          = all_data['numbuses']
  branches          = all_data['branches']
  numbranches       = all_data['numbranches']
  gens              = all_data['gens']
  IDtoCountmap      = all_data['IDtoCountmap']
  FeasibilityTol    = all_data['FeasibilityTol']
  threshold         = all_data['threshold']
  T                 = all_data['T']
  casename          = all_data['casename']
  casetype          = all_data['casetype']
  loadsfilename     = all_data['loadsfilename']
  
  all_data['alpha_threshold'] = 1e2
  
  ############################ LOAD SOLUTION ##################################

  if all_data['matpower_sol']:
    getsol_matpower(log,all_data)

  if all_data['ampl_sol']:
    getsol_ampl_mtp(log,all_data)
    #getsol_ampl(log,all_data)

  ################################ VARIABLES ##################################

  cvar    = {}
  svar    = {}
  Pvar_f  = {}
  Qvar_f  = {}
  Pvar_t  = {}
  Qvar_t  = {}
  Pinjvar = {}
  Qinjvar = {}
  GenPvar = {}
  GenQvar = {}
  GenTvar = {}

  # File with mtploads
  log.joint('the name of mtp file ' + loadsfilename + '\n')
  loads = open(loadsfilename,"r")
 
  if all_data['arpae']:
    Pd = getloads2(log,all_data,loads)
    all_data['Pd'] = Pd
    #all_data['Qd'] = Qd
  else:
    Pd = getloads(log,all_data,loads)
    all_data['Pd'] = Pd
    
  log.joint(" Loads obtained\n")
  
  rampfilename = '../../ampl_acopf/data/ramprates/' + casename + '_rampr_' + str(T) + '.txt'
  
  try :
    rampr = open(rampfilename,"r")
    log.joint(" Opening ramprates file\n")
    all_data['newrampr'] = 1
    rampru, ramprd       = getrampr(log,all_data,rampr)
    log.joint(" Ramp rates obtained\n")
  except:
    log.joint(" File with ramp rates does not exist, so we create it\n")
    log.joint(" File name " + rampfilename + "\n")    
    rampr               = open(rampfilename,"w+")
    all_data['newrampr'] = 0
    rampru               = {}
    ramprd               = {}    
  
  all_data['rampru'] = rampru
  all_data['ramprd'] = ramprd  
  
  for k in range(T):
    cvar[k]    = {}
    svar[k]    = {}
    Pvar_f[k]  = {}
    Qvar_f[k]  = {}
    Pvar_t[k]  = {}
    Qvar_t[k]  = {}
    Pinjvar[k] = {}
    Qinjvar[k] = {}
    GenPvar[k] = {}
    GenQvar[k] = {}
    GenTvar[k] = {}

  if all_data['newrampr'] == 0:
    for k in range(T):
      rampru[k] = {}
      ramprd[k] = {}

  log.joint(' %d time periods\n' %T)      
  log.joint(' creating variables...\n')

  varcount = 0

  #cbus, bus-injection, GenP, GenQ, and GenT variables
  for bus in buses.values():

    maxprod = bus.Vmax*bus.Vmax
    minprod = bus.Vmin*bus.Vmin
          
    ubound = maxprod
    lbound = minprod

    for k in range(T):
      Pubound, Plbound, Qubound, Qlbound = computebalbounds(log,all_data,bus,k)
      
      cvar[k][bus] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                     name = "c_" + str(bus.nodeID) + "_" 
                                     + str(bus.nodeID) + "_" + str(k))
      Pinjvar[k][bus] = themodel.addVar(obj = 0.0, lb = Plbound, ub = Pubound, 
                                        name = "IP_"+str(bus.nodeID) + "_" + str(k))
      Qinjvar[k][bus] = themodel.addVar(obj = 0.0, lb = Qlbound, ub = Qubound, 
                                        name = "IQ_"+str(bus.nodeID) + "_" + str(k))
    
      varcount += 3

    for genid in bus.genidsbycount:
      gen = gens[genid]

      lower = gen.Pmin*gen.status
      upper = gen.Pmax*gen.status

      if all_data['newrampr'] == 0:
        for k in range(T):
          rampru[k][gen] = 0.5
          ramprd[k][gen] = 0.5
          rampr.write("gen " + str(genid) + " bus " + str(gen.nodeID) + " k "
                      + str(k) + " rpu " + str(rampru[k][gen]) + " rpd "
                      + str(ramprd[k][gen]) + " GPmax " + str(upper) + "\n")
          GenPvar[k][gen] = themodel.addVar(obj = 0.0, lb = lower, ub = upper, 
                                            name = "GP_" + str(gen.count) + "_" 
                                            + str(gen.nodeID) + "_" + str(k))
          varcount += 1
      else:
        for k in range(T):
          GenPvar[k][gen] = themodel.addVar(obj = 0.0, lb = lower, ub = upper, 
                                            name = "GP_" + str(gen.count) + "_" 
                                            + str(gen.nodeID) + "_" + str(k))
          varcount += 1
      
      lower = gen.Qmin*gen.status
      upper = gen.Qmax*gen.status

      if bus.nodetype == 3:
        upper = GRB.INFINITY
        lower = -GRB.INFINITY

      for k in range(T):
        GenQvar[k][gen] = themodel.addVar(obj = 0.0, lb = lower, ub = upper, 
                                          name = "GQ_" + str(gen.count) + "_" 
                                          + str(gen.nodeID) + "_" + str(k))
        varcount += 1

        #if model has quadratic objective and we linearize
      # if( gen.costdegree == 2 and gen.costvector[0] != 0 and 
      #     (all_data['linear_objective'] or all_data['hybrid']) ):

      #   for k in range(T):
      #     GenTvar[k][gen] = themodel.addVar(obj = 0.0, lb = 0.0, ub = GRB.INFINITY,
      #                                       name = 't_g_' + str(gen.count) + '_' 
      #                                       + str(gen.nodeID))         
      #     varcount += 1


  if all_data['newrampr'] == 0:
    rampr.write("END\n")
    
  rampr.close()
  #exit(0) #
  
          
  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]
    maxprod     = buses[count_of_f].Vmax*buses[count_of_t].Vmax
    minprod     = buses[count_of_f].Vmin*buses[count_of_t].Vmin

    ubound      = maxprod
    lbound      = -maxprod
    maxanglerad = branch.maxangle_rad
    minanglerad = branch.minangle_rad

    # Cosine                                                                                                                                           
    if maxanglerad <= 0.5*math.pi:
      # In this case minangle <= 0                                                                                                                   
      if minanglerad >= -0.5*math.pi:
        lbound = minprod*min(math.cos(maxanglerad), math.cos(minanglerad))
      elif minanglerad >= -math.pi:
        lbound = maxprod*math.cos(minangle_rad)  # Which is negative                                                                               
      elif minanglerad >= -1.5*math.pi:
        lbound = -maxprod
      else:
        lbound = -maxprod

    elif maxanglerad <= math.pi:
      if minanglerad >= -0.5*math.pi:
        lbound = maxprod*math.cos(maxanglerad)
      elif minanglerad >= -math.pi:
        lbound = maxprod*min(math.cos(maxanglerad), math.cos(minanglerad))
      elif minanglerad >= -1.5*math.pi:
        lbound = -maxprod
      else:
        lbound = -maxprod

    elif maxanglerad <= 1.5*math.pi:
      lbound = -maxprod    

    elif maxanglerad <= 2*math.pi:
      lbound = -maxprod

    else:
      ubound = maxprod
      lbound = -maxprod

    for k in range(T):
      cvar[k][branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                        name = "c_" + str(branchcount) + "_" 
                                        + str(f) + "_" + str(t) + "_" + str(k))
      varcount += 1
      
    #s variables

    if maxanglerad <= math.pi/2:
      ubound = maxprod*math.sin(maxanglerad)

      if  minanglerad >= -0.5*math.pi:
        lbound = maxprod*math.sin(minanglerad)
      elif  minanglerad >= -math.pi:
        lbound = -maxprod
      elif  minanglerad >= -1.5*math.pi:
        ubound = maxprod*max( math.sin(maxanglerad), math.sin(minanglerad))
        lbound = -maxprod
      else:
        ubound = maxprod
        lbound = -maxprod

    elif maxanglerad <= math.pi:
      ubound = maxprod

      if minanglerad >= -0.5*math.pi:
        lbound = maxprod*math.sin(minanglerad)
      elif minanglerad >= -math.pi:
        lbound = -maxprod
      elif minanglerad >= -1.5*math.pi:
        lbound = -maxprod
      else:
        lbound = -maxprod

    elif maxanglerad <= 1.5*math.pi:
      ubound = maxprod

      if minanglerad >= -0.5*math.pi:
        lbound = maxprod*min(math.sin(maxanglerad), math.sin(minanglerad))
      elif minanglerad >= -math.pi:
        lbound = -maxprod
      elif minanglerad >= -1.5*math.pi:
        lbound = -maxprod
      else:
        lbound = -maxprod
    else:
      ubound = maxprod
      lbound = -maxprod

    for k in range(T):
      svar[k][branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                        name = "s_" + str(branchcount) + "_" 
                                        + str(f) + "_" + str(t) + "_" + str(k))

      varcount += 1
    

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]

    ubound = branch.limit 
    lbound = -branch.limit

    for k in range(T):
      Pvar_f[k][branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                          name = "P_" + str(branch.count) + "_" 
                                          + str(f) + "_" + str(t) + "_" + str(k))
      Pvar_t[k][branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                          name = "P_" + str(branch.count) + "_" 
                                          + str(t) + "_" + str(f) + "_" + str(k))
      Qvar_f[k][branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                          name = "Q_" + str(branch.count) + "_" 
                                          + str(f) + "_" + str(t) + "_" + str(k))
      Qvar_t[k][branch] = themodel.addVar(obj = 0.0, lb = lbound, ub = ubound, 
                                          name = "Q_" + str(branch.count) + "_" 
                                          + str(t) + "_" + str(f) + "_" + str(k))

      varcount +=4 

      
  #i2 variables
  if all_data['i2']:
    i2var_f = {}
    #i2var_t = {}

    alphadic = all_data['alphadic'] = {}
    
    for k in range(T):
      i2var_f[k] = {}
      

    for branch in branches.values():
      branchcount = branch.count
      f           = branch.f
      t           = branch.t
      count_of_f  = IDtoCountmap[f]
      count_of_t  = IDtoCountmap[t]
      ratio       = branch.ratio
      y           = branch.y
      g           = y.real
      b           = y.imag
      bshunt      = branch.bc

      alpha            = ( g*g + b*b + bshunt * (b + (bshunt/4)) ) / (ratio**4)
      alphadic[branch] = alpha

      if alpha < all_data['alpha_threshold']:
        bus_f        = buses[count_of_f]
        upperbound_f = branch.limit**2 / (bus_f.Vmin * bus_f.Vmin)
        #upperbound_t = branch.limit**2 / (bus_t.Vmin * bus_t.Vmin)
        
        for k in range(T):
          i2var_f[k][branch] = themodel.addVar(obj = 0.0, lb = 0, ub = upperbound_f ,
                                               name = "i2_" + str(branch.count) + "_"
                                               + str(f) + "_" + str(t) + "_" + str(k))

          #i2var_t[branch] = themodel.addVar(obj = 0.0, lb = 0, ub = upperbound_t,
          #name = "i2_" + str(branch.count) + "_"
          #+ str(t) + "_" + str(f))
    
          varcount += 1
      
  themodel.update()
  log.joint('   %d variables added\n' %varcount)

  all_data['themodel']    = themodel
  all_data['cvar']        = cvar
  all_data['svar']        = svar
  all_data['GenPvar']     = GenPvar
  all_data['GenTvar']     = GenTvar
  all_data['Pvar_f']      = Pvar_f
  all_data['Pvar_t']      = Pvar_t
  all_data['Qvar_f']      = Qvar_f
  all_data['Qvar_t']      = Qvar_t
  all_data['Pd']          = Pd

  
  if all_data['i2']:
    all_data['i2var_f']   = i2var_f
    #all_data['i2var_t']  = i2var_f

  ############################## OBJECTIVE ####################################

  log.joint(' creating objective...\n')
  varobjcount = 0

  # Constant term

  constexpr   = LinExpr()
  constobjval = 0
  
  for gen in gens.values():  
    if gen.status > 0:
      constobjval += gen.costvector[gen.costdegree]
  
  constvar   = themodel.addVar(lb = 1.0, ub = 1.0,name = "constant")
  constexpr += T * constobjval * constvar
  
  varobjcount +=1

  # Linear and Quadratic terms
  lincostexpr = LinExpr()
  qcostexpr   = QuadExpr()
  
  for gen in gens.values():
    lincoeff = gen.costvector[gen.costdegree-1]
    
    for k in range(T):
      lincostexpr += lincoeff * GenPvar[k][gen]
      varobjcount +=1
      
      if gen.costdegree == 2 and gen.costvector[0] != 0:
        qcostexpr += gen.costvector[0]*GenPvar[k][gen]*GenPvar[k][gen]
        varobjcount +=1

  log.joint('   %d terms in the objective\n' %varobjcount)
  
  themodel.setObjective(constexpr + lincostexpr + qcostexpr)
  
  # Old hybrid code
  #quad terms
  #objvar   = themodel.addVar(obj = 1.0, lb = -GRB.INFINITY, ub = GRB.INFINITY,
  #                           name = "objvar")
  #qcostvar = themodel.addVar(obj = 1.0, lb = 0, ub = GRB.INFINITY, 
  #                           name = "qcostvar") 
  #varobjcount +=2
  
  # if all_data['linear_objective'] == 0 or all_data['hybrid']:  
  #   objvar.setAttr("Obj",0)
  #   qcostexpr = QuadExpr()
  #   for gen in gens.values():
  #     for k in range(T):
  #       if gen.costdegree == 2 and gen.costvector[0] != 0:
  #         qcostexpr += gen.costvector[0]*GenPvar[k][gen]*GenPvar[k][gen]
  #   qcost = themodel.addConstr(qcostexpr <= qcostvar, name = "qcost")
  #   all_data['qcost'] = qcost
    
  # #linear terms
  # if all_data['linear_objective'] or all_data['hybrid']:
  #   lincostvar = themodel.addVar(obj = 0.0, lb = -GRB.INFINITY, 
  #                                ub = GRB.INFINITY, name = "lincostvar")
  #   lincost    = themodel.addConstr(lincostvar <= objvar, name= "lincost")
  #   all_data['lincost'] = lincost
  # else:
  #   lincostvar = themodel.addVar(obj = 1.0, lb = -GRB.INFINITY, 
  #                                ub = GRB.INFINITY, name = "lincostvar")

  # #varobjcount += 1

  # lincostexpr = LinExpr()
  # coeff       = [gen.costvector[gen.costdegree-1] for gen in gens.values()]
  
  # for k in range(T):
  #   variables    = [GenPvar[k][gen] for gen in gens.values()]
  #   lincostexpr += LinExpr(coeff, variables)
    
  # lincostdef  = themodel.addConstr(lincostexpr == lincostvar, 
  #                                  name= "lincostdef")

  # #sumTvars if using linear_objective
  # if all_data['linear_objective'] or all_data['hybrid']:
  #   qcostvar.setAttr("Obj",0.0)
  #   sumTvars = LinExpr()
  #   for gen in gens.values():
  #     if gen.costdegree == 2 and gen.costvector[0] != 0:
  #       sumTvars += GenTvar[gen]
  #   sumTconstr = themodel.addConstr(sumTvars <= objvar, name='obj_var_quad')
  #   all_data['sumTconstr'] = sumTconstr

  # if all_data['hybrid']: #we start with full objective
  #   objvar.setAttr("Obj",0)
  #   lincostvar.setAttr("Obj",1)
  #   qcostvar.setAttr("Obj",1)
  #   themodel.remove(sumTconstr)
  #   themodel.remove(lincost)
    
  # if all_data['linear_objective']:
  #   log.joint('  linear objective added\n')
  # elif all_data['hybrid']:
  #   log.joint('  hybrid algorithm, linear and quadratic objectives created\n')
  # else:
  #   log.joint('  quadratic objective added\n')

  themodel.update()
  
  #all_data['objvar']     = objvar
  #all_data['lincostvar'] = lincostvar
  #all_data['qcostvar']   = qcostvar

  ############################# CONSTRAINTS ###################################

  log.joint(' creating constraints...\n')

  constrcount = 0
  count       = 0

  #definition flow variables
  log.joint('  active power flow variables definition\n')

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]

    if branch.status == 0:
      log.joint(' branch ' + str(branch.count) + ' f ' + str(f) + ' t ' 
                + str(t) + ' is OFF\n')
      breakexit('check, reader does not include branches with status = 0')

    for k in range(T):
      #  Gff cff + Gft cft + Bft sft
      constrname = "Pdef_"+str(branch.count)+"_"+str(f)+"_"+str(t)+"_"+str(k)
      expr = LinExpr()
      expr += branch.Gff*cvar[k][buses[count_of_f]]
      expr += branch.Gft*cvar[k][branch]
      expr += branch.Bft*svar[k][branch]
    
      themodel.addConstr(expr == Pvar_f[k][branch], name = constrname)

      #  Gtt ctt + Gtf cft + Btf stf = Gtt ctt + Gtf cft - Btf sft
      constrname = "Pdef_"+str(branch.count)+"_"+str(t)+"_"+str(f)+"_"+str(k)
      expr = LinExpr()
      expr += branch.Gtt*cvar[k][buses[count_of_t]]
      expr += branch.Gtf*cvar[k][branch]
      expr += -branch.Btf*svar[k][branch]

      themodel.addConstr(expr == Pvar_t[k][branch], name = constrname)
    
      constrcount += 2
      count       += 2

  log.joint('   %d active power flow definition constraints added\n'%count)

  log.joint('  reactive power flow variables definition\n')
  count = 0

  for branch in branches.values():
    f = branch.f
    t = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]

    for k in range(T):
      # -Bff cff - Bft cft + Gft sft    
      constrname = "Qdef_"+str(branch.count)+"_"+str(f)+"_"+str(t)+"_"+str(k)
      expr = LinExpr()
      expr += -branch.Bff*cvar[k][buses[count_of_f]]
      expr += -branch.Bft*cvar[k][branch]
      expr += +branch.Gft*svar[k][branch]

      themodel.addConstr(expr == Qvar_f[k][branch], name = constrname)

      # -Btt ctt - Btf cft + Gtf stf = -Btt ctt - Btf cft - Gtf sft 
      constrname = "Qdef_"+str(branch.count)+"_"+str(t)+"_"+str(f)+"_"+str(k)
      expr = LinExpr()
      expr += -branch.Btt*cvar[k][buses[count_of_t]]
      expr += -branch.Btf*cvar[k][branch]
      expr += -branch.Gtf*svar[k][branch]

      themodel.addConstr(expr == Qvar_t[k][branch], name = constrname)

      constrcount += 2
      count       += 2
  log.joint('   %d reactive power flow definition constraints added\n'%count)

  #balance constraints
  log.joint('  active power injection constraints\n')
  count = 0

  for bus in buses.values():

    for k in range(T):
      constrname = "PBaldef"+str(bus.nodeID)+"_"+str(k)
      expr = LinExpr()
      for branchid in bus.frombranchids.values():
        expr += Pvar_f[k][branches[branchid]]
      for branchid in bus.tobranchids.values():
        expr += Pvar_t[k][branches[branchid]]

      if ( (bus.Gs != 0) and ( ( len(bus.frombranchids) != 0 ) 
                               or ( len(bus.tobranchids) != 0 ) ) ):
        expr += bus.Gs*cvar[k][bus]

      themodel.addConstr(expr == Pinjvar[k][bus], name = constrname)

      constrcount += 1
      count       += 1
  
  log.joint('   %d active power injection constraints added\n'%count)
  log.joint('  reactive power injection constraints\n')
  count = 0

  for bus in buses.values():
    for k in range(T):
      constrname = "QBaldef"+str(bus.nodeID)+"_"+str(k)
      expr = LinExpr()
      for branchid in bus.frombranchids.values():
        expr += Qvar_f[k][branches[branchid]]
      for branchid in bus.tobranchids.values():
        expr += Qvar_t[k][branches[branchid]]
 
      if ( (bus.Bs != 0) and ( ( len(bus.frombranchids) != 0 ) 
                               or ( len(bus.tobranchids) != 0 ) ) ):
        expr += (-bus.Bs)*cvar[k][bus]

      themodel.addConstr(expr == Qinjvar[k][bus], name = constrname)
    
      constrcount += 1
      count       += 1
      
  log.joint('   %d reactive power injection constraints added\n'%count)

  #definition bus-injection variables

  log.joint('  adding injection definition constraints...\n')

  count = 0

  for bus in buses.values():
    for k in range(T):
      constrname = "Bus_PInj_"+str(bus.nodeID)+"_"+str(k)
      expr = LinExpr()

      if len(bus.genidsbycount) > 0:
        for genid in bus.genidsbycount:
          gen = gens[genid]
          expr += GenPvar[k][gen]

      themodel.addConstr(Pinjvar[k][bus] == expr - Pd[k][bus], name = constrname)

      constrname = "Bus_QInj_"+str(bus.nodeID)+"_"+str(k)
      expr = LinExpr()

      if len(bus.genidsbycount) > 0:
        for genid in bus.genidsbycount:
          gen = gens[genid]
          expr += GenQvar[k][gen]

      if all_data['arpae2']:
        themodel.addConstr(Qinjvar[k][bus] == expr - all_data['Qd'][k][bus], name = constrname)
      else:
        themodel.addConstr(Qinjvar[k][bus] == expr - bus.Qd, name = constrname)

      constrcount += 2
      count       += 2

  log.joint('   %d power injection definitions added\n'%count)

  #ramp-up contraints

  log.joint('  adding ramp constraints...\n')

  count   = 0
  abs_gen = {}
  for k in range(T-1):
    abs_gen[k] = {}

  for gen in gens.values():

    genid  = gen.count
    nodeID = gen.nodeID

    for k in range(T-1):
      constrname_rup   = "rup_" + str(genid) + "_" + str(nodeID) + "_" + str(k) + "_" + str(k+1)
      constrname_rdown = "rdown_" + str(genid) + "_" + str(nodeID) + "_" + str(k) + "_" + str(k+1)      
      rpu        = rampru[k][gen]
      rpd        = ramprd[k][gen]

      #themodel.addConstr(GenPvar[k+1][gen] - (1 + rpu) * GenPvar[k][gen] <= 0, name = constrname_rup)
      #themodel.addConstr((1 - rpd) * GenPvar[k][gen] - GenPvar[k+1][gen] <= 0, name = constrname_rdown)
      #count += 2
      
      abs_gen[k][gen] = themodel.addVar(lb = 0, name = 'abs_gen_'+str(genid)+'_'+str(k))
      themodel.addConstr(GenPvar[k+1][gen] - GenPvar[k][gen] - rpu * abs_gen[k][gen] <= 0, name = constrname_rup)
      themodel.addConstr(GenPvar[k][gen] - rpd * abs_gen[k][gen] - GenPvar[k+1][gen] <= 0, name = constrname_rdown)

      themodel.addConstr(GenPvar[k][gen] - abs_gen[k][gen] <= 0, name = constrname_rup + '_1')
      themodel.addConstr(- GenPvar[k][gen] - abs_gen[k][gen] <= 0, name = constrname_rup + '_2')      

      count += 4
  
  #definition i2 variables
  if all_data['i2']:
    constrcount += i2_def(log,all_data)

  #active power loss-inequalities
  if all_data['loss_inequalities']:    
    constrcount += loss_inequalities(log,all_data)

  #reactive power loss-inequalities
  if all_data['qloss_inequalities']:    
    constrcount += qloss_inequalities(log,all_data)

  #jabr inequalities
  if all_data['jabr_inequalities']:
    constrcount += jabr_inequalities(log,all_data)

  #i2 inequalities
  if all_data['i2_inequalities']:
    constrcount += i2_inequalities(log,all_data)

  #limit constraints
  if all_data['limit_inequalities']:
    constrcount += limit_inequalities(log,all_data)
  
  log.joint('  %d constraints added\n'%constrcount)
    
  themodel.update()

  formulation_end = time.time()

  all_data['formulation_time'] = formulation_end - formulation_start
  all_data['numvars']          = themodel.NumVars
  all_data['numconstrs']       = themodel.NumConstrs
  
  log.joint(' formulation time: %g\n' % all_data['formulation_time'])
  log.joint(' numvars ' + str(all_data['numvars']) + ' numconstrs ' + str(all_data['numconstrs']) + '\n')
  
  # Write model to .lp file
  if all_data['writelps']:
    log.joint(' writing to lpfile ' + all_data['lpfilename'] + '\n')  
    themodel.write(all_data['lpfilename'])

    #breakexit('writelp without cuts, then remove break')

  #themodel.write("check.lp")
  #breakexit("c")
  
  ###################### INIT DATA STRUCTURES FOR CUTS ########################

  if all_data['jabrcuts']:
    jabr_cuts_info = all_data['jabr_cuts_info']
    for k in range(T):
      jabr_cuts_info[k] = {}
      for branch in branches.values():
        jabr_cuts_info[k][branch] = {}

  if all_data['i2cuts']:
    i2_cuts_info = all_data['i2_cuts_info']
    for k in range(T):
      i2_cuts_info[k] = {}
      for branch in branches.values():
        i2_cuts_info[k][branch] = {}

  if all_data['limitcuts']:
    limit_cuts_info = all_data['limit_cuts_info']
    for k in range(T):
      limit_cuts_info[k] = {}
      for branch in branches.values():
        limit_cuts_info[k][branch] = {}

    for branch in branches.values():
      limit_cuts_info[branch] = {}

  if all_data['linear_objective'] or all_data['hybrid']:
    dicGenPvalues = {}
    for gen in gens.values():
      dicGenPvalues[gen] = []
    all_data['dicGenPvalues'] = dicGenPvalues

  ######################## FIXING/WRITING AN AC SOLUTION ######################

  if all_data['fixflows']:
    fixflows(log,all_data)
    if all_data['fixcs'] == 0:
      return None

  if all_data['fixcs']:
    fixcs(log,all_data)
    return None

  if all_data['writeACsol']:
    writeACsol(log,all_data)
    return None

  ########################### SOLVER PARAMETERS ###############################

  themodel.Params.Method    = all_data['solver_method']
  themodel.Params.Crossover = all_data['crossover'] 
  themodel.Params.LogFile   = all_data['mylogfile']
  themodel.Params.TimeLimit = all_data['max_time']

  if all_data['solver_method'] == 2:
    themodel.Params.BarHomogeneous = 1
    themodel.Params.BarConvTol     = all_data['barconvtol']
    themodel.Params.FeasibilityTol = all_data['feastol']
    themodel.Params.OptimalityTol  = all_data['opttol']
    
  themodel.Params.NumericFocus = 1
  themodel.Params.OutPutFlag = 1
  
  ######################### READING AND LOADING CUTS ##########################

  if all_data['addcuts']:

    t0_cuts = time.time()

    if all_data['max_rounds'] > 1:
      add_cuts_ws(log,all_data)
    else:
      add_cuts(log,all_data)

    themodel.update()

    t1_cuts = time.time()

    all_data['addcuts_time'] = t1_cuts - t0_cuts

    log.joint(' pre-computed cuts added and model updated\n')

    log.joint(' reading and loading cuts time = '
              + str(all_data['addcuts_time']) + '\n')

    if all_data['writelps']:
      themodel.write(all_data['casename']+'_precomputed_cuts.lp')
      log.joint(' model with precomputed written to .lp file\n\n')
      #breakexit('precomputed cuts')

  
  ########################## CUTPLANE MAIN LOOP ###############################

  all_data['round']                  = 1
  all_data['runtime']                = time.time() - all_data['T0']
  all_data['round_time']             = time.time()
  all_data['cumulative_solver_time'] = 0
  all_data['ftol_counter']           = 0
  oldobj                             = 1
  gap                                = 1e20

  #  all_data['losstest'] 
  all_data['losstest_dic']        = {}
  all_data['losstest_loss']       = {}
  all_data['losstest_branchloss'] = {}

  while ((all_data['round'] <= all_data['max_rounds']) and 
         (all_data['runtime'] <= all_data['max_time']) and 
         (all_data['ftol_counter'] <= all_data['ftol_iterates'])):
    
    
    # #new tolerances
    # if ((all_data['ftol_counter'] == 4) or (all_data['max_rounds'] == 1) 
    #     or (all_data['runtime'] > 150)):

    #   themodel.Params.BarHomogeneous = 1
    #   themodel.Params.NumericFocus   = 1 #off and then on doesnt help
    #   themodel.Params.BarConvTol     = 1e-6
    #   themodel.Params.FeasibilityTol = 1e-6
    #   themodel.Params.OptimalityTol  = 1e-6


    ########################### HYBRID ALGORITHM ##############################

    if all_data['hybrid']:
      hybrid(log,all_data)
  
    ############################ SOLVING MODEL ################################

    cutplane_optimize(log,all_data)

    ########################### STORING SOLUTION ##############################

    log.joint(' storing current solution ...\n')

    all_data['Pfvalues']   = {}
    all_data['Qfvalues']   = {}
    all_data['Ptvalues']   = {}
    all_data['Qtvalues']   = {}
    all_data['cvalues']    = {}
    all_data['svalues']    = {}
    all_data['GenPvalues'] = {}
    all_data['GenQvalues'] = {}
    
    for k in range(T):
      Pfvalues   = {}
      Qfvalues   = {}
      Ptvalues   = {}
      Qtvalues   = {}
      cvalues    = {}
      svalues    = {}
      GenPvalues = {}
      GenQvalues = {}

      for bus in buses.values():
        cvalues[bus] = cvar[k][bus].X
      
      for branch in branches.values():
        Pfvalues[branch]   = Pvar_f[k][branch].X
        Qfvalues[branch]   = Qvar_f[k][branch].X
        Ptvalues[branch]   = Pvar_t[k][branch].X
        Qtvalues[branch]   = Qvar_t[k][branch].X
        cvalues[branch]    = cvar[k][branch].X
        svalues[branch]    = svar[k][branch].X

      for gen in gens.values():
        GenPvalues[gen] = GenPvar[k][gen].X
        GenQvalues[gen] = GenQvar[k][gen].X        

      all_data['Pfvalues'][k]   = Pfvalues
      all_data['Qfvalues'][k]   = Qfvalues
      all_data['Ptvalues'][k]   = Ptvalues
      all_data['Qtvalues'][k]   = Qtvalues
      all_data['cvalues'][k]    = cvalues
      all_data['svalues'][k]    = svalues
      all_data['GenPvalues'][k] = GenPvalues
      all_data['GenQvalues'][k] = GenQvalues      
      
    
    if all_data['i2']:
      all_data['i2fvalues'] = {}
      for k in range(T):
        i2fvalues = {}
        for branch in branches.values():   # CHANGE this, loop over selected br
          alpha = all_data['alphadic'][branch]
          if alpha < all_data['alpha_threshold']:
            i2fvalues[branch] = i2var_f[k][branch].X
        all_data['i2fvalues'][k] = i2fvalues
        
    if all_data['linear_objective'] or all_data['hybrid']:
      for gen in gens.values():
        if gen.costdegree == 2 and gen.costvector[0] != 0:
          gen_values = dicGenPvalues[gen]
          bisect.insort(gen_values,GenPvar[gen].x)

      all_data['dicGenPvalues'] = dicGenPvalues

    log.joint(' done storing values\n')
     
    ########################## CHECK OBJ IMPROVEMENT ##########################

    if ((all_data['objval'] - oldobj)/oldobj) < all_data['ftol']:
      all_data['ftol_counter'] += 1
    else:
      all_data['ftol_counter'] = 0

    oldobj              = all_data['objval']
    all_data['runtime'] = time.time() - all_data['T0']

    ########################### ROUND STATISTICS ##############################

    cutplane_stats(log,all_data)

    ########################### Flow Decomposition ############################

    # log.joint(' ploss ' + str(sum(plossvalues.values())) + '\n')

    # if all_data['writesol']:
    #   all_data['sol_Pfvalues']   = all_data['Pfvalues']
    #   all_data['sol_Ptvalues']   = all_data['Ptvalues']
    #   all_data['sol_GenPvalues'] = {}

    #   for gen in gens.values():
    #     all_data['sol_GenPvalues'][gen.count] = all_data['GenPvalues'][gen]
    
    #   log.joint(' sol_GenPvalues ' + str(all_data['sol_GenPvalues']) + '\n') 

    #   flowdecomp(log,all_data)

    ######################### SUMMARY EXPERIMENTS #############################

    all_data['runtime'] = time.time() - all_data['T0']

    log.joint("\n writing casename, opt stauts, obj and " +
              "runtime to summary_ws.log\n")

    summary_ws   = open("summary_ws.log","a+") #later add feasibility errors
    numcutsadded = ( all_data['ID_jabr_cuts'] + all_data['ID_i2_cuts']
                     + all_data['ID_limit_cuts'] )
    numcuts      = ( all_data['num_jabr_cuts'] + all_data['num_i2_cuts']
                     + all_data['num_limit_cuts'] )
    
    summary_ws.write(' case ' + all_data['casename'] + ' opt_status ' 
                     + str(all_data['optstatus']) + ' obj ' 
                     + str(all_data['objval']) + ' runtime ' 
                     + str(all_data['runtime']) + ' iterations ' 
                     + str(all_data['round']) + ' rndcuts '
                     + str((all_data['round']-1)) + ' numcutsadded '
                     + str(numcutsadded) + ' numcuts '
                     + str(numcuts) +  '\n')

    summary_ws.close()

    ############################ GET DUALS #################################

    if (themodel.status != GRB.status.NUMERIC):
      getduals(log,all_data)
    
    ############################ TERMINATION #################################



    if (all_data['round'] >= all_data['max_rounds']):

      if all_data['writesol']:
        writesol(log,all_data)
        writesol_allvars(log,all_data)

      if (themodel.status != GRB.status.NUMERIC):
        print_duals(log,all_data)
        print_duals_table(log,all_data)
      
      ###### TABLE
      objvalue_last = round(all_data['objval'],2)
      time_last     = round(time.time() - all_data['round_time'],2)
      status_last   = all_data['optstatus']
      numcuts_last  = numcuts
      dinfs_last    = round(all_data['dinfs_scaled'],10)
      numvars       = all_data['numvars']
      numconstrs    = all_data['numconstrs']
      ftime         = round(all_data['formulation_time'],2)
      if all_data['addcuts']:
        tablename     = 'table_cp_ws_' + all_data['datetime'] + '.txt'
        result_1st    = str(objvalue_1st) + " & " + str(time_1st) + " & " + str(numcuts_1st) + " & " + str(dinfs_1st)
      else:
        tablename     = 'table_cp_sc_' + all_data['datetime'] + '.txt'
        result_1st    = str(objvalue_1st) + " & " + str(time_1st) + " & " + str(dinfs_1st)        
      table         = open(tablename,"a+")
      result_last   = str(objvalue_last) + " & " + str(time_last) + " & " + str(numcuts_last) + " & " + str(dinfs_last)
      result =  casename + " & " + str(numvars) + " & " + str(numconstrs) + " & " + str(ftime) + " & " + result_1st + " & " + result_last + " & " + str(all_data['runtime']) + " & " + str(all_data['round'] - 1) + " &\\\ \n"
      table.write(result)
      table.close()
      ######

      if all_data['writelastLP']:
        log.joint(' writing down last lp...\n')        
        themodel.write(casename + '_' + str(T) + '_' + casetype + "_last.lp")
        
      summary_ws = open("summary_ws.log","a+") #later add feasibility errors, etc
      summary_ws.write(' rounds limit reached!\n\n')
      summary_ws.close()
      log.joint(' rounds limit reached!\n')
      log.joint(' bye\n')
      return None
          
    if all_data['runtime'] > all_data['max_time']:

      if all_data['writesol']:
        writesol(log,all_data)
        writesol_allvars(log,all_data)


      print_duals(log,all_data)
      print_duals_table(log,all_data)
      
      ###### TABLE
      objvalue_last = round(all_data['objval'],2)
      time_last     = round(time.time() - all_data['round_time'],2)
      status_last   = all_data['optstatus']
      numcuts_last  = numcuts
      dinfs_last    = round(all_data['dinfs_scaled'],10)
      numvars       = all_data['numvars']
      numconstrs    = all_data['numconstrs']
      ftime         = round(all_data['formulation_time'],2)
      if all_data['round'] == 1:
        objvalue_1st = obvalue_last
        time_1st     = time_last
        numcuts_1st  = numcuts_last
        dinfs_1st    = dinfs_last
      if all_data['addcuts']:
        tablename     = 'table_cp_ws_' + all_data['datetime'] + '.txt'
        result_1st    = str(objvalue_1st) + " & " + str(time_1st) + " & " + str(numcuts_1st) + " & " + str(dinfs_1st)
      else:
        tablename     = 'table_cp_sc_' + all_data['datetime'] + '.txt'
        result_1st    = str(objvalue_1st) + " & " + str(time_1st) + " & " + str(dinfs_1st)              
      table         = open(tablename,"a+")
      result_last   = str(objvalue_last) + " & " + str(time_last) + " & " + str(numcuts_last) + " & " + str(dinfs_last)
      result =  casename + " & " + str(numvars) + " & " + str(numconstrs) + " & " + str(ftime) + " & " + result_1st + " & " + result_last + " & " + str(all_data['runtime']) + " & " + str(all_data['round'] - 1) + " &\\\ \n"
      table.write(result)
      table.close()
      
      ######

      if all_data['writelastLP']:
        log.joint(' writing down last lp...\n')        
        themodel.write(casename + '_' + str(T) + '_' + casetype + "_last.lp")      
        
      summary_ws = open("summary_ws.log","a+") #later add feasibility errors, etc
      summary_ws.write(' time limit reached!\n\n')
      summary_ws.close()
      log.joint(' time limit reached!\n')
      log.joint(' bye\n')
      return None

    if (all_data['ftol_counter'] > all_data['ftol_iterates']):

      if all_data['writesol']:
        writesol(log,all_data)
        writesol_allvars(log,all_data)

        
      print_duals(log,all_data)
      print_duals_table(log,all_data)
      
      ###### TABLE
      objvalue_last = round(all_data['objval'],2)
      time_last     = round(time.time() - all_data['round_time'],2)
      status_last   = all_data['optstatus']
      numcuts_last  = numcuts
      dinfs_last    = round(all_data['dinfs_scaled'],10)
      numvars       = all_data['numvars']
      numconstrs    = all_data['numconstrs']
      ftime         = round(all_data['formulation_time'],2)
      if all_data['addcuts']:
        tablename     = 'table_cp_ws_' + all_data['datetime'] + '.txt'
        result_1st    = str(objvalue_1st) + " & " + str(time_1st) + " & " + str(numcuts_1st) + " & " + str(dinfs_1st)
      else:
        tablename     = 'table_cp_sc_' + all_data['datetime'] + '.txt'
        result_1st    = str(objvalue_1st) + " & " + str(time_1st) + " & " + str(dinfs_1st)              
      table         = open(tablename,"a+")
      result_last   = str(objvalue_last) + " & " + str(time_last) + " & " + str(numcuts_last) + " & " + str(dinfs_last)
      result =  casename + " & " + str(numvars) + " & " + str(numconstrs) + " & " + str(ftime) + " & " + result_1st + " & " + result_last + " & " + str(all_data['runtime']) + " & " + str(all_data['round'] - 1) + " &\\\ \n"
      table.write(result)
      table.close()      


      ######


      if all_data['writelastLP']:
        log.joint(' writing down last lp...\n')
        themodel.write(casename + '_' + str(T) + '_' + casetype + "_last.lp")      

        
      summary_ws = open("summary_ws.log","a+") 
      summary_ws.write(' poor consecutive obj improvement limit reached!\n\n')
      summary_ws.close()
      log.joint(' poor consecutive obj improvement limit reached\n')
      log.joint(' bye\n')
      return None

    ##########################################################################

    if all_data['max_rounds'] == 1:

      if all_data['writesol']:
        writesol(log,all_data)      
        writesol_allvars(log,all_data)
        
      log.joint(' bye!\n')
      sys.exit(0)

    ############################### CUTS ######################################

    # Cut computations and management
    cutplane_cuts(log,all_data)

    # Cut statistics
    cutplane_cutstats(log,all_data)
    
    themodel.update()

    log.joint(' model updated\n')
    log.joint('\n')

    ############################### WRITE CUTS ################################

    if all_data['writecuts']:
      write_cuts(log,all_data)

    ############################### WRITE LPS #################################
    
    if all_data['losstest']:
      themodel.write('perturbedjabr.lp')
      constr = themodel.getQConstrs()[all_data['rbranch']]
      constr.setAttr("QCRHS",0)
  
      themodel.update()
 
    ###########

    if all_data['writelps'] and ( all_data['round'] > 0 ):
      name = 'post_cuts' + '_' + str(all_data['round']) + '.lp'
      themodel.write(name)
      log.joint(' model with new cuts written to .lp file\n')


    if all_data['losstest']:
      valloss = 0
      for branch in branches.values():
        loss = all_data['Pfvalues'][branch] + all_data['Ptvalues'][branch]
        log.joint(' loss at branch ' + str(branch.count) + ' = ' 
                  + str(loss) + '\n')
        valloss += loss

      log.joint('ploss ' + str(sum(plossvalues.values())) + '\n')
      log.joint('value loss ' + str(valloss) + '\n')

    ################# losses test

    if all_data['losstest']:
      rbranch = all_data['branches'][all_data['rbranch']]

      rbranchloss = all_data['Pfvalues'][rbranch] + all_data['Ptvalues'][rbranch]
      all_data['losstest_dic'][all_data['round']] = (all_data['rbranch'],
                                                     all_data['objval'],
                                                     sum(plossvalues.values()),
                                                     rbranchloss)

      all_data['losstest_loss'][all_data['round']] = sum(plossvalues.values())
      all_data['losstest_branchloss'][all_data['round']] = rbranchloss

      avg        = sum(all_data['losstest_loss'].values()) / all_data['round']
      avg_branch = sum(all_data['losstest_branchloss'].values()) / all_data['round']
      
      maxloss = max(all_data['losstest_loss'].values())
      maxloss_branch = max(all_data['losstest_branchloss'].values())

      minloss = min(all_data['losstest_loss'].values())
      minloss_branch = min(all_data['losstest_branchloss'].values())

      log.joint(' round ' + str(all_data['round']) + '\n')
      log.joint(' random branch ' + str(all_data['rbranch']) + '\n')
      log.joint(' current losstest avg ' + str(avg) + '\n' )
      log.joint(' current losstest avg (branch)' + str(avg_branch) + '\n' )
      log.joint(' current max loss ' + str(maxloss) + '\n')
      log.joint(' current max loss (at the branch)' + str(maxloss_branch) + '\n')
      log.joint(' current min loss ' + str(minloss) + '\n')
      log.joint(' current min loss (at the branch)' + str(minloss_branch) + '\n')
      #log.joint(' losstest ' + str(all_data['losstest_loss']) + '\n')
      #log.joint(' losstest (branch) ' + str(all_data['losstest_branchloss']) + '\n')
      #breakexit('c')
      if all_data['round'] == 50:
        return None


    ###########################################################################        
    ###### TABLE                                                                                                                  

    if all_data['round'] == 1:
      objvalue_1st = round(all_data['objval'],2)
      time_1st     = round(time.time() - all_data['T0'],2)
      status_1st   = all_data['optstatus']
      numcuts_1st  = numcuts
      dinfs_1st    = round(all_data['dinfs_scaled'],10)
        
    ###########################################################################

    #Kappa retrieval
    #themodel.write('A.lp')
    #expr = LinExpr()
    #themodel.setObjective(expr)
    #themodel.update()
    #if all_data['round'] >= 1:
    #  themodel.printQuality()
    #  log.joint('condition number ' + str(themodel.KappaExact) + '\n')
    #breakexit('c')
    
                                              
    all_data['round']      += 1
    all_data['round_time']  = time.time()


    ###########################################################################

def fixflows(log,all_data):

  if (all_data['matpower_sol'] == 0) and (all_data['ampl_sol'] == 0):
    log.joint(' cannot fix flows since no solution has been loaded!\n')
    return None

  themodel    = all_data['themodel']
  branches    = all_data['branches']
  T           = all_data['T']
  tolerance   = all_data['fix_tolerance']
  Pvar_f      = all_data['Pvar_f']
  Pvar_t      = all_data['Pvar_t']
  Qvar_f      = all_data['Qvar_f']
  Qvar_t      = all_data['Qvar_t']

  sol_Pfvalues = all_data['sol_Pfvalues']
  sol_Ptvalues = all_data['sol_Ptvalues']
  sol_Qfvalues = all_data['sol_Qfvalues']
  sol_Qtvalues = all_data['sol_Qtvalues']

  for branch in branches.values():
    for k in range(T):
      sol_Pf = sol_Pfvalues[k][branch]
      sol_Pt = sol_Ptvalues[k][branch]
      sol_Qf = sol_Qfvalues[k][branch]
      sol_Qt = sol_Qtvalues[k][branch]

      #Pf
      ubound_Pf = sol_Pf + tolerance
      lbound_Pf = sol_Pf - tolerance

      Pvar_f[k][branch].setAttr("ub",ubound_Pf)
      Pvar_f[k][branch].setAttr("lb",lbound_Pf)

      #Pt
      ubound_Pt = sol_Pt + tolerance
      lbound_Pt = sol_Pt - tolerance

      Pvar_t[k][branch].setAttr("ub",ubound_Pt)
      Pvar_t[k][branch].setAttr("lb",lbound_Pt)

      #Qf
      ubound_Qf = sol_Qf + tolerance
      lbound_Qf = sol_Qf - tolerance

      Qvar_f[k][branch].setAttr("ub",ubound_Qf)
      Qvar_f[k][branch].setAttr("lb",lbound_Qf)

      #Qt
      ubound_Qt = sol_Qt + tolerance
      lbound_Qt = sol_Qt - tolerance

      Qvar_t[k][branch].setAttr("ub",ubound_Qt)
      Qvar_t[k][branch].setAttr("lb",lbound_Qt)

  themodel.update()
  themodel.write('fixflows.lp')
  log.joint('check fixflows.lp\n')  


def fixcs(log,all_data):

  if (all_data['matpower_sol'] == 0) and (all_data['ampl_sol'] == 0):
    log.joint(' cannot fix flows since no solution has been loaded!\n')
    return None

  T           = all_data['T']
  themodel    = all_data['themodel']
  tolerance   = all_data['fix_tolerance']
  buses       = all_data['buses']
  branches    = all_data['branches']
  cvar        = all_data['cvar']
  svar        = all_data['svar']
  sol_cvalues = all_data['sol_cvalues']
  sol_svalues = all_data['sol_svalues']

  for bus in buses.values():

    for k in range(T):
      sol_v2 = all_data['sol_cvalues'][k][bus]
      ubound_v2 = sol_v2 + tolerance
      lbound_v2 = sol_v2 - tolerance

      cvar[k][bus].setAttr("ub",ubound_v2)
      cvar[k][bus].setAttr("lb",lbound_v2)

  for branch in branches.values():
    for k in range(T):
    
      sol_c = all_data['sol_cvalues'][k][branch]
      sol_s = all_data['sol_svalues'][k][branch]

      #c
      ubound_c = sol_c + tolerance
      lbound_c = sol_c - tolerance

      cvar[k][branch].setAttr("ub",ubound_c)
      cvar[k][branch].setAttr("lb",lbound_c)

      #s
      ubound_s = sol_s + tolerance
      lbound_s = sol_s - tolerance

      svar[k][branch].setAttr("ub",ubound_s)
      svar[k][branch].setAttr("lb",lbound_s)

  themodel.update()
  themodel.write('fixCS.lp')
  log.joint('check fixCS.lp\n')


def computebalbounds(log, all_data, bus, k):

  #first let's get max/min generations

  loud = 0

  baseMVA = all_data['baseMVA']
  gens = all_data['gens']
  Pd   = all_data['Pd']

  Pubound = Plbound = 0
  Qubound = Qlbound = 0

  for gencounter in bus.genidsbycount:
    if gens[gencounter].status:
      Pubound += gens[gencounter].Pmax
      Plbound += gens[gencounter].Pmin
      Qubound += gens[gencounter].Qmax
      Qlbound += gens[gencounter].Qmin

      if bus.nodetype == 3:
        Qubound = + GRB.INFINITY
        Qlbound = - GRB.INFINITY
        
    if loud:
     #log.joint(" Pubound for " + str(bus.nodeID) + " " + str(Pubound) + " genc " + str(gencounter) + "\n")
     #log.joint(" Plbound for " + str(bus.nodeID) + " " + str(Plbound) + " genc " + str(gencounter) + "\n")
     log.joint(" Qubound for " + str(bus.nodeID) + " " + str(Qubound) + " genc " + str(gencounter) + "\n")
     log.joint(" Qlbound for " + str(bus.nodeID) + " " + str(Qlbound) + " genc " + str(gencounter) + "\n")


  #log.joint(' bus ' + str(bus) + ' busnode ' + str(bus.nodeID) +  ' t ' + str(k) + '\n')
  Pubound -= Pd[k][bus]
  Plbound -= Pd[k][bus]
  
  if all_data['arpae2']:
    Qubound = - all_data['Qd'][k][bus]
    Qlbound = - all_data['Qd'][k][bus]
  else:
    Qubound -= bus.Qd
    Qlbound -= bus.Qd

  if bus.nodetype == 4:
    Pubound = Plbound = Qubound = Qlbound = 0
    Pd[k][bus] = 0
    
  if loud:
     #log.joint(" Pubound for " + str(bus.nodeID) + " final " + str(Pubound) + "\n")
     #log.joint(" (Pd was %g)\n" %bus.Pd)
     #log.joint(" Plbound for " + str(bus.nodeID) + " final " + str(Plbound) + "\n")
     log.joint(" Qubound for " + str(bus.nodeID) + " final " + str(Qubound) + "\n")
     log.joint(" (Qd was %g)\n" %bus.Qd)
     log.joint(" Qlbound for " + str(bus.nodeID) + " final " + str(Qlbound) + "\n")
     breakexit(" ")
  
  return Pubound, Plbound, Qubound, Qlbound


def computeangles(log,all_data):
  
  cvalues  = all_data['cvalues']
  svalues  = all_data['svalues']
  branches = all_data['branches']

  for branch in branches.values():
    s = svalues[branch]
    c = cvalues[branch]

    log.joint(' branch ' + str(branch.count) + ' f ' + str(branch.f) + ' t ' + str(branch.t) + ' angle ' + str(math.atan(s/c)) + '\n')



def getsol_matpower(log,all_data):

  casename = all_data['casename']
  filename = 'mp_sols/solution_va_'+ casename +'.txt'

  log.joint(" reading file with matpower solution volt magnitudes and angles " + filename + "\n")

  try:
      thefile  = open(filename, "r")
      lines    = thefile.readlines()
      lenlines = len(lines)
      thefile.close()
  except:
      log.stateandquit("cannot open file", filename)
      sys.exit("failure")

  sol_vm       = {}
  sol_angle    = {}
  sol_cvalues  = {}
  sol_svalues  = {}

  buses        = all_data['buses']
  branches     = all_data['branches']
  IDtoCountmap = all_data['IDtoCountmap']


  linenum  = 1
  while linenum < lenlines:
    thisline = lines[linenum].split(',')
    buscount = int(thisline[0]) #bug, not bus_id -> matpower uses buscount
    bus      = buses[buscount]

    sol_vm[bus]        = float(thisline[1])
    sol_angle[bus]     = float(thisline[2]) * math.pi / 180
    sol_cvalues[bus]   = sol_vm[bus]**2

    linenum += 1

  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]
    bus_f       = buses[count_of_f]
    bus_t       = buses[count_of_t]


    vm_f        = sol_vm[bus_f]
    vm_t        = sol_vm[bus_t]

    sol_c      = vm_f * vm_t * math.cos(sol_angle[bus_f] - sol_angle[bus_t])
    sol_s      = vm_f * vm_t * math.sin(sol_angle[bus_f] - sol_angle[bus_t])

    sol_cvalues[branch] = sol_c
    sol_svalues[branch] = sol_s

  all_data['sol_vm']      = sol_vm
  all_data['sol_angle']   = sol_angle
  all_data['sol_cvalues'] = sol_cvalues
  all_data['sol_svalues'] = sol_svalues

  log.joint(' done loading volts, angles, and cs values\n')

  #########

  filename = 'mp_sols/solution_'+ casename +'.txt'
  log.joint(" reading file with matpower power flows " + filename + "\n")

  try:
      thefile  = open(filename, "r")
      lines    = thefile.readlines()
      lenlines = len(lines)
      thefile.close()
  except:
      log.stateandquit("cannot open file", filename)
      sys.exit("failure")

  sol_Pfvalues    = {}
  sol_Ptvalues    = {}
  sol_Qfvalues    = {}
  sol_Qtvalues    = {}

  linenum = 1

  while linenum < lenlines:
    thisline            = lines[linenum].split(',')
    branchcount         = int(thisline[0])
    branch              = branches[branchcount]
    f                   = branch.f
    t                   = branch.t

    if int(thisline[1]) != f or int(thisline[2]) != t:
      breakexit('check')

    sol_Pfvalues[branch]      = float(thisline[3])/all_data['baseMVA']
    sol_Ptvalues[branch]      = float(thisline[4])/all_data['baseMVA']
    sol_Qfvalues[branch]      = float(thisline[5])/all_data['baseMVA']
    sol_Qtvalues[branch]      = float(thisline[6])/all_data['baseMVA']

    linenum   += 1
  
  all_data['sol_Pfvalues'] = sol_Pfvalues
  all_data['sol_Ptvalues'] = sol_Ptvalues
  all_data['sol_Qfvalues'] = sol_Qfvalues
  all_data['sol_Qtvalues'] = sol_Qtvalues

  log.joint(' done loading power flows\n')
  log.joint(' matpower solution loaded\n')

def writeACsol(log,all_data):

  if (all_data['matpower_sol'] == 0) and (all_data['ampl_sol'] == 0):
    log.joint(' cannot write AC solution since no solution has been loaded!\n')
    return None

  branches     = all_data['branches']
  buses        = all_data['buses']
  IDtoCountmap = all_data['IDtoCountmap']
  sol_vm       = all_data['sol_vm']
  sol_angle    = all_data['sol_angle']
  sol_cvalues  = all_data['sol_cvalues']
  sol_svalues  = all_data['sol_svalues']
  sol_Pfvalues = all_data['sol_Pfvalues']
  sol_Ptvalues = all_data['sol_Ptvalues']
  sol_Qfvalues = all_data['sol_Qfvalues']
  sol_Qtvalues = all_data['sol_Qtvalues']
  tolerance    = all_data['tol_fix']

  filename   = 'ACsol_' + all_data['casename'] + '.lp'
  thefile    = open(filename,"w+")

  log.joint(' writing voltages to ' + filename + ' ...\n')

  for bus in buses.values():
    f        = bus.nodeID
    sol_v2   = sol_cvalues[bus]
    linelp_v = str(sol_v2 - tolerance) + ' <= c_' + str(f) + '_' + str(f) + ' <= ' + str(sol_v2 + tolerance) + '\n'
    thefile.write(linelp_v)

  log.joint(' writing cs values to ' + filename + ' ...\n')

  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]
    bus_f       = buses[count_of_f]
    bus_t       = buses[count_of_t]

    vm_f        = sol_vm[bus_f]
    vm_t        = sol_vm[bus_t]

    sol_c      = vm_f * vm_t * math.cos(sol_angle[bus_f] - sol_angle[bus_t])
    sol_s      = vm_f * vm_t * math.sin(sol_angle[bus_f] - sol_angle[bus_t])

    linelp_c = str(sol_c - tolerance) + ' <= c_' + str(f) + '_' + str(t) + ' <= ' + str(sol_c + tolerance) + '\n'
    thefile.write(linelp_c)

    linelp_s = str(sol_s - tolerance) + ' <= s_' + str(f) + '_' + str(t) + ' <= ' + str(sol_s + tolerance) + '\n'
    thefile.write(linelp_s)

  log.joint(' writing flows to ' + filename + ' ...\n')

  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    Pfval       = sol_Pfvalues[branch]
    Ptval       = sol_Ptvalues[branch]
    Qfval       = sol_Qfvalues[branch]
    Qtval       = sol_Qtvalues[branch]
  
    linelp_Pf = str(Pfval-tolerance) + ' <= P_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' <= ' + str(Pfval+tolerance) + '\n'
    thefile.write(linelp_Pf)
    
    linelp_Pt = str(Ptval-tolerance) + ' <= P_' + str(branchcount) + '_' + str(t) + '_' + str(f) + ' <= ' + str(Ptval+tolerance) + '\n'
    thefile.write(linelp_Pt)
    
    linelp_Qf = str(Qfval-tolerance) + ' <= Q_' + str(branchcount) + '_' + str(f) + '_' + str(t) + ' <= ' + str(Qfval+tolerance) + '\n'
    thefile.write(linelp_Qf)

    linelp_Qt = str(Qtval-tolerance) + ' <= Q_' + str(branchcount) + '_' + str(t) + '_' + str(f) + ' <= ' + str(Qtval+tolerance) + '\n'
    thefile.write(linelp_Qt)

  thefile.close()
  log.joint(' done writing AC solution to .lp file\n')
  

def computei2value(log,all_data,branch,mp_c,mp_s,mp_cbusf,mp_cbust):

  ratio  = branch.ratio
  y      = branch.y
  g      = y.real
  b      = y.imag
  bshunt = branch.bc
  angle  = branch.angle_rad

  #first f                                                                                                                                                                                           
  i2f  = 0
  i2f += (g*g + b*b)/(ratio*ratio) * ( (mp_cbusf/(ratio*ratio)) + mp_cbust - (2/ratio) * ( mp_c * math.cos(angle) + mp_s * math.sin(angle) ) )
  i2f += b*bshunt/(ratio**3) * ( (mp_cbusf/ratio) - (mp_c * math.cos(angle) + mp_s * math.sin(angle) ))
  i2f += g*bshunt/(ratio**3) * ( mp_s * math.cos(angle) - mp_c * math.sin(angle) )
  i2f += (bshunt*bshunt*mp_cbusf/(4*(ratio**4)) )

  return i2f
  
def loss_inequalities(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  Pvar_f         = all_data['Pvar_f']
  Pvar_t         = all_data['Pvar_t']
  FeasibilityTol = all_data['FeasibilityTol']

  if ( all_data['matpower_sol'] or all_data['ampl_sol'] ) and all_data['loss_validity']:
    FeasibilityTol = all_data['FeasibilityTol']
    log.joint('  adding and checking validity of loss inequalities wrt a solution\n')
    log.joint('   feasibility tolerance ' + str(FeasibilityTol) + '\n')
  else:
    log.joint('  loss inequalities\n')

  all_data['loss_cuts'] = {}
  counter_loss = 0

  for branch in branches.values():
    if branch.r < 0: ##!!!
      continue

    branchcount = branch.count
    f           = branch.f
    t           = branch.t

    if all_data['loss_validity']:
     sol_Pf     = all_data['sol_Pfvalues'][branch]
     sol_Pt     = all_data['sol_Ptvalues'][branch]
     violation  = - (sol_Pf + sol_Pt)

     if violation > FeasibilityTol:
       log.joint('   WARNING, the loss inequality associated to branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' is violated by the AC solution!\n')
       log.joint('   violation ' + str(violation) + '\n')
       log.joint('   values (AC solution) ' + ' Pf ' + str(sol_Pf) + ' Pt ' + str(sol_Pt) + '\n')
       breakexit('check!')
     else:
       log.joint('   AC solution satisfies loss inequality at branch ' + str(branchcount) + ' with slack ' + str(violation) + '\n')

    counter_loss                 += 1
    all_data['loss_cuts'][branch] = (0,0,FeasibilityTol)

    lossexp    = LinExpr()
    constrname = "loss_ineq_"+str(branch.count)+"_"+str(f)+"_"+str(t)
    lossexp   += Pvar_f[branch] + Pvar_t[branch]
    themodel.addConstr(lossexp >= 0, name = constrname)

  all_data['num_loss_cuts']         = counter_loss
  all_data['num_loss_cuts_dropped'] = 0
  all_data['dropped_loss']          = []

  log.joint('   %d loss inequalities added\n'%counter_loss) 
  
  return counter_loss
  
def qloss_inequalities(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  Pvar_f         = all_data['Pvar_f']
  Pvar_t         = all_data['Pvar_t']
  Qvar_f         = all_data['Qvar_f']
  Qvar_t         = all_data['Qvar_t']
  FeasibilityTol = all_data['FeasibilityTol']

  if ( all_data['matpower_sol'] or all_data['ampl_sol'] ) and all_data['qloss_validity']:
    FeasibilityTol = all_data['FeasibilityTol']
    log.joint('  adding and checking validity of loss inequalities wrt a solution\n')
    log.joint('   feasibility tolerance ' + str(FeasibilityTol) + '\n')
  else:
    log.joint('  qloss inequalities\n')

  all_data['qloss_cuts'] = {}
  qcounter_loss = 0

  for branch in branches.values():
    if branch.r <= 0 or (branch.bc != 0): ##!!!
      continue

    branchcount = branch.count
    f           = branch.f
    t           = branch.t

    if all_data['qloss_validity']:
     sol_Pf     = all_data['sol_Qfvalues'][branch]
     sol_Pt     = all_data['sol_Qtvalues'][branch]
     violation  = - (sol_Qf + sol_Qt)

     if violation > FeasibilityTol:
       log.joint('   WARNING, the qloss inequality associated to branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' is violated by the AC solution!\n')
       log.joint('   violation ' + str(violation) + '\n')
       log.joint('   values (AC solution) ' + ' Qf ' + str(sol_Qf) + ' Pt ' + str(sol_Qt) + '\n')
       breakexit('check!')
     else:
       log.joint('   AC solution satisfies loss inequality at branch ' + str(branchcount) + ' with slack ' + str(violation) + '\n')

    qcounter_loss                 += 1
    all_data['qloss_cuts'][branch] = (0,0,FeasibilityTol)

    qlossexp    = LinExpr()
    constrname = "qloss_ineq_"+str(branch.count)+"_"+str(f)+"_"+str(t)
    qlossexp   += (Qvar_f[branch] + Qvar_t[branch]) - (branch.x/branch.r)*(Pvar_f[branch] + Pvar_t[branch])
    themodel.addConstr(qlossexp == 0, name = constrname)

  all_data['num_qloss_cuts']         = qcounter_loss
  all_data['num_qloss_cuts_dropped'] = 0
  all_data['dropped_qloss']          = []

  log.joint('   %d qloss inequalities added\n'%qcounter_loss) 
  
  return qcounter_loss

def jabr_inequalities(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  buses          = all_data['buses']
  cvar           = all_data['cvar']
  svar           = all_data['svar']
  IDtoCountmap   = all_data['IDtoCountmap']
  FeasibilityTol = all_data['FeasibilityTol']
  T              = all_data['T']
  
  if ( all_data['matpower_sol'] or all_data['ampl_sol'] ) and all_data['jabr_validity']:
    log.joint('  adding and checking validity of Jabr inequalities wrt a solution\n')
    log.joint('   feasibility tolerance ' + str(FeasibilityTol) + '\n')
  else:
    log.joint('  Jabr inequalities\n')

  maxviolation = 0
  violated     = 0
  maxbranch    = -1
  maxbusf      = -1
  maxbust      = -1
  counter_jabr = 0

  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]

    if ( all_data['matpower_sol'] or all_data['ampl_sol'] ) and all_data['jabr_validity']:
      sol_c         = all_data['sol_cvalues'][branch]
      sol_s         = all_data['sol_svalues'][branch]
      sol_cbusf     = all_data['sol_cvalues'][buses[count_of_f]]
      sol_cbust     = all_data['sol_cvalues'][buses[count_of_t]]
      relviolation = violation = sol_c * sol_c + sol_s * sol_s - sol_cbusf * sol_cbust
      #relviolation = violation / ( ( coeff_Pft**2 + coeff_Qft**2 + coeff_cff**2 + coeff_i2ft**2 )**0.5 )
      if relviolation > maxviolation:
        maxviolation = relviolation
        maxbranch    = branch.count
        maxbusf      = f
        maxbust      = t

      if relviolation > FeasibilityTol:
        violated += 1
        log.joint('   WARNING, the Jabr-inequality associated to branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' is violated by the AC solution!\n')
        log.joint('   violation ' + str(violation) + '\n')
        log.joint('   relative violation ' + str(relviolation) + '\n')
        log.joint('   values (AC solution) ' + ' cft ' + str(sol_c) + ' sft ' + str(sol_s) + ' cff ' + str(sol_cbusf) + ' ctt ' + str(sol_cbust) + '\n' )
        #breakexit('check!')
      else:
        log.joint('   AC solution satisfies loss inequality at branch ' + str(branchcount) + ' with slack ' + str(relviolation) + '\n')

    for k in range(T):
      
      trigexp       = QuadExpr()
      constrname    = "jabr_"+str(branchcount)+"_"+str(f)+"_"+str(t)+"_"+str(k)
      trigexp      += cvar[k][branch]*cvar[k][branch] + svar[k][branch]*svar[k][branch] - cvar[k][buses[count_of_f]]*cvar[k][buses[count_of_t]]

      themodel.addConstr(trigexp <= 0, name = constrname)
      counter_jabr += 1
    
  if ( all_data['matpower_sol'] or all_data['ampl_sol'] ) and all_data['jabr_validity']:
    log.joint('  max violation of Jabr-inequalities by AC solution ' + str(maxviolation) + ' at branch ' + str(maxbranch) + ' f ' + str(maxbusf) + ' t ' + str(maxbust) + '\n')
    log.joint('  number of violated Jabr-inequalities ' + str(violated) + '\n')
    breakexit('  check Jabr violation')

  log.joint('   %d Jabr inequalities added\n'%counter_jabr) 

  return counter_jabr



def i2_def(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  buses          = all_data['buses']
  i2var_f        = all_data['i2var_f']
  cvar           = all_data['cvar']
  svar           = all_data['svar']
  IDtoCountmap   = all_data['IDtoCountmap']
  T              = all_data['T']
  alphadic       = all_data['alphadic']
  
  counter_i2def  = 0
  counter_i2con  = 0

  ######
  # alpha_max = -1e20
  # alpha_avg = 0
  # beta_max  = -1e20
  # beta_avg  = 0
  # gamma_max = -1e20
  # gamma_avg = 0
  # zeta_max  = -1e20
  # zeta_avg  = 0  

  # newbeta_max  = -1e20
  # newbeta_avg  = 0  
  # newgamma_max = -1e20
  # newgamma_avg = 0
  # newzeta_max  = -1e20
  # newzeta_avg  = 0
  # RHS_min      = 1e20
  # RHS_max      = -1e20
  # RHS_avg      = 0
  #######

  
  log.joint('  i2 variables definition and i2 linear inequalities\n')

  for branch in branches.values():
    expr_f      = LinExpr()
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]
    ratio       = branch.ratio
    y           = branch.y
    g           = y.real
    b           = y.imag
    bshunt      = branch.bc
    angle       = branch.angle_rad
    bus_f       = buses[count_of_f]
    bus_t       = buses[count_of_t]
    
    alpha = alphadic[branch]
    beta  = ( g*g + b*b ) / (ratio**2)
    gamma = ( math.cos(angle) * ( - 2 * (g*g + b*b) - b * bshunt ) + math.sin(angle) * ( - g * bshunt) ) / (ratio**3)
    zeta  = ( math.sin(angle) * ( - 2 * (g*g + b*b) - b * bshunt ) - math.cos(angle) * ( - g* bshunt) ) / (ratio**3)

    #log.joint(' branch ' + str(branch) + ' id ' + str(branch.count)
    #          + ' alpha ' + str(alpha) + '\n')
    
    if alpha < all_data['alpha_threshold']:
      for k in range(T):
        constrname_f = 'i2def_'+str(branch.count)+"_"+str(f) + "_" + str(t) + "_" + str(k)
        expr_f = alpha * cvar[k][bus_f] + beta * cvar[k][bus_t] + gamma * cvar[k][branch] + zeta * svar[k][branch]

        themodel.addConstr(expr_f == i2var_f[k][branch],name = constrname_f) 
        counter_i2def += 1
    else:
      upperbound_f = branch.limit**2 / (bus_f.Vmin * bus_f.Vmin)
      newalpha     = 1
      newbeta      = beta/alpha
      newgamma     = gamma/alpha
      newzeta      = zeta/alpha
      RHS          = upperbound_f / alpha

      ########
      # alpha_avg += alpha 
      # beta_avg  += beta
      # gamma_avg += gamma
      # zeta_avg  += zeta

      # if alpha >= alpha_max:
      #   alpha_max = alpha      
      # if beta  >= beta_max:
      #   beta_max  = beta
      # if newgamma >= gamma_max:
      #   gamma_max = gamma
      # if newzeta  >= zeta_max:
      #   zeta_max  = zeta
      
      # newbeta_avg  += newbeta
      # newgamma_avg += newgamma
      # newzeta_avg  += newzeta
      # RHS_avg      += RHS
      
      # if newbeta  >= newbeta_max:
      #   newbeta_max  = newbeta
      # if newgamma >= newgamma_max:
      #   newgamma_max = newgamma
      # if newzeta  >= newzeta_max:
      #   newzeta_max  = newzeta
      # if RHS >= RHS_max:
      #   RHS_max = RHS
      # if RHS <= RHS_min:
      #   RHS_min = RHS
      ########

      
      for k in range(T):
        uppconstrname_f = 'uppi2_'+str(branchcount)+"_"+str(f) + "_" + str(t) + "_" + str(k)
        lowconstrname_f = 'lowi2_'+str(branchcount)+"_"+str(f) + "_" + str(t) + "_" + str(k)

        expr_f = cvar[k][bus_f] + newbeta * cvar[k][bus_t] + newgamma * cvar[k][branch] + newzeta * svar[k][branch]
        themodel.addConstr(expr_f <= RHS, name = uppconstrname_f)
        themodel.addConstr(expr_f >= 0,   name = lowconstrname_f)
        counter_i2con += 2
     
  log.joint('   %d i2 definition constraints added\n'%counter_i2def) 
  log.joint('   %d i2 linear constraints added\n'%counter_i2con)

  ##### check coeffs
  # numbadi2s    = (1/(2 * T)) * counter_i2con

  # alpha_avg  = alpha_avg / numbadi2s
  # beta_avg  = beta_avg / numbadi2s
  # gamma_avg = gamma_avg / numbadi2s
  # zeta_avg  = zeta_avg / numbadi2s
  
  # newbeta_avg  = newbeta_avg / numbadi2s
  # newgamma_avg = newgamma_avg / numbadi2s
  # newzeta_avg  = newzeta_avg / numbadi2s
  # RHS_avg      = RHS_avg / numbadi2s

  # log.joint('  alpha_max/avg ' + str(alpha_max) + '/' + str(alpha_avg)
  #           + ' beta_max/avg ' + str(beta_max) + '/' + str(beta_avg)
  #           + ' gamma_max/avg ' + str(gamma_max) + '/' + str(gamma_avg)
  #           + ' zeta_max/avg ' + str(zeta_max) + '/' + str(zeta_avg)            
  #           + '\n')
  
  # log.joint('  newbeta_max/avg ' + str(newbeta_max) + '/' + str(newbeta_avg)
  #           + ' newgamma_max/avg ' + str(newgamma_max) + '/' + str(newgamma_avg)
  #           + ' newzeta_max/avg ' + str(newzeta_max) + '/' + str(newzeta_avg)            
  #           + '\n')
  # log.joint('  RHS_max/min/avg ' + str(RHS_max) + '/' + str(RHS_min) + '/' + str(RHS_avg) + '\n')  
  
  # breakexit('check coeffs')
  ######
  
  return counter_i2def + counter_i2con


def i2_inequalities(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  buses          = all_data['buses']
  cvar           = all_data['cvar']
  i2var_f        = all_data['i2var_f']
  Pvar_f         = all_data['Pvar_f']
  Qvar_f         = all_data['Qvar_f']
  IDtoCountmap   = all_data['IDtoCountmap']
  FeasibilityTol = all_data['FeasibilityTol']

  if ( all_data['matpower_sol'] or all_data['ampl_sol'] ) and all_data['i2_validity']:
    log.joint('  adding and checking validity of i2 inequalities wrt a solution\n')
    log.joint('   feasibility tolerance ' + str(FeasibilityTol) + '\n')
  else:
    log.joint('  i2 inequalities\n')

  maxviolation = 0
  violated     = 0
  maxbranch    = -1
  maxbusf      = -1
  maxbust      = -1
  maxi2f       = -1
  maxcff       = -1
  maxPf        = -1
  maxQf        = -1

  counter_i2   = 0

  for branch in branches.values():
    branchcount = branch.count
    f           = branch.f
    t           = branch.t
    count_of_f  = IDtoCountmap[f]
    count_of_t  = IDtoCountmap[t]

    if ( all_data['matpower_sol'] or all_data['ampl_sol'] ) and all_data['i2_validity']:
      sol_Pf        = all_data['sol_Pfvalues'][branch]
      sol_Qf        = all_data['sol_Qfvalues'][branch]
      sol_c         = all_data['sol_cvalues'][branch]
      sol_s         = all_data['sol_svalues'][branch]
      sol_cbusf     = all_data['sol_cvalues'][buses[count_of_f]]
      sol_cbust     = all_data['sol_cvalues'][buses[count_of_t]]
      sol_i2f       = computei2value(log,all_data,branch,sol_c,sol_s,sol_cbusf,sol_cbust)
      relviolation = violation = sol_Pf * sol_Pf + sol_Qf * sol_Qf - sol_cbusf * sol_i2f
      #relviolation = violation / ( ( coeff_Pft**2 + coeff_Qft**2 + coeff_cff**2 + coeff_i2ft**2 )**0.5 )
      if relviolation > maxviolation:
        maxviolation = relviolation
        maxbranch    = branch.count
        maxbusf      = f
        maxbust      = t
        maxi2f       = sol_i2f
        maxcff       = sol_cbusf
        maxPf        = sol_Pf
        maxQf        = sol_Qf

      if relviolation > FeasibilityTol:
        violated += 1
        log.joint('   WARNING, the i2 inequality associated to branch ' + str(branchcount) + ' f ' + str(f) + ' t ' + str(t) + ' is violated by the AC solution!\n')
        log.joint('   violation ' + str(violation) + '\n')
        log.joint('   relative violation ' + str(relviolation) + '\n')
        log.joint('   values (AC solution) ' + ' Pft ' + str(sol_Pf) + ' Qft ' + str(sol_Qf) + ' cff ' + str(sol_cbusf) + ' i2ft ' + str(sol_i2f) + '\n' )
        #breakexit('check!')
      else:
        log.joint('   AC solution satisfies i2 inequality at branch ' + str(branchcount) + ' with slack ' + str(relviolation) + '\n')

    counter_i2 += 1
    trigexp     = QuadExpr()
    constrname  = "i2_"+str(branchcount)+"_"+str(f)+"_"+str(t)
    trigexp    += Pvar_f[branch]**2 + Qvar_f[branch]**2 - cvar[buses[count_of_f]] * i2var_f[branch]
    themodel.addConstr(trigexp <= 0, name = constrname)

  if ( all_data['matpower_sol'] or all_data['ampl_sol'] ) and all_data['i2_validity']:
    log.joint('  max violation of i2 inequalities by AC solution ' + str(maxviolation) + ' at branch ' + str(maxbranch) + ' f ' + str(maxbusf) + ' t ' + str(maxbust) + '\n')
    log.joint('  values (AC solution) ' + ' Pft ' + str(maxPf) + ' Qft ' + str(maxQf) + ' cff ' + str(maxcff) + ' i2ft ' + str(maxi2f) + '\n' )
    log.joint('  number of violated i2 inequalities ' + str(violated) + '\n')
    breakexit('  check i2 violation')

  log.joint('   %d i2 inequalities added\n'%counter_i2) 

  return counter_i2  


def limit_inequalities(log,all_data):

  themodel       = all_data['themodel']
  branches       = all_data['branches']
  Pvar_f         = all_data['Pvar_f']
  Pvar_t         = all_data['Pvar_t']
  Qvar_f         = all_data['Qvar_f']
  Qvar_t         = all_data['Qvar_t']
  T              = all_data['T']
  
  counter_limit  = 0

  log.joint('  limit inequalities\n')

  for branch in branches.values():
    #if branch.constrainedflow:
    if True:
      branchcount = branch.count
      f           = branch.f
      t           = branch.t

      for k in range(T):
        constrname = "limit_f_"+str(branchcount)+"_"+str(f)+"_"+str(t)+"_"+str(k)
        limexp = QuadExpr()
        limexp += Pvar_f[k][branch]*Pvar_f[k][branch] + Qvar_f[k][branch]*Qvar_f[k][branch]
        themodel.addConstr(limexp <= branch.limit**2, name = constrname)

        constrname = "limit_t_"+str(branchcount)+"_"+str(t)+"_"+str(f)+"_"+str(k)
        limexp = QuadExpr()
        limexp += Pvar_t[k][branch]*Pvar_t[k][branch] + Qvar_t[k][branch]*Qvar_t[k][branch]
        themodel.addConstr(limexp <= branch.limit**2, name = constrname)

        counter_limit += 2

  log.joint('   %d limit inequalities added\n'%counter_limit) 

  return counter_limit


def hybrid(log,all_data):

  themodel   = all_data['themodel']
  gens       = all_data['gens']
  GenPvar    = all_data['GenPvar']
  objvar     = all_data['objvar']
  lincostvar = all_data['lincostvar']
  qcostvar   = all_data['qcostvar']
  qcost      = all_data['qcost']
  sumTconstr = all_data['sumTconstr']
  lincost    = all_data['lincost']

  QP                  = 50
  QP_to_LP            = 3
  no_objective_cuts   = 11

  log.joint(' running iteration ' + str(all_data['round']) + ' of hybrid algorithm\n')

  #QP
  if all_data['round'] < QP_to_LP or all_data['round'] >= QP:
    if all_data['linear_objective']:  #switching from linear to convex QP
      themodel.params.method       = 2
      themodel.Params.BarHomogeneous = 1
      themodel.Params.BarConvTol = 1e-06
      all_data['linear_objective'] = 0
      all_data['objective_cuts']   = 0
      objvar.setAttr("Obj",0)
      lincostvar.setAttr("Obj",1)
      qcostvar.setAttr("Obj",1)
      themodel.remove(lincost)
      themodel.remove(sumTconstr)
      all_data['qcost'] = qcost = themodel.addConstr(qcostexpr <= qcostvar, name = "qcost")

  #LP
  if QP_to_LP <= all_data['round'] and all_data['round'] < QP:
    if all_data['linear_objective'] == 0: #switching from convex QP to linear
      all_data['linear_objective'] = 1
      all_data['objective_cuts']   = 1
      themodel.params.method       = 1
      objvar.setAttr("Obj",1)
      lincostvar.setAttr("Obj",0)
      qcostvar.setAttr("Obj",0)
      themodel.remove(qcost)
      all_data['sumTconstr'] = sumTconstr = themodel.addConstr(sumTvars <= objvar, name= 'objvar_quad')
      all_data['lincost']    = lincost = themodel.addConstr(lincostvar <= objvar, name= 'objvar_linear')

    if all_data['round'] > no_objective_cuts:
      all_data['objective_cuts'] = 0


  themodel.update()
  hybridlpname = 'hybrid_' + all_data['casename'] + "_" + str(all_data['round']) + '.lp'
  log.joint(' writing down .lp file (hybrid)...\n')
  themodel.write(hybridlpname)

  #breakexit('check.lp')

def cutplane_stats(log,all_data):

  themodel = all_data['themodel']
  
  log.joint('\n ******************** round statistics **********************\n') 
    
  log.joint(' casename = %s\n' % all_data['casename'] )
  log.joint(' round = %g\n' % all_data['round'] )
  log.joint(' objective = %g\n' % all_data['objval'])
  log.joint(' solver status = ' + str(all_data['optstatus']) 
            + ' solver method = ' + str(themodel.params.method) + '\n')
  if themodel.params.method == 2:
    log.joint(' crossover ' + str(themodel.params.crossover) + '\n')
  log.joint(' BarConvTol = ' + str(themodel.params.BarConvTol) 
            + ' FeasTol = ' + str(themodel.params.FeasibilityTol) 
            + ' OptTol = ' + str(themodel.params.OptimalityTol) + '\n') 

  if all_data['optstatus'] == 2 or all_data['optstatus'] == 13:
    log.joint(' max (unscaled/scaled) dual constraint error =\n'
              + '  ' + str(themodel.DualResidual) + ' / '
              + str(themodel.DualSResidual) + '\n')
    #log.joint(' Dual constraint with maximum error (unscaled/scaled)'
    #          + str(themodel.getConstrs()[themodel.DualResidualIndex].constrname) + ' / '
    #          + str(themodel.getConstrs()[themodel.DualSResidualIndex].constrname) + '\n')
  #log.joint(' -- active power --\n')
  #log.joint(' total active power generation = %g\n' % all_data['total_active_gen'] )
  #log.joint(' active power losses = %g\n' % all_data['total_active_losses'] )

  #log.joint(' -- reactive power --\n')
  #log.joint(' total reactive power generation = %g\n' % all_data['total_reactive_gen'] )
  #log.joint(' reactive power net losses = %g\n' % all_data['total_reactive_losses'] )
  #log.joint(' reactive power gains = %g\n' % - all_data['total_reactive_gains'] )

  if all_data['addcuts'] and all_data['round'] == 1:
    log.joint(' -- precomputed cuts --\n')
    log.joint(' Jabr-envelope cuts = %d\n' 
               % all_data['addcuts_numjabrcuts'])
    log.joint(' i2-envelope cuts = %d\n' 
               % all_data['addcuts_numi2cuts'])
    log.joint(' Limit-envelope cuts = %d\n' 
               % all_data['addcuts_numlimitcuts'])
  else:
    log.joint(' -- cut parameters --\n')
    log.joint(' cut age limit = ' 
              + str(all_data['cut_age_limit']) + '\n')
    log.joint(' parallel-cuts threshold = ' 
              + str(all_data['threshold_dotprod']) + '\n')
    log.joint(' initial threshold = ' 
              + str(all_data['initial_threshold']) + '\n') #check what to do with this...
    log.joint(' threshold = ' 
              + str(all_data['tolerance']) + '\n') #check what to do with this...

    log.joint(' total number of cuts = ' + str(all_data['num_jabr_cuts']
                                               + all_data['num_i2_cuts']
                                               + all_data['num_limit_cuts'])
              + '\n')
    if all_data['jabrcuts']:
      log.joint(' -- Jabr-envelope cuts --\n')
      log.joint(' cuts = %d\n' 
                % all_data['num_jabr_cuts'])
      log.joint(' top percent of most violated cuts added = %g\n' 
                % (100*all_data['most_violated_fraction_jabr']) )
      log.joint(' max error (abs) = %g\n' 
                % all_data['max_error_jabr'])
      log.joint(' cut threshold = %g\n' 
                % all_data['threshold']) #check this, one threshold for all?

    if all_data['i2cuts']:
      log.joint(' -- i2-envelope cuts --\n')
      log.joint(' cuts = %d\n' 
                % all_data['num_i2_cuts'])
      log.joint(' top percent of most violated cuts added = %g\n' 
                % (100*all_data['most_violated_fraction_i2']))
      log.joint(' max error (abs) = %g\n'
                % all_data['max_error_i2'])
      log.joint(' cut threshold = %g\n' 
                % all_data['threshold_i2']) #check this as well

    if all_data['limitcuts']:
      log.joint(' -- Limit-envelope cuts --\n')
      log.joint(' cuts = %d\n' 
                % all_data['num_limit_cuts'])
      log.joint(' top percent of most violated added = %g\n' 
                % (100*all_data['most_violated_fraction_limit']))
      log.joint(' max error (abs) = %g\n' 
                % all_data['max_error_limit'])
      log.joint(' cut threshold = %g\n' 
                % all_data['threshold_limit'])

    if all_data['loss_inequalities']:
      log.joint(' -- loss inequalities --\n')
      log.joint(' inequalities = %d\n'
                % all_data['num_loss_cuts'])
      log.joint(' top percent of most violated cuts added = %g\n'
                % (100*all_data['most_violated_fraction_loss']))
      log.joint(' max error (abs) = %g\n'
                % all_data['max_error_loss'])
      

  if all_data['linear_objective'] or all_data['hybrid']:
    log.joint(' -- objective cuts --\n')
    log.joint(' objective-cuts = %d\n'
              % all_data['num_objective_cuts'])
    log.joint(' threshold = %g\n'
              % all_data['threshold_objcuts'])
    
  log.joint(' -- runtimes --\n')
  log.joint(' solver runtime (current round) = %g\n'
            % all_data['solvertime'])
  log.joint(' cumulative solver time = %g\n' 
            % all_data['cumulative_solver_time'])
  
  timenow = time.time()

  log.joint(' time so far (overall - formulation time) = %g\n'
            % (timenow - all_data['T0'] - all_data['formulation_time']))
  log.joint(' time so far (overall) = %g\n'
            % (timenow - all_data['T0']))
  log.joint(' max runtime = %g\n'
            % all_data['max_time'])
  if all_data['ftol_counter']:
    log.joint(' minimum obj improvement (ftol) = %g\n'
              % all_data['ftol'])
    log.joint(' consecutive rounds with poor obj improvement = %d\n' 
              % all_data['ftol_counter'])

  log.joint(' ************************************************************\n') 



def cutplane_cutstats(log,all_data):

  themodel = all_data['themodel']

  log.joint('\n *********************** cuts statistics ********************\n') 
  if all_data['jabrcuts']:
    log.joint(' -- Jabr-envelope cuts --\n')
    log.joint(' cuts = %d\n'
              % all_data['num_jabr_cuts'])
    if all_data['addcuts']:
      log.joint('  precomputed cuts = %d\n' 
              % all_data['addcuts_numjabrcuts'])
    log.joint('  added in current round = %d\n' 
              % all_data['num_jabr_cuts_added'])
    if all_data['dropjabr']:
      log.joint('  dropped in current round = %d\n' 
                % all_data['num_jabr_cuts_dropped'])
    log.joint(' added (overall) = %d\n'
              % all_data['ID_jabr_cuts'])
    if all_data['dropjabr']:
      log.joint(' dropped (overall) = %d\n'
                % all_data['total_jabr_dropped'])

  if all_data['i2cuts']:
    log.joint(' -- i2-envelope cuts --\n')
    log.joint(' cuts = %d\n'
              % all_data['num_i2_cuts'])
    if all_data['addcuts']:
      log.joint('  precomputed cuts = %d\n' 
                % all_data['addcuts_numi2cuts'])
    log.joint('  added in current round = %d\n' 
              % all_data['num_i2_cuts_added'])
    if all_data['dropi2']:
      log.joint('  dropped in current round = %d\n'
                % all_data['num_i2_cuts_dropped'] )
  
    log.joint(' added (overall) = %d\n'
              % all_data['ID_i2_cuts'])

    if all_data['dropi2']:
      log.joint(' dropped (overall) = %d\n'
                % all_data['total_i2_dropped'])

  if all_data['limitcuts']:
    log.joint(' -- Limit-envelope cuts --\n')
    log.joint(' cuts = %d\n'
              % all_data['num_limit_cuts'])
    if all_data['addcuts']:
      log.joint('  precomputed cuts = %d\n' 
                % all_data['addcuts_numlimitcuts'])
    log.joint('  added in current round = %d\n'
              % all_data['num_limit_cuts_added'])
    if all_data['droplimit']:
      log.joint('  dropped in current round = %d\n'
                % all_data['num_limit_cuts_dropped'] )
    log.joint(' added (overall) = %d\n'
               % all_data['ID_limit_cuts'])
    if all_data['droplimit']:
      log.joint(' dropped (overall) = %d\n'
                % all_data['total_limit_dropped'])

  if all_data['loss_inequalities']:
    log.joint(' -- loss inequalities --\n')
    log.joint(' Loss-inequalities = %d\n'
              % all_data['num_loss_cuts'])
    log.joint('  added in current round = %d\n'
              % all_data['num_loss_cuts_added'])
    if all_data['droploss']:
      log.joint('  dropped in current round = %d\n'
                % all_data['num_loss_cuts_dropped'])

  if all_data['linear_objective'] or all_data['hybrid']:
    log.joint(' -- objective cuts --\n')
    log.joint(' objective-cuts = %g\n' % all_data['num_objective_cuts'])
    log.joint(' objective-cuts threshold = %g\n' % all_data['threshold_objcuts'])    

  log.joint(' ---\n')
  log.joint(' total number of cuts = ' + str(all_data['num_jabr_cuts']
                                             + all_data['num_i2_cuts']
                                             + all_data['num_limit_cuts'])
            + '\n')
  log.joint(' ************************************************************\n\n') 


def cutplane_cuts(log,all_data):

  log.joint('\n starting cut procedure ...\n')

  t0_cuts = time.time()

  t0_jabr = time.time()

  if all_data['jabrcuts']:
    jabr_cuts(log,all_data)

    if all_data['NO_jabrs_violated']:
      if all_data['threshold'] > all_data['tolerance']:
        all_data['threshold'] *= 1e-01
        log.joint(' threshold updated to ' + str(all_data['threshold'])
                  + '\n' )

  t1_jabr = time.time()

  log.joint(' time spent on Jabr-cuts ' + str(t1_jabr - t0_jabr) + '\n')

  t0_i2 = time.time()

  if all_data['i2cuts']:
    i2_cuts(log,all_data)

    if all_data['NO_i2_cuts_violated']:
      if all_data['threshold_i2'] > all_data['tolerance']:
        all_data['threshold_i2'] *= 1e-01
        log.joint(' threshold updated to ' + str(all_data['threshold_i2'])
                  + '\n' )

  t1_i2 = time.time()
  log.joint(' time spent on i2-cuts ' + str(t1_i2 - t0_i2) + '\n')


  t0_lim = time.time()

  if all_data['limitcuts']:
    limit_cuts(log,all_data)
    if all_data['NO_limit_cuts_violated']:
      if all_data['threshold_limit'] > all_data['tolerance']:
        all_data['threshold_limit'] *= 1e-01
        log.joint(' threshold updated to ' + str(all_data['threshold_limit'])
                  + '\n' )

  t1_lim = time.time()
  log.joint(' time spent on lim-cuts ' + str(t1_lim - t0_lim) + '\n')

  
  if all_data['losscuts']:
    loss_cuts(log,all_data)
    if all_data['NO_loss_violated']:
      if all_data['threshold_limit'] > all_data['tolerance']:
        all_data['threshold_limit'] *= 1e-01
        log.joint(' threshold updated to ' + str(all_data['threshold_limit'])
                  + '\n' )

  if all_data['linear_objective']:
    if all_data['objective_cuts']:
      objective_cuts(log,all_data)

  t1_cuts = time.time()

  log.joint('\n time spent on cuts ' + str(t1_cuts - t0_cuts) + '\n')


def cutplane_cutmanagement(log,all_data):

  t0_cutmanag = time.time()

  if all_data['jabrcuts'] and all_data['dropjabr'] and (all_data['round'] >= all_data['cut_age_limit']):
    drop_jabr(log,all_data)

  if all_data['i2cuts'] and all_data['dropi2'] and (all_data['round'] >= all_data['cut_age_limit']):
    drop_i2(log,all_data)

  if all_data['limitcuts'] and all_data['droplimit'] and (all_data['round'] >= all_data['cut_age_limit']):
    drop_limit(log,all_data)

  if all_data['losscuts'] and all_data['droploss'] and (all_data['round'] >= all_data['cut_age_limit']):
    drop_loss(log,all_data)

  if all_data['loss_inequalities'] and all_data['droploss'] and (all_data['round'] >= all_data['cut_age_limit']):
    drop_loss(log,all_data)

  if all_data['cut_analysis']:
    cut_analysis(log,all_data)

  log.joint('\n')

  t1_cutmanag = time.time()

  log.joint(' time spent on cut management ' + str(t1_cutmanag - t0_cutmanag)
            + '\n')


def cutplane_optimize(log,all_data):

  themodel = all_data['themodel']

  log.joint(' solving model with method ' + str(themodel.params.method) + '\n')
  log.joint(' crossover ' + str(themodel.params.crossover) + '\n')

  # Experiments for arXiv paper
  
  if all_data['losstest']:
    all_data['rbranch'] = randombranch(log,all_data)

  #####
    
  t0_solve = time.time()
  themodel.optimize()
  t1_solve = time.time()

  if themodel.status == GRB.status.INF_OR_UNBD:
    log.joint(' -> LP infeasible or unbounded\n')
    log.joint(' turning presolve off and reoptimizing\n')
    themodel.params.presolve = 0

    t0_solve = time.time()
    themodel.optimize()
    t1_solve = time.time()

    #breakexit('Computing IIS (if infeasible)')
    #themodel.computeIIS()
    #themodel.write('model.ilp')

    all_data['runtime'] = t1_solve - all_data['T0']

    log.joint(' writing casename, opt status, and runtime to summary_ws.log\n')

    summary_ws = open("summary_ws.log","a+") 
    summary_ws.write(' case ' + all_data['casename'] + ' opt_status ' 
                     + str(themodel.status) + ' runtime ' 
                     + str(all_data['runtime']) + ' iterations ' 
                     + str(all_data['round']) + '\n')

    summary_ws.close()

    log.joint(' optimization status ' + str(themodel.status) + '\n')
    log.joint(' solver runtime current round = %g\n' % (t1_solve - t0_solve) )
    log.joint(' overall time = %g\n' % all_data['runtime'] )
    log.joint(' bye.\n')
    exit(0)

  elif themodel.status == GRB.status.INFEASIBLE:
    log.joint(' -> LP infeasible\n')

    all_data['runtime'] = time.time() - all_data['T0']

    #breakexit('Computing IIS')

    #themodel.computeIIS()
    #themodel.write('model.ilp')

    log.joint(' writing casename, opt status, and runtime to summary_ws.log\n')

    summary_ws = open("summary_ws.log","a+") #later add feasibility errors, etc                            
    summary_ws.write(' case ' + all_data['casename'] + ' opt_status ' 
                     + str(themodel.status) + ' runtime ' 
                     + str(all_data['runtime']) + ' iterations ' 
                     + str(all_data['round']) + '\n')

    summary_ws.close()

    log.joint(' solver runtime current round = %g\n' % (t1_solve - t0_solve) )
    log.joint(' overall time = %g\n' % all_data['runtime'] )
    log.joint(' bye.\n')
    exit(0)

    
  elif ( themodel.status != GRB.status.OPTIMAL and 
         themodel.status != GRB.status.SUBOPTIMAL and
         themodel.status != GRB.status.NUMERIC):

    log.joint(' -> solver terminated with status ' + str(themodel.status) + '\n')

    all_data['runtime'] = time.time() - all_data['T0']

    log.joint(' writing casename, opt status and runtime to summary_ws.log\n')

    summary_ws = open("summary_ws.log","a+") #later add feasibility errors, etc                            
    summary_ws.write(' case ' + all_data['casename'] + ' opt_status ' 
                     + str(themodel.status) + ' runtime ' 
                     + str(all_data['runtime']) + ' iterations ' 
                     + str(all_data['round']) + '\n')

    summary_ws.close()

    log.joint(' solver runtime current round = %g\n' % (t1_solve - t0_solve) )
    log.joint(' overall time = %g\n' % all_data['runtime'] )
    log.joint(' bye.\n')
    exit(0)

  all_data['objval']                  = themodel.ObjVal
  all_data['optstatus']               = themodel.status
  all_data['solvertime']              = t1_solve - t0_solve
  all_data['cumulative_solver_time'] += (t1_solve - t0_solve)
  if (themodel.status != GRB.status.NUMERIC):
    all_data['dinfs']                   = themodel.DualResidual
    all_data['dinfs_scaled']            = themodel.DualSResidual
  else:
    all_data['dinfs'] = -1
    all_data['dinfs_scaled'] = -1

def getsol_ampl(log,all_data):

  casename    = all_data['casename']
  filename    = 'sols/AMPLsol_'+ casename +'.txt'

  try:
    thefile   = open(filename, "r")
    lines     = thefile.readlines()
    lenlines  = len(lines)
    thefile.close()
  except:
    log.stateandquit("cannot open file " + datafilename)
    sys.exit("failure")

  branches     = all_data['branches']
  buses        = all_data['buses']
  IDtoCountmap = all_data['IDtoCountmap'] 

  sol_obj        = 0
  sol_vm         = {}
  sol_angle      = {}
  sol_cvalues    = {}
  sol_svalues    = {}
  sol_Pfvalues   = {}
  sol_Ptvalues   = {}
  sol_Qfvalues   = {}
  sol_Qtvalues   = {}
  sol_GenPvalues = {}

  linenum = 0
  log.joint(' reading file\n')
  while linenum < lenlines:
    thisline = lines[linenum].split()
    if thisline[0] == 'value':
      sol_obj              = float(lines[0].split()[1])
    elif thisline[0] == 'bus':
      buscount             = int(thisline[1])
      bus                  = buses[buscount]
      sol_vm[bus]          = float(thisline[3])
      sol_angle[bus]       = float(thisline[5]) 
      sol_cvalues[bus]     = sol_vm[bus]**2
    elif thisline[0] == 'branch':
      branchcount          = int(thisline[1])
      branch               = branches[branchcount]
      sol_Pfvalues[branch] = float(thisline[7])
      sol_Ptvalues[branch] = float(thisline[9])
      sol_Qfvalues[branch] = float(thisline[11])
      sol_Qtvalues[branch] = float(thisline[13])
    elif thisline[0] == 'genid':
      genid      = int(thisline[1])
      gen_nodeID = int(thisline[3])
      sol_GenPvalues[genid] = float(thisline[5])

    linenum += 1

  for branch in branches.values():
    f          = branch.f
    t          = branch.t
    count_of_f = IDtoCountmap[f]
    count_of_t = IDtoCountmap[t]
    bus_f      = buses[count_of_f]
    bus_t      = buses[count_of_t]
    vm_f       = sol_vm[bus_f]
    vm_t       = sol_vm[bus_t]

    sol_cvalues[branch] = vm_f * vm_t * math.cos(sol_angle[bus_f] - sol_angle[bus_t])
    sol_svalues[branch] = vm_f * vm_t * math.sin(sol_angle[bus_f] - sol_angle[bus_t])

  all_data['sol_vm']       = sol_vm
  all_data['sol_angle']    = sol_angle
  all_data['sol_cvalues']  = sol_cvalues
  all_data['sol_svalues']  = sol_svalues
  all_data['sol_Pfvalues'] = sol_Pfvalues
  all_data['sol_Ptvalues'] = sol_Ptvalues
  all_data['sol_Qfvalues'] = sol_Qfvalues
  all_data['sol_Qtvalues'] = sol_Qtvalues
  all_data['sol_GenPvalues'] = sol_GenPvalues

  log.joint(' knitro solution loaded\n')

  #log.joint(' sol_GenPvalues ' + str(sol_GenPvalues) + '\n')

  #flowdecomp(log,all_data)
  #breakexit('c')

def getsol_ampl_mtp(log,all_data):

  T           = all_data['T']
  casename    = all_data['casename']
  casetype    = all_data['casetype']


  filename      = 'ACsol_' + all_data['casename'] + '_' + str(T) + '_' + casetype + '.txt'

  try:
    thefile   = open(filename, "r")
    lines     = thefile.readlines()
    lenlines  = len(lines)
    thefile.close()
  except:
    log.stateandquit("cannot open file " + filename)
    sys.exit("failure")

  branches     = all_data['branches']
  buses        = all_data['buses']
  IDtoCountmap = all_data['IDtoCountmap'] 
  gens         = all_data['gens']
  tolerance    = 1e-05

  
  sol_obj        = 0
  sol_vm         = {}
  sol_angle      = {}
  sol_cvalues    = {}
  sol_svalues    = {}
  sol_Pfvalues   = {}
  sol_Ptvalues   = {}
  sol_Qfvalues   = {}
  sol_Qtvalues   = {}
  sol_GenPvalues = {}
  sol_GenQvalues = {}
  
  for k in range(T):
      sol_vm[k]         = {}
      sol_angle[k]      = {}
      sol_cvalues[k]    = {}
      sol_svalues[k]    = {}
      sol_Pfvalues[k]   = {}
      sol_Ptvalues[k]   = {}
      sol_Qfvalues[k]   = {}
      sol_Qtvalues[k]   = {}
      sol_GenPvalues[k] = {}
      sol_GenQvalues[k] = {}
      
  linenum = 0
  log.joint(' reading file\n')
  while linenum < lenlines:
    thisline = lines[linenum].split()
    if thisline[0] == 'objvalue':
      sol_obj              = float(thisline[1])
    elif thisline[0] == 'bus':
      #log.joint(' thisline ' + str(thisline) + '\n')
      buscount             = int(thisline[1])
      k                    = int(thisline[7])
      bus                  = buses[buscount]
      sol_vm[k][bus]       = float(thisline[3])
      sol_angle[k][bus]    = float(thisline[5]) 
      sol_cvalues[k][bus]  = float(thisline[3])**2
    elif thisline[0] == 'branch':
      #log.joint(' thisline ' + str(thisline) + '\n')      
      branchcount             = int(thisline[1])
      k                       = int(thisline[19])
      branch                  = branches[branchcount]
      sol_Pfvalues[k][branch] = float(thisline[7])
      sol_Ptvalues[k][branch] = float(thisline[9])
      sol_Qfvalues[k][branch] = float(thisline[11])
      sol_Qtvalues[k][branch] = float(thisline[13])
      sol_cvalues[k][branch]  = float(thisline[15])
      sol_svalues[k][branch]  = float(thisline[17])      
    elif thisline[0] == 'genid':
      genid      = int(thisline[1])
      k          = int(thisline[9])
      gen_nodeID = int(thisline[3])
      sol_GenPvalues[k][genid] = float(thisline[5])
      sol_GenQvalues[k][genid] = float(thisline[7])      

    linenum += 1

  all_data['sol_obj']        = sol_obj
  all_data['sol_vm']         = sol_vm
  all_data['sol_angle']      = sol_angle
  all_data['sol_cvalues']    = sol_cvalues
  all_data['sol_svalues']    = sol_svalues
  all_data['sol_Pfvalues']   = sol_Pfvalues
  all_data['sol_Ptvalues']   = sol_Ptvalues
  all_data['sol_Qfvalues']   = sol_Qfvalues
  all_data['sol_Qtvalues']   = sol_Qtvalues
  all_data['sol_GenPvalues'] = sol_GenPvalues
  all_data['sol_GenQvalues'] = sol_GenQvalues
  
  log.joint(' knitroampl multi-time period solution loaded\n')
  

def writesol(log,all_data):
    
  branches     = all_data['branches']
  buses        = all_data['buses']
  gens         = all_data['gens']
  IDtoCountmap = all_data['IDtoCountmap']
  cvalues      = all_data['cvalues']
  svalues      = all_data['svalues']
  GenPvalues   = all_data['GenPvalues']
  GenQvalues   = all_data['GenQvalues']
  Pfvalues     = all_data['Pfvalues']
  Ptvalues     = all_data['Ptvalues']
  Qfvalues     = all_data['Qfvalues']
  Qtvalues     = all_data['Qtvalues']
  i2fvalues    = all_data['i2fvalues']
  T            = all_data['T']
  alphadic     = all_data['alphadic']
  casetype     = all_data['casetype']
  #i2tvalues     = {}    

  filename = all_data['sols'] + 'CPsol_' + all_data['casename'] + '_' + str(T) + '_' + casetype + '.txt'
  thefile  = open(filename,'w+')
  
  log.joint(' writing solution to ' + filename + '\n')
  
  machinename    = "cool1"
  now            = time.time()
  solver_version = 'Gurobi 10.0.1'
  opsystem       = 'Fedora 34 (Workstation Edition)'
  processor      = 'Intel(R) Xeon(R) Linux64 CPU E5-2687W v3 3.10GHz'
  cores          = '20 physical cores, 40 logical processors'
  ram            = '256 GB RAM'
    
    
  thefile.write('/CUTPLANEsolution : ' + all_data['casename'] + '\n')
  thefile.write('/Date : ' + str(time.strftime('%m-%d-%Y %H:%M:%S %Z', time.localtime(now))) + '\n')
  thefile.write('/MachineName : ' + machinename + '\n')
  thefile.write('/Processor : ' + processor + '\n')
  thefile.write('/OS : ' + opsystem + '\n')
  thefile.write('/Cores : ' + cores + '\n')
  thefile.write('/RAM : ' + ram + '\n')
  thefile.write('/Solver : ' + solver_version + '\n')
  thefile.write('objvalue ' + str(all_data['objval']) + '\n')
  thefile.write('round ' + str(all_data['round']) + '\n')
  thefile.write('time-periods ' + str(all_data['T']) + '\n')
  
  thefile.write('voltages:\n')

  for buscount in buses.keys():
    for k in range(T):
      bus   = buses[buscount]
      v2val = cvalues[k][bus]
      vval  = (v2val)**(0.5)
      line  = 'bus ' + str(buscount) + ' M ' + str(vval) + ' k ' + str(k) + '\n'
    thefile.write(line)

  thefile.write('power flows and cs variables:\n')

  for branchid in branches.keys():
    for k in range(T):
      branch     = branches[branchid]
      f          = branch.f
      t          = branch.t
      count_of_f = IDtoCountmap[f]
      count_of_t = IDtoCountmap[t]
      Pfval      = Pfvalues[k][branch]
      Ptval      = Ptvalues[k][branch]
      Qfval      = Qfvalues[k][branch]
      Qtval      = Qtvalues[k][branch]
      cftval     = cvalues[k][branch]
      sftval     = svalues[k][branch]
      alpha      = alphadic[branch]

      if alpha < all_data['alpha_threshold']:
        i2fval     = i2fvalues[k][branch]
        line = 'branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' Pft ' + str(Pfval) + ' Ptf ' + str(Ptval) + ' Qft ' + str(Qfval) + ' Qtf ' + str(Qtval) + ' cft ' + str(cftval) + ' sft ' + str(sftval) + ' i2ft ' + str(i2fval) + ' k ' + str(k) + '\n'
      else:
        line = 'branch ' + str(branchid) + ' f ' + str(f) + ' t ' + str(t) + ' Pft ' + str(Pfval) + ' Ptf ' + str(Ptval) + ' Qft ' + str(Qfval) + ' Qtf ' + str(Qtval) + ' cft ' + str(cftval) + ' sft ' + str(sftval) + ' k ' + str(k) + '\n'
      thefile.write(line)

  thefile.write('generation:\n')

  for genid in gens.keys():
    for k in range(T):
      gen     = gens[genid] 
      nodeID  = gen.nodeID
      GenPval = GenPvalues[k][gen]
      GenQval = GenQvalues[k][gen]
      line_gen = 'genid ' + str(genid) + ' bus ' + str(nodeID) + ' GP ' + str(GenPval) + ' GQ ' + str(GenQval) + ' k ' + str(k) + '\n'
      thefile.write(line_gen)

  thefile.close()

  log.joint(' done writing CUTPLANE solution to .txt file\n\n')



def writesol_allvars(log,all_data):

  branches     = all_data['branches']
  buses        = all_data['buses']
  gens         = all_data['gens']
  IDtoCountmap = all_data['IDtoCountmap']
  cvalues      = all_data['cvalues']
  svalues      = all_data['svalues']
  GenPvalues   = all_data['GenPvalues']
  GenQvalues   = all_data['GenQvalues']
  Pfvalues     = all_data['Pfvalues']
  Ptvalues     = all_data['Ptvalues']
  Qfvalues     = all_data['Qfvalues']
  Qtvalues     = all_data['Qtvalues']
  i2fvalues    = all_data['i2fvalues']
  T            = all_data['T']
  alphadic     = all_data['alphadic']
  Pd           = all_data['Pd']
  casetype     = all_data['casetype']
  #i2tvalues     = {}    


  filename    = all_data['sols'] + 'CPsol_' + all_data['casename'] + '_' + str(T) + '_' + casetype + '.sol'  
  thefilevars = open(filename,'w+')

  log.joint(' writing solution to ' + filename + '\n')

  machinename    = "cool1"
  now            = time.time()
  solver_version = 'Gurobi 10.0.1'
  opsystem       = 'Fedora 34 (Workstation Edition)'
  processor      = 'Intel(R) Xeon(R) Linux64 CPU E5-2687W v3 3.10GHz'
  cores          = '20 physical cores, 40 logical processors'
  ram            = '256 GB RAM'
    
    
  thefilevars.write('/CUTPLANEsolution : ' + all_data['casename'] + '\n')
  thefilevars.write('/Date : ' + str(time.strftime('%m-%d-%Y %H:%M:%S %Z', time.localtime(now))) + '\n')
  thefilevars.write('/MachineName : ' + machinename + '\n')
  thefilevars.write('/Processor : ' + processor + '\n')
  thefilevars.write('/OS : ' + opsystem + '\n')
  thefilevars.write('/Cores : ' + cores + '\n')
  thefilevars.write('/RAM : ' + ram + '\n')
  thefilevars.write('/Solver : ' + solver_version + '\n')
  thefilevars.write('/Objvalue ' + str(all_data['objval']) + '\n')
  thefilevars.write('/Time-periods ' + str(all_data['T']) + '\n')

  
  for buscount in buses.keys():
    for k in range(T):
      bus        = buses[buscount]
      f          = bus.nodeID
      v2value    = cvalues[k][bus]        

      v2name     = 'c_' + str(f) + '_' + str(f) + '_' + str(k)
      v2line     = v2name + ' = ' + str(v2value) + '\n'
      thefilevars.write(v2line)

      IPvalue = - Pd[k][bus]
      if all_data['arpae2']:
        IQvalue = - all_data['Qd'][k][bus]
      else:
        IQvalue = - bus.Qd
      IPname  = 'IP_' + str(bus.nodeID) + '_' + str(k)
      IQname  = 'IQ_' + str(bus.nodeID) + '_' + str(k)

      for gencounter in bus.genidsbycount:
        if gens[gencounter].status:
          gen      = gens[gencounter]
          IPvalue += GenPvalues[k][gen]
          IQvalue += GenQvalues[k][gen]

      IPline = IPname + ' = ' + str(IPvalue) + '\n'
      thefilevars.write(IPline)

      IQline = IQname + ' = ' + str(IQvalue) + '\n'
      thefilevars.write(IQline)


  for branchid in branches.keys():
    for k in range(T):
      branch     = branches[branchid]
      f          = branch.f
      t          = branch.t
      count_of_f = IDtoCountmap[f]
      count_of_t = IDtoCountmap[t]
      Pfval      = Pfvalues[k][branch]
      Ptval      = Ptvalues[k][branch]
      Qfval      = Qfvalues[k][branch]
      Qtval      = Qtvalues[k][branch]
      cftval     = cvalues[k][branch]
      sftval     = svalues[k][branch]
      alpha      = alphadic[branch]
      
      Pfname  = 'P_' + str(branchid) + '_' + str(f) + '_' + str(t) + '_' + str(k)
      Ptname  = 'P_' + str(branchid) + '_' + str(t) + '_' + str(f) + '_' + str(k)
      Qfname  = 'Q_' + str(branchid) + '_' + str(f) + '_' + str(t) + '_' + str(k)
      Qtname  = 'Q_' + str(branchid) + '_' + str(t) + '_' + str(f) + '_' + str(k)
      cftname = 'c_' + str(branchid) + '_' + str(f) + '_' + str(t) + '_' + str(k)
      sftname = 's_' + str(branchid) + '_' + str(f) + '_' + str(t) + '_' + str(k)

      Pfline  = Pfname + ' = ' + str(Pfval) + '\n'
      thefilevars.write(Pfline)

      Ptline  = Ptname + ' = ' + str(Ptval) + '\n'
      thefilevars.write(Ptline)

      Qfline  = Qfname + ' = ' + str(Qfval) + '\n'
      thefilevars.write(Qfline)

      Qtline  = Qtname + ' = ' + str(Qtval) + '\n'
      thefilevars.write(Qtline)

      cftline  = cftname + ' = ' + str(cftval) + '\n'
      thefilevars.write(cftline)

      sftline  = sftname + ' = ' + str(sftval) + '\n'
      thefilevars.write(sftline)

      if alpha < all_data['alpha_threshold']:
        i2fval     = i2fvalues[k][branch]
        i2fname = 'i2_'+ str(branchid) + '_' + str(f) + '_' + str(t) + '_' + str(k)      
        i2fline  = i2fname + ' = ' + str(i2fval) + '\n'
        thefilevars.write(i2fline)                        


  for genid in gens.keys():
    for k in range(T):
      gen     = gens[genid] 
      nodeID  = gen.nodeID
      GenPval = GenPvalues[k][gen]
      GenQval = GenQvalues[k][gen]
      GPname  = "GP_" + str(genid) + "_" + str(nodeID) + "_" + str(k)
      GQname  = "GQ_" + str(genid) + "_" + str(nodeID) + "_" + str(k)

      GPline  = GPname + ' = ' + str(GenPval) + '\n'
      thefilevars.write(GPline)

      GQline  = GQname + ' = ' + str(GenQval) + '\n'
      thefilevars.write(GQline)


  log.joint(' done writing CUTPLANE allvars solution to .sol file\n\n')

  thefilevars.close()

def getloads(log,all_data,loads):

  lines = loads.readlines()
  lenlines = len(lines) - 1 # END                                                                                        
  
  buses        = all_data['buses']
  IDtoCountmap = all_data['IDtoCountmap']
  T            = all_data['T']

  Pd      = {}
  for k in range(T):
    Pd[k] = {}
    
  linenum = 0

  log.joint(' reading file with loads\n')
  while linenum < lenlines:
    thisline = lines[linenum].split()
    buscount         = int(thisline[1])
    bus              = buses[buscount]
    k                = int(thisline[5])
    load             = float(thisline[7])
    Pd[k][bus]       = load
    linenum         += 1

  return  Pd


def getloads2(log,all_data,loads):

  lines = loads.readlines()
  lenlines = len(lines) - 1 # END                                                                                        
  buses        = all_data['buses']
  IDtoCountmap = all_data['IDtoCountmap']
  T            = all_data['T']

  Pd      = {}
  #Qd      = {}
  for k in range(T):
    Pd[k] = {}
    #Qd[k] = {}
    for bus in buses.values():
      Pd[k][bus] = 0
      #Qd[k][bus] = 0
    
  linenum = 0

  log.joint(' reading file with loads\n')
  while linenum < lenlines:
    thisline = lines[linenum].split()
    bus_id           = int(thisline[1])
    buscount         = IDtoCountmap[bus_id]
    bus              = buses[buscount]
    k                = int(thisline[3])
    load_P           = float(thisline[5])
    #load_Q           = float(thisline[7])    
    Pd[k][bus]       = load_P
    #Qd[k][bus]       = load_Q    
    linenum         += 1

  all_data['Pd'] = Pd
  #all_data['Qd'] = Qd
    
  #return  Pd, Qd
  return  Pd



def getrampr(log,all_data,loads):

  lines = loads.readlines()
  lenlines = len(lines) - 1 # END                                                                                        
  
  gens         = all_data['gens']
  IDtoCountmap = all_data['IDtoCountmap']
  T            = all_data['T']

  rampru      = {}
  ramprd      = {}
  for k in range(T):
    rampru[k] = {}
    ramprd[k] = {}    
    
  linenum = 0

  log.joint(' reading file with ramprates\n')
  while linenum < lenlines:
    thisline = lines[linenum].split()
    gencount         = int(thisline[1])
    gen              = gens[gencount]
    k                = int(thisline[5])
    rpru             = float(thisline[7])
    rprd             = float(thisline[9])
    rampru[k][gen]   = rpru
    ramprd[k][gen]   = rprd
    linenum         += 1

  return  rampru, ramprd


def getduals(log,all_data):


  buses    = all_data['buses']
  themodel = all_data['themodel']
  T        = all_data['T']
  rnd      = all_data['round']
  duals    = {}
  
  for k in range(T):
    for bus in buses.values():
      nodeID     = bus.nodeID
      constrname = "PBaldef"+str(bus.nodeID)+"_"+str(k)
      constr     = themodel.getConstrByName(constrname)
      dual       = constr.Pi
      duals[constrname] = dual 
  

  all_data['duals'][rnd] = duals
  
  if rnd > 1:
    sqdiff = 0
    duals_prev = all_data['duals'][rnd-1]
    for constrname in duals.keys():
      sqdiff += (duals_prev[constrname] - duals[constrname])**2

    dual_diff = round(math.sqrt(sqdiff),4)
    all_data['dual_diff'][rnd] = dual_diff
    log.joint(' dual diff at round  ' + str(rnd) + ' = ' + str(dual_diff) + '\n')
    #breakexit(' check dual diff')
    

def print_duals(log,all_data):

  rnd = all_data['round']
  
  log.joint('dual diffs:\n')
  for t in range(2,rnd+1):
    log.joint('(' + str(t) + ',' + str(all_data['dual_diff'][t]) + ')')

  log.joint('\n')

def print_duals_table(log,all_data):

  rnd = all_data['round']
  
  log.joint('(table) dual diffs:\n')
  for t in range(2,rnd+1):
    log.joint(str(all_data['dual_diff'][t]) + ' & ') 

  log.joint('\n')    
