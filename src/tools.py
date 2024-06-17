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

def getangles(log,all_data):

    svalues  = all_data['svalues']
    cvalues  = all_data['cvalues']
    refbus   = all_data['refbus']
    buses    = all_data['buses']
    branches = all_data['branches']
    IDtoCountmap = all_data['IDtoCountmap']

    angles           = {}
    stack            = []
    unexplored_buses = list(buses.keys())
    anglediff        = {}

    angle = 0
    angles[refbus] = 0
    unexplored_buses.remove(refbus)
    stack.append(refbus)

    loud = False

    if loud:
        log.joint(' refbus ' + str(refbus) + '\n')
        log.joint(' stack ' + str(stack) + '\n')
        log.joint(' unexplored buses ' + str(unexplored_buses) + '\n')
        log.joint(' --- \n')

    while unexplored_buses:
        busid          = stack[0]
        bus            = buses[busid]
        from_branchids = list(bus.frombranchids.values())
        to_branchids   = list(bus.tobranchids.values())
        angle          = angles[busid]
        stack.remove(busid)

        if loud:
            log.joint(' current bus ' + str(busid) + ' angle ' + str(angle) 
                      + '\n')
            log.joint(' from branchids ' + str(from_branchids) + '\n')
            log.joint(' to branchids ' + str(to_branchids) + '\n')
            log.joint(' stack ' + str(stack) + '\n')
            log.joint(' unexplored buses ' + str(unexplored_buses) + '\n')

        for branchid in from_branchids:
            branch = branches[branchid] 
            t      = branch.t
            id_t   = IDtoCountmap[t]
            if id_t in unexplored_buses:
                s = svalues[branch]
                c = cvalues[branch]
                theta = math.atan2(s,c)
                anglediff[branchid] = theta
                angles[id_t] = - theta + angle
                if loud:
                    log.joint(' new bus found, busid ' + str(id_t) + ' t ' 
                              + str(t) + '\n')
                    log.joint(' theta_f_t ' + str(theta) + ' theta_f ' 
                              + str(angle) + '\n')
                    log.joint(' theta_t ' + str(angles[id_t]) + '\n')
                
                unexplored_buses.remove(id_t)
                stack.append(id_t)
                if loud:
                    log.joint(' stack updated ' + str(stack) + '\n')
                    log.joint(' unexplored buses updated ' + str(unexplored_buses) 
                              + '\n')

        for branchid in to_branchids:
            branch = branches[branchid] 
            f      = branch.f
            id_f   = IDtoCountmap[f]
            if id_f in unexplored_buses:
                s = svalues[branch]
                c = cvalues[branch]
                theta = math.atan2(s,c)
                anglediff[branchid] = theta
                angles[id_f] = theta + angle
                
                if loud:
                    log.joint(' new bus found, busid ' + str(id_f) + ' f ' 
                              + str(f) + '\n')
                    log.joint(' theta_f_t ' + str(theta) + ' theta_t ' 
                              + str(angle) + '\n')
                    log.joint(' theta_f ' + str(angles[id_f]) + '\n')

                unexplored_buses.remove(id_f)
                stack.append(id_f)
                if loud:
                    log.joint(' stack updated ' + str(stack) + '\n')
                    log.joint(' unexplored buses updated ' + str(unexplored_buses) 
                              + '\n')
        if loud:
            log.joint(' we are done with busid ' + str(busid) + '\n')
            log.joint('\n')
        #stack.remove(busid)
    
    if loud:
        log.joint(' angles ' + str(angles) + '\n')
        log.joint(' anglediff ' + str(anglediff) + '\n')
        log.joint(' branches ' + str(branches.keys()) + '\n')
        log.joint(' buses ' + str(buses.keys()) + '\n')
        log.joint(' unexplored buses ' + str(unexplored_buses) + '\n')
        log.joint(' refbus ' + str(refbus) + '\n')



    all_data['angles'] = angles
                


    

        





