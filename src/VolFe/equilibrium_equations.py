# equilibrium_equations.py

import pandas as pd
import numpy as np
import warnings as w

import VolFe.differential_equations as de
import VolFe.melt_gas as mg
import VolFe.model_dependent_variables as mdv
import VolFe.calculations as c


def set_system(melt_wf,models):
    wt_H = melt_wf["HT"]
    wt_C = melt_wf["CT"]
    wt_S = melt_wf["ST"]
    wt_X = melt_wf["XT"]
    if wt_H > 0. and wt_C == 0. and wt_S == 0. and wt_X == 0.:
        sys = "HOFe"
    elif wt_H == 0. and wt_C > 0. and wt_S == 0. and wt_X == 0.:
        sys = "COFe"
    elif wt_H == 0. and wt_C == 0. and wt_S > 0. and wt_X == 0.:
        sys = "SOFe"  
    elif wt_H > 0. and wt_C > 0. and wt_S == 0. and wt_X == 0.:
        sys = "CHOFe"
    elif wt_H > 0. and wt_C == 0. and wt_S > 0. and wt_X == 0.:
        sys = "SHOFe"
    elif wt_H == 0. and wt_C > 0. and wt_S > 0. and wt_X == 0.:
        sys = "SCOFe"
    elif wt_H == 0. and wt_C > 0. and wt_S == 0. and wt_X > 0.:
        sys = "COXFe"
    elif wt_H > 0. and wt_C == 0. and wt_S == 0. and wt_X > 0.:
        sys = "HOXFe"
    elif wt_H > 0. and wt_C > 0. and wt_S > 0. and wt_X == 0.:
        sys = "SCHOFe"
    elif wt_H > 0. and wt_C > 0. and wt_S == 0. and wt_X > 0.:
        sys = "CHOXFe"
    elif wt_H > 0. and wt_C > 0. and wt_S > 0. and wt_X > 0.:
        sys = "SCHOXFe"
    return sys

def initial_guesses(run,PT,melt_wf,setup,models,system): ### CHECK ###
    starting_P = models.loc["starting_P", "option"]
    #xenia = models.loc["xenia","option"]
    solve_species = models.loc["solve_species","option"]
    
    if starting_P == "set":
        guessx = setup.loc[run,"xg_O2"]
    else:
        guessx = mg.xg_O2(PT,melt_wf,models)
    
    if models.loc["COH_species","option"] == "H2O-CO2 only":
        guessx = mg.xg_CO2(PT,melt_wf,models)
    
    if system in ["COFe","HOFe","SOFe","CHOFe","COXFe","SHOFe","SCOFe","CHOXFe","HOXFe"]:
        guessw = 0.
        
    if system in ["COFe","HOFe","SOFe"]:
        guessy = 0.
        guessz = 0.
    
    if system in ["COXFe","CHOFe","SHOFe","SCOFe","HOXFe"]:
        guessz = 0.
    
    if system in ["CHOFe","COXFe"]:
        if starting_P == "set":
            guessy = setup.loc[run,"xg_CO"]
        elif models.loc["COH_species","option"] == "H2O-CO2 only":
            guessy = 0.
        else:
            guessy = mg.xg_CO(PT,melt_wf,models) 
    elif system == "SHOFe":
        if starting_P == "set": 
            guessy = setup.loc[run,"xg_S2"]
        else:
            guessy = mg.xg_S2(PT,melt_wf,models)
    elif system == "HOXFe":
        if starting_P == "set": 
            guessy = setup.loc[run,"xg_H2"]
        else:
            guessy = mg.xg_H2(PT,melt_wf,models)
    elif system == "SCOFe":
        guessy = mg.xg_S2(PT,melt_wf,models) 
    
    elif system == "CHOXFe":
        if starting_P == "set":
            guessy = setup.loc[run,"xg_CO"]
            guessz = setup.loc[run,"xg_X"]
        else:
            guessy = mg.xg_CO(PT,melt_wf,models) 
            guessz = mg.xg_X(PT,melt_wf,models)
    elif system == "SCHOFe" or system == "SCHOXFe":
        if solve_species in ["OCS"]:
            if starting_P == "set":
                guessy = setup.loc[run,"xg_CO"]
                guessz = setup.loc[run,"xg_S2"]
                guessw = 0.
            else:
                guessy = mg.xg_CO(PT,melt_wf,models) 
                guessz = mg.xg_S2(PT,melt_wf,models)
                guessw = mg.xg_H2(PT,melt_wf,models)
        elif solve_species == "OHS":
            if starting_P == "set":
                guessy = setup.loc[run,"xg_H2"]
                guessz = setup.loc[run,"xg_S2"]
                guessw = 0.
            else:
                guessy = mg.xg_H2(PT,melt_wf,models) 
                guessz = mg.xg_S2(PT,melt_wf,models)
                guessw = mg.xg_CO(PT,melt_wf,models)
        elif solve_species == "OCH":
            if starting_P == "set":
                guessy = setup.loc[run,"xg_CO"]
                guessz = setup.loc[run,"xg_H2"]
                guessw = 0.
            else:
                guessy = mg.xg_CO(PT,melt_wf,models) 
                guessz = mg.xg_H2(PT,melt_wf,models)
                guessw = mg.xg_S2(PT,melt_wf,models)
        if system == "SCHOXFe":
            if starting_P == "set":
                guessw = setup.loc[run,"xg_X"]
            else:
                guessw = mg.xg_X(run,PT,melt_wf,setup,models) 
    guesses = {"guessx":guessx, "guessy":guessy, "guessz":guessz, "guessw":guessw}
    return guesses

def mg_equilibrium(PT,melt_wf,bulk_wf,models,nr_step,nr_tol,guesses): ### CHECK ###
    system = set_system(melt_wf,models)
    solve_species = models.loc["solve_species","option"]

    if system in ["COFe","SOFe","SCOFe","COXFe"]: # no H
        wt_H_, xg_H2_, xg_H2O_, xg_CH4_, xg_H2S_, xm_H2O_, wm_H2O_, wm_H2_, wm_CH4_, wm_H2S_, wt_g_H = 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,""
    if system in ["COFe","HOFe","CHOFe","CHOXFe","COXFe","HOXFe"]: # no S
        wt_S_, xg_S2_, xg_SO2_, xg_H2S_, xg_OCS_, wm_S_, wm_SO3_, wm_ST_, wm_H2S_, S62, S6T, wt_g_S = 0.,0.,0.,0.,0.,0.,0.,0.,0.,"","",""
    if system in ["HOFe","SOFe","SHOFe","HOXFe"]: # no C
        wt_C_, xg_CO_, xg_CO2_, xg_CH4_, xg_OCS_, xm_CO2_, wm_CO2mol_, wm_CO2carb_, wm_CO2_, wm_CO_, wm_CH4_, wt_g_C = 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,""
    if system in ["COFe","HOFe","SOFe","SHOFe","CHOFe","SCOFe","CHOXFe","SCHOFe"]: # no X
        wt_X_, xg_X_, wm_X_, wt_g_X = 0.,0.,0.,""
    if system in ["COFe","SOFe","HOFe","HOFe_xenia"]: # one component
        guessy,guessz,guessw = "","",""
        solve_species = "O"
    if system in ["SHOFe","CHOFe","SCOFe","COXFe","HOXFe"]: # two components
        guessz,guessw = "",""
    if system in ["CHOXFe","SCHOFe"]: # three components
        guessw = ""
    
    if system in ["HOFe","SHOFe","CHOFe","CHOXFe","SCHOFe","HOXFe"]:
        if models.loc["Hspeciation","option"] == "ideal":
            print("not currently possible")
        if models.loc["Hspeciation","option"] == "regular":
            print("not currently possible")
    
    if system == "COFe": 
        xg_O2_,A,B,C = eq_COFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses) # COFe system
        xg_CO2_, xg_CO_, xm_CO2_, Xm_t, Xg_t, Fe32, Fe3T, wm_CO2_, wm_CO_ = A
        diff,wt_g_O,wt_g_C = B
        wt_g, wt_O_, wt_C_ = C
        guessx = xg_O2_
    elif system == "HOFe": 
        xg_O2_,A,B,C = eq_HOFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses) # HOFe system
        xg_H2O_, xg_H2_, xm_H2O_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2O_, wm_H2_ = A
        diff,wt_g_O,wt_g_H = B
        wt_g, wt_O_, wt_H_ = C
        guessx = xg_O2_
    elif system == "HOFe_xenia": 
        xg_O2_,A,B,C = eq_HOFe_xenia(PT,bulk_wf,models,nr_step,nr_tol,guesses) # HOFe system
        xg_H2O_, xg_H2_, xm_H2Omol_, xm_OH_, xm_H2O_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2O_, wm_H2_ = A
        diff,wt_g_O,wt_g_H = B
        wt_g, wt_O_, wt_H_ = C
        guessx = xg_O2_
    elif system == "SOFe":
        xg_O2_,A,B,C = eq_SOFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses) # SOFe system
        xg_SO2_, xg_S2_, Xg_t, Fe32, Fe3T, wm_S_, wm_SO3_, S62, S6T, wm_ST_ = A
        diff,wt_g_O,wt_g_S = B
        wt_g, wt_O_, wt_S_ = C
        guessx = xg_O2_
        Xm_t,= "",
    elif system == "CHOFe": 
        if models.loc["COH_species","option"] == "no_H2_CO_CH4_melt":
            xg_O2_, xg_CO_, A,B,C = eq_CHOFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses) # CHOFe system
            xg_CO2_, xg_H2O_, xg_H2_, xg_CH4_, xm_H2O_, xm_CO2_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2O_, wm_CO2_ = A
            mbCO, mbCH, wt_g_O, wt_g_C, wt_g_H = B
            wt_g, wt_O_, wt_C_, wt_H_ = C
            guessx, guessy = xg_O2_,xg_CO_
            solve_species = "OC"
        elif models.loc["COH_species","option"] == "H2O-CO2 only":
            xg_CO2_, A,B,C = eq_CH(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses) # CHOFe system
            xg_CO2_, xg_H2O_, xm_H2O_, xm_CO2_, Xm_t, Xg_t, wm_H2O_, wm_CO2_ = A
            mbCH, wt_g_C, wt_g_H = B
            wt_g, wt_O_, wt_C_, wt_H_ = C
            guessx = xg_CO2_
            xg_O2_, xg_H2_, xg_CO_, xg_CH4_, wm_H2_, wm_CO_, wm_CH4_, wt_g_O = 0.,0.,0.,0.,0.,0.,0.,""
            Fe3T = melt_wf["Fe3FeT"]
            Fe32 = mg.overtotal2ratio(Fe3T)
        elif models.loc["COH_species","option"] == "yes_H2_CO_CH4_melt":
            xg_O2_, xg_CO_, A,B,C = eq_CHOFe_2(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses) # CHOFe system
            xg_CO2_, xg_H2O_, xg_H2_, xg_CH4_, xm_H2O_, xm_CO2_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_ = A
            mbCO, mbCH, wt_g_O, wt_g_C, wt_g_H = B
            wt_g, wt_O_, wt_C_, wt_H_ = C
            guessx, guessy = xg_O2_, xg_CO_
            solve_species = "OC"
    elif system == "SHOFe": 
        if models.loc["COH_species","option"] == "no_H2_CO_CH4_melt" and models.loc["H2S_m","option"] == "False":
            xg_O2_, xg_S2_, A,B,C, = eq_SHOFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses) # SHOFe system
            xg_SO2_, xg_H2O_, xg_H2_, xg_H2S_, xm_H2O_, wm_S_, wm_SO3_, Xm_t, Xg_t, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_ST_ = A
            mbSO, mbSH, wt_g_O, wt_g_S, wt_g_H = B
            wt_g, wt_O_, wt_S_, wt_H_ = C
            wm_H2_, wm_H2S_ = 0.,0.
        elif models.loc["COH_species","option"] == "yes_H2_CO_CH4_melt" and models.loc["H2S_m","option"] == "True":
            A,B,C,D = eq_SHOFe_2(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses)
            xg_O2_, xg_H2_, xg_S2_, xg_H2O_, xg_SO2_, xg_H2S_, Xg_t, xm_H2O_, wm_S_, wm_SO3_, Xm_t, Xm_t_ox, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_ST_, wm_H2S_, wm_H2_ = B
            mbX, mbZ, wt_g_O, wt_g_H, wt_g_S = C
            wt_g, wt_O_, wt_H_, wt_S_ = D
        elif models.loc["COH_species","option"] == "yes_H2_CO_CH4_melt" and models.loc["H2S_m","option"] == "False":
            print("not possible")
        elif models.loc["COH_species","option"] == "no_H2_CO_CH4_melt" and models.loc["H2S_m","option"] == "True": 
            print("not possible")
        guessx, guessy = xg_O2_, xg_S2_
        solve_species = "OS"
    elif system == "SCOFe":
        xg_O2_, xg_CO_, A,B,C, = eq_SCOFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses)
        xg_SO2_, xg_CO2_, xg_CO_, xg_OCS_, xg_S2_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Xg_t, Fe32, Fe3T, S62, S6T, wm_CO2_, wm_ST_, wm_CO_ = A
        mbCO, mbCS, wt_g_O, wt_g_S, wt_g_C = B
        wt_g, wt_O_, wt_S_, wt_C_= C
        guessx, guessy = xg_O2_, xg_CO_
    elif system == "COXFe":
        xg_O2_, xg_CO_, A,B,C, = eq_COXFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses)
        xg_X_, xg_CO2_, xg_CO_, xm_CO2_, wm_X_, Xm_t, Xg_t, Fe32, Fe3T, wm_CO2_, wm_CO_ = A
        mbCO, mbCX, wt_g_O, wt_g_X, wt_g_C = B
        wt_g, wt_O_, wt_X_, wt_C_= C
        guessx, guessy = xg_O2_, xg_CO_
    elif system == "HOXFe":
        xg_O2_, xg_H2_, A,B,C, = eq_HOXFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guessx,guessy)
        xg_X_, xg_H2O_, xg_H2_, xm_H2O_, wm_X_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2O_, wm_H2_ = A
        mbHO, mbHX, wt_g_O, wt_g_X, wt_g_H = B
        wt_g, wt_O_, wt_X_, wt_H_= C
        guessx, guessy = xg_O2_, xg_H2_
    elif system == "SCHOFe":
        if models.loc["COH_species","option"] == "H2O-CO2 only":
            print("change COH_species option to yes_H2_CO_CH4_melt or no_H2_CO_CH4_melt")
        elif models.loc["COH_species","option"] == "no_H2_CO_CH4_melt" and models.loc["H2S_m","option"] == "False":
            A,B,C,D = eq_SCHOFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses,solve_species) # SHOFe system
            if solve_species == "OCS":
                xg_O2_, xg_CO_, xg_S2_ = A
                xg_CO2_, xg_H2O_, xg_H2_, xg_CH4_, xg_SO2_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_ = B
                guessx, guessy, guessz, guessw = xg_O2_, xg_CO_, xg_S2_, xg_H2_
            elif solve_species == "OHS":
                xg_O2_, xg_H2_, xg_S2_ = A
                xg_CO_, xg_CO2_, xg_H2O_, xg_CH4_, xg_SO2_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_ = B
                guessx, guessy, guessz, guessw = xg_O2_, xg_H2_, xg_S2_, xg_CO_
            elif solve_species == "OCH":
                xg_O2_, xg_CO_, xg_H2_ = A
                xg_CO2_, xg_H2O_, xg_CH4_, xg_S2_, xg_SO2_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_ = B
                guessx, guessy, guessz, guessw = xg_O2_, xg_CO_, xg_H2_, xg_S2_
            wm_H2_, wm_CO_, wm_CH4_, wm_H2S_ = 0.,0.,0.,0.
        elif models.loc["COH_species","option"] == "yes_H2_CO_CH4_melt" and models.loc["H2S_m","option"] == "True":
            if models.loc["Hspeciation","option"] in ["none",'none+regular','none+ideal']:
                A,B,C,D, models, solve_species = eq_SCHOFe_2(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses,solve_species)
            elif models.loc["Hspeciation","option"] == "linear":
                A,B,C,D = eq_SCHOFe_3(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses,solve_species)
            xg_O2_, xg_H2_, xg_S2_, xg_H2O_, xg_CO_, xg_CO2_, xg_SO2_, xg_CH4_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Xm_t_ox, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_, wm_H2S_, wm_H2_, wm_CH4_, wm_CO_ = B
            if solve_species == "OCS":
                guessx, guessy, guessz, guessw = xg_O2_, xg_CO_, xg_S2_, xg_H2_
            elif solve_species == "OHS":
                guessx, guessy, guessz, guessw = xg_O2_, xg_H2_, xg_S2_, xg_CO_ 
            elif solve_species == "OCH":
                guessx, guessy, guessz, guessw = xg_O2_, xg_CO_, xg_H2_, xg_S2_
        elif models.loc["COH_species","option"] == "yes_H2_CO_CH4_melt" and models.loc["H2S_m","option"] == "False":
            print("not possible")
        mbX, mbY, mbZ, wt_g_O, wt_g_C, wt_g_H, wt_g_S = C
        wt_g, wt_O_, wt_C_, wt_H_, wt_S_ = D
    elif system == "CHOXFe":
        if models.loc["COH_species","option"] == "no_H2_CO_CH4_melt":
            print("Not working yet")
        elif models.loc["COH_species","option"] == "yes_H2_CO_CH4_melt" :
            A,B,C,D = eq_CHOXFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses,solve_species) # SHOFe system
            xg_O2_, xg_CO_, xg_X_ = A
            guessx, guessy, guessz = xg_O2_, xg_CO_, xg_X_
            xg_O2__, xg_H2_, xg_X__, xg_H2O_, xg_CO__, xg_CO2_, xg_CH4_, Xg_t, xm_H2O_, xm_CO2_, wm_X_, Xm_t, Xm_t_ox, Fe32, Fe3T, wm_H2O_, wm_CO2_, wm_H2_, wm_CH4_, wm_CO_ = B
        mbX, mbY, mbZ, wt_g_O, wt_g_C, wt_g_H, wt_g_X = C
        wt_g, wt_O_, wt_C_, wt_H_, wt_X_ = D
    elif system == "SCHOXFe":
        if models.loc["COH_species","option"] == "no_H2_CO_CH4_melt" and models.loc["H2S_m","option"] == "False":
            print("cannot do")
        elif models.loc["COH_species","option"] == "yes_H2_CO_CH4_melt" and models.loc["H2S_m","option"] == "True":
            xg_O2_,xg_CO_,xg_S2_,xg_X_,xg_H2_,xg_H2O_,xg_CO2_,xg_OCS_,xg_SO2_,xg_H2S_,xg_CH4_,Xg_t,wm_CO2_,wm_CH4_,wm_CO_,wm_H2O_,wm_H2_,wm_SO3_,wm_S_,wm_H2S_,wm_ST_,wm_X_,xm_H2O_,xm_CO2_,Xm_t,Fe32,Fe3T,S62,S6T,wt_g = eq_SCHOXFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses,solve_species) # SCHOXFe system
            Xm_t_ox,mbX, mbY, mbZ, wt_g_O, wt_g_C, wt_g_H, wt_g_S, mbW, wt_g_X,wt_O_, wt_C_, wt_H_, wt_S_, wt_X_ = "","","","","","","","","","","","","","",""
            
    melt_wf["CO2"] = wm_CO2_
    melt_wf["H2OT"] = wm_H2O_
    if system in ["HOFe","HOFe","CHOFe","CHOXFe","SHOFe","SCHOFe","SCHOXFe"]: # contains H
        wm_H2Omol_, wm_OH_ = mg.wm_H2Omol_OH(PT,melt_wf,models)
    else:
        wm_H2Omol_, wm_OH_ = 0., 0.
    if system in ["COFe","COXFe","CHOFe","CHOXFe","SCHOFe","SCHOXFe"]: # contains C
        wm_CO2carb_, wm_CO2mol_ = mg.wm_CO32_CO2mol(PT,melt_wf,models)
    else:
        wm_CO2carb_, wm_CO2mol_ = 0., 0.
    
    wm_S6p_ = (wm_SO3_/mdv.species.loc["SO3","M"])*mdv.species.loc["S","M"]

    xg = {"xg_O2":xg_O2_, "xg_H2":xg_H2_, "xg_S2":xg_S2_, "xg_H2O":xg_H2O_, "xg_CO":xg_CO_, "xg_CO2":xg_CO2_, "xg_SO2":xg_SO2_, "xg_CH4":xg_CH4_, "xg_H2S":xg_H2S_, "xg_OCS":xg_OCS_, "xg_X":xg_X_, "Xg_t":Xg_t}
    melt = {"xm_H2O":xm_H2O_, "xm_CO2":xm_CO2_, "Xm_t":Xm_t, "wm_H2O":wm_H2O_, "wm_H2Omol":wm_H2Omol_, "wm_OH":wm_OH_,"wm_H2":wm_H2_, "wm_CO2":wm_CO2_, "wm_CO2carb":wm_CO2carb_, "wm_CO2mol":wm_CO2mol_, "wm_CO":wm_CO_, "wm_CH4":wm_CH4_, "wm_S":wm_S_,"wm_S2m":wm_S_, "wm_SO3":wm_SO3_,"wm_S6p":wm_S6p_, "wm_H2S":wm_H2S_, "wm_ST":wm_ST_, "wm_X":wm_X_, "Fe32":Fe32, "Fe3T":Fe3T, "S62":S62, "S6T":S6T}
    melt_and_gas = {"wt_g_O":wt_g_O, "wt_g_C":wt_g_C, "wt_g_H":wt_g_H, "wt_g_S":wt_g_S, "wt_g_X":wt_g_X, "wt_g":wt_g, "wt_O":wt_O_, "wt_C":wt_C_, "wt_H":wt_H_, "wt_S":wt_S_, "wt_X":wt_X_, 'wt_Fe':bulk_wf['Fe']}

    # chech mass balance
    mass_balance = c.check_mass_balance(xg,melt,melt_and_gas)

    guesses["guessx"] = guessx
    guesses["guessy"] = guessy
    guesses["guessz"] = guessz
    guesses["guessw"] = guessw

    return xg, melt, melt_and_gas, guesses, models, solve_species, mass_balance


#################################################
### specitation of H and C at given P and fO2 ###
#################################################

def eq_C_melt(PT,melt_wf,models): # equilibrium partitioning of C in the melt in CO system
    
    wt_C = melt_wf['CT'] # weight fraction
    K1 = mdv.KCOg(PT,models) 
    K2 = mdv.C_CO3(PT,melt_wf,models) # mole fraction
    K3 = (mdv.C_CO(PT,melt_wf,models))/1000000. # weight fraction
    M_C = mdv.species.loc['C','M']
    M_CO = mdv.species.loc['CO','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_m_ = mg.M_m_SO(melt_wf)
    fO2 = mdv.f_O2(PT,melt_wf,models)
    a = K3*M_C*(M_CO2 - M_m_)
    b = K3*M_m_*M_C + K1*K2*fO2**0.5*M_CO*(M_C - M_CO2*wt_C + M_m_*wt_C)
    c = -wt_C*K1*K2*fO2**0.5*M_CO*M_m_
    xm_CO2_ = (-b + (b**2 - 4.*a*c)**0.5)/(2.*a) # mole fraction CO32-
    wm_CO2_ = (xm_CO2_*M_CO2)/((xm_CO2_*M_CO2) + ((1.-xm_CO2_)*M_m_)) # weight fraction CO32- as CO2
    wm_CO_ = ((wt_C - (wm_CO2_/M_CO2)*M_C)/M_C)*M_CO # weight fraction CO
    check = M_C*((wm_CO2_/M_CO2)+(wm_CO_/M_CO))
    result = {"xm_CO2":xm_CO2_, "wm_CO2":wm_CO2_, "wm_CO":wm_CO_}
    return result

def eq_H_melt(PT,melt_wf,models,nr_step,nr_tol): # equilibrium partitioning of H in the melt in HO system
    if models.loc["Hspeciation","option"] != "none":
        print("not possible yet")
    wt_H = melt_wf['HT'] # weight fraction
    K1 = mdv.KHOg(PT,models) 
    K2 = mdv.C_H2O(PT,melt_wf,models) # mole fraction
    K3 = mdv.C_H2(PT,melt_wf,models) # ppm
    M_H = mdv.species.loc['H','M']
    M_H2 = mdv.species.loc['H2','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_m_ = mg.M_m_SO(melt_wf)
    fO2 = mdv.f_O2(PT,melt_wf,models)
    X = ((K1*K2)/K3)*fO2**0.5 # xH2O**2/wH2 with H2O in mole fraction and H2 in ppm
    Y = X/(X+1.) # xH2O**2/(xH2O**2+wH2) with H2O in mole fration and H2 in ppm 
    a = Y*M_H2O - Y*M_m_ - M_H2O + M_m_
    b = (Y - 1.)*M_m_
    c = 1000000.*Y*(wt_H*M_H2O - wt_H*M_m_ - 2.*M_H)
    d = 1000000.*Y*wt_H*M_m_
    constants = [a, b, c, d]
    def f(x, constants):
        return a*x**3 + b*x**2 + c*x + d
    def df(x, constants):
        return 3.*a*x**2 + 2.*b*x + c
    def dx(x, constants):
        f_ = f(x, constants)
        result =(abs(0-f_))
        return result
    x0 = mg.xm_H2OT_so(melt_wf)    
    delta1 = dx(x0, constants)
    while delta1 > nr_tol:
        f_ = f(x0, constants)
        df_ = df(x0, constants)
        x0 = x0 - nr_step*(f_/df_)
        delta1 = dx(x0, constants)
    xm_H2O_ = x0        
    wm_H2O_ = (xm_H2O_*M_H2O)/(xm_H2O_*M_H2O + (1.-xm_H2O_)*M_m_) # weight fraction H2O
    wm_H2_ = wt_H - ((wm_H2O_/M_H2O)*(2*M_H)) # weight fraction H2
    result = {"xm_H2O":xm_H2O_, "wm_H2O":wm_H2O_, "wm_H2":wm_H2_}
    return result

def eq_S_melt(PT,melt_wf,models):
    S6p_ST = mg.S6ST(PT,melt_wf,models)
    S2m_ST = 1. - S6p_ST
    wm_S2m_ = S2m_ST*melt_wf['ST']
    wm_S6p_ = S6p_ST*melt_wf['ST']
    result = {"wm_S2m":wm_S2m_, "wm_S6p":wm_S6p_}
    return result
    
def eq_CH_melt(PT,melt_wf,models,nr_step,nr_tol,guesses):
    P = PT["P"]
    wt_C = melt_wf['CT']
    wt_H = melt_wf['HT']
    fO2 = mdv.f_O2(PT,melt_wf,models)
    guessx = guesses["guessx"]
    guessy = guesses["guessy"]
    
    # equilibrium constants
    K1_ = mdv.KHOg(PT,models)
    K2_ = mdv.KCOg(PT,models)
    K3_ = mdv.KCOHg(PT,models)
    K4_ = mdv.C_H2O(PT,melt_wf,models) # mole fraction
    K5_ = mdv.C_CO3(PT,melt_wf,models) # mole fraction
    K6_ = mdv.C_H2(PT,melt_wf,models) # ppm
    K7_ = mdv.C_CO(PT,melt_wf,models) # ppm
    K8_ = mdv.C_CH4(PT,melt_wf,models) # ppm
   
    # molecular masses
    M_H = mdv.species.loc['H','M']
    M_C = mdv.species.loc['C','M']
    M_CO = mdv.species.loc['CO','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_H2 = mdv.species.loc['H2','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_CH4 = mdv.species.loc['CH4','M']
    M_m_ = mg.M_m_SO(melt_wf)
    
    constants = [wt_C, wt_H, K1_, K2_, K3_, K4_, K5_, K6_, K7_, K8_, M_C, M_H, M_CO, M_CO2, M_H2, M_H2O, M_CH4, M_m_, fO2]
    
    def mg_CH(xm_CO2_,xm_H2O_):
        Xm_t = xm_CO2_*M_CO2 + xm_H2O_*M_H2O + (1.0-xm_CO2_-xm_H2O_)*M_m_
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t # weight fraction
        wm_CO2_ = (xm_CO2_*M_CO2)/Xm_t # weight fraction
        fH2O = (xm_H2O_**2.)/K4_
        fCO2 = xm_CO2_/K5_
        fH2 = fH2O/(K1_*fO2**0.5)
        wm_H2_ = (fH2*K6_)/1000000. # weight fraction
        fCO = fCO2/(K2_*fO2**0.5)
        wm_CO_ = (fCO*K7_)/1000000. # weight fraction
        fCH4 = (fCO2*fH2O**2.)/(K3_*fO2**2.)
        wm_CH4_ = (fCH4*K8_)/1000000. # weight fraction
        return Xm_t, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_
    
    def f_CH(xm_CO2_,xm_H2O_):
        Xm_t, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_ = mg_CH(xm_CO2_,xm_H2O_)
        wt_m_C = M_C*((wm_CO2_/M_CO2) + (wm_CO_/M_CO) + (wm_CH4_/M_CH4))
        wt_m_H = M_H*(2.*(wm_H2O_/M_H2O) + 2.*(wm_H2_/M_H2) + 4.*(wm_CH4_/M_CH4))
        mbC = (wt_C - wt_m_C)
        mbH = (wt_H - wt_m_H)
        return mbC, mbH, wt_m_C, wt_m_H, 0.
    
    def df_CH(xm_CO2_,xm_H2O_,constants):
        dmbC_C = -M_C*(xm_CO2_*(-M_CO2 + M_m_)/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0))**2 + 1/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0)) + 1.0e-6*K8_*fO2**(-2.0)*(xm_H2O_**2.0/K4_)**2.0/(K3_*K5_*M_CH4) + 1.0e-6*K7_*fO2**(-0.5)/(K2_*K5_*M_CO))
        dmbC_H = -M_C*(xm_CO2_*(-M_H2O + M_m_)/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0))**2 + 4.0e-6*K8_*fO2**(-2.0)*xm_CO2_*xm_H2O_**(-1.0)*(xm_H2O_**2.0/K4_)**2.0/(K3_*K5_*M_CH4))
        dmbH_C = -M_H*(2.0*xm_H2O_*(-M_CO2 + M_m_)/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0))**2 + 4.0e-6*K8_*fO2**(-2.0)*(xm_H2O_**2.0/K4_)**2.0/(K3_*K5_*M_CH4))
        dmbH_H = -M_H*(2.0*xm_H2O_*(-M_H2O + M_m_)/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0))**2 + 2.0/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0)) + 1.6e-5*K8_*fO2**(-2.0)*xm_CO2_*xm_H2O_**(-1.0)*(xm_H2O_**2.0/K4_)**2.0/(K3_*K5_*M_CH4) + 4.0e-6*K6_*fO2**(-0.5)*xm_H2O_**1.0/(K1_*K4_*M_H2))        
        return dmbC_C, dmbC_H, dmbH_C, dmbH_H
    
    xm_CO2_,xm_H2O_ = jac_newton(guessx,guessy,constants,f_CH,df_CH,nr_step,nr_tol)
    results1 = mg_CH(xm_CO2_,xm_H2O_)
    results2 = f_CH(xm_CO2_,xm_H2O_)
    return xm_CO2_,xm_H2O_, results1, results2

def eq_HS_melt(PT,melt_wf,models,nr_step,nr_tol): # not sure this is right?
    wt_S = melt_wf['ST']
    wt_H = melt_wf['HT']
    fO2 = mdv.f_O2(PT,melt_wf,models)
    
    # equilibrium constants
    K1_ = mdv.KHOSg(PT,models)
    K2_ = mdv.C_H2O(PT,melt_wf,models) # mole fraction
    K3_ = mdv.C_S(PT,melt_wf,models)/1000000. # weight fraction
    K4_ = mdv.C_SO4(PT,melt_wf,models)/1000000. # weight fraction
    K5_ = mdv.C_H2S(PT,melt_wf,models)/1000000. # weight fraction
    K6_ = mdv.KHOg(PT,models) 
    K7_ = mdv.C_H2(PT,melt_wf,models)/1000000. # weight fraction
   
    # molecular masses
    M_H = mdv.species.loc['H','M']
    M_S = mdv.species.loc['S','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_H2S = mdv.species.loc['H2S','M']
    M_m_ = mg.M_m_SO(melt_wf)
    
    #constants = [wt_C, wt_H, K1_, K2_, K3_, K4_, K5_, K6_, K7_, K8_, M_C, M_H, M_CO, M_CO2, M_H2, M_H2O, M_CH4, M_m_, fO2]
    
    def dx(x):
        f_ = f(x)
        result =(abs(0-f_))
        return result

    def mg_HS(xm_H2O_):
        Xm_t = xm_H2O_*M_H2O + (1.0-xm_H2O_)*M_m_
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t # weight fraction
        fH2O = (xm_H2O_**2.)/K2_
        fH2 = fH2O/(K6_*fO2**0.5)
        if models.loc['COH_species','option'] == 'yes_H2_CO_CH4_melt':
            wm_H2_ = fH2*K7_ # weight fraction
        else:
            wm_H2_ = 0.
        wm_S2m_ = wt_S/(1 + (K4_*fO2**2)/K3_ + (K1_*K5_*M_S*xm_H2O_**2)/(K2_*K3_*M_H2S)) # weight fraction
        if models.loc['H2S_m','option'] == 'True':
            wm_H2S_ = (K1_*K5_*xm_H2O_**2*wm_S2m_)/(K2_*K3_) # weight fraction
        else:
            wm_H2S_ = 0.
        wm_S6p_ = (K4_*fO2**2*wm_S2m_)/K3_ # weight fraction
        return Xm_t, wm_H2O_, wm_H2_, wm_S2m_, wm_S6p_, wm_H2S_
    
    def f(xm_H2O_):
        Xm_t, wm_H2O_, wm_H2_, wm_S2m_, wm_S6p_, wm_H2S_ = mg_HS(xm_H2O_)
        result = wm_H2_ + ((2.*M_H*wm_H2S_)/M_H2S) + (2.*M_H*xm_H2O_)/(M_H2O*xm_H2O_ + (1 - xm_H2O_)*M_m_) - wt_H
        return result
    
    def df(xm_H2O_):
        result = -4*K1_**2*K5_**2*M_H*M_S*wt_S*xm_H2O_**3/(K2_**2*K3_**2*M_H2S**2*(K1_*K5_*M_S*xm_H2O_**2/(K2_*K3_*M_H2S) + 1 + K4_*fO2**2/K3_)**2) + 4*K1_*K5_*M_H*wt_S*xm_H2O_/(K2_*K3_*M_H2S*(K1_*K5_*M_S*xm_H2O_**2/(K2_*K3_*M_H2S) + 1 + K4_*fO2**2/K3_)) + 2*M_H*xm_H2O_*(-M_H2O + M_m_)/(M_H2O*xm_H2O_ + M_m_*(1 - xm_H2O_))**2 + 2*M_H/(M_H2O*xm_H2O_ + M_m_*(1 - xm_H2O_)) + 2*K7_*fO2**(-0.5)*xm_H2O_/(K2_*K6_) 
        return result
    
    x0 = mg.xm_H2OT_so(melt_wf)    
    delta1 = dx(x0)
    while delta1 > nr_tol:
        f_ = f(x0)
        df_ = df(x0)
        x0 = x0 - nr_step*(f_/df_)
        delta1 = dx(x0)
    xm_H2O_ = x0     
    Xm_t, wm_H2O_, wm_H2_, wm_S2m_, wm_S6p_, wm_H2S_ = mg_HS(xm_H2O_)
    return xm_H2O_, wm_H2O_, wm_H2_, wm_S2m_, wm_S6p_, wm_H2S_ 

def eq_CHS_melt(PT,melt_wf,models,nr_step,nr_tol,guesses):
    wt_S = melt_wf['ST']
    wt_H = melt_wf['HT']
    wt_C = melt_wf['CT']
    fO2 = mdv.f_O2(PT,melt_wf,models)
    guessx = guesses["guessx"]
    guessy = guesses["guessy"]
    guessz = guesses["guessz"]
    
    # equilibrium constants
    K1_ = mdv.KHOSg(PT,models)
    K2_ = mdv.C_H2O(PT,melt_wf,models) # mole fraction
    K3_ = mdv.C_S(PT,melt_wf,models)/1000000. # weight fraction
    K4_ = mdv.C_SO4(PT,melt_wf,models)/1000000. # weight fraction
    K5_ = mdv.C_H2S(PT,melt_wf,models)/1000000. # weight fraction
    K6_ = mdv.KHOg(PT,models) 
    K7_ = mdv.C_H2(PT,melt_wf,models)/1000000. # weight fraction
    K8_ = mdv.KHOg(PT,models)
    K9_ = mdv.KCOg(PT,models)
    K10_ = mdv.KCOHg(PT,models)
    K11_ = mdv.C_CO3(PT,melt_wf,models) # mole fraction
    K12_ = mdv.C_CO(PT,melt_wf,models)/1000000. # weight fraction
    K13_ = mdv.C_CH4(PT,melt_wf,models)/1000000. # weight fraction
   
    # molecular masses
    M_H = mdv.species.loc['H','M']
    M_S = mdv.species.loc['S','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_H2S = mdv.species.loc['H2S','M']
    M_m_ = mg.M_m_SO(melt_wf)
    M_C = mdv.species.loc['C','M']
    M_CO = mdv.species.loc['CO','M']
    M_H2 = mdv.species.loc['H2','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_CH4 = mdv.species.loc['CH4','M']

    constants = [wt_C, wt_H, wt_S, K1_, K2_, K3_, K4_, K5_, K6_, K7_, K8_, K9_, K10_, K11_, K12_, K13_, M_C, M_H, M_S, M_CO, M_CO2, M_H2, M_H2O, M_CH4, M_m_, fO2]
    
    def dx(x):
        f_ = f(x)
        result =(abs(0-f_))
        return result

    def mg_CHS(xm_CO2_,xm_H2O_,wm_S2m_):
        Xm_t = xm_CO2_*M_CO2 + xm_H2O_*M_H2O + (1.0-xm_CO2_-xm_H2O_)*M_m_
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t # weight fraction
        wm_CO2_ = (xm_CO2_*M_CO2)/Xm_t # weight fraction
        if models.loc["Hspeciation","option"] in ["none",'none+regular','none+ideal']:
            fH2O = (xm_H2O_**2.)/K2_
        elif models.loc["Hspeciation","option"] == "linear":
            fH2O = xm_H2O_/K2_
        fCO2 = xm_CO2_/K11_
        if models.loc["COH_species","option"] == "yes_H2_CO_CH4_melt":
            fH2 = fH2O/(K6_*fO2**0.5)
            wm_H2_ = fH2*K7_ # weight fraction
            fCO = fCO2/(K9_*fO2**0.5)
            wm_CO_ = fCO*K12_ # weight fraction
            fCH4 = (fCO2*fH2O**2.)/(K10_*fO2**2.)
            wm_CH4_ = fCH4*K13_ # weight fraction
        else:
            wm_H2_, wm_CO_,wm_CH4 = 0.,0.,0.
        fS2 = (fO2*wm_S2m_**2.)/K3_**2.
        if models.loc["H2S_m","option"] == "True":
            fH2S = (K1_*fS2**0.5*fH2O)/fO2**0.5
            wm_H2S_ = fH2S*K5_ # weight fraction
        else:
            wm_H2S_ = 0.
        wm_S6p_ = K4_*fS2**0.5*fO2**1.5 # weight fraction
        return Xm_t, wm_H2O_, wm_H2_, wm_CO2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_
    
    def f_CHS(xm_CO2_,xm_H2O_,wm_S2m_):
        Xm_t, wm_H2O_, wm_H2_, wm_CO2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_ = mg_CHS(xm_CO2_,xm_H2O_,wm_S2m_)
        wt_m_C = M_C*((wm_CO2_/M_CO2) + (wm_CO_/M_CO) + (wm_CH4_/M_CH4))
        wt_m_H = M_H*(2.*(wm_H2O_/M_H2O) + 2.*(wm_H2_/M_H2) + 4.*(wm_CH4_/M_CH4) + 2.*(wm_H2S_/M_H2S))
        wt_m_S = M_S*((wm_S2m_/M_S) + (wm_S6p_/M_S) + (wm_H2S_/M_H2S))
        wt_m_O = "na"
        mbC = (wt_C - wt_m_C)
        mbH = (wt_H - wt_m_H)
        mbS = (wt_S - wt_m_S)
        return mbC, mbH, mbS, wt_m_C, wt_m_H, wt_m_S, wt_m_O

    def df_CHS(xm_CO2_,xm_H2O_,wm_S2m_,constants):
        if models.loc["Hspeciation","option"] in ["none",'none+regular','none+ideal']:
            dmbC_C = -M_C*(xm_CO2_*(-M_CO2 + M_m_)/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0))**2 + 1/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0)) + K12_*fO2**(-0.5)/(K11_*K9_*M_CO) + K13_*fO2**(-2.0)*(xm_H2O_**2.0/K2_)**2.0/(K10_*K11_*M_CH4))
            dmbC_H = -M_C*(xm_CO2_*(-M_H2O + M_m_)/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0))**2 + 4.0*K13_*fO2**(-2.0)*xm_CO2_*xm_H2O_**(-1.0)*(xm_H2O_**2.0/K2_)**2.0/(K10_*K11_*M_CH4))
            dmbC_S = 0.
            dmbH_C = -M_H*(2.0*xm_H2O_*(-M_CO2 + M_m_)/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0))**2 + 4.0*K13_*fO2**(-2.0)*(xm_H2O_**2.0/K2_)**2.0/(K10_*K11_*M_CH4))
            dmbH_H = -M_H*(4.0*K1_*K5_*fO2**(-0.5)*xm_H2O_**1.0*(K3_**(-2.0)*fO2*wm_S2m_**2.0)**0.5/(K2_*M_H2S) + 2.0*xm_H2O_*(-M_H2O + M_m_)/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0))**2 + 2.0/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0)) + 4.0*K7_*fO2**(-0.5)*xm_H2O_**1.0/(K2_*K6_*M_H2) + 16.0*K13_*fO2**(-2.0)*xm_CO2_*xm_H2O_**(-1.0)*(xm_H2O_**2.0/K2_)**2.0/(K10_*K11_*M_CH4))
            dmbH_S = -2.0*K1_*K5_*M_H*fO2**(-0.5)*wm_S2m_**(-1.0)*xm_H2O_**2.0*(K3_**(-2.0)*fO2*wm_S2m_**2.0)**0.5/(K2_*M_H2S)
            dmbS_C = 0.
            dmbS_H = -2.0*K1_*K5_*M_S*fO2**(-0.5)*xm_H2O_**1.0*(K3_**(-2.0)*fO2*wm_S2m_**2.0)**0.5/(K2_*M_H2S)
            dmbS_S = -M_S*(1.0*K1_*K5_*fO2**(-0.5)*wm_S2m_**(-1.0)*xm_H2O_**2.0*(K3_**(-2.0)*fO2*wm_S2m_**2.0)**0.5/(K2_*M_H2S) + 1.0*K4_*fO2**1.5*wm_S2m_**(-1.0)*(K3_**(-2.0)*fO2*wm_S2m_**2.0)**0.5/M_S + 1/M_S)
        elif models.loc["Hspeciation","option"] == "linear":
            dmbC_C = -M_C*(xm_CO2_*(-M_CO2 + M_m_)/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0))**2 + 1/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0)) + K12_*fO2**(-0.5)/(K11_*K9_*M_CO) + K13_*fO2**(-2.0)*(xm_H2O_/K2_)**2.0/(K10_*K11_*M_CH4))
            dmbC_H = -M_C*(xm_CO2_*(-M_H2O + M_m_)/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0))**2 + 2.0*K13_*fO2**(-2.0)*xm_CO2_*(xm_H2O_/K2_)**2.0/(K10_*K11_*M_CH4*xm_H2O_))
            dmbC_S = 0.
            dmbH_C = -M_H*(2.0*xm_H2O_*(-M_CO2 + M_m_)/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0))**2 + 4.0*K13_*fO2**(-2.0)*(xm_H2O_/K2_)**2.0/(K10_*K11_*M_CH4))
            dmbH_H = -M_H*(2.0*K1_*K5_*fO2**(-0.5)*(K3_**(-2.0)*fO2*wm_S2m_**2.0)**0.5/(K2_*M_H2S) + 2.0*xm_H2O_*(-M_H2O + M_m_)/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0))**2 + 2.0/(M_CO2*xm_CO2_ + M_H2O*xm_H2O_ + M_m_*(-xm_CO2_ - xm_H2O_ + 1.0)) + 2.0*K7_*fO2**(-0.5)/(K2_*K6_*M_H2) + 8.0*K13_*fO2**(-2.0)*xm_CO2_*(xm_H2O_/K2_)**2.0/(K10_*K11_*M_CH4*xm_H2O_))
            dmbH_S = -2.0*K1_*K5_*M_H*fO2**(-0.5)*wm_S2m_**(-1.0)*xm_H2O_*(K3_**(-2.0)*fO2*wm_S2m_**2.0)**0.5/(K2_*M_H2S)
            dmbS_C = 0.
            dmbS_H = -K1_*K5_*M_S*fO2**(-0.5)*(K3_**(-2.0)*fO2*wm_S2m_**2.0)**0.5/(K2_*M_H2S)
            dmbS_S = -M_S*(1.0*K1_*K5_*fO2**(-0.5)*wm_S2m_**(-1.0)*xm_H2O_*(K3_**(-2.0)*fO2*wm_S2m_**2.0)**0.5/(K2_*M_H2S) + 1.0*K4_*fO2**1.5*wm_S2m_**(-1.0)*(K3_**(-2.0)*fO2*wm_S2m_**2.0)**0.5/M_S + 1/M_S)
        return dmbC_C, dmbC_H, dmbC_S, dmbH_C, dmbH_H, dmbH_S, dmbS_C, dmbS_H, dmbS_S
    
    xm_CO2_,xm_H2O_,wm_S2m_ = jac_newton3(guessx,guessy,guessz,constants,f_CHS,df_CHS,nr_step,nr_tol)
    results1 = mg_CHS(xm_CO2_,xm_H2O_,wm_S2m_)
    results2 = f_CHS(xm_CO2_,xm_H2O_,wm_S2m_)
    return xm_CO2_,xm_H2O_,wm_S2m_, results1, results2

def melt_speciation(PT,melt_wf,models,nr_step,nr_tol):
    system = set_system(melt_wf,models)
    wt_C = melt_wf['CT']
    wt_H = melt_wf['HT']
    wt_S = melt_wf['ST']
    M_H = mdv.species.loc['H','M']
    M_S = mdv.species.loc['S','M']
    M_C = mdv.species.loc['C','M']
    M_CO = mdv.species.loc['CO','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_H2 = mdv.species.loc['H2','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_CH4 = mdv.species.loc['CH4','M']
    M_H2S = mdv.species.loc['H2S','M']
    
    if system in ["COFe","COXFe","SOFe","SCOFe"]: # no H
        wm_CH4_, wm_H2S_, xm_H2O_, wm_H2O_, wm_H2_, wm_H2Omol_, wm_OH_  = 0.,0.,0.,0.,0.,0.,0.
    if system in ["HOFe","HOXFe","SOFe","SHOFe"]: # no C
        wm_CH4_, xm_CO2_, wm_CO2_, wm_CO_, wm_CO2carb_, wm_CO2mol_ = 0.,0.,0.,0.,0.,0.
    if system in ["COFe","COXFe","HOFe","HOXFe","CHOFe","CHOXFe"]: # no S
        wm_H2S_, wm_S2m_, wm_S6p_ = 0.,0.,0.
    if system in ["COFe","SOFe","SCOFe","HOFe","CHOFe","SHOFe"]: # no X
        wm_X_ = 0.
    else:
        wm_X_ = melt_wf['XT']
    
    if models.loc['COH_species','option'] == 'no_H2_CO_CH4_melt':
        wm_CO_, wm_CH4_, wm_H2_ = 0., 0., 0.
    if models.loc['H2S_m','option'] == 'False':
        wm_H2S_ = 0.
    
    if system in ["HOFe","HOXFe"]:
        if models.loc['COH_species','option'] == 'yes_H2_CO_CH4_melt':
            result = eq_H_melt(PT,melt_wf,models,nr_step,nr_tol)
            xm_H2O_, wm_H2O_, wm_H2_ = result["xm_H2O"], result["wm_H2O"], result["wm_H2"]
        else:
            wm_H2O_, xm_H2O_ = melt_wf["H2OT"], mg.xm_H2O_so(melt_wf)
    if system in ["COFe","COXFe",'SCOFe']:
        if models.loc['COH_species','option'] == 'yes_H2_CO_CH4_melt':    
            result = eq_C_melt(PT,melt_wf,models)
            xm_CO2_, wm_CO2_, wm_CO_ = result['xm_CO2'], result['wm_CO2'], result['wm_CO']
        else:
            wm_CO2_, xm_CO2_ = melt_wf["CO2"], mg.xm_CO2_so(melt_wf)
    if system in ["SOFe",'SCOFe']:
        result = eq_S_melt(PT,melt_wf,models)
        wm_S2m_, wm_S6p_ = result['wm_S2m'],result['wm_S6p']
    if system in ["CHOFe","CHOXFe"]:
        if models.loc['COH_species','option'] == 'yes_H2_CO_CH4_melt':    
            guesses = {"guessx":mg.xm_CO2_so(melt_wf), "guessy":mg.xm_H2OT_so(melt_wf)}
            xm_CO2_,xm_H2O_, A, B = eq_CH_melt(PT,melt_wf,models,nr_step,nr_tol,guesses)
            Xm_t, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_ = A
        else:
            wm_H2O_, xm_H2O_ = melt_wf["H2OT"], mg.xm_H2O_so(melt_wf)
            wm_CO2_, xm_CO2_ = melt_wf["CO2"], mg.xm_CO2_so(melt_wf)
    if system == "SHOFe":
        xm_H2O_, wm_H2O_, wm_H2_, wm_S2m_, wm_S6p_, wm_H2S_ = eq_HS_melt(PT,melt_wf,models,nr_step,nr_tol)
    if system in ["SCHOFe","SCHOXFe"]:
        guesses = {"guessx":mg.xm_CO2_so(melt_wf),"guessy":mg.xm_H2OT_so(melt_wf),"guessz":wt_S}
        xm_CO2_,xm_H2O_,wm_S2m_, A, B = eq_CHS_melt(PT,melt_wf,models,nr_step,nr_tol,guesses)
        Xm_t, wm_H2O_, wm_H2_, wm_CO2_, wm_CO_, wm_CH4_, wm_S2m_, wm_S6p_, wm_H2S_ = A
    
    if system in ["HOFe","HOFe","CHOFe","CHOXFe","SHOFe","SCHOFe","SCHOXFe"]: # contains H
        wm_H2Omol_, wm_OH_ = mg.wm_H2Omol_OH(PT,melt_wf,models)
    if system in ["COFe","COXFe","SCOFe","CHOFe","CHOXFe","SCHOFe","SCHOXFe"]: # contains C
        wm_CO2carb_, wm_CO2mol_ = mg.wm_CO32_CO2mol(PT,melt_wf,models)
    
    wm_SO3_ = (wm_S6p_*mdv.species.loc["SO3","M"])/mdv.species.loc["S","M"]

    Fe3FeT = melt_wf["Fe3FeT"]

    conc = {"xm_H2O":xm_H2O_, "wm_H2O":wm_H2O_, "wm_H2Omol":wm_H2Omol_, "wm_OH":wm_OH_, "xm_CO2":xm_CO2_, "wm_CO2":wm_CO2_, "wm_CO2carb":wm_CO2carb_,"wm_CO2mol":wm_CO2mol_,"wm_H2":wm_H2_, "wm_CO":wm_CO_, "wm_CH4":wm_CH4_, "wm_H2S":wm_H2S_, "wm_S2m":wm_S2m_, "wm_S6p":wm_S6p_, "wm_SO3":wm_SO3_,"wm_ST": wt_S,'Fe3FeT':Fe3FeT,"wm_X":wm_X_,"Fe3T":Fe3FeT}
    return conc

def eq_SOFe_melt(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses): # equilibrium between S and Fe in the melt
    P = PT["P"]
    wt_O = bulk_wf['O']
    wt_S = bulk_wf['S']
    wt_Fe = bulk_wf['Fe']
    wt_H = bulk_wf['H']
    wt_C = bulk_wf['C']
    guess = guesses["guessx"]
    
    # equilibrium constants
    K2_ = mdv.C_S(PT,melt_wf,models)/1000000.
    K3_ = mdv.C_SO4(PT,melt_wf,models)/1000000.
    KD1_, KD2, y = mdv.FefO2_KC91_EqA_terms(PT,melt_wf,models)
    
    # molecular masses
    M_S = mdv.species.loc['S','M']
    M_O = mdv.species.loc['O','M']
    M_H = mdv.species.loc['H','M']
    M_C = mdv.species.loc['C','M']
    M_Fe = mdv.species.loc['Fe','M']
    M_O2 = mdv.species.loc['O2','M']
    M_FeO15 = mdv.species.loc['FeO1.5','M']
    M_FeO = mdv.species.loc['FeO','M']
    M_m_ = mg.M_m_ox(melt_wf,models)
    
    constants = [P, wt_O, wt_S, wt_Fe, wt_H, wt_C, K2_, K3_, KD1_, KD2, y, M_S, M_O, M_Fe, M_O2, M_FeO, M_FeO15, M_m_]
     
    def mg_SOFe(fO2):
        Fe32 = ((KD1_*fO2**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*(fO2**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*(fO2**(0.5*y)))
        Fe3T = Fe32/(1+Fe32)
        S62 = (K3_/K2_)*fO2**2.
        S6T = S62/(1+S62)
        return Fe32, Fe3T, S62, S6T
    
    def mb_SOFe(fO2):
        return "", '', ''
    
    def f_SOFe(fO2):
        Fe32, Fe3T, S62, S6T = mg_SOFe(fO2)
        mb = wt_O - (((wt_H/M_H)/2.)*M_O) - ((wt_C/M_C)*2.*M_O) - (M_O*((wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))) - (((S6T*wt_S)/M_S)*3.*M_O)
        return mb,0.,0.
    
    def df_SOFe(fO2,constants):
        result = -wt_Fe*(-0.75*KD1_**(2.0*y)*KD2*fO2**(0.5*y)*y*(1.0 - 2.0*y)*(KD1_*fO2**0.25 + 2.0*KD1_**(2.0*y)*KD2*fO2**(0.5*y)*y)/(fO2*(KD1_**(2.0*y)*KD2*fO2**(0.5*y)*(1.0 - 2.0*y) + 1.0)**2) + 1.5*(0.25*KD1_*fO2**(-0.75) + 1.0*KD1_**(2.0*y)*KD2*fO2**(0.5*y)*y**2/fO2)/(KD1_**(2.0*y)*KD2*fO2**(0.5*y)*(1.0 - 2.0*y) + 1.0))/(M_Fe*((KD1_*fO2**0.25 + 2.0*KD1_**(2.0*y)*KD2*fO2**(0.5*y)*y)/(KD1_**(2.0*y)*KD2*fO2**(0.5*y)*(1.0 - 2.0*y) + 1.0) + 1.0)) - wt_Fe*(1.5*(KD1_*fO2**0.25 + 2.0*KD1_**(2.0*y)*KD2*fO2**(0.5*y)*y)/(KD1_**(2.0*y)*KD2*fO2**(0.5*y)*(1.0 - 2.0*y) + 1.0) + 1.0)*(0.5*KD1_**(2.0*y)*KD2*fO2**(0.5*y)*y*(1.0 - 2.0*y)*(KD1_*fO2**0.25 + 2.0*KD1_**(2.0*y)*KD2*fO2**(0.5*y)*y)/(fO2*(KD1_**(2.0*y)*KD2*fO2**(0.5*y)*(1.0 - 2.0*y) + 1.0)**2) - (0.25*KD1_*fO2**(-0.75) + 1.0*KD1_**(2.0*y)*KD2*fO2**(0.5*y)*y**2/fO2)/(KD1_**(2.0*y)*KD2*fO2**(0.5*y)*(1.0 - 2.0*y) + 1.0))/(M_Fe*((KD1_*fO2**0.25 + 2.0*KD1_**(2.0*y)*KD2*fO2**(0.5*y)*y)/(KD1_**(2.0*y)*KD2*fO2**(0.5*y)*(1.0 - 2.0*y) + 1.0) + 1.0)**2) - 6.0*K3_*M_O*fO2**1.0*wt_S/(K2_*M_S*(1 + K3_*fO2**2.0/K2_)) + 6.0*K3_**2*M_O*fO2**3.0*wt_S/(K2_**2*M_S*(1 + K3_*fO2**2.0/K2_)**2)
        return result
    
    fO2_ = newton_raphson(guessx,constants,nr_tol,nr_step,f_SOFe,df_SOFe)
    result1 = mg_SOFe(fO2_)
    result2 = f_SOFe(fO2_)
    result3 = mb_SOFe(fO2_)
    #print(wt_S, wt_O, xg_O2_, result1, result2, result3)
    return fO2_, result1, result2, result3
    
##############
### solver ###
##############

def test_f(x,y):
    eq1 = x**2 + y**3 + 7*x*y
    eq2 = x + 2*y**2 + 3*x*y
    return eq1, eq2,0,0,0

def test_df(x,y):
    eq1_x = 2*x + 7*y
    eq1_y = 3*y**2 + 7*x
    eq2_x = 1 + 3*y
    eq2_y = 4*y + 3*x
    return eq1_x, eq1_y, eq2_x, eq2_y

def test3_f(x,y,z):
    eq1 = x**2 + y**3 + 7*x*y + x*z
    eq2 = x + 2*y**2 + 3*x*y + x*z**3
    eq3 = x + y + z + 3*z**2
    return eq1, eq2, eq3,0,0,0,0

def test3_df(x,y,z):
    eq1_x = 2*x + 7*y + z
    eq1_y = 3*y**2 + 7*x
    eq1_z = x
    eq2_x = 1 + 3*y + z**3
    eq2_y = 4*y + 3*x
    eq2_z = 3*x*z**2
    eq3_x = 1
    eq3_y = 1
    eq3_z = 1 + 3*z
    return eq1_x, eq1_y, eq1_z, eq2_x, eq2_y, eq2_z, eq3_x, eq3_y, eq3_z

def newton_raphson(x0,constants,e1,step,eqs,deriv):
    # create results table
    results = pd.DataFrame([["guessx","diff","step"]])  
    results.to_csv('results_newtraph.csv', index=False, header=False)
    
    def dx(x,eqs):
        f_,wtg1,wtg2 = eqs(x)
        result =(abs(0-f_))
        return result
    
    delta1 = dx(x0,eqs)
    
    results1 = pd.DataFrame([[x0,delta1,step]])
    n = 0.
    results = pd.concat([results, results1], ignore_index=True)

    while delta1 > e1:
        while x0 < 0.:
            n = n+1.
            results1 = pd.DataFrame([[x0,delta1,step]]) 
            results = pd.concat([results, results1], ignore_index=True)
            if(i % 50==0):
                results.to_csv('results_newtraph.csv', index=False, header=False)     
            step = step/10.
            x0 = x0 - step*(f_/df_)
        n = n+1.
        f_,wtg1,wtg2 = eqs(x0)
        df_ = deriv(x0,constants)
        x0 = x0 - step*(f_/df_)
        delta1 = dx(x0,eqs)
        results1 = pd.DataFrame([[x0,delta1,step]]) 
        results = pd.concat([results, results1], ignore_index=True)
        if(n % 50==0):
            results.to_csv('results_newtraph.csv', index=False, header=False)     
    return x0        

#jac_newton(1,1,test_f,test_df,1)

def jac_newton(x0,y0,constants,eqs,deriv,step,tol,maxiter=1000):

    # create results table
    results = pd.DataFrame([["guessx","guessy","diff1","diff2","step"]])  
    diff1, diff2, wtg1,wtg2,wtg3 = eqs(x0,y0)
    results1 = pd.DataFrame([[x0,y0,diff1,diff2,step]]) 
    results = pd.concat([results, results1], ignore_index=True)
    n = 0.

    def F(eqs,x,y):
        a = eqs(x,y)
        return np.array([a[0],a[1]])
    
    def x2jac(step,deriv,eqs,guessx,guessy):
        eq1_x, eq1_y, eq2_x, eq2_y = deriv
        Func = F(eqs,guessx,guessy)
        J = np.array([[eq1_x, eq1_y],[eq2_x, eq2_y]])
        det = J[0][0]*J[1][1] - J[0][1]*J[1][0]
        inv_J = (1/det)*np.array(([J[-1][-1],-(J[0][-1])],[-(J[-1][0]),J[0][0]]), dtype=object)
        #inv_J = np.linalg.inv(J)
        new_guess = np.array([guessx, guessy], dtype=object) - step*np.dot(inv_J, Func)
        #new_guess = np.array([guessx,guessy]) - step*np.dot(inv_J,Func)
        return new_guess[0], new_guess[-1], J
    
    x00, y00 = x0, y0
    step0 = step

    for iter in range(maxiter):
        n = n+1.
        deriv_ = deriv(x0,y0,constants)
        guessx, guessy, J = x2jac(step,deriv_,eqs,x0, y0)
        while guessx < 0.0 or guessy < 0.0:
            step = step/10.
            guessx, guessy, J = x2jac(step,deriv_,eqs,x0,y0)
        #print(guessx,guessy)
        #print(guessx,guessy,diff1,diff2,wtg1,wtg2,wtg3)
        #if (abs(x0 - guessx)/x0)*100 < tol and (abs(y0-guessy)/y0)*100 < tol: tol = 1e-5
        diff1, diff2, wtg1,wtg2,wtg3 = eqs(guessx,guessy)
        if abs(diff1) < tol and abs(diff2) < tol:
            return guessx, guessy
        #elif np.isnan(float(guessx)) or np.isnan(float(guessy)):
            #print("nan encountered")
        x0 = guessx
        y0 = guessy
        results1 = pd.DataFrame([[guessx, guessy,diff1,diff2,step]])
        results = pd.concat([results, results1], ignore_index=True)
        if(n % 50==0):
            results.to_csv('results_jacnewton2.csv', index=False, header=False) 

    #step = step0/10.
    #x0, y0 = x00, y00
    #for iter in range(maxiter):
    #    deriv_ = deriv(x0,y0,constants)
    #    guessx, guessy, J = x2jac(step,deriv_,eqs,x0,y0)
    #    while guessx < 0.0 or guessy < 0.0:
    #        step = step/10.
    #        guessx, guessy, J = x2jac(step,deriv_,eqs,x0,y0)
    #    diff1, diff2,wtg1,wtg2,wtg3 = eqs(guessx,guessy)
    #    if abs(diff1) < tol and abs(diff2) < tol:
    #        return guessx, guessy
    #    #elif np.isnan(float(guessx)) or np.isnan(float(guessy)):
    #        #print("nan encountered")
    #    x0 = guessx
    #    y0 = guessy
    #    results1 = pd.DataFrame([[guessx, guessy,diff1,diff2,step]])
    #    results = pd.concat([results, results1], ignore_index=True)
    #    results.to_csv('results_jacnewton2.csv', index=False, header=False)

    for iter in range(9):
        step = step0 - (step0/10.)
        x0, y0 = x00, y00
        for iter in range(maxiter):
            deriv_ = deriv(x0,y0,constants)
            guessx, guessy, J = x2jac(step,deriv_,eqs,x0,y0)
            while guessx < 0.0 or guessy < 0.0:
                step = step/10.
                guessx, guessy, J = x2jac(step,deriv_,eqs,x0,y0)
            diff1, diff2,wtg1,wtg2,wtg3 = eqs(guessx,guessy)
            if abs(diff1) < tol and abs(diff2) < tol:
                return guessx, guessy
            #elif np.isnan(float(guessx)) or np.isnan(float(guessy)):
                #print("nan encountered")
            x0 = guessx
            y0 = guessy
            results1 = pd.DataFrame([[guessx, guessy,diff1,diff2,step]])
            results = pd.concat([results, results1], ignore_index=True)
            results.to_csv('results_jacnewton2.csv', index=False, header=False)

    guessx,guessy = 1.,1.
    return guessx,guessy

def jac_newton3(x0,y0,z0,constants,eqs,deriv,step,tol,maxiter=50):
# create results table
    results = pd.DataFrame([["guessx","guessy","guessz","diff1","diff2","diff3","step"]])  
    diff1, diff2, diff3, wtg1,wtg2,wtg3,wtg4 = eqs(x0,y0,z0)
    results1 = pd.DataFrame([[x0,y0,z0,diff1,diff2,diff3,step]]) 
    results = pd.concat([results, results1], ignore_index=True)
    n = 0.

    def F(eqs,x,y,z):
        a = eqs(x,y,z)
        return np.array([a[0],a[1],a[2]])
    
    def x3jac(step,deriv,eqs,guessx,guessy,guessz,constants):
        eq1_x, eq1_y, eq1_z, eq2_x, eq2_y, eq2_z, eq3_x, eq3_y, eq3_z = deriv
        Func = F(eqs,guessx,guessy,guessz)
        J = np.array([[eq1_x, eq1_y, eq1_z],[eq2_x, eq2_y, eq2_z],[eq3_x, eq3_y, eq3_z]])
        m1, m2, m3, m4, m5, m6, m7, m8, m9 = J.ravel()
        determinant = m1*m5*m9 + m4*m8*m3 + m7*m2*m6 - m1*m6*m8 - m3*m5*m7 - m2*m4*m9
        inv_jac = np.array([[m5*m9-m6*m8, m3*m8-m2*m9, m2*m6-m3*m5],
                        [m6*m7-m4*m9, m1*m9-m3*m7, m3*m4-m1*m6],
                        [m4*m8-m5*m7, m2*m7-m1*m8, m1*m5-m2*m4]])/determinant
        dot = np.dot(inv_jac, Func)  # To get the arrays as 2 3x1 columns
        #inv_J = np.linalg.inv(J)
        new_guess = np.array([guessx, guessy, guessz]) - step*dot  # Returns a 3x1 array of guessx, guessy, guessz
        #new_guess = np.array([guessx,guessy]) - step*np.dot(inv_J,Func)
        return new_guess[0], new_guess[1], new_guess[2], J  # guessx, guessy, guessz, Jacobian
    
    x00, y00, z00 = x0, y0, z0
    step0 = step
        
    for iter in range(maxiter):
        n = n+1.
        deriv_ = deriv(x0,y0,z0,constants)
        guessx, guessy, guessz, J = x3jac(step,deriv_,eqs,x0,y0,z0,constants)
        if guessx < 0.0 or guessy < 0.0 or guessz < 0.0 or guessx > 1.0 or guessy > 1.0 or guessz > 1.0 and step:
            step = step/10.
            guessx, guessy, guessz, J = x3jac(step,deriv_,eqs,x0,y0,z0,constants)
        if guessx < 0.0 or guessy < 0.0 or guessz < 0.0 or guessx > 1.0 or guessy > 1.0 or guessz > 1.0 and step:
            break
        diff1, diff2, diff3, wtg1,wtg2,wtg3,wtg4 = eqs(guessx,guessy,guessz)
        if abs(diff1) < tol and abs(diff2) < tol and abs(diff3) < tol:
            return guessx, guessy, guessz
        #elif np.isnan(float(guessx)) or np.isnan(float(guessy)) or np.isnan(float(guessz)):
            #print("nan encountered")
        x0 = guessx
        y0 = guessy
        z0 = guessz
        results1 = pd.DataFrame([[guessx, guessy,guessz,diff1,diff2,diff3,step]])
        results = pd.concat([results, results1], ignore_index=True)
        if(n % 50==0):
            results.to_csv('results_jacnewton3.csv', index=False, header=False)

    #step = step0/10.
    #x0, y0, z0 = x00, y00, z00
    #for iter in range(maxiter):
    #    deriv_ = deriv(x0,y0,z0,constants)
    #    guessx, guessy, guessz, J = x3jac(step,deriv_,eqs,x0,y0,z0,constants)
    #    while guessx < 0.0 or guessy < 0.0 or guessz < 0.0:
    #        step = step/10.
    #        guessx, guessy, guessz, J = x3jac(step,deriv_,eqs,x0,y0,z0,constants)
    #    diff1, diff2, diff3, wtg1,wtg2,wtg3,wtg4 = eqs(guessx,guessy,guessz)
    #    if abs(diff1) < tol and abs(diff2) < tol and abs(diff3) < tol:
    #        return guessx, guessy, guessz
    #    #elif np.isnan(float(guessx)) or np.isnan(float(guessy)) or np.isnan(float(guessz)):
    #        #print("nan encountered")
    #    x0 = guessx
    #    y0 = guessy
    #    z0 = guessz
    #    results1 = pd.DataFrame([[guessx, guessy,guessz,diff1,diff2,diff3,step]])
    #    results = pd.concat([results, results1], ignore_index=True)
    #    results.to_csv('results_jacnewton3.csv', index=False, header=False)

    x0, y0, z0 = x00, y00, z00
    step = step0 - (step0/10.)
    for iter in range(9):
        for iter in range(maxiter):
            deriv_ = deriv(x0,y0,z0,constants)
            guessx, guessy, guessz, J = x3jac(step,deriv_,eqs,x0,y0,z0,constants)
            while guessx < 0.0 or guessy < 0.0 or guessz < 0.0 or guessx > 1.0 or guessy > 1.0 or guessz > 1.0:
                x0, y0, z0 = x00, y00, z00
                break
                #step = step/10.
                #guessx, guessy, guessz, J = x3jac(step,deriv_,eqs,x0,y0,z0,constants)
            diff1, diff2, diff3, wtg1,wtg2,wtg3,wtg4 = eqs(guessx,guessy,guessz)
            if abs(diff1) < tol and abs(diff2) < tol and abs(diff3) < tol:
                return guessx, guessy, guessz
            #elif np.isnan(float(guessx)) or np.isnan(float(guessy)) or np.isnan(float(guessz)):
                #print("nan encountered")
            x0 = guessx
            y0 = guessy
            z0 = guessz
            results1 = pd.DataFrame([[guessx, guessy,guessz,diff1,diff2,diff3,step]])
            results = pd.concat([results, results1], ignore_index=True)
            results.to_csv('results_jacnewton3.csv', index=False, header=False)
        step = step - (step0/10.)
    
    guessx,guessy,guessz = 1.,1.,1.
    return guessx,guessy,guessz
             
############        
### COFe ###
############

def eq_COFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses):
    P = PT["P"]
    wt_O = bulk_wf['O']
    wt_C = bulk_wf['C']
    wt_Fe = bulk_wf['Fe']
    guessx = guesses["guessx"]
    
    # equilibrium constants
    K1_ = mdv.KCOg(PT,models)
    K2_ = mdv.C_CO3(PT,melt_wf,models)
    K3_ = (mdv.C_CO(PT,melt_wf,models))/1000000.
    KD1_, KD2, y = mdv.FefO2_KC91_EqA_terms(PT,melt_wf,models)
    
    # fugacity coefficients
    y_CO2_ = mdv.y_CO2(PT,models)
    y_O2_ = mdv.y_O2(PT,models)
    y_CO_ = mdv.y_CO(PT,models)
    
    # molecular masses
    M_C = mdv.species.loc['C','M']
    M_O = mdv.species.loc['O','M']
    M_Fe = mdv.species.loc['Fe','M']
    M_O2 = mdv.species.loc['O2','M']
    M_CO = mdv.species.loc['CO','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_FeO = mdv.species.loc['FeO','M']
    M_FeO15 = mdv.species.loc['FeO1.5','M']
    M_m_ = mg.M_m_SO(melt_wf)
    
    constants = [P, wt_O, wt_C, wt_Fe, K1_, K2_, K3_, KD1_, KD2, y, y_CO2_, y_O2_, y_CO_, M_C, M_O, M_Fe, M_O2, M_CO, M_CO2, M_FeO, M_FeO15, M_m_]
     
    def mg_COFe(xg_O2_):
        xg_CO2_ = (1.0-xg_O2_)/(1.0 + (y_CO2_/(K1_*y_CO_*(y_O2_*xg_O2_*P)**0.5)))
        xg_CO_ = (y_CO2_*xg_CO2_)/(K1_*y_CO_*(y_O2_*xg_O2_*P)**0.5)
        xm_CO2_ = K2_*y_CO2_*xg_CO2_*P
        Xm_t = xm_CO2_*M_CO2 + (1.0-xm_CO2_)*M_m_
        Xg_t = xg_CO2_*M_CO2 + xg_CO_*M_CO + xg_O2_*M_O2
        Fe32 = ((KD1_*(y_O2_*xg_O2_*P)**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y)))
        Fe3T = Fe32/(1+Fe32)
        wm_CO2_ = (xm_CO2_*M_CO2)/Xm_t
        if models.loc["COH_species","option"] == "yes_H2_CO_CH4_melt":
            wm_CO_ = K3_*y_CO_*xg_CO_*P
        else:
            wm_CO_ = 0.
        return xg_CO2_, xg_CO_, xm_CO2_, Xm_t, Xg_t, Fe32, Fe3T, wm_CO2_, wm_CO_
    
    def mb_COFe(xg_O2_):
        xg_CO2_, xg_CO_, xm_CO2_, Xm_t, Xg_t, Fe32, Fe3T, wm_CO2_, wm_CO_ = mg_COFe(xg_O2_)
        diff,wt_g_O,wt_g_C = f_COFe(xg_O2_)
        wt_g = (wt_g_O+wt_g_C)/2.
        wt_C_ = M_C*((wt_g*(((xg_CO2_+xg_CO_)/Xg_t) - (xm_CO2_/Xm_t) - (wm_CO_/M_CO)))+(xm_CO2_/Xm_t) + (wm_CO_/M_CO))
        wt_O_ = M_O*((wt_g*(((2.0*xg_CO2_ + xg_CO_ + 2.0*xg_O2_)/Xg_t)-(2.0*xm_CO2_/Xm_t) - (wm_CO_/M_CO)))+(2.0*xm_CO2_/Xm_t) + (wm_CO_/M_CO) + (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))
        return wt_g, wt_O_, wt_C_
    
    def f_COFe(xg_O2_):
        xg_CO2_, xg_CO_, xm_CO2_, Xm_t, Xg_t, Fe32, Fe3T, wm_CO2_, wm_CO_ = mg_COFe(xg_O2_)
        wt_g_C = ((wt_C/M_C) - (xm_CO2_/Xm_t) - (wm_CO_/M_CO))/(((xg_CO2_+xg_CO_)/Xg_t) - (xm_CO2_/Xm_t) - (wm_CO_/M_CO))
        wt_g_O = ((wt_O/M_O) - (2.0*xm_CO2_/Xm_t) - (wm_CO_/M_CO) - (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))/(((2.0*xg_CO2_ + xg_CO_ + 2.0*xg_O2_)/Xg_t)-(2.0*xm_CO2_/Xm_t) - (wm_CO_/M_CO))
        diff = wt_g_C - wt_g_O
        return diff,wt_g_O,wt_g_C
    
    def df_COFe(xg_O2_,constants):
        result = de.COFe_O2(xg_O2_,constants)
        return result
    
    xg_O2_ = newton_raphson(guessx,constants,nr_tol,nr_step,f_COFe,df_COFe)
    result1 = mg_COFe(xg_O2_)
    result2 = f_COFe(xg_O2_)
    result3 = mb_COFe(xg_O2_)
    #print(wm_CO2,xg_O2,xg_CO,xg_CO2,F_,wt_g_C,wt_g_O)
    return xg_O2_, result1, result2, result3


############
### HOFe ###
############

def eq_HOFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses):
    P = PT["P"]
    wt_O = bulk_wf['O']
    wt_H = bulk_wf['H']
    wt_Fe = bulk_wf['Fe']
    guessx = guesses["guessx"]
    
    # equilibrium constants
    K1_ = mdv.KHOg(PT,models)
    K2_ = mdv.C_H2O(PT,melt_wf,models)
    K3_ = mdv.C_H2(PT,melt_wf,models)/1000000.
    KD1_, KD2, y = mdv.FefO2_KC91_EqA_terms(PT,melt_wf,models)
   
    # fugacity coefficients
    y_H2O_ = mdv.y_H2O(PT,models)
    y_O2_ = mdv.y_O2(PT,models)
    y_H2_ = mdv.y_H2(PT,models)
    
    # molecular masses
    M_H = mdv.species.loc['H','M']
    M_O = mdv.species.loc['O','M']
    M_Fe = mdv.species.loc['Fe','M']
    M_O2 = mdv.species.loc['O2','M']
    M_H2 = mdv.species.loc['H2','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_FeO = mdv.species.loc['FeO','M']
    M_FeO15 = mdv.species.loc['FeO1.5','M']
    M_m_ = mg.M_m_SO(melt_wf)
    
    constants = [P, wt_O, wt_H, wt_Fe, K1_, K2_, K3_, KD1_, KD2, y, y_H2O_, y_O2_, y_H2_, M_H, M_O, M_Fe, M_O2, M_H2, M_H2O, M_FeO, M_FeO15, M_m_]
     
    def mg_HOFe(xg_O2_):
        xg_H2O_ = (1.0-xg_O2_)/(1.0 + (y_H2O_/(K1_*y_H2_*y_O2_**0.5*xg_O2_**0.5*P**0.5)))
        xg_H2_ = (y_H2O_*xg_H2O_)/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)
        xm_H2O_ = (K2_*y_H2O_*xg_H2O_*P)**0.5
        Xm_t = xm_H2O_*M_H2O + (1.0-xm_H2O_)*M_m_
        Xg_t = xg_H2O_*M_H2O + xg_H2_*M_H2 + xg_O2_*M_O2
        Fe32 = ((KD1_*(y_O2_*xg_O2_*P)**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y)))
        Fe3T = Fe32/(1.0+Fe32)
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t
        if models.loc["COH_species","option"] == "yes_H2_CO_CH4_melt":
            wm_H2_ = K3_*y_H2_*xg_H2_*P
        else:
            wm_H2_ = 0.
        return xg_H2O_, xg_H2_, xm_H2O_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2O_, wm_H2_
    
    def mb_HOFe(xg_O2_):
        xg_H2O_, xg_H2_, xm_H2O_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2O_, wm_H2_ = mg_HOFe(xg_O2_)
        diff,wt_g_O,wt_g_H = f_HOFe(xg_O2_)
        wt_g = (wt_g_O+wt_g_H)/2
        wt_H_ = 2.0*M_H*((wt_g*(((xg_H2O_+xg_H2_)/Xg_t) - (xm_H2O_/Xm_t) - (wm_H2_/M_H2)))+(xm_H2O_/Xm_t) + (wm_H2_/M_H2))
        wt_O_ = M_O*((wt_g*(((xg_H2O_ + 2.0*xg_O2_)/Xg_t)-(xm_H2O_/Xm_t)))+(xm_H2O_/Xm_t) + (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))
        return wt_g, wt_O_, wt_H_
    
    def f_HOFe(xg_O2_):
        xg_H2O_, xg_H2_, xm_H2O_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2O_,wm_H2_ = mg_HOFe(xg_O2_)
        wt_g_H = ((wt_H/(2.0*M_H)) - (xm_H2O_/Xm_t) - (wm_H2_/M_H2))/(((xg_H2O_+xg_H2_)/Xg_t) - (xm_H2O_/Xm_t) - (wm_H2_/M_H2))
        wt_g_O = ((wt_O/M_O) - (xm_H2O_/Xm_t) - (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))/(((xg_H2O_ + 2.0*xg_O2_)/Xg_t)-(xm_H2O_/Xm_t))
        diff = wt_g_H - wt_g_O
        return diff,wt_g_O,wt_g_H
    
    def df_HOFe(xg_O2_,constants):
        if models.loc["COH_species","option"] == "no_H2_CO_CH4_melt":
            result = de.HOFe_O2(xg_O2_,constants) 
        else:
            result = de.HOFe_O2_2(xg_O2_,constants) 
        return result
    
    xg_O2_ = newton_raphson(guessx,constants,nr_tol,nr_step,f_HOFe,df_HOFe)
    result1 = mg_HOFe(xg_O2_)
    result2 = f_HOFe(xg_O2_)
    result3 = mb_HOFe(xg_O2_)
    return xg_O2_, result1, result2, result3


####################
### HOFe - Xenia ###
####################

def eq_HOFe_xenia(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses):
    P = PT["P"]
    wt_O = bulk_wf['O']
    wt_H = bulk_wf['H']
    wt_Fe = bulk_wf['Fe']
    guessx = guesses["guessx"]
    
    # equilibrium constants
    K1_ = 1.0/(6.9617e-10) # H2O = H2 + 0.5O2
    K2_ = exp(-13.869+3890.0/T_K)# H2Og = H2Omol Zhang 1999 eq14 @1bar
    K3_ = 0.368 # H2O + O = 2OH
    K4a_ = mdv.C_H2_a()/100. # H2g = H2m, converted to weight fraction
    K4b_ = mdv.C_H2_b() # power for H2g = H2m
    K5_ = mdv.PT_KCterm(PT) # PT terms in K&C91
   
    # fugacity coefficients
    y_H2O_ = 1.0 #y_H2O(P,T_K,ideal_gas)
    y_O2_ = 1.0 #y_O2(P,T_K,ideal_gas)
    y_H2_ = 1.0 #y_H2(P,T_K,ideal_gas)
    
    # molecular masses
    M_H = mdv.species.loc['H','M']
    M_O = mdv.species.loc['O','M']
    M_Fe = mdv.species.loc['Fe','M']
    M_O2 = mdv.species.loc['O2','M']
    M_H2 = mdv.species.loc['H2','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_FeO = mdv.species.loc['FeO','M']
    M_FeO15 = mdv.species.loc['FeO1.5','M']
    M_m_ = mg.M_m(melt_wf)
    
    constants = [P, wt_O, wt_H, wt_Fe, K1_, K2_, K3_, K4a_, K4b_, K5_, y_H2O_, y_O2_, y_H2_, M_H, M_O, M_Fe, M_O2, M_H2, M_H2O, M_FeO, M_FeO15, M_m_]
     
    def mg_HOFe(xg_O2_):
        xg_H2O_ = (1.0-xg_O2_)/(1.0 + (y_H2O_/(K1_*y_H2_*y_O2_**0.5*xg_O2_**0.5*P**0.5)))
        xg_H2_ = (y_H2O_*xg_H2O_)/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)
        xm_H2Omol_ = K2_*y_H2O_*xg_H2O_*P
        a = 4.0
        b = 2.0*xm_H2Omol_*(K3_ - 4.0)
        c = xm_H2Omol_*(xm_H2Omol_*(4.0 - K3_) - K3_)
        xm_H2OT_ = (-b + (b**2.0 - 4*a*c)**0.5)/(2.0*a)
        xm_OH_ = 2.0*(xm_H2OT_ - xm_H2Omol_)
        wm_H2_ = K4a_*(y_H2_*xg_H2_*P)**K4b_
        Xm_t = xm_H2OT_*M_H2O + (1.0-xm_H2OT_)*M_m_
        Xg_t = xg_H2O_*M_H2O + xg_H2_*M_H2 + xg_O2_*M_O2
        Fe32 = 2.0*exp(0.196*log(y_O2_*xg_O2_*P)+K5_)
        Fe3T = Fe32/(1.0+Fe32)
        wm_H2OT_ = (xm_H2OT_*M_H2O)/Xm_t
        return xg_H2O_, xg_H2_, xm_H2Omol_, xm_OH_, xm_H2OT_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2OT_, wm_H2_
    
    def mb_HOFe(xg_O2_):
        xg_H2O_, xg_H2_, xm_H2Omol_, xm_OH_, xm_H2OT_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2OT_, wm_H2_ = mg_HOFe(xg_O2_)
        diff,wt_g_O,wt_g_H = f_HOFe(xg_O2_)
        wt_g = (wt_g_O+wt_g_H)/2
        wt_H_ = 2.0*M_H*((wt_g*(((xg_H2O_+xg_H2_)/Xg_t) - (xm_H2OT_/Xm_t) - (wm_H2_/M_H2)))+(xm_H2OT_/Xm_t) + (wm_H2_/M_H2))
        wt_O_ = M_O*((wt_g*(((xg_H2O_ + 2.0*xg_O2_)/Xg_t)-(xm_H2OT_/Xm_t)))+(xm_H2OT_/Xm_t) + (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))
        return wt_g, wt_O_, wt_H_
    
    def f_HOFe(xg_O2_):
        xg_H2O_, xg_H2_, xm_H2Omol_, xm_OH_, xm_H2OT_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2OT_, wm_H2_ = mg_HOFe(xg_O2_)
        wt_g_H = ((wt_H/(2.0*M_H)) - (xm_H2OT_/Xm_t) - (wm_H2_/M_H2))/(((xg_H2O_+xg_H2_)/Xg_t) - (xm_H2OT_/Xm_t) - (wm_H2_/M_H2))
        wt_g_O = ((wt_O/M_O) - (xm_H2OT_/Xm_t) - (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))/(((xg_H2O_ + 2.0*xg_O2_)/Xg_t)-(xm_H2OT_/Xm_t))
        diff = wt_g_H - wt_g_O
        return diff,wt_g_O,wt_g_H
    
    def df_HOFe(xg_O2_,constants):
        result = de.HOFe_O2_xenia(xg_O2_,constants) 
        return result
    
    xg_O2_ = newton_raphson(guessx,nr_tol,nr_step,f_HOFe,df_HOFe)
    result1 = mg_HOFe(xg_O2_)
    result2 = f_HOFe(xg_O2_)
    result3 = mb_HOFe(xg_O2_)
    #print(wt_H, wt_O, xg_O2_, result1, result2, result3)
    return xg_O2_, result1, result2, result3


############
### SOFe ###
############

def eq_SOFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses):
    P = PT["P"]
    wt_O = bulk_wf['O']
    wt_S = bulk_wf['S']
    wt_Fe = bulk_wf['Fe']
    guessx = guesses["guessx"]
    
    # equilibrium constants
    K1_ = mdv.KOSg(PT,models)
    K2_ = mdv.C_S(PT,melt_wf,models)/1000000.
    K3_ = mdv.C_SO4(PT,melt_wf,models)/1000000.
    KD1_, KD2, y = mdv.FefO2_KC91_EqA_terms(PT,melt_wf,models)
    
    # fugacity coefficients
    y_SO2_ = mdv.y_SO2(PT,models)
    y_O2_ = mdv.y_O2(PT,models)
    y_S2_ = mdv.y_S2(PT,models)
    
    # molecular masses
    M_S = mdv.species.loc['S','M']
    M_O = mdv.species.loc['O','M']
    M_Fe = mdv.species.loc['Fe','M']
    M_O2 = mdv.species.loc['O2','M']
    M_S2 = mdv.species.loc['S2','M']
    M_SO2 = mdv.species.loc['SO2','M']
    M_SO3 = mdv.species.loc['SO3','M']        
    M_FeO = mdv.species.loc['FeO','M']
    M_FeO15 = mdv.species.loc['FeO1.5','M']
    M_m_ = mg.M_m_ox(melt_wf,models)
    
    constants = [P, wt_O, wt_S, wt_Fe, K1_, K2_, K3_, KD1_, KD2, y, y_SO2_, y_O2_, y_S2_, M_S, M_O, M_Fe, M_O2, M_S2, M_SO2, M_SO3, M_FeO, M_FeO15, M_m_]
     
    def mg_SOFe(xg_O2_):
        a = (y_SO2_**2.)/(K1_**2.*y_O2_**2.*xg_O2_**2.*y_S2_*P)
        b = 1.
        c = xg_O2_ - 1.
        xg_SO2_ = (-b + (b**2.-(4.*a*c))**0.5)/(2.*a)
        xg_S2_ = (((y_SO2_*xg_SO2_)/(K1_*y_O2_*xg_O2_))**2.)/(y_S2_*P)
        wm_S_ = K2_*(((y_S2_*xg_S2_)/(y_O2_*xg_O2_))**0.5)
        wm_SO3_ = M_SO3*((((y_S2_*xg_S2_)**0.5*(y_O2_*xg_O2_)**1.5*P**2.0)*K3_)/M_S)
        wm_ST_ = wm_S_ + ((M_S*wm_SO3_)/M_SO3)        
        Xg_t = xg_SO2_*M_SO2 + xg_S2_*M_S2 + xg_O2_*M_O2
        Fe32 = ((KD1_*(y_O2_*xg_O2_*P)**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y)))
        Fe3T = Fe32/(1+Fe32)
        S62 = (wm_SO3_/M_SO3)/(wm_S_/M_S)
        S6T = S62/(1+S62)
        return xg_SO2_, xg_S2_, Xg_t, Fe32, Fe3T, wm_S_, wm_SO3_, S62, S6T, wm_ST_
    
    def mb_SOFe(xg_O2_):
        xg_SO2_, xg_S2_, Xg_t, Fe32, Fe3T, wm_S_, wm_SO3_, S62, S6T, wm_ST_ = mg_SOFe(xg_O2_)
        diff,wt_g_O,wt_g_S = f_SOFe(xg_O2_)
        wt_g = (wt_g_O+wt_g_S)/2.
        wt_S_ = M_S*((wt_g*(((xg_SO2_+2.0*xg_S2_)/Xg_t) - (wm_S_/M_S) - (wm_SO3_/M_SO3)))+((wm_S_/M_S) + (wm_SO3_/M_SO3)))
        wt_O_ = M_O*((wt_g*(((2.0*xg_SO2_ + 2.0*xg_O2_)/Xg_t)-(3.0*wm_SO3_/M_SO3)))+(3.0*wm_SO3_/M_SO3) + (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))
        return wt_g, wt_O_, wt_S_
    
    def f_SOFe(xg_O2_):
        xg_SO2_, xg_S2_, Xg_t, Fe32, Fe3T, wm_S_, wm_SO3_, S62, S6T, wm_ST_ = mg_SOFe(xg_O2_)
        wt_g_S = ((wt_S/M_S) - (wm_S_/M_S) - (wm_SO3_/M_SO3))/(((xg_SO2_+2.0*xg_S2_)/Xg_t) - (wm_S_/M_S) - (wm_SO3_/M_SO3))
        wt_g_O = ((wt_O/M_O) - (3.0*wm_SO3_/M_SO3) - (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))/(((2.0*xg_SO2_ + 2.0*xg_O2_)/Xg_t) - (3.0*wm_SO3_/M_SO3))
        diff = wt_g_S - wt_g_O
        return diff,wt_g_O,wt_g_S
    
    def df_SOFe(xg_O2_,constants):
        result = de.SOFe_O2(xg_O2_,constants)
        return result
    
    xg_O2_ = newton_raphson(guessx,constants,nr_tol,nr_step,f_SOFe,df_SOFe)
    result1 = mg_SOFe(xg_O2_)
    result2 = f_SOFe(xg_O2_)
    result3 = mb_SOFe(xg_O2_)
    #print(wt_S, wt_O, xg_O2_, result1, result2, result3)
    return xg_O2_, result1, result2, result3


#############
### SHOFe ###
#############

# no H2S or H2
def eq_SHOFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses):
    P = PT["P"]
    wt_O = bulk_wf['O']
    wt_S = bulk_wf['S']
    wt_H = bulk_wf['H']
    wt_Fe = bulk_wf['Fe']
    guessx = guesses["guessx"]
    guessy = guesses["guessy"]
    
    # fugacity coefficients
    y_S2_ = mdv.y_S2(PT,models)
    y_O2_ = mdv.y_O2(PT,models)
    y_SO2_ = mdv.y_SO2(PT,models)
    y_H2S_ = mdv.y_H2S(PT,models)
    y_H2O_ = mdv.y_H2O(PT,models)
    y_H2_ = mdv.y_H2(PT,models)
    
    # molecular masses
    M_H = mdv.species.loc['H','M']
    M_S = mdv.species.loc['S','M']
    M_O = mdv.species.loc['O','M']
    M_Fe = mdv.species.loc['Fe','M']
    M_O2 = mdv.species.loc['O2','M']
    M_H2 = mdv.species.loc['H2','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_S2 = mdv.species.loc['S2','M']
    M_SO2 = mdv.species.loc['SO2','M']
    M_SO3 = mdv.species.loc['SO3','M']
    M_H2S = mdv.species.loc['H2S','M']
    M_FeO = mdv.species.loc['FeO','M']
    M_FeO15 = mdv.species.loc['FeO1.5','M']
    M_m_ = mg.M_m_SO(melt_wf)

    # equilibrium constants
    K1_ = mdv.KHOg(PT,models)
    K2_ = mdv.KOSg(PT,models)
    K3_ = mdv.KHOSg(PT,models)
    K4_ = mdv.C_H2O(PT,melt_wf,models)
    K5_ = mdv.C_S(PT,melt_wf,models)/1000000.0
    K6_ = (mdv.C_SO4(PT,melt_wf,models)/1000000.0)
    KD1_, KD2, y = mdv.FefO2_KC91_EqA_terms(PT,melt_wf,models)
    
    constants = [P, wt_O, wt_S, wt_H, wt_Fe, K1_, K2_, K3_, K4_, K5_, K6_, KD1_, KD2, y, y_S2_, y_O2_, y_SO2_, y_H2S_, y_H2O_, y_H2_, M_H, M_S, M_O, M_Fe, M_O2, M_H2, M_H2O, M_S2, M_SO2, M_SO3, M_H2S, M_FeO, M_FeO15, M_m_]
     
    def mg_SHOFe(xg_O2_,xg_S2_):
        xg_SO2_ = (K2_*y_O2_*xg_O2_*(y_S2_*xg_S2_*P)**0.5)/y_SO2_
        xg_H2O_ = (1.0 - xg_O2_ - xg_S2_ - xg_SO2_)/(1.0 + ((K3_*(y_S2_*xg_S2_)**0.5*y_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5)) + (y_H2O_/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)))
        xg_H2_ = (y_H2O_*xg_H2O_)/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)
        xg_H2S_ = (K3_*(y_S2_*xg_S2_)**0.5*y_H2O_*xg_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5)
        xm_H2O_ = (K4_*y_H2O_*xg_H2O_*P)**0.5
        wm_S_ = K5_*((y_S2_*xg_S2_)/(y_O2_*xg_O2_))**0.5
        wm_SO3_ = ((K6_*(y_S2_*xg_S2_)**0.5*(y_O2_*xg_O2_)**1.5*P**2.0)/M_S)*M_SO3
        wm_ST_ = wm_S_ + ((M_S*wm_SO3_)/M_SO3)
        Xm_t = xm_H2O_*M_H2O + (1.0-xm_H2O_)*M_m_
        Xg_t = xg_SO2_*M_SO2 + xg_S2_*M_S2 + xg_O2_*M_O2 + xg_H2O_*M_H2O + xg_H2_*M_H2 + xg_H2S_*M_H2S
        Fe32 = ((KD1_*(y_O2_*xg_O2_*P)**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y)))
        Fe3T = Fe32/(1.0+Fe32)
        S62 = (wm_SO3_/M_SO3)/(wm_S_/M_S)
        S6T = S62/(1.0+S62)
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t
        return xg_SO2_, xg_H2O_, xg_H2_, xg_H2S_, xm_H2O_, wm_S_, wm_SO3_, Xm_t, Xg_t, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_ST_
    
    def mb_SHOFe(xg_O2_, xg_S2_):
        xg_SO2_, xg_H2O_, xg_H2_, xg_H2S_, xm_H2O_, wm_S_, wm_SO3_, Xm_t, Xg_t, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_ST_ = mg_SHOFe(xg_O2_, xg_S2_)
        mbSO, mbSH, wt_g_O, wt_g_S, wt_g_H = f_SHOFe(xg_O2_,xg_S2_)
        wt_g = (wt_g_O+wt_g_H+wt_g_S)/3.0
        wt_H_ = 2.0*M_H*((wt_g*(((xg_H2O_+xg_H2_+xg_H2S_)/Xg_t) - (xm_H2O_/Xm_t)))+(xm_H2O_/Xm_t))
        wt_O_ = M_O*((wt_g*(((2.0*xg_SO2_ + 2.0*xg_O2_ + xg_H2O_)/Xg_t) - (xm_H2O_/Xm_t) - (3.0*wm_SO3_/M_SO3)))+(xm_H2O_/Xm_t) + (3.0*wm_SO3_/M_SO3) + (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))
        wt_S_ = M_S*((wt_g*(((xg_SO2_+2.0*xg_S2_+xg_H2S_)/Xg_t) - (wm_S_/M_S) - (wm_SO3_/M_SO3))) + (wm_S_/M_S) + (wm_SO3_/M_SO3))
        return wt_g, wt_O_, wt_S_, wt_H_
    
    def f_SHOFe(xg_O2_,xg_S2_):
        xg_SO2_, xg_H2O_, xg_H2_, xg_H2S_, xm_H2O_, wm_S_, wm_SO3_, Xm_t, Xg_t, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_ST_ = mg_SHOFe(xg_O2_, xg_S2_)
        wt_g_S = ((wt_S/M_S) - (wm_S_/M_S) - (wm_SO3_/M_SO3))/(((xg_SO2_+2.0*xg_S2_+xg_H2S_)/Xg_t) - (wm_S_/M_S) - (wm_SO3_/M_SO3))
        wt_g_O = ((wt_O/M_O) - (xm_H2O_/Xm_t) - (3.0*wm_SO3_/M_SO3) - (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))/(((2.0*xg_SO2_ + 2.0*xg_O2_ + xg_H2O_)/Xg_t) - (xm_H2O_/Xm_t) - (3.0*wm_SO3_/M_SO3))
        wt_g_H = ((wt_H/(2.0*M_H)) - (xm_H2O_/Xm_t))/(((xg_H2O_+xg_H2_+xg_H2S_)/Xg_t) - (xm_H2O_/Xm_t))
        mbSO = (wt_g_S - wt_g_O)
        mbSH = (wt_g_S - wt_g_H)
        return mbSO, mbSH, wt_g_O, wt_g_S, wt_g_H
    
    def df_SHOFe(xg_O2_,xg_S2_,constants):
        dmbSO_O = de.SHOFe_SO_O2(xg_O2_,xg_S2_,constants) 
        dmbSO_S = de.SHOFe_SO_S2(xg_O2_,xg_S2_,constants)
        dmbSH_O = de.SHOFe_SH_O2(xg_O2_,xg_S2_,constants) 
        dmbSH_S = de.SHOFe_SH_S2(xg_O2_,xg_S2_,constants)
        return dmbSO_O, dmbSO_S, dmbSH_O, dmbSH_S
    
    xg_O2_, xg_S2_ = jac_newton(guessx,guessy,constants,f_SHOFe,df_SHOFe,nr_step,nr_tol)
    results1 = mg_SHOFe(xg_O2_,xg_S2_)
    results2 = f_SHOFe(xg_O2_,xg_S2_)
    results3 = mb_SHOFe(xg_O2_,xg_S2_)
    #print(wt_O, wt_S, wt_H, xg_O2, xg_S2, results1, results2, results3)
    return xg_O2_, xg_S2_, results1, results2, results3

# includes H2S and H2
def eq_SHOFe_2(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses):
    P = PT["P"]
    wt_O = bulk_wf['O']
    wt_S = bulk_wf['S']
    wt_H = bulk_wf['H']
    wt_Fe = bulk_wf['Fe']
    guessx = guesses["guessx"]
    guessy = guesses["guessy"]
        
    # equilibrium constants
    K1_ = mdv.KHOg(PT,models)
    K4_ = mdv.C_H2O(PT,melt_wf,models)
    K6_ = mdv.KOSg(PT,models)
    K7_ = mdv.KHOSg(PT,models)
    K8_ = mdv.C_S(PT,melt_wf,models)/1000000.0
    K9_ = mdv.C_SO4(PT,melt_wf,models)/1000000.0
    K13_ = mdv.C_H2(PT,melt_wf,models)/1000000.0
    K14_ = mdv.C_H2S(PT,melt_wf,models)/1000000.0
    KD1_, KD2, y = mdv.FefO2_KC91_EqA_terms(PT,melt_wf,models)
   
    # fugacity coefficients
    y_O2_ = mdv.y_O2(PT,models)
    y_H2O_ = mdv.y_H2O(PT,models)
    y_H2_ = mdv.y_H2(PT,models)
    y_S2_ = mdv.y_S2(PT,models)
    y_SO2_ = mdv.y_SO2(PT,models)
    y_H2S_ = mdv.y_H2S(PT,models)
    
    # molecular masses
    M_H = mdv.species.loc['H','M']
    M_O = mdv.species.loc['O','M']
    M_S = mdv.species.loc['S','M']
    M_Fe = mdv.species.loc['Fe','M']
    M_O2 = mdv.species.loc['O2','M']
    M_OH = mdv.species.loc['OH','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_H2 = mdv.species.loc['H2','M']
    M_S2 = mdv.species.loc['S2','M']
    M_SO2 = mdv.species.loc['SO2','M']
    M_SO3 = mdv.species.loc['SO3','M']
    M_H2S = mdv.species.loc['H2S','M']
    M_FeO = mdv.species.loc['FeO','M']
    M_FeO15 = mdv.species.loc['FeO1.5','M']
    M_m_ = mg.M_m_SO(melt_wf)
    M_m_ox = mg.M_m_ox(melt_wf,models)
    
    constants = [P, wt_O, wt_H, wt_S, wt_Fe, K1_, K4_, K6_, K7_, K8_, K9_, K13_, K14_, KD1_, KD2, y, y_H2_, y_O2_, y_H2O_, y_SO2_, y_S2_, y_H2S_, M_H, M_O, M_S, M_Fe, M_O2, M_OH, M_H2, M_H2O, M_S2, M_SO2, M_SO3, M_H2S, M_FeO, M_FeO15, M_m_, M_m_ox]

    def mg_SHOFe(xg_O2_,xg_S2_):
        xg_O2__ = xg_O2_
        xg_S2__ = xg_S2_
        xg_SO2_ = (K6_*y_O2_*xg_O2_*(y_S2_*xg_S2_*P)**0.5)/y_SO2_
        xg_H2O_ = (1. - xg_O2_ - xg_S2_ - xg_SO2_)/(1. + ((K7_*(y_S2_*xg_S2_)**0.5*y_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5)) + (y_H2O_/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)))
        xg_H2_ = (y_H2O_*xg_H2O_)/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)
        xg_H2S_ = (K7_*(y_S2_*xg_S2_)**0.5*y_H2O_*xg_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5)
        Xg_t = xg_SO2_*M_SO2 + xg_S2_*M_S2 + xg_O2_*M_O2 + xg_H2O_*M_H2O + xg_H2_*M_H2 + xg_H2S_*M_H2S
        xm_H2O_ = (K4_*y_H2O_*xg_H2O_*P)**0.5
        Xm_t = xm_H2O_*M_H2O + (1.0-xm_H2O_)*M_m_
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t
        wm_H2_ = K13_*y_H2_*xg_H2_*P
        wm_H2S_ = K14_*y_H2S_*xg_H2S_*P
        wm_S_ = K8_*((y_S2_*xg_S2_)/(y_O2_*xg_O2_))**0.5
        wm_SO3_ = ((K9_*(y_S2_*xg_S2_)**0.5*(y_O2_*xg_O2_)**1.5*P**2.0)/M_S)*M_SO3
        wm_ST_ = wm_S_ + ((M_S*wm_SO3_)/M_SO3) + ((M_S*wm_H2S_)/M_H2S)
        Fe32 = ((KD1_*(y_O2_*xg_O2_*P)**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y)))
        Fe3T = Fe32/(1.0+Fe32)
        S62 = (wm_SO3_/M_SO3)/(wm_S_/M_S)
        S6T = S62/(1.0+S62)
        Xm_t_ox = ""
        return xg_O2__, xg_H2_, xg_S2__, xg_H2O_, xg_SO2_, xg_H2S_, Xg_t, xm_H2O_, wm_S_, wm_SO3_, Xm_t, Xm_t_ox, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_ST_, wm_H2S_, wm_H2_
    
    def mb_SHOFe(xg_O2_,xg_S2_): 
        xg_O2__, xg_H2_, xg_S2_, xg_H2O_, xg_SO2_, xg_H2S_, Xg_t, xm_H2O_, wm_S_, wm_SO3_, Xm_t, Xm_t_ox, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_ST_, wm_H2S_, wm_H2_ = mg_SHOFe(xg_O2_, xg_S2_)
        mba, mbc, wt_g_O, wt_g_H, wt_g_S = f_SHOFe(xg_O2_,xg_S2_)
        wt_g = (wt_g_O+wt_g_H+wt_g_S)/3.
        wt_H_ = 2.*M_H*((wt_g*(((xg_H2O_+xg_H2_+xg_H2S_)/Xg_t) - (wm_H2O_/M_H2O) - (wm_H2_/M_H2) - (wm_H2S_/M_H2S))) + (wm_H2O_/M_H2O) + (wm_H2_/M_H2) + (wm_H2S_/M_H2S))
        wt_O_ = M_O*((wt_g*(((2.*xg_O2_ + xg_H2O_ + 2.*xg_SO2_)/Xg_t) - (wm_H2O_/M_H2O) - (3.*wm_SO3_/M_SO3))) + (wm_H2O_/M_H2O) + (3.*wm_SO3_/M_SO3) + ((wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0))))
        wt_S_ = M_S*((wt_g*(((xg_SO2_+2.0*xg_S2_+xg_H2S_)/Xg_t) - (wm_S_/M_S) - (wm_SO3_/M_SO3) - (wm_H2S_/M_H2S))) + (wm_S_/M_S) + (wm_SO3_/M_SO3) + (wm_H2S_/M_H2S))
        return wt_g, wt_O_, wt_H_, wt_S_
    
    def f_SHOFe(xg_O2_,xg_S2_):
        xg_O2__, xg_H2_, xg_S2_, xg_H2O_, xg_SO2_, xg_H2S_, Xg_t, xm_H2O_, wm_S_, wm_SO3_, Xm_t, Xm_t_ox, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_ST_, wm_H2S_, wm_H2_ = mg_SHOFe(xg_O2_, xg_S2_)
        wt_g_O = ((wt_O/M_O) - (wm_H2O_/M_H2O) - (3.*wm_SO3_/M_SO3) - (wt_Fe/M_Fe)*((1.5*Fe32+1.)/(Fe32+1.)))/(((2.*xg_O2_ + xg_H2O_ + 2.0*xg_SO2_)/Xg_t) - (wm_H2O_/M_H2O) - (3.0*wm_SO3_/M_SO3))
        wt_g_H = ((wt_H/(2.*M_H)) - (wm_H2O_/M_H2O) - (wm_H2_/M_H2) - (wm_H2S_/M_H2S))/(((xg_H2O_+xg_H2_+xg_H2S_)/Xg_t) - (wm_H2O_/M_H2O) - (wm_H2_/M_H2) - (wm_H2S_/M_H2S))
        wt_g_S = ((wt_S/M_S) - (wm_S_/M_S) - (wm_SO3_/M_SO3) - (wm_H2S_/M_H2S))/(((xg_SO2_+2.0*xg_S2_+xg_H2S_)/Xg_t) - (wm_S_/M_S) - (wm_SO3_/M_SO3) - (wm_H2S_/M_H2S))
        mba = (wt_g_H - wt_g_O) # mbHO
        mbc = (wt_g_H - wt_g_S) # mbHS
        return mba, mbc, wt_g_O, wt_g_H, wt_g_S
    
    def df_SHOFe(xg_O2_,xg_S2_,constants): 
        dmba_O = de.SHOFe2_HO_O2(xg_O2_,xg_S2_,constants) # dmbHO_O
        dmba_S = de.SHOFe2_HO_S2(xg_O2_,xg_S2_,constants) # dmbHO_S
        dmbc_O = de.SHOFe2_HS_O2(xg_O2_,xg_S2_,constants) # dmbHS_O
        dmbc_S = de.SHOFe2_HS_S2(xg_O2_,xg_S2_,constants) # dmbHS_S
        return dmba_O, dmba_S, dmbc_O, dmbc_S

    xg_O2_, xg_S2_ = jac_newton(guessx,guessy,constants,f_SHOFe,df_SHOFe,nr_step,nr_tol)
    results1 = xg_O2_, xg_S2_
    results2 = mg_SHOFe(xg_O2_,xg_S2_)
    results3 = f_SHOFe(xg_O2_,xg_S2_)
    results4 = mb_SHOFe(xg_O2_,xg_S2_)
    return results1, results2, results3, results4

#############
### SCOFe ###
#############

### not finished ###
def eq_SCOFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses):
    P = PT["P"]
    wt_O = bulk_wf['O']
    wt_S = bulk_wf['S']
    wt_C = bulk_wf['C']
    wt_Fe = bulk_wf['Fe']
    guessx = guesses["guessx"]
    guessy = guesses["guessy"]
           
    # fugacity coefficients
    y_S2_ = mdv.y_S2(PT,models)
    y_O2_ = mdv.y_O2(PT,models)
    y_SO2_ = mdv.y_SO2(PT,models)
    y_CO_ = mdv.y_CO(PT,models)
    y_CO2_ = mdv.y_CO2(PT,models)
    y_OCS_ = mdv.y_OCS(PT,models)
    
    # molecular masses
    M_C = mdv.species.loc['C','M']
    M_S = mdv.species.loc['S','M']
    M_O = mdv.species.loc['O','M']
    M_Fe = mdv.species.loc['Fe','M']
    M_O2 = mdv.species.loc['O2','M']
    M_CO = mdv.species.loc['CO','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_S2 = mdv.species.loc['S2','M']
    M_SO2 = mdv.species.loc['SO2','M']
    M_SO3 = mdv.species.loc['SO3','M']
    M_OCS = mdv.species.loc['OCS','M']
    M_FeO = mdv.species.loc['FeO','M']
    M_FeO15 = mdv.species.loc['FeO1.5','M']
    M_m_ = mg.M_m_SO(melt_wf)

    # equilibrium constants
    K1_ = mdv.KCOg(PT,models)
    K2_ = mdv.KOSg(PT,models)
    K3_ = mdv.KOCSg(PT,models) 
    K4_ = mdv.C_CO3(PT,melt_wf,models)
    K5_ = mdv.C_S(PT,melt_wf,models)/1000000.0
    K6_ = mdv.C_SO4(PT,melt_wf,models)/1000000.0
    K7_ = (mdv.C_CO(PT,melt_wf,models)/1000000.0)
    KD1_, KD2, y = mdv.FefO2_KC91_EqA_terms(PT,melt_wf,models)
    
    constants = [P, wt_O, wt_S, wt_C, wt_Fe, K1_, K2_, K3_, K4_, K5_, K6_, K7_, KD1_, KD2, y, y_S2_, y_O2_, y_SO2_, y_CO_, y_CO2_, y_OCS_, M_C, M_S, M_O, M_Fe, M_O2, M_CO, M_CO2, M_S2, M_SO2, M_SO3, M_OCS, M_FeO, M_FeO15, M_m_]
     
    def mg_SCOFe(xg_O2_,xg_S2_):
        xg_SO2_ = (K2_*y_O2_*xg_O2_*(y_S2_*xg_S2_*P)**0.5)/y_SO2_
        xg_CO2_ = (1.-xg_O2_-xg_S2_-xg_SO2_)/(1.+(y_CO2_/(K1_*y_CO_*(y_O2_*xg_O2_*P)**0.5)+((y_SO2_*xg_SO2_*y_CO2_)/(K1_**3.*K3_*y_OCS_*(y_O2_*xg_O2_)**1.5*P**0.5))))
        xg_CO_ = (y_CO2_*xg_CO2_)/(y_CO_*K1_*(y_O2_*xg_O2_*P)**0.5)
        xg_OCS_ = ((y_CO_*xg_CO_)**3.*y_SO2_*xg_SO2_*P)/(K3_*(y_CO2_*xg_CO2_)**2.*y_OCS_)
        wm_S_ = K5_*((y_S2_*xg_S2_)/(y_O2_*xg_O2_))**0.5
        wm_SO3_ = ((K6_*(y_S2_*xg_S2_)**0.5*(y_O2_*xg_O2_)**1.5*P**2.)/M_S)*M_SO3
        wm_ST_ = wm_S_ + ((M_S*wm_SO3_)/M_SO3)
        xm_CO2_ = K4_*y_CO2_*xg_CO2_*P
        Xm_t = xm_CO2_*M_CO2 + (1.-xm_CO2_)*M_m_
        Xg_t = xg_CO2_*M_CO2 + xg_CO_*M_CO + xg_O2_*M_O2 + xg_SO2_*M_SO2 + xg_S2_*M_S2 + xg_OCS_*M_OCS
        Fe32 = ((KD1_*(y_O2_*xg_O2_*P)**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y)))
        Fe3T = Fe32/(1.0+Fe32)
        S62 = (wm_SO3_/M_SO3)/(wm_S_/M_S)
        S6T = S62/(1.0+S62)
        wm_CO2_ = (xm_CO2_*M_CO2)/Xm_t
        if models.loc["COH_species","option"] == "no_H2_CO_CH4_melt":
            wm_CO_ = 0.
        elif models.loc["COH_species","option"] == "yes_H2_CO_CH4_melt":
            wm_CO_ = K7_*y_CO_*xg_CO_*P
        return xg_SO2_, xg_CO2_, xg_CO_, xg_OCS_, xg_S2_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Xg_t, Fe32, Fe3T, S62, S6T, wm_CO2_, wm_ST_, wm_CO_
    
    def mb_SCOFe(xg_O2_, xg_S2_):
        xg_SO2_, xg_CO2_, xg_CO_, xg_OCS_, xg_S2_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Xg_t, Fe32, Fe3T, S62, S6T, wm_CO2_, wm_ST_, wm_CO_ = mg_SCOFe(xg_O2_, xg_S2_)
        mbSO, mbSH, wt_g_O, wt_g_S, wt_g_C = f_SCOFe(xg_O2_,xg_S2_)
        wt_g = (wt_g_O+wt_g_C+wt_g_S)/3.
        wt_C_ = M_C*((wt_g*(((xg_CO2_+xg_CO_+xg_OCS_)/Xg_t) - (xm_CO2_/Xm_t) - (wm_CO_/M_CO)))+(xm_CO2_/Xm_t) + (wm_CO_/M_CO))
        wt_O_ = M_O*((wt_g*(((2.*xg_SO2_ + 2.*xg_O2_ + 2.*xg_CO2_ + xg_OCS_ + xg_CO_)/Xg_t) - 2.*(xm_CO2_/Xm_t) - (3.*wm_SO3_/M_SO3) - (wm_CO_/M_CO)))+(2.*xm_CO2_/Xm_t) + (3.*wm_SO3_/M_SO3) + (wm_CO_/M_CO) + (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))
        wt_S_ = M_S*((wt_g*(((xg_SO2_+2.*xg_S2_+xg_OCS_)/Xg_t) - (wm_S_/M_S) - (wm_SO3_/M_SO3))) + (wm_S_/M_S) + (wm_SO3_/M_SO3))
        return wt_g, wt_O_, wt_S_, wt_C_
    
    def f_SCOFe(xg_O2_,xg_S2_):
        xg_SO2_, xg_CO2_, xg_CO_, xg_OCS_, xg_S2_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Xg_t, Fe32, Fe3T, S62, S6T, wm_CO2_, wm_ST_, wm_CO_ = mg_SCOFe(xg_O2_, xg_S2_)
        wt_g_S = ((wt_S/M_S) - (wm_S_/M_S) - (wm_SO3_/M_SO3))/(((xg_SO2_+2.*xg_S2_+xg_OCS_)/Xg_t) - (wm_S_/M_S) - (wm_SO3_/M_SO3))
        wt_g_O = ((wt_O/M_O) - 2.*(xm_CO2_/Xm_t) - (3.*wm_SO3_/M_SO3) - (wm_CO_/M_CO) - (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))/(((2.*xg_SO2_ + 2.*xg_O2_ + 2.*xg_CO2_ + xg_CO_ + xg_OCS_)/Xg_t) - 2.*(xm_CO2_/Xm_t) - (3.*wm_SO3_/M_SO3) - (wm_CO_/M_CO))
        wt_g_C = ((wt_C/M_C) - (xm_CO2_/Xm_t) - (wm_CO_/M_CO))/(((xg_CO2_+xg_CO_+xg_OCS_)/Xg_t) - (xm_CO2_/Xm_t) - (wm_CO_/M_CO))
        mbSO = (wt_g_S - wt_g_O)
        mbSC = (wt_g_S - wt_g_C)
        return mbSO, mbSC, wt_g_O, wt_g_S, wt_g_C
    
    def df_SCOFe(xg_O2_,xg_S2_,constants):
        dmbSO_O = de.SCOFe_SO_O2(xg_O2_,xg_S2_,constants)
        dmbSO_S = de.SCOFe_SO_S2(xg_O2_,xg_S2_,constants)
        dmbSC_O = de.SCOFe_SC_O2(xg_O2_,xg_S2_,constants) 
        dmbSC_S = de.SCOFe_SC_S2(xg_O2_,xg_S2_,constants)
        return dmbSO_O, dmbSO_S, dmbSC_O, dmbSC_S
    
    xg_O2_, xg_S2_ = jac_newton(guessx,guessy,constants,f_SCOFe,df_SCOFe,nr_step,nr_tol)
    results1 = mg_SCOFe(xg_O2_,xg_S2_)
    results2 = f_SCOFe(xg_O2_,xg_S2_)
    results3 = mb_SCOFe(xg_O2_,xg_S2_)
    return xg_O2_, xg_S2_, results1, results2, results3

#############
### COXFe ###
#############

def eq_COXFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses):
    species_X = models.loc["species X","option"]
    P = PT["P"]
    wt_O = bulk_wf['O']
    wt_X = bulk_wf['X']
    wt_C = bulk_wf['C']
    wt_Fe = bulk_wf['Fe']
    guessx = guesses["guessx"]
    guessy = guesses["guessy"]
           
    # fugacity coefficients
    y_X_ = mdv.y_X(PT,models)
    y_O2_ = mdv.y_O2(PT,models)
    y_CO_ = mdv.y_CO(PT,models)
    y_CO2_ = mdv.y_CO2(PT,models)
    
    # molecular masses
    M_C = mdv.species.loc['C','M']
    M_X = mdv.species.loc[species_X,'M'] 
    M_O = mdv.species.loc['O','M']
    M_Fe = mdv.species.loc['Fe','M']
    M_O2 = mdv.species.loc['O2','M']
    M_CO = mdv.species.loc['CO','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_FeO = mdv.species.loc['FeO','M']
    M_FeO15 = mdv.species.loc['FeO1.5','M']
    M_m_ = mg.M_m_SO(melt_wf)

    # equilibrium constants
    K1_ = mdv.KCOg(PT,models)
    K2_ = mdv.C_CO3(PT,melt_wf,models) # was 4
    K3_ = (mdv.C_CO(PT,melt_wf,models)/1000000.0) # was 7
    K4_ = (mdv.C_X(PT,melt_wf,models)/1000000.0)
    KD1_, KD2, y = mdv.FefO2_KC91_EqA_terms(PT,melt_wf,models)
    
    constants = [P, wt_O, wt_X, wt_C, wt_Fe, K1_, K2_, K3_, K4_, KD1_, KD2, y, y_X_, y_O2_, y_CO_, y_CO2_, M_C, M_X, M_O, M_Fe, M_O2, M_CO, M_CO2, M_FeO, M_FeO15, M_m_]
     
    def mg_COXFe(xg_O2_,xg_CO_):
        xg_CO2_ = (K1_*y_CO_*xg_CO_*(y_O2_*xg_O2_*P)**0.5)/y_CO2_
        xg_X_ = 1. - xg_O2_ - xg_CO_ - xg_CO2_
        wm_X_ = K4_*xg_X_*y_X_*P
        xm_CO2_ = K2_*y_CO2_*xg_CO2_*P
        Xm_t = xm_CO2_*M_CO2 + (1.-xm_CO2_)*M_m_
        Xg_t = xg_CO2_*M_CO2 + xg_CO_*M_CO + xg_O2_*M_O2 + xg_X_*M_X
        Fe32 = ((KD1_*(y_O2_*xg_O2_*P)**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y)))
        Fe3T = Fe32/(1.0+Fe32)
        wm_CO2_ = (xm_CO2_*M_CO2)/Xm_t
        if models.loc["COH_species","option"] == "no_H2_CO_CH4_melt":
            wm_CO_ = 0.
        elif models.loc["COH_species","option"] == "yes_H2_CO_CH4_melt":
            wm_CO_ = K3_*y_CO_*xg_CO_*P
        return xg_X_, xg_CO2_, xg_CO_, xm_CO2_, wm_X_, Xm_t, Xg_t, Fe32, Fe3T, wm_CO2_, wm_CO_
    
    def mb_COXFe(xg_O2_, xg_CO_):
        xg_X_, xg_CO2_, xg_CO_, xm_CO2_, wm_X_, Xm_t, Xg_t, Fe32, Fe3T, wm_CO2_, wm_CO_ = mg_COXFe(xg_O2_, xg_CO_)
        mbXO, mbXC, wt_g_O, wt_g_X, wt_g_C = f_COXFe(xg_O2_,xg_CO_)
        wt_g = (wt_g_O+wt_g_C+wt_g_X)/3.
        wt_C_ = M_C*((wt_g*(((xg_CO2_+xg_CO_)/Xg_t) - (xm_CO2_/Xm_t) - (wm_CO_/M_CO)))+(xm_CO2_/Xm_t) + (wm_CO_/M_CO))
        wt_O_ = M_O*((wt_g*(((2.*xg_O2_ + 2.*xg_CO2_ + xg_CO_)/Xg_t) - 2.*(xm_CO2_/Xm_t) - (wm_CO_/M_CO)))+(2.*xm_CO2_/Xm_t) + (wm_CO_/M_CO) + (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))
        wt_X_ = M_X*((wt_g*(((xg_X_)/Xg_t) - (wm_X_/M_X))) + (wm_X_/M_X))
        return wt_g, wt_O_, wt_X_, wt_C_
    
    def f_COXFe(xg_O2_,xg_CO_):
        xg_X_, xg_CO2_, xg_CO_, xm_CO2_, wm_X_, Xm_t, Xg_t, Fe32, Fe3T, wm_CO2_, wm_CO_ = mg_COXFe(xg_O2_, xg_CO_)
        wt_g_X = ((wt_X/M_X) - (wm_X_/M_X))/(((xg_X_)/Xg_t) - (wm_X_/M_X))
        wt_g_O = ((wt_O/M_O) - 2.*(xm_CO2_/Xm_t) - (wm_CO_/M_CO) - (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))/(((2.*xg_O2_ + 2.*xg_CO2_ + xg_CO_)/Xg_t) - 2.*(xm_CO2_/Xm_t) - (wm_CO_/M_CO))
        wt_g_C = ((wt_C/M_C) - (xm_CO2_/Xm_t) - (wm_CO_/M_CO))/(((xg_CO2_+xg_CO_)/Xg_t) - (xm_CO2_/Xm_t) - (wm_CO_/M_CO))
        mbXO = (wt_g_X - wt_g_O)
        mbXC = (wt_g_X - wt_g_C)
        return mbXO, mbXC, wt_g_O, wt_g_X, wt_g_C
    
    def df_COXFe(xg_O2_,xg_CO_,constants):
        dmbXO_O = de.COXFe_XO_O2(xg_O2_,xg_CO_,constants)
        dmbXO_C = de.COXFe_XO_CO(xg_O2_,xg_CO_,constants)
        dmbXC_O = de.COXFe_XC_O2(xg_O2_,xg_CO_,constants) 
        dmbXC_C = de.COXFe_XC_CO(xg_O2_,xg_CO_,constants)
        return dmbXO_O, dmbXO_C, dmbXC_O, dmbXC_C
    
    xg_O2_, xg_CO_ = jac_newton(guessx,guessy,constants,f_COXFe,df_COXFe,nr_step,nr_tol)
    results1 = mg_COXFe(xg_O2_,xg_CO_)
    results2 = f_COXFe(xg_O2_,xg_CO_)
    results3 = mb_COXFe(xg_O2_,xg_CO_)
    return xg_O2_, xg_CO_, results1, results2, results3

#############
### HOXFe ###
#############

def eq_HOXFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses):
    species_X = models.loc["species X","option"]
    P = PT["P"]
    wt_O = bulk_wf['O']
    wt_X = bulk_wf['X']
    wt_H = bulk_wf['H']
    wt_Fe = bulk_wf['Fe']
    guessx = guesses["guessx"]
    guessy = guesses["guessy"]
           
    # fugacity coefficients
    y_X_ = mdv.y_X(PT,models)
    y_O2_ = mdv.y_O2(PT,models)
    y_H2_ = mdv.y_H2(PT,models)
    y_H2O_ = mdv.y_H2O(PT,models)
    
    # molecular masses
    M_H = mdv.species.loc['H','M']
    M_X = mdv.species.loc[species_X,'M'] 
    M_O = mdv.species.loc['O','M']
    M_Fe = mdv.species.loc['Fe','M']
    M_O2 = mdv.species.loc['O2','M']
    M_H2 = mdv.species.loc['H2','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_FeO = mdv.species.loc['FeO','M']
    M_FeO15 = mdv.species.loc['FeO1.5','M']
    M_m_ = mg.M_m_SO(melt_wf)

    # equilibrium constants
    K1_ = mdv.KHOg(PT,models)
    K2_ = mdv.C_H2O(PT,melt_wf,models)
    K3_ = (mdv.C_H2(PT,melt_wf,models)/1000000.0)
    K4_ = (mdv.C_X(PT,melt_wf,models)/1000000.0)
    KD1_, KD2, y = mdv.FefO2_KC91_EqA_terms(PT,melt_wf,models)
    
    constants = [P, wt_O, wt_X, wt_H, wt_Fe, K1_, K2_, K3_, K4_, KD1_, KD2, y, y_X_, y_O2_, y_H2_, y_H2O_, M_H, M_X, M_O, M_Fe, M_O2, M_H2, M_H2O, M_FeO, M_FeO15, M_m_]
     
    def mg_HOXFe(xg_O2_,xg_H2_):
        xg_H2O_ = (K1_*y_H2_*xg_H2_*(y_O2_*xg_O2_*P)**0.5)/y_H2O_
        xg_X_ = 1. - xg_O2_ - xg_H2_ - xg_H2O_
        wm_X_ = K4_*xg_X_*y_X_*P
        xm_H2O_ = (K2_*y_H2O_*xg_H2O_*P)**0.5
        Xm_t = xm_H2O_*M_H2O + (1.-xm_H2O_)*M_m_
        Xg_t = xg_H2O_*M_H2O + xg_H2_*M_H2 + xg_O2_*M_O2 + xg_X_*M_X
        Fe32 = ((KD1_*(y_O2_*xg_O2_*P)**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y)))
        Fe3T = Fe32/(1.0+Fe32)
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t
        if models.loc["COH_species","option"] == "no_H2_CO_CH4_melt":
            wm_H2_ = 0.
        elif models.loc["COH_species","option"] == "yes_H2_CO_CH4_melt":
            wm_H2_ = K3_*y_H2_*xg_H2_*P
        return xg_X_, xg_H2O_, xg_H2_, xm_H2O_, wm_X_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2O_, wm_H2_
    
    def mb_HOXFe(xg_O2_, xg_H2_):
        xg_X_, xg_H2O_, xg_H2_, xm_H2O_, wm_X_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2O_, wm_H2_ = mg_HOXFe(xg_O2_, xg_H2_)
        mbXO, mbXH, wt_g_O, wt_g_X, wt_g_H = f_HOXFe(xg_O2_,xg_H2_)
        wt_g = (wt_g_O+wt_g_H+wt_g_X)/3.
        wt_H_ = 2.0*M_H*((wt_g*(((xg_H2O_+xg_H2_)/Xg_t) - (xm_H2O_/Xm_t) - (wm_H2_/M_H2)))+(xm_H2O_/Xm_t) + (wm_H2_/M_H2))
        wt_O_ = M_O*((wt_g*(((2.0*xg_O2_ + xg_H2O_)/Xg_t) - (xm_H2O_/Xm_t)))+(xm_H2O_/Xm_t) + (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))
        wt_X_ = M_X*((wt_g*(((xg_X_)/Xg_t) - (wm_X_/M_X))) + (wm_X_/M_X))
        return wt_g, wt_O_, wt_X_, wt_H_
    
    def f_HOXFe(xg_O2_,xg_H2_):
        xg_X_, xg_H2O_, xg_H2_, xm_H2O_, wm_X_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2O_, wm_H2_ = mg_HOXFe(xg_O2_, xg_H2_)
        wt_g_X = ((wt_X/M_X) - (wm_X_/M_X))/(((xg_X_)/Xg_t) - (wm_X_/M_X))
        wt_g_O = ((wt_O/M_O) - (xm_H2O_/Xm_t) - (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))/(((2.0*xg_O2_ + xg_H2O_)/Xg_t) - (xm_H2O_/Xm_t))
        wt_g_H = ((wt_H/(2.0*M_H)) - (xm_H2O_/Xm_t)  - (wm_H2_/M_H2))/(((xg_H2O_+xg_H2_)/Xg_t) - (xm_H2O_/Xm_t) - (wm_H2_/M_H2))
        mbXO = (wt_g_X - wt_g_O)
        mbXH = (wt_g_X - wt_g_H)
        return mbXO, mbXH, wt_g_O, wt_g_X, wt_g_H
    
    def df_HOXFe(xg_O2_,xg_H2_,constants):
        dmbXO_O = de.HOXFe_XO_O2(xg_O2_,xg_H2_,constants)
        dmbXO_H = de.HOXFe_XO_H2(xg_O2_,xg_H2_,constants)
        dmbXH_O = de.HOXFe_XH_O2(xg_O2_,xg_H2_,constants) 
        dmbXH_H = de.HOXFe_XH_H2(xg_O2_,xg_H2_,constants)
        return dmbXO_O, dmbXO_H, dmbXH_O, dmbXH_H
    
    xg_O2_, xg_H2_ = jac_newton(guessx,guessy,constants,f_HOXFe,df_HOXFe,nr_step,nr_tol)
    results1 = mg_HOXFe(xg_O2_,xg_H2_)
    results2 = f_HOXFe(xg_O2_,xg_H2_)
    results3 = mb_HOXFe(xg_O2_,xg_H2_)
    return xg_O2_, xg_H2_, results1, results2, results3

#############
### CHOFe ###
#############

### insolubles not included

def eq_CHOFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses):
    P = PT["P"]
    wt_O = bulk_wf['O']
    wt_C = bulk_wf['C']
    wt_H = bulk_wf['H']
    wt_Fe = bulk_wf['Fe']
    guessx = guesses["guessx"]
    guessy = guesses["guessy"]
    
    # equilibrium constants
    K1_ = mdv.KHOg(PT,models)
    K2_ = mdv.KCOg(PT,models)
    K3_ = mdv.KCOHg(PT,models)
    K4_ = mdv.C_H2O(PT,melt_wf,models)
    K5_ = mdv.C_CO3(PT,melt_wf,models)
    KD1_, KD2, y = mdv.FefO2_KC91_EqA_terms(PT,melt_wf,models)
   
    # fugacity coefficients
    y_CO_ = mdv.y_CO(PT,models)
    y_O2_ = mdv.y_O2(PT,models)
    y_CO2_ = mdv.y_CO2(PT,models)
    y_CH4_ = mdv.y_CH4(PT,models)
    y_H2O_ = mdv.y_H2O(PT,models)
    y_H2_ = mdv.y_H2(PT,models)
    
    # molecular masses
    M_H = mdv.species.loc['H','M']
    M_C = mdv.species.loc['C','M']
    M_O = mdv.species.loc['O','M']
    M_Fe = mdv.species.loc['Fe','M']
    M_O2 = mdv.species.loc['O2','M']
    M_CO = mdv.species.loc['CO','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_H2 = mdv.species.loc['H2','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_CH4 = mdv.species.loc['CH4','M']
    M_FeO = mdv.species.loc['FeO','M']
    M_FeO15 = mdv.species.loc['FeO1.5','M']
    M_m_ = mg.M_m_SO(melt_wf)
    
    constants = [P, wt_O, wt_C, wt_H, wt_Fe, K1_, K2_, K3_, K4_, K5_, KD1_, KD2, y, y_H2_, y_O2_, y_H2O_, y_CO_, y_CO2_, y_CH4_, M_C, M_H, M_O, M_Fe, M_O2, M_CO, M_CO2, M_H2, M_H2O, M_CH4, M_FeO, M_FeO15, M_m_]
     
    def mg_CHOFe(xg_O2_,xg_CO_):
        xg_CO2_ = (K2_*y_CO_*xg_CO_*(y_O2_*xg_O2_*P)**0.5)/y_CO2_
        a = (y_CO2_*xg_CO2_*y_H2O_**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
        b = 1.0 + (y_H2O_/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5))
        c = xg_CO2_ + xg_CO_ + xg_O2_ - 1.0
        xg_H2O_ = (-b + (b**2.0-(4.0*a*c))**0.5)/(2.0*a)
        xg_H2_ = (y_H2O_*xg_H2O_)/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)
        xg_CH4_ = (y_CO2_*xg_CO2_*(y_H2O_*xg_H2O_)**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
        xm_H2O_ = (K4_*y_H2O_*xg_H2O_*P)**0.5
        xm_CO2_ = K5_*y_CO2_*xg_CO2_*P
        Xm_t = xm_CO2_*M_CO2 + xm_H2O_*M_H2O + (1.0-xm_CO2_-xm_H2O_)*M_m_
        Xg_t = xg_CO2_*M_CO2 + xg_CO_*M_CO + xg_O2_*M_O2 + xg_H2O_*M_H2O + xg_H2_*M_H2 + xg_CH4_*M_CH4
        Fe32 = ((KD1_*(y_O2_*xg_O2_*P)**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y)))
        Fe3T = Fe32/(1.0+Fe32)
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t
        wm_CO2_ = (xm_CO2_*M_CO2)/Xm_t
        return xg_CO2_, xg_H2O_, xg_H2_, xg_CH4_, xm_H2O_, xm_CO2_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2O_, wm_CO2_
    
    def mb_CHOFe(xg_O2_, xg_CO_):
        xg_CO2_, xg_H2O_, xg_H2_, xg_CH4_, xm_H2O_, xm_CO2_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2O_, wm_CO2_ = mg_CHOFe(xg_O2_, xg_CO_)
        mbCO, mbCH, wt_g_O, wt_g_C, wt_g_H = f_CHOFe(xg_O2_,xg_CO_)
        wt_g = (wt_g_O+wt_g_H+wt_g_C)/3.0
        wt_H_ = 2.0*M_H*((wt_g*(((xg_H2O_+xg_H2_+2.0*xg_CH4_)/Xg_t) - (xm_H2O_/Xm_t)))+(xm_H2O_/Xm_t))
        wt_O_ = M_O*((wt_g*(((2.0*xg_CO2_ + xg_CO_ + 2.0*xg_O2_ + xg_H2O_)/Xg_t) - (xm_H2O_/Xm_t) - ((2.0*xm_CO2_)/Xm_t)))+(xm_H2O_/Xm_t) + ((2.0*xm_CO2_)/Xm_t) + (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))
        wt_C_ = M_C*((wt_g*(((xg_CO2_+xg_CO_+xg_CH4_)/Xg_t) - (xm_CO2_/Xm_t))) + (xm_CO2_/Xm_t))
        return wt_g, wt_O_, wt_C_, wt_H_
    
    def f_CHOFe(xg_O2_,xg_CO_):
        xg_CO2_, xg_H2O_, xg_H2_, xg_CH4_, xm_H2O_, xm_CO2_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2O_, wm_CO2_ = mg_CHOFe(xg_O2_, xg_CO_)
        wt_g_C = ((wt_C/M_C) - (xm_CO2_/Xm_t))/(((xg_CO2_+xg_CO_+xg_CH4_)/Xg_t) - (xm_CO2_/Xm_t))
        wt_g_O = ((wt_O/M_O) - (xm_H2O_/Xm_t) - ((2.0*xm_CO2_)/Xm_t) - (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))/(((2.0*xg_CO2_ + xg_CO_ + 2.0*xg_O2_ + xg_H2O_)/Xg_t) - (xm_H2O_/Xm_t) - ((2.0*xm_CO2_)/Xm_t))
        wt_g_H = ((wt_H/(2.0*M_H)) - (xm_H2O_/Xm_t))/(((xg_H2O_+xg_H2_+2.0*xg_CH4_)/Xg_t) - (xm_H2O_/Xm_t))
        mbCO = (wt_g_C - wt_g_O)
        mbCH = (wt_g_C - wt_g_H)
        return mbCO, mbCH, wt_g_O, wt_g_C, wt_g_H
    
    def df_CHOFe(xg_O2_,xg_CO_,constants):
        dmbCO_O = de.CHOFe_CO_O2(xg_O2_,xg_CO_,constants)
        dmbCO_C = de.CHOFe_CO_CO(xg_O2_,xg_CO_,constants)
        dmbCH_O = de.CHOFe_CH_O2(xg_O2_,xg_CO_,constants)
        dmbCH_C = de.CHOFe_CH_CO(xg_O2_,xg_CO_,constants)
        return dmbCO_O, dmbCO_C, dmbCH_O, dmbCH_C
    
    xg_O2_, xg_CO_ = jac_newton(guessx,guessy,constants,f_CHOFe,df_CHOFe,nr_step,nr_tol)
    results1 = mg_CHOFe(xg_O2_,xg_CO_)
    results2 = f_CHOFe(xg_O2_,xg_CO_)
    results3 = mb_CHOFe(xg_O2_,xg_CO_)
    return xg_O2_, xg_CO_, results1, results2, results3

### insolubles included
def eq_CHOFe_2(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses):
    P = PT["P"]
    wt_O = bulk_wf['O']
    wt_C = bulk_wf['C']
    wt_H = bulk_wf['H']
    wt_Fe = bulk_wf['Fe']
    guessx = guesses["guessx"]
    guessy = guesses["guessy"]
    
    # equilibrium constants
    K1_ = mdv.KHOg(PT,models)
    K2_ = mdv.KCOg(PT,models)
    K3_ = mdv.KCOHg(PT,models)
    K4_ = mdv.C_H2O(PT,melt_wf,models)
    K5_ = mdv.C_CO3(PT,melt_wf,models)
    K6_ = mdv.C_H2(PT,melt_wf,models)/1000000.
    K7_ = mdv.C_CO(PT,melt_wf,models)/1000000.
    K8_ = mdv.C_CH4(PT,melt_wf,models)/1000000.
    KD1_, KD2, y = mdv.FefO2_KC91_EqA_terms(PT,melt_wf,models)
   
    # fugacity coefficients
    y_CO_ = mdv.y_CO(PT,models)
    y_O2_ = mdv.y_O2(PT,models)
    y_CO2_ = mdv.y_CO2(PT,models)
    y_CH4_ = mdv.y_CH4(PT,models)
    y_H2O_ = mdv.y_H2O(PT,models)
    y_H2_ = mdv.y_H2(PT,models)
    
    # molecular masses
    M_H = mdv.species.loc['H','M']
    M_C = mdv.species.loc['C','M']
    M_O = mdv.species.loc['O','M']
    M_Fe = mdv.species.loc['Fe','M']
    M_O2 = mdv.species.loc['O2','M']
    M_CO = mdv.species.loc['CO','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_H2 = mdv.species.loc['H2','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_CH4 = mdv.species.loc['CH4','M']
    M_FeO = mdv.species.loc['FeO','M']
    M_FeO15 = mdv.species.loc['FeO1.5','M']
    M_m_ = mg.M_m_SO(melt_wf)
    
    constants = [P, wt_O, wt_C, wt_H, wt_Fe, K1_, K2_, K3_, K4_, K5_, K6_, K7_, K8_, KD1_, KD2, y, y_H2_, y_O2_, y_H2O_, y_CO_, y_CO2_, y_CH4_, M_C, M_H, M_O, M_Fe, M_O2, M_CO, M_CO2, M_H2, M_H2O, M_CH4, M_FeO, M_FeO15, M_m_]
     
    def mg_CHOFe(xg_O2_,xg_CO_):
        xg_CO2_ = (K2_*y_CO_*xg_CO_*(y_O2_*xg_O2_*P)**0.5)/y_CO2_
        a = (y_CO2_*xg_CO2_*y_H2O_**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
        b = 1.0 + (y_H2O_/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5))
        c = xg_CO2_ + xg_CO_ + xg_O2_ - 1.0
        xg_H2O_ = (-b + (b**2.0-(4.0*a*c))**0.5)/(2.0*a)
        xg_H2_ = (y_H2O_*xg_H2O_)/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)
        xg_CH4_ = (y_CO2_*xg_CO2_*(y_H2O_*xg_H2O_)**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
        xm_H2O_ = (K4_*y_H2O_*xg_H2O_*P)**0.5
        xm_CO2_ = K5_*y_CO2_*xg_CO2_*P
        Xm_t = xm_CO2_*M_CO2 + xm_H2O_*M_H2O + (1.0-xm_CO2_-xm_H2O_)*M_m_
        Xg_t = xg_CO2_*M_CO2 + xg_CO_*M_CO + xg_O2_*M_O2 + xg_H2O_*M_H2O + xg_H2_*M_H2 + xg_CH4_*M_CH4
        Fe32 = ((KD1_*(y_O2_*xg_O2_*P)**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y)))
        Fe3T = Fe32/(1.0+Fe32)
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t
        wm_CO2_ = (xm_CO2_*M_CO2)/Xm_t
        wm_H2_ = K6_*y_H2_*xg_H2_*P
        wm_CO_ = K7_*y_CO_*xg_CO_*P
        wm_CH4_ = K8_*y_CH4_*xg_CH4_*P
        return xg_CO2_, xg_H2O_, xg_H2_, xg_CH4_, xm_H2O_, xm_CO2_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_
    
    def mb_CHOFe(xg_O2_, xg_CO_):
        xg_CO2_, xg_H2O_, xg_H2_, xg_CH4_, xm_H2O_, xm_CO2_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_ = mg_CHOFe(xg_O2_, xg_CO_)
        mbCO, mbCH, wt_g_O, wt_g_C, wt_g_H = f_CHOFe(xg_O2_,xg_CO_)
        wt_g = (wt_g_O+wt_g_H+wt_g_C)/3.0
        wt_H_ = 2.0*M_H*((wt_g*(((xg_H2O_+xg_H2_+2.0*xg_CH4_)/Xg_t) - (xm_H2O_/Xm_t) - (wm_H2_/M_H2) - ((2.*wm_CH4_)/M_CH4)))+(xm_H2O_/Xm_t) + (wm_H2_/M_H2) + ((2.*wm_CH4_)/M_CH4))
        wt_O_ = M_O*((wt_g*(((2.0*xg_CO2_ + xg_CO_ + 2.0*xg_O2_ + xg_H2O_)/Xg_t) - (xm_H2O_/Xm_t) - ((2.0*xm_CO2_)/Xm_t) - (wm_CO_/M_CO)))+(xm_H2O_/Xm_t) + ((2.0*xm_CO2_)/Xm_t) + (wm_CO_/M_CO) + (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))
        wt_C_ = M_C*((wt_g*(((xg_CO2_+xg_CO_+xg_CH4_)/Xg_t) - (xm_CO2_/Xm_t) - (wm_CO_/M_CO) - (wm_CH4_/M_CH4))) + (xm_CO2_/Xm_t) + (wm_CO_/M_CO) + (wm_CH4_/M_CH4))
        return wt_g, wt_O_, wt_C_, wt_H_
    
    def f_CHOFe(xg_O2_,xg_CO_):
        xg_CO2_, xg_H2O_, xg_H2_, xg_CH4_, xm_H2O_, xm_CO2_, Xm_t, Xg_t, Fe32, Fe3T, wm_H2O_, wm_CO2_, wm_H2_, wm_CO_, wm_CH4_ = mg_CHOFe(xg_O2_, xg_CO_)
        wt_g_C = ((wt_C/M_C) - (xm_CO2_/Xm_t) - (wm_CO_/M_CO) - (wm_CH4_/M_CH4))/(((xg_CO2_+xg_CO_+xg_CH4_)/Xg_t) - (xm_CO2_/Xm_t) - (wm_CO_/M_CO) - (wm_CH4_/M_CH4))
        wt_g_O = ((wt_O/M_O) - (xm_H2O_/Xm_t) - ((2.0*xm_CO2_)/Xm_t) - (wm_CO_/M_CO) - (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))/(((2.0*xg_CO2_ + xg_CO_ + 2.0*xg_O2_ + xg_H2O_)/Xg_t) - (xm_H2O_/Xm_t) - ((2.0*xm_CO2_)/Xm_t) -(wm_CO_/M_CO))
        wt_g_H = ((wt_H/(2.0*M_H)) - (xm_H2O_/Xm_t)  - (wm_H2_/M_H2) - ((2.*wm_CH4_)/M_CH4))/(((xg_H2O_+xg_H2_+2.0*xg_CH4_)/Xg_t) - (xm_H2O_/Xm_t) - (wm_H2_/M_H2) - ((2.*wm_CH4_)/M_CH4))
        mbCO = (wt_g_C - wt_g_O)
        mbCH = (wt_g_C - wt_g_H)
        return mbCO, mbCH, wt_g_O, wt_g_C, wt_g_H
    
    def df_CHOFe(xg_O2_,xg_CO_,constants):
        dmbCO_O = de.CHOFe_CO_O2_2(xg_O2_,xg_CO_,constants)
        dmbCO_C = de.CHOFe_CO_CO_2(xg_O2_,xg_CO_,constants)
        dmbCH_O = de.CHOFe_CH_O2_2(xg_O2_,xg_CO_,constants)
        dmbCH_C = de.CHOFe_CH_CO_2(xg_O2_,xg_CO_,constants)
        return dmbCO_O, dmbCO_C, dmbCH_O, dmbCH_C
    
    xg_O2_, xg_CO_ = jac_newton(guessx,guessy,constants,f_CHOFe,df_CHOFe,nr_step,nr_tol)
    results1 = mg_CHOFe(xg_O2_,xg_CO_)
    results2 = f_CHOFe(xg_O2_,xg_CO_)
    results3 = mb_CHOFe(xg_O2_,xg_CO_)
    return xg_O2_, xg_CO_, results1, results2, results3


#########################
### CH - H2O-CO2 only ###
#########################

def eq_CH(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses): # H2O
    P = PT["P"]
    wt_C = bulk_wf['C']
    wt_H = bulk_wf['H']
    guessx = guesses["guessx"]
    
    # equilibrium constants
    K4_ = mdv.C_H2O(PT,melt_wf,models)
    K5_ = mdv.C_CO3(PT,melt_wf,models)
   
    # fugacity coefficients
    y_CO2_ = mdv.y_CO2(PT,models)
    y_H2O_ = mdv.y_H2O(PT,models)
    
    # molecular masses
    M_H = mdv.species.loc['H','M']
    M_C = mdv.species.loc['C','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_m_ = mg.M_m_SO(melt_wf)
    
    constants = [P, wt_C, wt_H, K4_, K5_, y_H2O_, y_CO2_, M_H, M_C, M_CO2, M_H2O, M_m_]
     
    def mg_CH(xg_CO2_):
        xg_H2O_ = 1. - xg_CO2_
        xm_H2O_ = (K4_*y_H2O_*xg_H2O_*P)**0.5
        xm_CO2_ = K5_*y_CO2_*xg_CO2_*P
        Xm_t = xm_CO2_*M_CO2 + xm_H2O_*M_H2O + (1.0-xm_CO2_-xm_H2O_)*M_m_
        Xg_t = xg_CO2_*M_CO2 + xg_H2O_*M_H2O
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t
        wm_CO2_ = (xm_CO2_*M_CO2)/Xm_t
        return xg_CO2_, xg_H2O_, xm_H2O_, xm_CO2_, Xm_t, Xg_t, wm_H2O_, wm_CO2_
    
    def mb_CH(xg_CO2_):
        xg_CO2_, xg_H2O_, xm_H2O_, xm_CO2_, Xm_t, Xg_t, wm_H2O_, wm_CO2_ = mg_CH(xg_CO2_)
        mbCH, wt_g_C, wt_g_H = f_CH(xg_CO2_)
        wt_g = (wt_g_H+wt_g_C)/2.
        wt_H_ = 2.0*M_H*((wt_g*((xg_H2O_/Xg_t) - (xm_H2O_/Xm_t)))+(xm_H2O_/Xm_t))
        wt_C_ = M_C*((wt_g*((xg_CO2_/Xg_t) - (xm_CO2_/Xm_t))) + (xm_CO2_/Xm_t))
        return wt_g, "", wt_C_, wt_H_
    
    def f_CH(xg_CO2_):
        xg_CO2_, xg_H2O_, xm_H2O_, xm_CO2_, Xm_t, Xg_t, wm_H2O_, wm_CO2_ = mg_CH(xg_CO2_)
        wt_g_C = ((wt_C/M_C) - (xm_CO2_/Xm_t))/((xg_CO2_/Xg_t) - (xm_CO2_/Xm_t))
        wt_g_H = ((wt_H/(2.0*M_H)) - (xm_H2O_/Xm_t))/((xg_H2O_/Xg_t) - (xm_H2O_/Xm_t))
        mbCH = (wt_g_C - wt_g_H)
        return mbCH, wt_g_C, wt_g_H
    
    def df_CH(xg_CO2_,constants):
        dmbCH_C = de.CH_CO2(xg_CO2_,constants)
        return dmbCH_C
    
    xg_CO2_ = newton_raphson(guessx,constants,nr_tol,nr_step,f_CH,df_CH)
    result1 = mg_CH(xg_CO2_)
    result2 = f_CH(xg_CO2_)
    result3 = mb_CH(xg_CO2_)
    return xg_CO2_, result1, result2, result3

##############                        
### CHOXFe ###
##############
                        
def eq_CHOXFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses,solve_species):
    species_X = models.loc["species X","option"]
    P = PT["P"]
    wt_O = bulk_wf['O']
    wt_X = bulk_wf['X']
    wt_H = bulk_wf['H']
    wt_C = bulk_wf['C']
    wt_Fe = bulk_wf['Fe']
    guessx = guesses["guessx"]
    guessy = guesses["guessy"]
    guessz = guesses["guessz"]
        
    # equilibrium constants
    K1_ = mdv.KHOg(PT,models)
    K2_ = mdv.KCOg(PT,models)
    K3_ = mdv.KCOHg(PT,models)
    K4_ = mdv.C_H2O(PT,melt_wf,models)
    K5_ = mdv.C_CO3(PT,melt_wf,models)
    K6_ = mdv.C_X(PT,melt_wf,models)/1000000.0
    K11_ = mdv.C_CO(PT,melt_wf,models)/1000000.0
    K12_ = mdv.C_CH4(PT,melt_wf,models)/1000000.0
    K13_ = mdv.C_H2(PT,melt_wf,models)/1000000.0
    KD1_, KD2, y = mdv.FefO2_KC91_EqA_terms(PT,melt_wf,models)
   
    # fugacity coefficients
    y_CO_ = mdv.y_CO(PT,models)
    y_O2_ = mdv.y_O2(PT,models)
    y_CO2_ = mdv.y_CO2(PT,models)
    y_CH4_ = mdv.y_CH4(PT,models)
    y_H2O_ = mdv.y_H2O(PT,models)
    y_H2_ = mdv.y_H2(PT,models)
    y_X_ = mdv.y_X(PT,models)
    
    # molecular masses
    M_H = mdv.species.loc['H','M']
    M_C = mdv.species.loc['C','M']
    M_O = mdv.species.loc['O','M']
    M_Fe = mdv.species.loc['Fe','M']
    M_O2 = mdv.species.loc['O2','M']
    M_CO = mdv.species.loc['CO','M']
    M_OH = mdv.species.loc['OH','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_H2 = mdv.species.loc['H2','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_CH4 = mdv.species.loc['CH4','M']
    M_FeO = mdv.species.loc['FeO','M']
    M_FeO15 = mdv.species.loc['FeO1.5','M']
    M_X = mdv.species.loc[species_X,'M']                 
    M_m_ = mg.M_m_SO(melt_wf)
    M_m_ox = mg.M_m_ox(melt_wf,models)
    
    constants = [P, wt_O, wt_C, wt_H, wt_X, wt_Fe, K1_, K2_, K3_, K4_, K5_, K6_, K11_, K12_, K13_, KD1_, KD2, y, y_H2_, y_O2_, y_H2O_, y_CO_, y_CO2_, y_CH4_, y_X_, M_C, M_H, M_O, M_X, M_Fe, M_O2, M_CO, M_CO2, M_OH, M_H2, M_H2O, M_CH4, M_FeO, M_FeO15, M_m_, M_m_ox]
       
    def mg_CHOXFe(xg_O2_,xg_CO_,xg_X_):
        xg_O2__ = xg_O2_
        xg_CO__ = xg_CO_
        xg_X__ = xg_X_                  
        xg_CO2_ = (K2_*y_CO_*xg_CO_*(y_O2_*xg_O2_*P)**0.5)/y_CO2_
        a = (y_CO2_*xg_CO2_*y_H2O_**2.)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.)
        b = 1. + (y_H2O_/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5))
        c = xg_CO2_ + xg_CO_ + xg_O2_ + xg_X_ - 1.
        xg_H2O_ = (-b + (b**2.0-(4.0*a*c))**0.5)/(2.0*a)
        xg_H2_ = (y_H2O_*xg_H2O_)/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)
        xg_CH4_ = (y_CO2_*xg_CO2_*(y_H2O_*xg_H2O_)**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
        Xg_t = xg_CO2_*M_CO2 + xg_CO_*M_CO + xg_O2_*M_O2 + xg_H2O_*M_H2O + xg_H2_*M_H2 + xg_CH4_*M_CH4 + xg_X_*M_X
        # melt composition
        xm_H2O_ = (K4_*y_H2O_*xg_H2O_*P)**0.5
        xm_CO2_ = K5_*y_CO2_*xg_CO2_*P
        Xm_t = xm_CO2_*M_CO2 + xm_H2O_*M_H2O + (1.0-xm_CO2_-xm_H2O_)*M_m_
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t
        wm_CO2_ = (xm_CO2_*M_CO2)/Xm_t
        wm_CO_ = K11_*y_CO_*xg_CO_*P
        wm_CH4_ = K12_*y_CH4_*xg_CH4_*P
        wm_H2_ = K13_*y_H2_*xg_H2_*P
        wm_X_ = K6_*y_X_*xg_X_*P                
        Fe32 = ((KD1_*(y_O2_*xg_O2_*P)**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y)))
        Fe3T = Fe32/(1.0+Fe32)
        Xm_t_ox = ""
        return xg_O2__, xg_H2_, xg_X__, xg_H2O_, xg_CO__, xg_CO2_, xg_CH4_, Xg_t, xm_H2O_, xm_CO2_, wm_X_, Xm_t, Xm_t_ox, Fe32, Fe3T, wm_H2O_, wm_CO2_, wm_H2_, wm_CH4_, wm_CO_    
    
    def mb_CHOXFe(xg_O2_,xg_CO_,xg_X_):
        xg_O2__, xg_H2_, xg_X__, xg_H2O_, xg_CO__, xg_CO2_, xg_CH4_, Xg_t, xm_H2O_, xm_CO2_, wm_X_, Xm_t, Xm_t_ox, Fe32, Fe3T, wm_H2O_, wm_CO2_, wm_H2_, wm_CH4_, wm_CO_ = mg_CHOXFe(xg_O2_, xg_CO_, xg_X_)
        mba, mbb, mbc, wt_g_O, wt_g_C, wt_g_H, wt_g_X = f_CHOXFe(xg_O2_, xg_CO_, xg_X_)
        wt_g = (wt_g_O+wt_g_H+wt_g_C+wt_g_X)/4.
        wt_H_ = 2.*M_H*((wt_g*(((xg_H2O_+xg_H2_+2.*xg_CH4_)/Xg_t) - (wm_H2O_/M_H2O) - (wm_H2_/M_H2) - (2.*wm_CH4_/M_CH4))) + (wm_H2O_/M_H2O) + (wm_H2_/M_H2) + (2.*wm_CH4_/M_CH4))
        wt_O_ = M_O*((wt_g*(((2.0*xg_CO2_ + xg_CO_ + 2.0*xg_O2_ + xg_H2O_)/Xg_t) - (wm_H2O_/M_H2O) - ((2.0*wm_CO2_)/M_CO2) - (wm_CO_/M_CO))) + (wm_H2O_/M_H2O) + ((2.0*wm_CO2_)/M_CO2) + (wm_CO_/M_CO) + (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))
        wt_C_ = M_C*((wt_g*(((xg_CO2_+xg_CO_+xg_CH4_)/Xg_t) - (wm_CO2_/M_CO2) - (wm_CH4_/M_CH4) - (wm_CO_/M_CO))) + (wm_CO2_/M_CO2) + (wm_CH4_/M_CH4) + (wm_CO_/M_CO))
        wt_X_ = M_X*((wt_g*(((xg_X_)/Xg_t) - (wm_X_/M_X))) + (wm_X_/M_X))
        return wt_g, wt_O_, wt_C_, wt_H_, wt_X_
    
    def f_CHOXFe(xg_O2_,xg_CO_,xg_X_):
        xg_O2__, xg_H2_, xg_X__, xg_H2O_, xg_CO__, xg_CO2_, xg_CH4_, Xg_t, xm_H2O_, xm_CO2_, wm_X_, Xm_t, Xm_t_ox, Fe32, Fe3T, wm_H2O_, wm_CO2_, wm_H2_, wm_CH4_, wm_CO_ = mg_CHOXFe(xg_O2_,xg_CO_,xg_X_)
        wt_g_C = ((wt_C/M_C) - (wm_CO2_/M_CO2) - (wm_CH4_/M_CH4) - (wm_CO_/M_CO))/(((xg_CO2_+xg_CO_+xg_CH4_)/Xg_t) - (wm_CO2_/M_CO2) - (wm_CH4_/M_CH4) - (wm_CO_/M_CO))
        wt_g_O = ((wt_O/M_O) - (wm_H2O_/M_H2O) - ((2.0*wm_CO2_)/M_CO2) - (wm_CO_/M_CO) - (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))/(((2.0*xg_CO2_ + xg_CO_ + 2.0*xg_O2_ + xg_H2O_)/Xg_t) - (wm_H2O_/M_H2O) - ((2.0*wm_CO2_)/M_CO2) - (wm_CO_/M_CO))
        wt_g_H = ((wt_H/(2.0*M_H)) - (wm_H2O_/M_H2O) - (wm_H2_/M_H2) - (2.*wm_CH4_/M_CH4))/(((xg_H2O_+xg_H2_+2.0*xg_CH4_)/Xg_t) - (wm_H2O_/M_H2O) - (wm_H2_/M_H2) - (2.*wm_CH4_/M_CH4))
        wt_g_X = ((wt_X/M_X) - (wm_X_/M_X))/(((xg_X_)/Xg_t) - (wm_X_/M_X))
        mba = (wt_g_C - wt_g_O) # mbCO
        mbb = (wt_g_C - wt_g_H) # mbCH
        mbc = (wt_g_C - wt_g_X) # mbCS
        return mba, mbb, mbc, wt_g_O, wt_g_C, wt_g_H, wt_g_X
    
    def df_CHOXFe(xg_O2_,xg_CO_,xg_X_,constants): 
        dmba_O = de.CHOXFe_CO_O2(xg_O2_,xg_CO_,xg_X_,constants) # dmbCO_O
        dmba_A = de.CHOXFe_CO_CO(xg_O2_,xg_CO_,xg_X_,constants) # dmbCO_CO
        dmba_B = de.CHOXFe_CO_X(xg_O2_,xg_CO_,xg_X_,constants) # dmbCO_S 
        dmbb_O = de.CHOXFe_CH_O2(xg_O2_,xg_CO_,xg_X_,constants) # dmbCH_O 
        dmbb_A = de.CHOXFe_CH_CO(xg_O2_,xg_CO_,xg_X_,constants) # dmbCH_CO
        dmbb_B = de.CHOXFe_CH_X(xg_O2_,xg_CO_,xg_X_,constants) # dmbCH_S
        dmbc_O = de.CHOXFe_CX_O2(xg_O2_,xg_CO_,xg_X_,constants) # dmbCS_O 
        dmbc_A = de.CHOXFe_CX_CO(xg_O2_,xg_CO_,xg_X_,constants) # dmbCS_CO
        dmbc_B = de.CHOXFe_CX_X(xg_O2_,xg_CO_,xg_X_,constants) # dmbCS_S
        return dmba_O, dmba_A, dmba_B, dmbb_O, dmbb_A, dmbb_B, dmbc_O, dmbc_A, dmbc_B

    xg_O2_, xg_CO_, xg_X_ = jac_newton3(guessx,guessy,guessz,constants,f_CHOXFe,df_CHOXFe,nr_step,nr_tol)
    results1 = xg_O2_, xg_CO_, xg_X_
    results2 = mg_CHOXFe(xg_O2_,xg_CO_,xg_X_)
    results3 = f_CHOXFe(xg_O2_,xg_CO_,xg_X_)
    results4 = mb_CHOXFe(xg_O2_,xg_CO_,xg_X_)
    return results1, results2, results3, results4
                        
                        
##############
### SCHOFe ###
##############

# all gas species - only oxidised melt species - no H2S #

def eq_SCHOFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses,solve_species):
    P = PT["P"]
    wt_O = bulk_wf['O']
    wt_S = bulk_wf['S']
    wt_H = bulk_wf['H']
    wt_C = bulk_wf['C']
    wt_Fe = bulk_wf['Fe']
    guessx = guesses["guessx"]
    guessy = guesses["guessy"]
    guessz = guesses["guessz"]
        
    # equilibrium constants
    K1_ = mdv.KHOg(PT,models)
    K2_ = mdv.KCOg(PT,models)
    K3_ = mdv.KCOHg(PT,models)
    K4_ = mdv.C_H2O(PT,melt_wf,models)
    K5_ = mdv.C_CO3(PT,melt_wf,models)
    K6_ = mdv.KOSg(PT,models)
    K7_ = mdv.KHOSg(PT,models)
    K8_ = mdv.C_S(PT,melt_wf,models)/1000000.0
    K9_ = (mdv.C_SO4(PT,melt_wf,models)/1000000.0)
    K10_ = mdv.KOCSg(PT,models)
    KD1_, KD2, y = mdv.FefO2_KC91_EqA_terms(PT,melt_wf,models)
   
    # fugacity coefficients
    y_CO_ = mdv.y_CO(PT,models)
    y_O2_ = mdv.y_O2(PT,models)
    y_CO2_ = mdv.y_CO2(PT,models)
    y_CH4_ = mdv.y_CH4(PT,models)
    y_H2O_ = mdv.y_H2O(PT,models)
    y_H2_ = mdv.y_H2(PT,models)
    y_S2_ = mdv.y_S2(PT,models)
    y_SO2_ = mdv.y_SO2(PT,models)
    y_H2S_ = mdv.y_H2S(PT,models)
    y_OCS_ = mdv.y_OCS(PT,models)
    
    # molecular masses
    M_H = mdv.species.loc['H','M']
    M_C = mdv.species.loc['C','M']
    M_O = mdv.species.loc['O','M']
    M_S = mdv.species.loc['S','M']
    M_Fe = mdv.species.loc['Fe','M']
    M_O2 = mdv.species.loc['O2','M']
    M_CO = mdv.species.loc['CO','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_H2 = mdv.species.loc['H2','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_CH4 = mdv.species.loc['CH4','M']
    M_S2 = mdv.species.loc['S2','M']
    M_SO2 = mdv.species.loc['SO2','M']
    M_SO3 = mdv.species.loc['SO3','M']
    M_H2S = mdv.species.loc['H2S','M']
    M_FeO = mdv.species.loc['FeO','M']
    M_FeO15 = mdv.species.loc['FeO1.5','M']
    M_OCS = mdv.species.loc['OCS','M']
    M_m_ = mg.M_m_SO(melt_wf)
    
    constants = [P, wt_O, wt_C, wt_H, wt_S, wt_Fe, K1_, K2_, K3_, K4_, K5_, K6_, K7_, K8_, K9_, K10_, KD1_, KD2, y, y_H2_, y_O2_, y_H2O_, y_CO_, y_CO2_, y_CH4_, y_SO2_, y_S2_, y_H2S_, y_OCS_, M_C, M_H, M_O, M_S, M_Fe, M_O2, M_CO, M_CO2, M_H2, M_H2O, M_CH4, M_S2, M_SO2, M_SO3, M_H2S, M_OCS, M_FeO, M_FeO15, M_m_]
    
    def mg_SCHOFe_OCS(xg_O2_,xg_CO_,xg_S2_):
        xg_CO2_ = (K2_*y_CO_*xg_CO_*(y_O2_*xg_O2_*P)**0.5)/y_CO2_
        xg_SO2_ = (K6_*y_O2_*xg_O2_*(y_S2_*xg_S2_*P)**0.5)/y_SO2_
        xg_OCS_ = ((xg_CO_*y_CO_)**3.0*xg_SO2_*y_SO2_*P)/(y_OCS_*(xg_CO2_*y_CO2_)**2.0*K10_)
        a = (y_CO2_*xg_CO2_*y_H2O_**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
        b = 1.0 + (y_H2O_/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)) + ((K7_*(y_S2_*xg_S2_)**0.5*y_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5))
        c = xg_CO2_ + xg_CO_ + xg_O2_ + xg_S2_ + xg_SO2_ + xg_OCS_ - 1.0
        xg_H2O_ = (-b + (b**2.0-(4.0*a*c))**0.5)/(2.0*a)
        xg_H2_ = (y_H2O_*xg_H2O_)/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)
        xg_CH4_ = (y_CO2_*xg_CO2_*(y_H2O_*xg_H2O_)**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
        xg_H2S_ = (K7_*(y_S2_*xg_S2_)**0.5*y_H2O_*xg_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5)
        Xg_t = xg_CO2_*M_CO2 + xg_CO_*M_CO + xg_O2_*M_O2 + xg_H2O_*M_H2O + xg_H2_*M_H2 + xg_CH4_*M_CH4 + xg_SO2_*M_SO2 + xg_S2_*M_S2 + xg_H2S_*M_H2S + xg_OCS_*M_OCS
        xm_H2O_ = (K4_*y_H2O_*xg_H2O_*P)**0.5
        xm_CO2_ = K5_*y_CO2_*xg_CO2_*P
        wm_S_ = K8_*((y_S2_*xg_S2_)/(y_O2_*xg_O2_))**0.5
        wm_SO3_ = ((K9_*(y_S2_*xg_S2_)**0.5*(y_O2_*xg_O2_)**1.5*P**2.0)/M_S)*M_SO3
        wm_ST_ = wm_S_ + ((M_S*wm_SO3_)/M_SO3)
        Xm_t = xm_CO2_*M_CO2 + xm_H2O_*M_H2O + (1.0-xm_CO2_-xm_H2O_)*M_m_
        Fe32 = ((KD1_*(y_O2_*xg_O2_*P)**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y)))
        Fe3T = Fe32/(1.0+Fe32)
        S62 = (wm_SO3_/M_SO3)/(wm_S_/M_S)
        S6T = S62/(1.0+S62)
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t
        wm_CO2_ = (xm_CO2_*M_CO2)/Xm_t
        return xg_CO2_, xg_H2O_, xg_H2_, xg_CH4_, xg_SO2_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_

    def mg_SCHOFe_OHS(xg_O2_,xg_H2_,xg_S2_):
        xg_H2O_ = (xg_H2_*K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)/y_H2O_
        xg_SO2_ = (K6_*y_O2_*xg_O2_*(y_S2_*xg_S2_*P)**0.5)/y_SO2_
        xg_H2S_ = (K7_*(y_S2_*xg_S2_)**0.5*y_H2O_*xg_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5)
        xg_CO2_ = (1.0 - xg_S2_ - xg_O2_ - xg_H2_ - xg_H2O_ - xg_H2S_ - xg_SO2_)/(1.0 + (y_CO2_/(K2_*y_CO_*(y_O2_*xg_O2_*P)**0.5)) + ((y_CO2_*y_SO2_*xg_SO2_)/(K2_**3.0*K10_*y_OCS_*(y_O2_*xg_O2_)**1.5*P**0.5)) + ((y_CO2_*(y_H2O_*xg_H2O_)**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)))
        xg_CO_ = (y_CO2_*xg_CO2_)/(K2_*y_CO_*(y_O2_*xg_O2_*P)**0.5)
        xg_CH4_ = (y_CO2_*xg_CO2_*(y_H2O_*xg_H2O_)**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
        xg_OCS_ = ((xg_CO_*y_CO_)**3.0*xg_SO2_*y_SO2_*P)/(y_OCS_*(xg_CO2_*y_CO2_)**2.0*K10_)
        Xg_t = xg_CO2_*M_CO2 + xg_CO_*M_CO + xg_O2_*M_O2 + xg_H2O_*M_H2O + xg_H2_*M_H2 + xg_CH4_*M_CH4 + xg_SO2_*M_SO2 + xg_S2_*M_S2 + xg_H2S_*M_H2S + xg_OCS_*M_OCS
        xm_H2O_ = (K4_*y_H2O_*xg_H2O_*P)**0.5
        xm_CO2_ = K5_*y_CO2_*xg_CO2_*P
        wm_S_ = K8_*((y_S2_*xg_S2_)/(y_O2_*xg_O2_))**0.5
        wm_SO3_ = ((K9_*(y_S2_*xg_S2_)**0.5*(y_O2_*xg_O2_)**1.5*P**2.0)/M_S)*M_SO3
        wm_ST_ = wm_S_ + ((M_S*wm_SO3_)/M_SO3)
        Xm_t = xm_CO2_*M_CO2 + xm_H2O_*M_H2O + (1.0-xm_CO2_-xm_H2O_)*M_m_
        Fe32 = ((KD1_*(y_O2_*xg_O2_*P)**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y)))
        Fe3T = Fe32/(1.0+Fe32)
        S62 = (wm_SO3_/M_SO3)/(wm_S_/M_S)
        S6T = S62/(1.0+S62)
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t
        wm_CO2_ = (xm_CO2_*M_CO2)/Xm_t
        return xg_CO_, xg_CO2_, xg_H2O_, xg_CH4_, xg_SO2_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_    

    def mg_SCHOFe_OCH(xg_O2_,xg_CO_,xg_H2_):
        xg_H2O_ = (xg_H2_*K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)/y_H2O_
        xg_CO2_ = (K2_*y_CO_*xg_CO_*(y_O2_*xg_O2_*P)**0.5)/y_CO2_
        xg_CH4_ = (y_CO2_*xg_CO2_*(y_H2O_*xg_H2O_)**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
        a = 1.0
        b = ((K6_*(y_S2_*P)**0.5*y_O2_*xg_O2_)/y_SO2_) + (((K7_*y_H2O_*xg_H2O_)/y_H2S_)*(y_S2_/(y_O2_*xg_O2_))**0.5) + (((y_CO_*xg_CO_)**3.0*P**1.5*K6_*y_S2_**0.5*y_O2_*xg_O2_)/(K10_*y_OCS_*(y_CO2_*xg_CO2_)**2.0))
        c = xg_O2_ + xg_CO_ + xg_CO2_ + xg_H2_ + xg_H2O_ + xg_CH4_ - 1.0
        xg_S2_ = ((-b + (b**2.0-(4.0*a*c))**0.5)/(2.0*a))**2.0
        xg_SO2_ = (K6_*y_O2_*xg_O2_*(y_S2_*xg_S2_*P)**0.5)/y_SO2_
        xg_H2S_ = (K7_*(y_S2_*xg_S2_)**0.5*y_H2O_*xg_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5)
        xg_OCS_ = ((xg_CO_*y_CO_)**3.0*xg_SO2_*y_SO2_*P)/(y_OCS_*(xg_CO2_*y_CO2_)**2.0*K10_)
        Xg_t = xg_CO2_*M_CO2 + xg_CO_*M_CO + xg_O2_*M_O2 + xg_H2O_*M_H2O + xg_H2_*M_H2 + xg_CH4_*M_CH4 + xg_SO2_*M_SO2 + xg_S2_*M_S2 + xg_H2S_*M_H2S + xg_OCS_*M_OCS
        xm_H2O_ = (K4_*y_H2O_*xg_H2O_*P)**0.5
        xm_CO2_ = K5_*y_CO2_*xg_CO2_*P
        wm_S_ = K8_*((y_S2_*xg_S2_)/(y_O2_*xg_O2_))**0.5
        wm_SO3_ = ((K9_*(y_S2_*xg_S2_)**0.5*(y_O2_*xg_O2_)**1.5*P**2.0)/M_S)*M_SO3
        wm_ST_ = wm_S_ + ((M_S*wm_SO3_)/M_SO3)
        Xm_t = xm_CO2_*M_CO2 + xm_H2O_*M_H2O + (1.0-xm_CO2_-xm_H2O_)*M_m_
        Fe32 = ((KD1_*(y_O2_*xg_O2_*P)**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y)))
        Fe3T = Fe32/(1.0+Fe32)
        S62 = (wm_SO3_/M_SO3)/(wm_S_/M_S)
        S6T = S62/(1.0+S62)
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t
        wm_CO2_ = (xm_CO2_*M_CO2)/Xm_t
        return xg_CO2_, xg_H2O_, xg_CH4_, xg_S2_, xg_SO2_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_    
        
    def mb_SCHOFe(xg_O2_,xg_CO_,xg_H2_,xg_S2_):
        if solve_species == "OCS":
            xg_CO2_, xg_H2O_, xg_H2_, xg_CH4_, xg_SO2_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_ = mg_SCHOFe_OCS(xg_O2_, xg_CO_, xg_S2_)
            mbCO, mbCH, mbCS, wt_g_O, wt_g_C, wt_g_H, wt_g_S = f_SCHOFe_OCS(xg_O2_,xg_CO_,xg_S2_)
        elif solve_species == "OHS":
            xg_CO_, xg_CO2_, xg_H2O_, xg_CH4_, xg_SO2_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_ = mg_SCHOFe_OHS(xg_O2_, xg_H2_, xg_S2_)
            mbHO, mbHC, mbHS, wt_g_O, wt_g_C, wt_g_H, wt_g_S = f_SCHOFe_OHS(xg_O2_,xg_H2_,xg_S2_)
        elif solve_species == "OCH":
            xg_CO2_, xg_H2O_, xg_CH4_, xg_S2_, xg_SO2_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_ = mg_SCHOFe_OCH(xg_O2_, xg_CO_, xg_H2_)
            mbCO, mbCH, mbCS, wt_g_O, wt_g_C, wt_g_H, wt_g_S = f_SCHOFe_OCH(xg_O2_,xg_CO_,xg_H2_)
        wt_g = (wt_g_O+wt_g_H+wt_g_C+wt_g_S)/4.0
        wt_H_ = 2.0*M_H*((wt_g*(((xg_H2O_+xg_H2_+2.0*xg_CH4_+xg_H2S_)/Xg_t) - (xm_H2O_/Xm_t)))+(xm_H2O_/Xm_t))
        wt_O_ = M_O*((wt_g*(((2.0*xg_CO2_ + xg_CO_ + 2.0*xg_O2_ + xg_H2O_ + 2.0*xg_SO2_ + xg_OCS_)/Xg_t) - (xm_H2O_/Xm_t) - ((2.0*xm_CO2_)/Xm_t) - (3.0*wm_SO3_/M_SO3)))+(xm_H2O_/Xm_t) + ((2.0*xm_CO2_)/Xm_t) + (3.0*wm_SO3_/M_SO3) + (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))
        wt_C_ = M_C*((wt_g*(((xg_CO2_+xg_CO_+xg_CH4_+xg_OCS_)/Xg_t) - (xm_CO2_/Xm_t))) + (xm_CO2_/Xm_t))
        wt_S_ = M_S*((wt_g*(((xg_SO2_+2.0*xg_S2_+xg_H2S_+xg_OCS_)/Xg_t) - (wm_S_/M_S) - (wm_SO3_/M_SO3))) + (wm_S_/M_S) + (wm_SO3_/M_SO3))
        return wt_g, wt_O_, wt_C_, wt_H_, wt_S_
    
    def f_SCHOFe_OCS(xg_O2_,xg_CO_,xg_S2_):
        xg_CO2_, xg_H2O_, xg_H2_, xg_CH4_, xg_SO2_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_ = mg_SCHOFe_OCS(xg_O2_, xg_CO_, xg_S2_)
        wt_g_C = ((wt_C/M_C) - (xm_CO2_/Xm_t))/(((xg_CO2_+xg_CO_+xg_CH4_+xg_OCS_)/Xg_t) - (xm_CO2_/Xm_t))
        wt_g_O = ((wt_O/M_O) - (xm_H2O_/Xm_t) - ((2.0*xm_CO2_)/Xm_t) - (3.0*wm_SO3_/M_SO3) - (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))/(((2.0*xg_CO2_ + xg_CO_ + 2.0*xg_O2_ + xg_H2O_ + 2.0*xg_SO2_ + xg_OCS_)/Xg_t) - (xm_H2O_/Xm_t) - ((2.0*xm_CO2_)/Xm_t) - (3.0*wm_SO3_/M_SO3))
        wt_g_H = ((wt_H/(2.0*M_H)) - (xm_H2O_/Xm_t))/(((xg_H2O_+xg_H2_+2.0*xg_CH4_+xg_H2S_)/Xg_t) - (xm_H2O_/Xm_t))
        wt_g_S = ((wt_S/M_S) - (wm_S_/M_S) - (wm_SO3_/M_SO3))/(((xg_SO2_+2.0*xg_S2_+xg_H2S_+xg_OCS_)/Xg_t) - (wm_S_/M_S) - (wm_SO3_/M_SO3))
        mbCO = (wt_g_C - wt_g_O)
        mbCH = (wt_g_C - wt_g_H)
        mbCS = (wt_g_C - wt_g_S)
        return mbCO, mbCH, mbCS, wt_g_O, wt_g_C, wt_g_H, wt_g_S

    def f_SCHOFe_OHS(xg_O2_,xg_H2_,xg_S2_):
        xg_CO_, xg_CO2_, xg_H2O_, xg_CH4_, xg_SO2_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_ = mg_SCHOFe_OHS(xg_O2_, xg_H2_, xg_S2_)
        wt_g_C = ((wt_C/M_C) - (xm_CO2_/Xm_t))/(((xg_CO2_+xg_CO_+xg_CH4_+xg_OCS_)/Xg_t) - (xm_CO2_/Xm_t))
        wt_g_O = ((wt_O/M_O) - (xm_H2O_/Xm_t) - ((2.0*xm_CO2_)/Xm_t) - (3.0*wm_SO3_/M_SO3) - (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))/(((2.0*xg_CO2_ + xg_CO_ + 2.0*xg_O2_ + xg_H2O_ + 2.0*xg_SO2_ + xg_OCS_)/Xg_t) - (xm_H2O_/Xm_t) - ((2.0*xm_CO2_)/Xm_t) - (3.0*wm_SO3_/M_SO3))
        wt_g_H = ((wt_H/(2.0*M_H)) - (xm_H2O_/Xm_t))/(((xg_H2O_+xg_H2_+2.0*xg_CH4_+xg_H2S_)/Xg_t) - (xm_H2O_/Xm_t))
        wt_g_S = ((wt_S/M_S) - (wm_S_/M_S) - (wm_SO3_/M_SO3))/(((xg_SO2_+2.0*xg_S2_+xg_H2S_+xg_OCS_)/Xg_t) - (wm_S_/M_S) - (wm_SO3_/M_SO3))
        mbHO = (wt_g_H - wt_g_O)
        mbHC = (wt_g_H - wt_g_C)
        mbHS = (wt_g_H - wt_g_S)
        return mbHO, mbHC, mbHS, wt_g_O, wt_g_C, wt_g_H, wt_g_S 

    def f_SCHOFe_OCH(xg_O2_,xg_CO_,xg_H2_):
        xg_CO2_, xg_H2O_, xg_CH4_, xg_S2_, xg_SO2_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_ = mg_SCHOFe_OCH(xg_O2_, xg_CO_, xg_H2_)
        wt_g_C = ((wt_C/M_C) - (xm_CO2_/Xm_t))/(((xg_CO2_+xg_CO_+xg_CH4_+xg_OCS_)/Xg_t) - (xm_CO2_/Xm_t))
        wt_g_O = ((wt_O/M_O) - (xm_H2O_/Xm_t) - ((2.0*xm_CO2_)/Xm_t) - (3.0*wm_SO3_/M_SO3) - (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))/(((2.0*xg_CO2_ + xg_CO_ + 2.0*xg_O2_ + xg_H2O_ + 2.0*xg_SO2_ + xg_OCS_)/Xg_t) - (xm_H2O_/Xm_t) - ((2.0*xm_CO2_)/Xm_t) - (3.0*wm_SO3_/M_SO3))
        wt_g_H = ((wt_H/(2.0*M_H)) - (xm_H2O_/Xm_t))/(((xg_H2O_+xg_H2_+2.0*xg_CH4_+xg_H2S_)/Xg_t) - (xm_H2O_/Xm_t))
        wt_g_S = ((wt_S/M_S) - (wm_S_/M_S) - (wm_SO3_/M_SO3))/(((xg_SO2_+2.0*xg_S2_+xg_H2S_+xg_OCS_)/Xg_t) - (wm_S_/M_S) - (wm_SO3_/M_SO3))
        mbCO = (wt_g_C - wt_g_O)
        mbCH = (wt_g_C - wt_g_H)
        mbCS = (wt_g_C - wt_g_S)
        return mbCO, mbCH, mbCS, wt_g_O, wt_g_C, wt_g_H, wt_g_S
    
    def df_SCHOFe_OCS(xg_O2_,xg_CO_,xg_S2_,constants): 
        dmbCO_O = de.SCHOFe_OCS_CO_O2(xg_O2_,xg_CO_,xg_S2_,constants)
        dmbCO_CO = de.SCHOFe_OCS_CO_CO(xg_O2_,xg_CO_,xg_S2_,constants) 
        dmbCO_S = de.SCHOFe_OCS_CO_S2(xg_O2_,xg_CO_,xg_S2_,constants)  
        dmbCH_O = de.SCHOFe_OCS_CH_O2(xg_O2_,xg_CO_,xg_S2_,constants)
        dmbCH_CO = de.SCHOFe_OCS_CH_CO(xg_O2_,xg_CO_,xg_S2_,constants)
        dmbCH_S = de.SCHOFe_OCS_CH_S2(xg_O2_,xg_CO_,xg_S2_,constants) 
        dmbCS_O = de.SCHOFe_OCS_CS_O2(xg_O2_,xg_CO_,xg_S2_,constants) 
        dmbCS_CO = de.SCHOFe_OCS_CS_CO(xg_O2_,xg_CO_,xg_S2_,constants) 
        dmbCS_S = de.SCHOFe_OCS_CS_S2(xg_O2_,xg_CO_,xg_S2_,constants)
        return dmbCO_O, dmbCO_CO, dmbCO_S, dmbCH_O, dmbCH_CO, dmbCH_S, dmbCS_O, dmbCS_CO, dmbCS_S
    
    def df_SCHOFe_OHS(xg_O2_,xg_H2_,xg_S2_,constants):     
        dmbHO_O = de.SCHOFe_OHS_HO_O2(xg_O2_,xg_H2_,xg_S2_,constants)
        dmbHO_H = de.SCHOFe_OHS_HO_H2(xg_O2_,xg_H2_,xg_S2_,constants)
        dmbHO_S = de.SCHOFe_OHS_HO_S2(xg_O2_,xg_H2_,xg_S2_,constants)
        dmbHC_O = de.SCHOFe_OHS_HC_O2(xg_O2_,xg_H2_,xg_S2_,constants)
        dmbHC_H = de.SCHOFe_OHS_HC_H2(xg_O2_,xg_H2_,xg_S2_,constants)
        dmbHC_S = de.SCHOFe_OHS_HC_S2(xg_O2_,xg_H2_,xg_S2_,constants)
        dmbHS_O = de.SCHOFe_OHS_HS_O2(xg_O2_,xg_H2_,xg_S2_,constants)
        dmbHS_H = de.SCHOFe_OHS_HS_H2(xg_O2_,xg_H2_,xg_S2_,constants)
        dmbHS_S = de.SCHOFe_OHS_HS_S2(xg_O2_,xg_H2_,xg_S2_,constants)
        return dmbHO_O, dmbHO_H, dmbHO_S, dmbHC_O, dmbHC_H, dmbHC_S, dmbHS_O, dmbHS_H, dmbHS_S
           
    def df_SCHOFe_OCH(xg_O2_,xg_CO_,xg_H2_,constants):
        dmbCO_O = de.SCHOFe_OCH_CO_O2(xg_O2_,xg_CO_,xg_H2_,constants)
        dmbCO_CO = de.SCHOFe_OCH_CO_CO(xg_O2_,xg_CO_,xg_H2_,constants)
        dmbCO_H = de.SCHOFe_OCH_CO_H2(xg_O2_,xg_CO_,xg_H2_,constants)
        dmbCH_O = de.SCHOFe_OCH_CH_O2(xg_O2_,xg_CO_,xg_H2_,constants)
        dmbCH_CO = de.SCHOFe_OCH_CH_CO(xg_O2_,xg_CO_,xg_H2_,constants)
        dmbCH_H = de.SCHOFe_OCH_CH_H2(xg_O2_,xg_CO_,xg_H2_,constants)
        dmbCS_O = de.SCHOFe_OCH_CS_O2(xg_O2_,xg_CO_,xg_H2_,constants)
        dmbCS_CO = de.SCHOFe_OCH_CS_CO(xg_O2_,xg_CO_,xg_H2_,constants)
        dmbCS_H = de.SCHOFe_OCH_CS_H2(xg_O2_,xg_CO_,xg_H2_,constants)
        return dmbCO_O, dmbCO_CO, dmbCO_H, dmbCH_O, dmbCH_CO, dmbCH_H, dmbCS_O, dmbCS_CO, dmbCS_H
    
    if solve_species == "OCS":
        xg_O2_, xg_CO_, xg_S2_ = jac_newton3(guessx,guessy,guessz,constants,f_SCHOFe_OCS,df_SCHOFe_OCS,nr_step,nr_tol)
        results1 = xg_O2_, xg_CO_, xg_S2_
        results2 = mg_SCHOFe_OCS(xg_O2_,xg_CO_,xg_S2_)
        results3 = f_SCHOFe_OCS(xg_O2_,xg_CO_,xg_S2_)
        xg_H2_ = results2[2]
    elif solve_species == "OHS":
        xg_O2_, xg_H2_, xg_S2_ = jac_newton3(guessx,guessy,guessz,constants,f_SCHOFe_OHS,df_SCHOFe_OHS,nr_step,nr_tol)
        results1 = xg_O2_, xg_H2_, xg_S2_
        results2 = mg_SCHOFe_OHS(xg_O2_,xg_H2_,xg_S2_)
        results3 = f_SCHOFe_OHS(xg_O2_,xg_H2_,xg_S2_)
        xg_CO_ = results2[0]
    elif solve_species == "OCH":
        xg_O2_, xg_CO_, xg_H2_ = jac_newton3(guessx,guessy,guessz,constants,f_SCHOFe_OCH,df_SCHOFe_OCH,nr_step,nr_tol)
        results1 = xg_O2_, xg_CO_, xg_H2_
        results2 = mg_SCHOFe_OCH(xg_O2_,xg_CO_,xg_H2_)
        results3 = f_SCHOFe_OCH(xg_O2_,xg_CO_,xg_H2_)
        xg_S2_ = results2[3]    
    results4 = mb_SCHOFe(xg_O2_,xg_CO_,xg_H2_,xg_S2_)
    return results1, results2, results3, results4


# all gas species - oxidised and reduced melt species - H2S #

def eq_SCHOFe_2(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses,solve_species):
    P = PT["P"]
    wt_O = bulk_wf['O']
    wt_S = bulk_wf['S']
    wt_H = bulk_wf['H']
    wt_C = bulk_wf['C']
    wt_Fe = bulk_wf['Fe']
    guessx = guesses["guessx"]
    guessy = guesses["guessy"]
    guessz = guesses["guessz"]
    guessw = guesses["guessw"]
        
    # equilibrium constants
    K1_ = mdv.KHOg(PT,models)
    K2_ = mdv.KCOg(PT,models)
    K3_ = mdv.KCOHg(PT,models)
    K4_ = mdv.C_H2O(PT,melt_wf,models)
    K5_ = mdv.C_CO3(PT,melt_wf,models)
    K6_ = mdv.KOSg(PT,models)
    K7_ = mdv.KHOSg(PT,models)
    K8_ = mdv.C_S(PT,melt_wf,models)/1000000.0
    K9_ = mdv.C_SO4(PT,melt_wf,models)/1000000.0
    K10_ = mdv.KOCSg(PT,models)
    K11_ = mdv.C_CO(PT,melt_wf,models)/1000000.0
    K12_ = mdv.C_CH4(PT,melt_wf,models)/1000000.0
    K13_ = mdv.C_H2(PT,melt_wf,models)/1000000.0
    K14_ = mdv.C_H2S(PT,melt_wf,models)/1000000.0
    KD1_, KD2, y = mdv.FefO2_KC91_EqA_terms(PT,melt_wf,models)
   
    # fugacity coefficients
    y_CO_ = mdv.y_CO(PT,models)
    y_O2_ = mdv.y_O2(PT,models)
    y_CO2_ = mdv.y_CO2(PT,models)
    y_CH4_ = mdv.y_CH4(PT,models)
    y_H2O_ = mdv.y_H2O(PT,models)
    y_H2_ = mdv.y_H2(PT,models)
    y_S2_ = mdv.y_S2(PT,models)
    y_SO2_ = mdv.y_SO2(PT,models)
    y_H2S_ = mdv.y_H2S(PT,models)
    y_OCS_ = mdv.y_OCS(PT,models)
    
    # molecular masses
    M_H = mdv.species.loc['H','M']
    M_C = mdv.species.loc['C','M']
    M_O = mdv.species.loc['O','M']
    M_S = mdv.species.loc['S','M']
    M_Fe = mdv.species.loc['Fe','M']
    M_O2 = mdv.species.loc['O2','M']
    M_CO = mdv.species.loc['CO','M']
    M_OH = mdv.species.loc['OH','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_H2 = mdv.species.loc['H2','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_CH4 = mdv.species.loc['CH4','M']
    M_S2 = mdv.species.loc['S2','M']
    M_SO2 = mdv.species.loc['SO2','M']
    M_SO3 = mdv.species.loc['SO3','M']
    M_H2S = mdv.species.loc['H2S','M']
    M_FeO = mdv.species.loc['FeO','M']
    M_FeO15 = mdv.species.loc['FeO1.5','M']
    M_OCS = mdv.species.loc['OCS','M']
    M_m_ = mg.M_m_SO(melt_wf)
    M_m_ox = mg.M_m_ox(melt_wf,models)
    
    constants = [P, wt_O, wt_C, wt_H, wt_S, wt_Fe, K1_, K2_, K3_, K4_, K5_, K6_, K7_, K8_, K9_, K10_, K11_, K12_, K13_, K14_, KD1_, KD2, y, y_H2_, y_O2_, y_H2O_, y_CO_, y_CO2_, y_CH4_, y_SO2_, y_S2_, y_H2S_, y_OCS_, M_C, M_H, M_O, M_S, M_Fe, M_O2, M_CO, M_CO2, M_OH, M_H2, M_H2O, M_CH4, M_S2, M_SO2, M_SO3, M_H2S, M_OCS, M_FeO, M_FeO15, M_m_, M_m_ox]
       
    def mg_SCHOFe(xg_O2_,xg_A,xg_B):
        xg_O2__ = xg_O2_
        if solve_species == "OCS": # A = CO, B = S2
            xg_CO_ = xg_A # mole fractions in the gas
            xg_S2_ = xg_B                 
            xg_CO2_ = (K2_*y_CO_*xg_CO_*(y_O2_*xg_O2_*P)**0.5)/y_CO2_
            xg_SO2_ = (K6_*y_O2_*xg_O2_*(y_S2_*xg_S2_*P)**0.5)/y_SO2_
            xg_OCS_ = ((xg_CO_*y_CO_)**3.0*xg_SO2_*y_SO2_*P)/(y_OCS_*(xg_CO2_*y_CO2_)**2.0*K10_)
            a = (y_CO2_*xg_CO2_*y_H2O_**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
            b = 1.0 + (y_H2O_/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)) + ((K7_*(y_S2_*xg_S2_)**0.5*y_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5))
            c = xg_CO2_ + xg_CO_ + xg_O2_ + xg_S2_ + xg_SO2_ + xg_OCS_ - 1.0
            xg_H2O_ = (-b + (b**2.0-(4.0*a*c))**0.5)/(2.0*a)
            xg_H2_ = (y_H2O_*xg_H2O_)/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)
            xg_CH4_ = (y_CO2_*xg_CO2_*(y_H2O_*xg_H2O_)**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
            xg_H2S_ = (K7_*(y_S2_*xg_S2_)**0.5*y_H2O_*xg_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5)
        elif solve_species == "OHS": # A = H2, B = S2
            xg_H2_ = xg_A # mole fractions in the gas
            xg_S2_ = xg_B
            xg_H2O_ = (xg_H2_*K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)/y_H2O_
            xg_H2S_ = (K7_*(y_S2_*xg_S2_)**0.5*y_H2O_*xg_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5)
            xg_SO2_ = (K6_*y_O2_*xg_O2_*(y_S2_*xg_S2_*P)**0.5)/y_SO2_
            xg_CO2_ = (1.0 - xg_S2_ - xg_O2_ - xg_H2_ - xg_H2O_ - xg_H2S_ - xg_SO2_)/(1.0 + (y_CO2_/(K2_*y_CO_*(y_O2_*xg_O2_*P)**0.5)) + ((y_CO2_*y_SO2_*xg_SO2_)/(K2_**3.0*K10_*y_OCS_*(y_O2_*xg_O2_)**1.5*P**0.5)) + ((y_CO2_*(y_H2O_*xg_H2O_)**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)))
            xg_CO_ = (y_CO2_*xg_CO2_)/(K2_*y_CO_*(y_O2_*xg_O2_*P)**0.5)
            xg_CH4_ = (y_CO2_*xg_CO2_*(y_H2O_*xg_H2O_)**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
            xg_OCS_ = ((xg_CO_*y_CO_)**3.0*xg_SO2_*y_SO2_*P)/(y_OCS_*(xg_CO2_*y_CO2_)**2.0*K10_)        
        elif solve_species == "OCH": # A = CO, B = H2
            xg_CO_ = xg_A # mole fractions in the gas
            xg_H2_ = xg_B
            xg_H2O_ = (xg_H2_*K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)/y_H2O_
            xg_CO2_ = (K2_*y_CO_*xg_CO_*(y_O2_*xg_O2_*P)**0.5)/y_CO2_
            xg_CH4_ = (y_CO2_*xg_CO2_*(y_H2O_*xg_H2O_)**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
            a = 1.0
            b = ((K6_*(y_S2_*P)**0.5*y_O2_*xg_O2_)/y_SO2_) + (((K7_*y_H2O_*xg_H2O_)/y_H2S_)*(y_S2_/(y_O2_*xg_O2_))**0.5) + (((y_CO_*xg_CO_)**3.0*P**1.5*K6_*y_S2_**0.5*y_O2_*xg_O2_)/(K10_*y_OCS_*(y_CO2_*xg_CO2_)**2.0))
            c = xg_O2_ + xg_CO_ + xg_CO2_ + xg_H2_ + xg_H2O_ + xg_CH4_ - 1.0
            xg_S2_ = ((-b + (b**2.0-(4.0*a*c))**0.5)/(2.0*a))**2.0
            xg_SO2_ = (K6_*y_O2_*xg_O2_*(y_S2_*xg_S2_*P)**0.5)/y_SO2_
            xg_H2S_ = (K7_*(y_S2_*xg_S2_)**0.5*y_H2O_*xg_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5)
            xg_OCS_ = ((xg_CO_*y_CO_)**3.0*xg_SO2_*y_SO2_*P)/(y_OCS_*(xg_CO2_*y_CO2_)**2.0*K10_)        
        Xg_t = xg_CO2_*M_CO2 + xg_CO_*M_CO + xg_O2_*M_O2 + xg_H2O_*M_H2O + xg_H2_*M_H2 + xg_CH4_*M_CH4 + xg_SO2_*M_SO2 + xg_S2_*M_S2 + xg_H2S_*M_H2S + xg_OCS_*M_OCS
        # melt composition
        xm_H2O_ = (K4_*y_H2O_*xg_H2O_*P)**0.5
        xm_CO2_ = K5_*y_CO2_*xg_CO2_*P
        Xm_t = xm_CO2_*M_CO2 + xm_H2O_*M_H2O + (1.0-xm_CO2_-xm_H2O_)*M_m_
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t
        wm_CO2_ = (xm_CO2_*M_CO2)/Xm_t
        wm_CO_ = K11_*y_CO_*xg_CO_*P
        wm_CH4_ = K12_*y_CH4_*xg_CH4_*P
        wm_H2_ = K13_*y_H2_*xg_H2_*P
        wm_H2S_ = K14_*y_H2S_*xg_H2S_*P
        wm_S_ = K8_*((y_S2_*xg_S2_)/(y_O2_*xg_O2_))**0.5
        wm_SO3_ = ((K9_*(y_S2_*xg_S2_)**0.5*(y_O2_*xg_O2_)**1.5*P**2.0)/M_S)*M_SO3
        wm_ST_ = wm_S_ + ((M_S*wm_SO3_)/M_SO3) + ((M_S*wm_H2S_)/M_H2S)
        Fe32 = ((KD1_*(y_O2_*xg_O2_*P)**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y)))
        Fe3T = Fe32/(1.0+Fe32)
        S62 = (wm_SO3_/M_SO3)/(wm_S_/M_S)
        S6T = S62/(1.0+S62)
        Xm_t_ox = ""
        return xg_O2__, xg_H2_, xg_S2_, xg_H2O_, xg_CO_, xg_CO2_, xg_SO2_, xg_CH4_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Xm_t_ox, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_, wm_H2S_, wm_H2_, wm_CH4_, wm_CO_    
    
    def mb_SCHOFe(xg_O2_,xg_CO_,xg_H2_,xg_S2_): # A = CO, B = H2, C = S2
        if solve_species == "OCS":
            xg_A = xg_CO_
            xg_B = xg_S2_
        elif solve_species == "OHS":
            xg_A = xg_H2_
            xg_B = xg_S2_
        elif solve_species == "OCH":
            xg_A = xg_CO_
            xg_B = xg_H2_
        xg_O2__, xg_H2_, xg_S2_, xg_H2O_, xg_CO_, xg_CO2_, xg_SO2_, xg_CH4_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Xm_t_ox, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_, wm_H2S_, wm_H2_, wm_CH4_, wm_CO_ = mg_SCHOFe(xg_O2_, xg_A, xg_B)
        mba, mbb, mbc, wt_g_O, wt_g_C, wt_g_H, wt_g_S = f_SCHOFe(xg_O2_,xg_A,xg_B)
        wt_g = (wt_g_O+wt_g_H+wt_g_C+wt_g_S)/4.
        wt_H_ = 2.*M_H*((wt_g*(((xg_H2O_+xg_H2_+2.*xg_CH4_+xg_H2S_)/Xg_t) - (wm_H2O_/M_H2O) - (wm_H2_/M_H2) - (wm_H2S_/M_H2S) - (2.*wm_CH4_/M_CH4))) + (wm_H2O_/M_H2O) + (wm_H2_/M_H2) + (wm_H2S_/M_H2S) + (2.*wm_CH4_/M_CH4))
        wt_O_ = M_O*((wt_g*(((2.0*xg_CO2_ + xg_CO_ + 2.0*xg_O2_ + xg_H2O_ + 2.0*xg_SO2_ + xg_OCS_)/Xg_t) - (wm_H2O_/M_H2O) - ((2.0*wm_CO2_)/M_CO2) - (3.0*wm_SO3_/M_SO3) - (wm_CO_/M_CO))) + (wm_H2O_/M_H2O) + ((2.0*wm_CO2_)/M_CO2) + (3.0*wm_SO3_/M_SO3) + (wm_CO_/M_CO) + (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))
        wt_C_ = M_C*((wt_g*(((xg_CO2_+xg_CO_+xg_CH4_+xg_OCS_)/Xg_t) - (wm_CO2_/M_CO2) - (wm_CH4_/M_CH4) - (wm_CO_/M_CO))) + (wm_CO2_/M_CO2) + (wm_CH4_/M_CH4) + (wm_CO_/M_CO))
        wt_S_ = M_S*((wt_g*(((xg_SO2_+2.0*xg_S2_+xg_H2S_+xg_OCS_)/Xg_t) - (wm_S_/M_S) - (wm_SO3_/M_SO3) - (wm_H2S_/M_H2S))) + (wm_S_/M_S) + (wm_SO3_/M_SO3) + (wm_H2S_/M_H2S))
        return wt_g, wt_O_, wt_C_, wt_H_, wt_S_
    
    def f_SCHOFe(xg_O2_,xg_A,xg_B):
        xg_O2__, xg_H2_, xg_S2_, xg_H2O_, xg_CO_, xg_CO2_, xg_SO2_, xg_CH4_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Xm_t_ox, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_, wm_H2S_, wm_H2_, wm_CH4_, wm_CO_ = mg_SCHOFe(xg_O2_, xg_A, xg_B)
        wt_g_C = ((wt_C/M_C) - (wm_CO2_/M_CO2) - (wm_CH4_/M_CH4) - (wm_CO_/M_CO))/(((xg_CO2_+xg_CO_+xg_CH4_+xg_OCS_)/Xg_t) - (wm_CO2_/M_CO2) - (wm_CH4_/M_CH4) - (wm_CO_/M_CO))
        wt_g_O = ((wt_O/M_O) - (wm_H2O_/M_H2O) - ((2.0*wm_CO2_)/M_CO2) - (3.0*wm_SO3_/M_SO3) - (wm_CO_/M_CO) - (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))/(((2.0*xg_CO2_ + xg_CO_ + 2.0*xg_O2_ + xg_H2O_ + 2.0*xg_SO2_ + xg_OCS_ )/Xg_t) - (wm_H2O_/M_H2O) - ((2.0*wm_CO2_)/M_CO2) - (3.0*wm_SO3_/M_SO3) - (wm_CO_/M_CO))
        wt_g_H = ((wt_H/(2.0*M_H)) - (wm_H2O_/M_H2O) - (wm_H2_/M_H2) - (2.*wm_CH4_/M_CH4) - (wm_H2S_/M_H2S))/(((xg_H2O_+xg_H2_+2.0*xg_CH4_+xg_H2S_)/Xg_t) - (wm_H2O_/M_H2O) - (wm_H2_/M_H2) - (2.*wm_CH4_/M_CH4) - (wm_H2S_/M_H2S))
        wt_g_S = ((wt_S/M_S) - (wm_S_/M_S) - (wm_SO3_/M_SO3) - (wm_H2S_/M_H2S))/(((xg_SO2_+2.0*xg_S2_+xg_H2S_+xg_OCS_)/Xg_t) - (wm_S_/M_S) - (wm_SO3_/M_SO3) - (wm_H2S_/M_H2S))
        if solve_species == "OCS":
            mba = (wt_g_C - wt_g_O) # mbCO
            mbb = (wt_g_C - wt_g_H) # mbCH
            mbc = (wt_g_C - wt_g_S) # mbCS
        elif solve_species == "OHS":
            mba = (wt_g_H - wt_g_O) # mbHO
            mbb = (wt_g_H - wt_g_C) # mbHC
            mbc = (wt_g_H - wt_g_S) # mbHS
        elif solve_species == "OCH":
            mba = (wt_g_C - wt_g_O) # mbCO
            mbb = (wt_g_C - wt_g_H) # mbCH
            mbc = (wt_g_C - wt_g_S) # mbCS
        return mba, mbb, mbc, wt_g_O, wt_g_C, wt_g_H, wt_g_S
    
    def df_SCHOFe(xg_O2_,xg_A,xg_B,constants): 
        if solve_species == "OCS": # A = CO, B = S2 and a = C-O, b = C-H, c = C-S
            dmba_O = de.SCHOFe2_OCS_CO_O2(xg_O2_,xg_A,xg_B,constants) # dmbCO_O
            dmba_A = de.SCHOFe2_OCS_CO_CO(xg_O2_,xg_A,xg_B,constants) # dmbCO_CO
            dmba_B = de.SCHOFe2_OCS_CO_S2(xg_O2_,xg_A,xg_B,constants) # dmbCO_S 
            dmbb_O = de.SCHOFe2_OCS_CH_O2(xg_O2_,xg_A,xg_B,constants) # dmbCH_O 
            dmbb_A = de.SCHOFe2_OCS_CH_CO(xg_O2_,xg_A,xg_B,constants) # dmbCH_CO
            dmbb_B = de.SCHOFe2_OCS_CH_S2(xg_O2_,xg_A,xg_B,constants) # dmbCH_S
            dmbc_O = de.SCHOFe2_OCS_CS_O2(xg_O2_,xg_A,xg_B,constants) # dmbCS_O 
            dmbc_A = de.SCHOFe2_OCS_CS_CO(xg_O2_,xg_A,xg_B,constants) # dmbCS_CO
            dmbc_B = de.SCHOFe2_OCS_CS_S2(xg_O2_,xg_A,xg_B,constants) # dmbCS_S
        elif solve_species == "OHS": # A = H2, B = S2 and a = H-O, b = H-C, c = H-S TO DO
            dmba_O = de.SCHOFe2_OHS_HO_O2(xg_O2_,xg_A,xg_B,constants) # dmbHO_O
            dmba_A = de.SCHOFe2_OHS_HO_H2(xg_O2_,xg_A,xg_B,constants) # dmbHO_H
            dmba_B = de.SCHOFe2_OHS_HO_S2(xg_O2_,xg_A,xg_B,constants) # dmbHO_S
            dmbb_O = de.SCHOFe2_OHS_HC_O2(xg_O2_,xg_A,xg_B,constants) # dmbHC_O
            dmbb_A = de.SCHOFe2_OHS_HC_H2(xg_O2_,xg_A,xg_B,constants) # dmbHC_H
            dmbb_B = de.SCHOFe2_OHS_HC_S2(xg_O2_,xg_A,xg_B,constants) # dmbHC_S
            dmbc_O = de.SCHOFe2_OHS_HS_O2(xg_O2_,xg_A,xg_B,constants) # dmbHS_O
            dmbc_A = de.SCHOFe2_OHS_HS_H2(xg_O2_,xg_A,xg_B,constants) # dmbHS_H
            dmbc_B = de.SCHOFe2_OHS_HS_S2(xg_O2_,xg_A,xg_B,constants) # dmbHS_S
        elif solve_species == "OCH": # A = CO, B = H2 and a = C-O, b = C-H, c = C-S TO DO
            dmba_O = de.SCHOFe2_OCH_CO_O2(xg_O2_,xg_A,xg_B,constants) # dmbCO_O
            dmba_A = de.SCHOFe2_OCH_CO_CO(xg_O2_,xg_A,xg_B,constants) # dmbCO_CO
            dmba_B = de.SCHOFe2_OCH_CO_H2(xg_O2_,xg_A,xg_B,constants) # dmbCO_H
            dmbb_O = de.SCHOFe2_OCH_CH_O2(xg_O2_,xg_A,xg_B,constants) # dmbCH_O
            dmbb_A = de.SCHOFe2_OCH_CH_CO(xg_O2_,xg_A,xg_B,constants) # dmbCH_CO
            dmbb_B = de.SCHOFe2_OCH_CH_H2(xg_O2_,xg_A,xg_B,constants) # dmbCH_H
            dmbc_O = de.SCHOFe2_OCH_CS_O2(xg_O2_,xg_A,xg_B,constants) # dmbCS_O
            dmbc_A = de.SCHOFe2_OCH_CS_CO(xg_O2_,xg_A,xg_B,constants) # dmbCS_CO
            dmbc_B = de.SCHOFe2_OCH_CS_H2(xg_O2_,xg_A,xg_B,constants) # dmbCS_H
        return dmba_O, dmba_A, dmba_B, dmbb_O, dmbb_A, dmbb_B, dmbc_O, dmbc_A, dmbc_B

    if solve_species == "OCS":
        xg_O2_, xg_CO_, xg_S2_ = jac_newton3(guessx,guessy,guessz,constants,f_SCHOFe,df_SCHOFe,nr_step,nr_tol)
    elif solve_species == "OHS":
        xg_O2_, xg_H2_, xg_S2_ = jac_newton3(guessx,guessy,guessz,constants,f_SCHOFe,df_SCHOFe,nr_step,nr_tol)
    elif solve_species == "OCH":
        xg_O2_, xg_CO_, xg_H2_ = jac_newton3(guessx,guessy,guessz,constants,f_SCHOFe,df_SCHOFe,nr_step,nr_tol)

    if xg_O2_ == 1.: # switch solve species once
        if solve_species == "OCS":
            print(PT["P"],": Switching solve species from OCS to OCH (first time)")
            solve_species = "OCH"
            models.loc["solve_species","option"] = "OCH"
            guessw_hold = guessw
            guessw = guessz # xgS2 is guessw
            guessz = guessw_hold # xgH2 is guessz
            xg_O2_, xg_CO_, xg_H2_ = jac_newton3(guessx,guessy,guessz,constants,f_SCHOFe,df_SCHOFe,nr_step,nr_tol)
        elif solve_species == "OHS":
            print(PT["P"],": Switching solve species from OHS to OCS (first time)")
            solve_species = "OCS"
            models.loc["solve_species","option"] = "OCS"
            guessw_hold = guessw
            guessw = guessy # xgH2 is guessw
            guessy = guessw_hold # xgCO is guessy
            xg_O2_, xg_CO_, xg_S2_ = jac_newton3(guessx,guessy,guessz,constants,f_SCHOFe,df_SCHOFe,nr_step,nr_tol)
        elif solve_species == "OCH":
            print(PT["P"],": Switching solve species from OCH to OHS (first time)")
            solve_species = "OHS"
            models.loc["solve_species","option"] = "OHS"
            guessw_hold = guessw
            guessz_hold = guessz
            guessw = guessy # xgCO is guessw
            guessy = guessz_hold # xgH2 is guessy
            guessz = guessw_hold # xgS2 is guess z
            xg_O2_, xg_H2_, xg_S2_ = jac_newton3(guessx,guessy,guessz,constants,f_SCHOFe,df_SCHOFe,nr_step,nr_tol)
    
    if xg_O2_ == 1.: # switch solve species second time
        if solve_species == "OCS":
            print(PT["P"],": Switching solve species from OCS to OCH (second time)")
            solve_species = "OCH"
            models.loc["solve_species","option"] = "OCH"
            guessw_hold = guessw
            guessw = guessz # xgS2 is guessw
            guessz = guessw_hold # xgH2 is guessz
            xg_O2_, xg_CO_, xg_H2_ = jac_newton3(guessx,guessy,guessz,constants,f_SCHOFe,df_SCHOFe,nr_step,nr_tol)
        elif solve_species == "OHS":
            print(PT["P"],": Switching solve species from OHS to OCS (second time)")
            solve_species = "OCS"
            models.loc["solve_species","option"] = "OCS"
            guessw_hold = guessw
            guessw = guessy # xgH2 is guessw
            guessy = guessw_hold # xgCO is guessy
            xg_O2_, xg_CO_, xg_S2_ = jac_newton3(guessx,guessy,guessz,constants,f_SCHOFe,df_SCHOFe,nr_step,nr_tol)
        elif solve_species == "OCH":
            print(PT["P"],": Switching solve species from OCH to OHS (second time)")
            solve_species = "OHS"
            models.loc["solve_species","option"] = "OHS"
            guessw_hold = guessw
            guessz_hold = guessz
            guessw = guessy # xgCO is guessw
            guessy = guessz_hold # xgH2 is guessy
            guessz = guessw_hold # xgS2 is guess z
            xg_O2_, xg_H2_, xg_S2_ = jac_newton3(guessx,guessy,guessz,constants,f_SCHOFe,df_SCHOFe,nr_step,nr_tol)

    if solve_species == "OCS":
        results1 = xg_O2_, xg_CO_, xg_S2_
        results2 = mg_SCHOFe(xg_O2_,xg_CO_,xg_S2_)
        results3 = f_SCHOFe(xg_O2_,xg_CO_,xg_S2_)
        xg_H2_ = results2[2]
    elif solve_species == "OHS":
        results1 = xg_O2_, xg_H2_, xg_S2_
        results2 = mg_SCHOFe(xg_O2_,xg_H2_,xg_S2_)
        results3 = f_SCHOFe(xg_O2_,xg_H2_,xg_S2_)
        xg_CO_ = results2[0]
    elif solve_species == "OCH":
        results1 = xg_O2_, xg_CO_, xg_H2_
        results2 = mg_SCHOFe(xg_O2_,xg_CO_,xg_H2_)
        results3 = f_SCHOFe(xg_O2_,xg_CO_,xg_H2_)
        xg_S2_ = results2[3]    
    results4 = mb_SCHOFe(xg_O2_,xg_CO_,xg_H2_,xg_S2_)
    return results1, results2, results3, results4, models, solve_species

# all gas species - oxidised and reduced melt species - H2S - H2O is linear with xmH2O #
def eq_SCHOFe_3(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses,solve_species):
    P = PT["P"]
    wt_O = bulk_wf['O']
    wt_S = bulk_wf['S']
    wt_H = bulk_wf['H']
    wt_C = bulk_wf['C']
    wt_Fe = bulk_wf['Fe']
    guessx = guesses["guessx"]
    guessy = guesses["guessy"]
    guessz = guesses["guessz"]
        
    # equilibrium constants
    K1_ = mdv.KHOg(PT,models)
    K2_ = mdv.KCOg(PT,models)
    K3_ = mdv.KCOHg(PT,models)
    K4_ = mdv.C_H2O(PT,melt_wf,models)
    K5_ = mdv.C_CO3(PT,melt_wf,models)
    K6_ = mdv.KOSg(PT,models)
    K7_ = mdv.KHOSg(PT,models)
    K8_ = mdv.C_S(PT,melt_wf,models)/1000000.0
    K9_ = mdv.C_SO4(PT,melt_wf,models)/1000000.0
    K10_ = mdv.KOCSg(PT,models)
    K11_ = mdv.C_CO(PT,melt_wf,models)/1000000.0
    K12_ = mdv.C_CH4(PT,melt_wf,models)/1000000.0
    K13_ = mdv.C_H2(PT,melt_wf,models)/1000000.0
    K14_ = mdv.C_H2S(PT,melt_wf,models)/1000000.0
    KD1_, KD2, y = mdv.FefO2_KC91_EqA_terms(PT,melt_wf,models)
   
    # fugacity coefficients
    y_CO_ = mdv.y_CO(PT,models)
    y_O2_ = mdv.y_O2(PT,models)
    y_CO2_ = mdv.y_CO2(PT,models)
    y_CH4_ = mdv.y_CH4(PT,models)
    y_H2O_ = mdv.y_H2O(PT,models)
    y_H2_ = mdv.y_H2(PT,models)
    y_S2_ = mdv.y_S2(PT,models)
    y_SO2_ = mdv.y_SO2(PT,models)
    y_H2S_ = mdv.y_H2S(PT,models)
    y_OCS_ = mdv.y_OCS(PT,models)
    
    # molecular masses
    M_H = mdv.species.loc['H','M']
    M_C = mdv.species.loc['C','M']
    M_O = mdv.species.loc['O','M']
    M_S = mdv.species.loc['S','M']
    M_Fe = mdv.species.loc['Fe','M']
    M_O2 = mdv.species.loc['O2','M']
    M_CO = mdv.species.loc['CO','M']
    M_OH = mdv.species.loc['OH','M']
    M_H2O = mdv.species.loc['H2O','M']
    M_H2 = mdv.species.loc['H2','M']
    M_CO2 = mdv.species.loc['CO2','M']
    M_CH4 = mdv.species.loc['CH4','M']
    M_S2 = mdv.species.loc['S2','M']
    M_SO2 = mdv.species.loc['SO2','M']
    M_SO3 = mdv.species.loc['SO3','M']
    M_H2S = mdv.species.loc['H2S','M']
    M_FeO = mdv.species.loc['FeO','M']
    M_FeO15 = mdv.species.loc['FeO1.5','M']
    M_OCS = mdv.species.loc['OCS','M']
    M_m_ = mg.M_m_SO(melt_wf)
    M_m_ox = mg.M_m_ox(melt_wf,models)
    
    constants = [P, wt_O, wt_C, wt_H, wt_S, wt_Fe, K1_, K2_, K3_, K4_, K5_, K6_, K7_, K8_, K9_, K10_, K11_, K12_, K13_, K14_, KD1_, KD2, y, y_H2_, y_O2_, y_H2O_, y_CO_, y_CO2_, y_CH4_, y_SO2_, y_S2_, y_H2S_, y_OCS_, M_C, M_H, M_O, M_S, M_Fe, M_O2, M_CO, M_CO2, M_OH, M_H2, M_H2O, M_CH4, M_S2, M_SO2, M_SO3, M_H2S, M_OCS, M_FeO, M_FeO15, M_m_, M_m_ox]
       
    def mg_SCHOFe(xg_O2_,xg_A,xg_B):
        xg_O2__ = xg_O2_
        if solve_species == "OCS": # A = CO, B = S2
            xg_CO_ = xg_A # mole fractions in the gas
            xg_S2_ = xg_B                 
            xg_CO2_ = (K2_*y_CO_*xg_CO_*(y_O2_*xg_O2_*P)**0.5)/y_CO2_
            xg_SO2_ = (K6_*y_O2_*xg_O2_*(y_S2_*xg_S2_*P)**0.5)/y_SO2_
            xg_OCS_ = ((xg_CO_*y_CO_)**3.0*xg_SO2_*y_SO2_*P)/(y_OCS_*(xg_CO2_*y_CO2_)**2.0*K10_)
            a = (y_CO2_*xg_CO2_*y_H2O_**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
            b = 1.0 + (y_H2O_/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)) + ((K7_*(y_S2_*xg_S2_)**0.5*y_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5))
            c = xg_CO2_ + xg_CO_ + xg_O2_ + xg_S2_ + xg_SO2_ + xg_OCS_ - 1.0
            xg_H2O_ = (-b + (b**2.0-(4.0*a*c))**0.5)/(2.0*a)
            xg_H2_ = (y_H2O_*xg_H2O_)/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)
            xg_CH4_ = (y_CO2_*xg_CO2_*(y_H2O_*xg_H2O_)**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
            xg_H2S_ = (K7_*(y_S2_*xg_S2_)**0.5*y_H2O_*xg_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5)
        #elif solve_species == "OHS": # A = H2, B = S2
            #xg_H2_ = xg_A # mole fractions in the gas
            #xg_S2_ = xg_B
            #xg_H2O_ = (xg_H2_*K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)/y_H2O_
            #xg_H2S_ = (K7_*(y_S2_*xg_S2_)**0.5*y_H2O_*xg_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5)
            #xg_SO2_ = (K6_*y_O2_*xg_O2_*(y_S2_*xg_S2_*P)**0.5)/y_SO2_
            #xg_CO2_ = (1.0 - xg_S2_ - xg_O2_ - xg_H2_ - xg_H2O_ - xg_H2S_ - xg_SO2_)/(1.0 + (y_CO2_/(K2_*y_CO_*(y_O2_*xg_O2_*P)**0.5)) + ((y_CO2_*y_SO2_*xg_SO2_)/(K2_**3.0*K10_*y_OCS_*(y_O2_*xg_O2_)**1.5*P**0.5)) + ((y_CO2_*(y_H2O_*xg_H2O_)**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)))
            #xg_CO_ = (y_CO2_*xg_CO2_)/(K2_*y_CO_*(y_O2_*xg_O2_*P)**0.5)
            #xg_CH4_ = (y_CO2_*xg_CO2_*(y_H2O_*xg_H2O_)**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
            #xg_OCS_ = ((xg_CO_*y_CO_)**3.0*xg_SO2_*y_SO2_*P)/(y_OCS_*(xg_CO2_*y_CO2_)**2.0*K10_)        
        #elif solve_species == "OCH": # A = CO, B = H2
            #xg_CO_ = xg_A # mole fractions in the gas
            #xg_H2_ = xg_B
            #xg_H2O_ = (xg_H2_*K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)/y_H2O_
            #xg_CO2_ = (K2_*y_CO_*xg_CO_*(y_O2_*xg_O2_*P)**0.5)/y_CO2_
            #xg_CH4_ = (y_CO2_*xg_CO2_*(y_H2O_*xg_H2O_)**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
            #a = 1.0
            #b = ((K6_*(y_S2_*P)**0.5*y_O2_*xg_O2_)/y_SO2_) + (((K7_*y_H2O_*xg_H2O_)/y_H2S_)*(y_S2_/(y_O2_*xg_O2_))**0.5) + (((y_CO_*xg_CO_)**3.0*P**1.5*K6_*y_S2_**0.5*y_O2_*xg_O2_)/(K10_*y_OCS_*(y_CO2_*xg_CO2_)**2.0))
            #c = xg_O2_ + xg_CO_ + xg_CO2_ + xg_H2_ + xg_H2O_ + xg_CH4_ - 1.0
            #xg_S2_ = ((-b + (b**2.0-(4.0*a*c))**0.5)/(2.0*a))**2.0
            #xg_SO2_ = (K6_*y_O2_*xg_O2_*(y_S2_*xg_S2_*P)**0.5)/y_SO2_
            #xg_H2S_ = (K7_*(y_S2_*xg_S2_)**0.5*y_H2O_*xg_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5)
            #xg_OCS_ = ((xg_CO_*y_CO_)**3.0*xg_SO2_*y_SO2_*P)/(y_OCS_*(xg_CO2_*y_CO2_)**2.0*K10_)        
        Xg_t = xg_CO2_*M_CO2 + xg_CO_*M_CO + xg_O2_*M_O2 + xg_H2O_*M_H2O + xg_H2_*M_H2 + xg_CH4_*M_CH4 + xg_SO2_*M_SO2 + xg_S2_*M_S2 + xg_H2S_*M_H2S + xg_OCS_*M_OCS
        # melt composition
        xm_H2O_ = K4_*y_H2O_*xg_H2O_*P
        xm_CO2_ = K5_*y_CO2_*xg_CO2_*P
        Xm_t = xm_CO2_*M_CO2 + xm_H2O_*M_H2O + (1.0-xm_CO2_-xm_H2O_)*M_m_
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t
        wm_CO2_ = (xm_CO2_*M_CO2)/Xm_t
        wm_CO_ = K11_*y_CO_*xg_CO_*P
        wm_CH4_ = K12_*y_CH4_*xg_CH4_*P
        wm_H2_ = K13_*y_H2_*xg_H2_*P
        wm_H2S_ = K14_*y_H2S_*xg_H2S_*P
        wm_S_ = K8_*((y_S2_*xg_S2_)/(y_O2_*xg_O2_))**0.5
        wm_SO3_ = ((K9_*(y_S2_*xg_S2_)**0.5*(y_O2_*xg_O2_)**1.5*P**2.0)/M_S)*M_SO3
        wm_ST_ = wm_S_ + ((M_S*wm_SO3_)/M_SO3) + ((M_S*wm_H2S_)/M_H2S)
        Fe32 = ((KD1_*(y_O2_*xg_O2_*P)**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y)))
        Fe3T = Fe32/(1.0+Fe32)
        S62 = (wm_SO3_/M_SO3)/(wm_S_/M_S)
        S6T = S62/(1.0+S62)
        Xm_t_ox = ""
        return xg_O2__, xg_H2_, xg_S2_, xg_H2O_, xg_CO_, xg_CO2_, xg_SO2_, xg_CH4_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Xm_t_ox, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_, wm_H2S_, wm_H2_, wm_CH4_, wm_CO_    
    
    def mb_SCHOFe(xg_O2_,xg_CO_,xg_H2_,xg_S2_): # A = CO, B = H2, C = S2
        if solve_species == "OCS":
            xg_A = xg_CO_
            xg_B = xg_S2_
        #elif solve_species == "OHS":
            #xg_A = xg_H2_
            #xg_B = xg_S2_
        #elif solve_species == "OCH":
            #xg_A = xg_CO_
            #ßxg_B = xg_H2_
        xg_O2__, xg_H2_, xg_S2_, xg_H2O_, xg_CO_, xg_CO2_, xg_SO2_, xg_CH4_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Xm_t_ox, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_, wm_H2S_, wm_H2_, wm_CH4_, wm_CO_ = mg_SCHOFe(xg_O2_, xg_A, xg_B)
        mba, mbb, mbc, wt_g_O, wt_g_C, wt_g_H, wt_g_S = f_SCHOFe(xg_O2_,xg_A,xg_B)
        wt_g = (wt_g_O+wt_g_H+wt_g_C+wt_g_S)/4.
        wt_H_ = 2.*M_H*((wt_g*(((xg_H2O_+xg_H2_+2.*xg_CH4_+xg_H2S_)/Xg_t) - (wm_H2O_/M_H2O) - (wm_H2_/M_H2) - (wm_H2S_/M_H2S) - (2.*wm_CH4_/M_CH4))) + (wm_H2O_/M_H2O) + (wm_H2_/M_H2) + (wm_H2S_/M_H2S) + (2.*wm_CH4_/M_CH4))
        wt_O_ = M_O*((wt_g*(((2.0*xg_CO2_ + xg_CO_ + 2.0*xg_O2_ + xg_H2O_ + 2.0*xg_SO2_ + xg_OCS_)/Xg_t) - (wm_H2O_/M_H2O) - ((2.0*wm_CO2_)/M_CO2) - (3.0*wm_SO3_/M_SO3) - (wm_CO_/M_CO))) + (wm_H2O_/M_H2O) + ((2.0*wm_CO2_)/M_CO2) + (3.0*wm_SO3_/M_SO3) + (wm_CO_/M_CO) + (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))
        wt_C_ = M_C*((wt_g*(((xg_CO2_+xg_CO_+xg_CH4_+xg_OCS_)/Xg_t) - (wm_CO2_/M_CO2) - (wm_CH4_/M_CH4) - (wm_CO_/M_CO))) + (wm_CO2_/M_CO2) + (wm_CH4_/M_CH4) + (wm_CO_/M_CO))
        wt_S_ = M_S*((wt_g*(((xg_SO2_+2.0*xg_S2_+xg_H2S_+xg_OCS_)/Xg_t) - (wm_S_/M_S) - (wm_SO3_/M_SO3) - (wm_H2S_/M_H2S))) + (wm_S_/M_S) + (wm_SO3_/M_SO3) + (wm_H2S_/M_H2S))
        return wt_g, wt_O_, wt_C_, wt_H_, wt_S_
    
    def f_SCHOFe(xg_O2_,xg_A,xg_B):
        xg_O2__, xg_H2_, xg_S2_, xg_H2O_, xg_CO_, xg_CO2_, xg_SO2_, xg_CH4_, xg_H2S_, xg_OCS_, Xg_t, xm_H2O_, xm_CO2_, wm_S_, wm_SO3_, Xm_t, Xm_t_ox, Fe32, Fe3T, S62, S6T, wm_H2O_, wm_CO2_, wm_ST_, wm_H2S_, wm_H2_, wm_CH4_, wm_CO_ = mg_SCHOFe(xg_O2_, xg_A, xg_B)
        wt_g_C = ((wt_C/M_C) - (wm_CO2_/M_CO2) - (wm_CH4_/M_CH4) - (wm_CO_/M_CO))/(((xg_CO2_+xg_CO_+xg_CH4_+xg_OCS_)/Xg_t) - (wm_CO2_/M_CO2) - (wm_CH4_/M_CH4) - (wm_CO_/M_CO))
        wt_g_O = ((wt_O/M_O) - (wm_H2O_/M_H2O) - ((2.0*wm_CO2_)/M_CO2) - (3.0*wm_SO3_/M_SO3) - (wm_CO_/M_CO) - (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))/(((2.0*xg_CO2_ + xg_CO_ + 2.0*xg_O2_ + xg_H2O_ + 2.0*xg_SO2_ + xg_OCS_ )/Xg_t) - (wm_H2O_/M_H2O) - ((2.0*wm_CO2_)/M_CO2) - (3.0*wm_SO3_/M_SO3) - (wm_CO_/M_CO))
        wt_g_H = ((wt_H/(2.0*M_H)) - (wm_H2O_/M_H2O) - (wm_H2_/M_H2) - (2.*wm_CH4_/M_CH4) - (wm_H2S_/M_H2S))/(((xg_H2O_+xg_H2_+2.0*xg_CH4_+xg_H2S_)/Xg_t) - (wm_H2O_/M_H2O) - (wm_H2_/M_H2) - (2.*wm_CH4_/M_CH4) - (wm_H2S_/M_H2S))
        wt_g_S = ((wt_S/M_S) - (wm_S_/M_S) - (wm_SO3_/M_SO3) - (wm_H2S_/M_H2S))/(((xg_SO2_+2.0*xg_S2_+xg_H2S_+xg_OCS_)/Xg_t) - (wm_S_/M_S) - (wm_SO3_/M_SO3) - (wm_H2S_/M_H2S))
        if solve_species == "OCS":
            mba = (wt_g_C - wt_g_O) # mbCO
            mbb = (wt_g_C - wt_g_H) # mbCH
            mbc = (wt_g_C - wt_g_S) # mbCS
        #elif solve_species == "OHS":
            #mba = (wt_g_H - wt_g_O) # mbHO
            #mbb = (wt_g_H - wt_g_C) # mbHC
            #mbc = (wt_g_H - wt_g_S) # mbHS
        #elif solve_species == "OCH":
            #mba = (wt_g_C - wt_g_O) # mbCO
            #mbb = (wt_g_C - wt_g_H) # mbCH
            #mbc = (wt_g_C - wt_g_S) # mbCS
        return mba, mbb, mbc, wt_g_O, wt_g_C, wt_g_H, wt_g_S
    
    def df_SCHOFe(xg_O2_,xg_A,xg_B,constants): 
        if solve_species == "OCS": # A = CO, B = S2 and a = C-O, b = C-H, c = C-S
            dmba_O = de.SCHOFe3_OCS_CO_O2(xg_O2_,xg_A,xg_B,constants) # dmbCO_O
            dmba_A = de.SCHOFe3_OCS_CO_CO(xg_O2_,xg_A,xg_B,constants) # dmbCO_CO
            dmba_B = de.SCHOFe3_OCS_CO_S2(xg_O2_,xg_A,xg_B,constants) # dmbCO_S 
            dmbb_O = de.SCHOFe3_OCS_CH_O2(xg_O2_,xg_A,xg_B,constants) # dmbCH_O 
            dmbb_A = de.SCHOFe3_OCS_CH_CO(xg_O2_,xg_A,xg_B,constants) # dmbCH_CO
            dmbb_B = de.SCHOFe3_OCS_CH_S2(xg_O2_,xg_A,xg_B,constants) # dmbCH_S
            dmbc_O = de.SCHOFe3_OCS_CS_O2(xg_O2_,xg_A,xg_B,constants) # dmbCS_O 
            dmbc_A = de.SCHOFe3_OCS_CS_CO(xg_O2_,xg_A,xg_B,constants) # dmbCS_CO
            dmbc_B = de.SCHOFe3_OCS_CS_S2(xg_O2_,xg_A,xg_B,constants) # dmbCS_S
        #elif solve_species == "OHS": # A = H2, B = S2 and a = H-O, b = H-C, c = H-S TO DO
            #dmba_O = de.SCHOFe3_OHS_HO_O2(xg_O2_,xg_A,xg_B,constants) # dmbHO_O
            #dmba_A = de.SCHOFe3_OHS_HO_H2(xg_O2_,xg_A,xg_B,constants) # dmbHO_H
            #dmba_B = de.SCHOFe3_OHS_HO_S2(xg_O2_,xg_A,xg_B,constants) # dmbHO_S
            #dmbb_O = de.SCHOFe3_OHS_HC_O2(xg_O2_,xg_A,xg_B,constants) # dmbHC_O
            #dmbb_A = de.SCHOFe3_OHS_HC_H2(xg_O2_,xg_A,xg_B,constants) # dmbHC_H
            #dmbb_B = de.SCHOFe3_OHS_HC_S2(xg_O2_,xg_A,xg_B,constants) # dmbHC_S
            #dmbc_O = de.SCHOFe3_OHS_HS_O2(xg_O2_,xg_A,xg_B,constants) # dmbHS_O
            #dmbc_A = de.SCHOFe3_OHS_HS_H2(xg_O2_,xg_A,xg_B,constants) # dmbHS_H
            #dmbc_B = de.SCHOFe3_OHS_HS_S2(xg_O2_,xg_A,xg_B,constants) # dmbHS_S
        #elif solve_species == "OCH": # A = CO, B = H2 and a = C-O, b = C-H, c = C-S TO DO
            #dmba_O = de.SCHOFe3_OCH_CO_O2(xg_O2_,xg_A,xg_B,constants) # dmbCO_O
            #dmba_A = de.SCHOFe3_OCH_CO_CO(xg_O2_,xg_A,xg_B,constants) # dmbCO_CO
            #dmba_B = de.SCHOFe3_OCH_CO_H2(xg_O2_,xg_A,xg_B,constants) # dmbCO_H
            #dmbb_O = de.SCHOFe3_OCH_CH_O2(xg_O2_,xg_A,xg_B,constants) # dmbCH_O
            #dmbb_A = de.SCHOFe3_OCH_CH_CO(xg_O2_,xg_A,xg_B,constants) # dmbCH_CO
            #dmbb_B = de.SCHOFe3_OCH_CH_H2(xg_O2_,xg_A,xg_B,constants) # dmbCH_H
            #dmbc_O = de.SCHOFe3_OCH_CS_O2(xg_O2_,xg_A,xg_B,constants) # dmbCS_O
            #dmbc_A = de.SCHOFe3_OCH_CS_CO(xg_O2_,xg_A,xg_B,constants) # dmbCS_CO
            #dmbc_B = de.SCHOFe3_OCH_CS_H2(xg_O2_,xg_A,xg_B,constants) # dmbCS_H
        return dmba_O, dmba_A, dmba_B, dmbb_O, dmbb_A, dmbb_B, dmbc_O, dmbc_A, dmbc_B

    if solve_species == "OCS":
        xg_O2_, xg_CO_, xg_S2_ = jac_newton3(guessx,guessy,guessz,constants,f_SCHOFe,df_SCHOFe,nr_step,nr_tol)
        results1 = xg_O2_, xg_CO_, xg_S2_
        results2 = mg_SCHOFe(xg_O2_,xg_CO_,xg_S2_)
        results3 = f_SCHOFe(xg_O2_,xg_CO_,xg_S2_)
        xg_H2_ = results2[2]
    #elif solve_species == "OHS":
        #xg_O2_, xg_H2_, xg_S2_ = jac_newton3(guessx,guessy,guessz,constants,f_SCHOFe,df_SCHOFe,nr_step,nr_tol)
        #results1 = xg_O2_, xg_H2_, xg_S2_
        #results2 = mg_SCHOFe(xg_O2_,xg_H2_,xg_S2_)
        #results3 = f_SCHOFe(xg_O2_,xg_H2_,xg_S2_)
        #xg_CO_ = results2[0]
    #elif solve_species == "OCH":
        #xg_O2_, xg_CO_, xg_H2_ = jac_newton3(guessx,guessy,guessz,constants,f_SCHOFe,df_SCHOFe,nr_step,nr_tol)
        #results1 = xg_O2_, xg_CO_, xg_H2_
        #results2 = mg_SCHOFe(xg_O2_,xg_CO_,xg_H2_)
        #results3 = f_SCHOFe(xg_O2_,xg_CO_,xg_H2_)
        #xg_S2_ = results2[3]    
    results4 = mb_SCHOFe(xg_O2_,xg_CO_,xg_H2_,xg_S2_)
    return results1, results2, results3, results4

###############
### EDITING ###
###############

def eq_SCHOXFe(PT,bulk_wf,melt_wf,models,nr_step,nr_tol,guesses,solve_species):
    species_X = models.loc["species X","option"]
    
    P = float(PT["P"])
    wt_O = float(bulk_wf['O'])
    wt_S = float(bulk_wf['S'])
    wt_H = float(bulk_wf['H'])
    wt_C = float(bulk_wf['C'])
    wt_X = float(bulk_wf['X'])
    wt_Fe = float(bulk_wf['Fe'])
    guessx = guesses["guessx"]
    guessy = guesses["guessy"]
    guessz = guesses["guessz"]
    guessz = guesses["guessz"]
        
    # equilibrium constants
    K1_ = float(mdv.KHOg(PT,models))
    K2_ = float(mdv.KCOg(PT,models))
    K3_ = float(mdv.KCOHg(PT,models))
    K4_ = float(mdv.C_H2O(rT,melt_wf,models))
    K5_ = float(mdv.C_CO3(PT,melt_wf,models))
    K6_ = float(mdv.KOSg(PT,models))
    K7_ = float(mdv.KHOSg(PT,models))
    K8_ = float(mdv.C_S(PT,melt_wf,models)/1000000.0)
    K9_ = float(mdv.C_SO4(PT,melt_wf,models)/1000000.0)
    K10_ = float(mdv.KOCSg(PT,models))
    K11_ = float(mdv.C_CO(PT,melt_wf,models)/1000000.0)
    K12_ = float(mdv.C_CH4(PT,melt_wf,models)/1000000.0)
    K13_ = float(mdv.C_H2(PT,melt_wf,models)/1000000.0)
    K14_ = float(mdv.C_H2S(PT,melt_wf,models)/1000000.0)
    K15_ = float(mdv.C_X(PT,melt_wf,models)/1000000.0)
    KD1_, KD2, y = mdv.FefO2_KC91_EqA_terms(PT,melt_wf,models)
    KD1_=float(KD1_)
    KD2=float(KD2)
    y=float(y)
   
    # fugacity coefficients
    y_CO_ = float(mdv.y_CO(PT,models))
    y_O2_ = float(mdv.y_O2(PT,models))
    y_CO2_ = float(mdv.y_CO2(PT,models))
    y_CH4_ = float(mdv.y_CH4(PT,models))
    y_H2O_ = float(mdv.y_H2O(PT,models))
    y_H2_ = float(mdv.y_H2(PT,models))
    y_S2_ = float(mdv.y_S2(PT,models))
    y_SO2_ = float(mdv.y_SO2(PT,models))
    y_H2S_ = float(mdv.y_H2S(PT,models))
    y_OCS_ = float(mdv.y_OCS(PT,models))
    y_X_ = float(mdv.y_X(PT,models))
    
    # molecular masses
    M_H = float(mdv.species.loc['H','M'])
    M_C = float(mdv.species.loc['C','M'])
    M_O = float(mdv.species.loc['O','M'])
    M_S = float(mdv.species.loc['S','M'])
    M_Fe = float(mdv.species.loc['Fe','M'])
    M_O2 = float(mdv.species.loc['O2','M'])
    M_CO = float(mdv.species.loc['CO','M'])
    M_OH = float(mdv.species.loc['OH','M'])
    M_H2O = float(mdv.species.loc['H2O','M'])
    M_H2 = float(mdv.species.loc['H2','M'])
    M_CO2 = float(mdv.species.loc['CO2','M'])
    M_CH4 = float(mdv.species.loc['CH4','M'])
    M_S2 = float(mdv.species.loc['S2','M'])
    M_SO2 = float(mdv.species.loc['SO2','M'])
    M_SO3 = float(mdv.species.loc['SO3','M'])
    M_H2S = float(mdv.species.loc['H2S','M'])
    M_FeO = float(mdv.species.loc['FeO','M'])
    M_FeO15 = float(mdv.species.loc['FeO1.5','M'])
    M_OCS = float(mdv.species.loc['OCS','M'])
    M_X = float(mdv.species.loc[species_X,'M'])
    M_m_ = float(mg.M_m_SO(melt_wf))
    M_m_ox = float(mg.M_m_ox(melt_wf,models))
    
    def system(guesses):
        xg_O2_ = guesses[0]
        xg_CO_ = guesses[1]
        xg_S2_ = guesses[2]
        xg_X_ = guesses[3]
        xg_CO2_ = (K2_*y_CO_*xg_CO_*(y_O2_*xg_O2_*P)**0.5)/y_CO2_
        xg_SO2_ = (K6_*y_O2_*xg_O2_*(y_S2_*xg_S2_*P)**0.5)/y_SO2_
        xg_OCS_ = ((xg_CO_*y_CO_)**3.0*xg_SO2_*y_SO2_*P)/(y_OCS_*(xg_CO2_*y_CO2_)**2.0*K10_)
        a = (y_CO2_*xg_CO2_*y_H2O_**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
        b = 1.0 + (y_H2O_/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)) + ((K7_*(y_S2_*xg_S2_)**0.5*y_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5))
        c = xg_CO2_ + xg_CO_ + xg_O2_ + xg_S2_ + xg_SO2_ + xg_OCS_ - 1.0
        xg_H2O_ = (-b + (b**2.0-(4.0*a*c))**0.5)/(2.0*a)
        xg_H2_ = (y_H2O_*xg_H2O_)/(K1_*y_H2_*(y_O2_*xg_O2_*P)**0.5)
        xg_CH4_ = (y_CO2_*xg_CO2_*(y_H2O_*xg_H2O_)**2.0)/(K3_*y_CH4_*(y_O2_*xg_O2_)**2.0)
        xg_H2S_ = (K7_*(y_S2_*xg_S2_)**0.5*y_H2O_*xg_H2O_)/(y_H2S_*(y_O2_*xg_O2_)**0.5)
        Xg_t = xg_CO2_*M_CO2 + xg_CO_*M_CO + xg_O2_*M_O2 + xg_H2O_*M_H2O + xg_H2_*M_H2 + xg_CH4_*M_CH4 + xg_SO2_*M_SO2 + xg_S2_*M_S2 + xg_H2S_*M_H2S + xg_OCS_*M_OCS + xg_X_*M_X
        # melt composition
        xm_H2O_ = (K4_*y_H2O_*xg_H2O_*P)**0.5
        xm_CO2_ = K5_*y_CO2_*xg_CO2_*P
        Xm_t = xm_CO2_*M_CO2 + xm_H2O_*M_H2O + (1.0-xm_CO2_-xm_H2O_)*M_m_
        wm_H2O_ = (xm_H2O_*M_H2O)/Xm_t
        wm_CO2_ = (xm_CO2_*M_CO2)/Xm_t
        wm_CO_ = K11_*y_CO_*xg_CO_*P
        wm_CH4_ = K12_*y_CH4_*xg_CH4_*P
        wm_H2_ = K13_*y_H2_*xg_H2_*P
        wm_H2S_ = K14_*y_H2S_*xg_H2S_*P
        wm_S_ = K8_*((y_S2_*xg_S2_)/(y_O2_*xg_O2_))**0.5
        wm_SO3_ = ((K9_*(y_S2_*xg_S2_)**0.5*(y_O2_*xg_O2_)**1.5*P**2.0)/M_S)*M_SO3
        wm_ST_ = wm_S_ + ((M_S*wm_SO3_)/M_SO3) + ((M_S*wm_H2S_)/M_H2S)
        wm_X_ = K15_*y_X_*xg_X_*P
        Fe32 = ((KD1_*(y_O2_*xg_O2_*P)**0.25)+(2.0*y*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y))))/(1.0 + (1.0 - 2.0*y)*KD2*(KD1_**(2.0*y))*((y_O2_*xg_O2_*P)**(0.5*y)))
        Fe3T = Fe32/(1.0+Fe32)
        S62 = (wm_SO3_/M_SO3)/(wm_S_/M_S)
        S6T = S62/(1.0+S62)
        return xg_O2_,xg_CO_,xg_S2_,xg_X_,xg_H2_,xg_H2O_,xg_CO2_,xg_OCS_,xg_SO2_,xg_H2S_,xg_CH4_,Xg_t,wm_CO2_,wm_CH4_,wm_CO_,wm_H2O_,wm_H2_,wm_SO3_,wm_S_,wm_H2S_,wm_ST_,wm_X_,xm_H2O_,xm_CO2_,Xm_t,Fe32,Fe3T,S62,S6T
    
    def solvers(guesses):
        xg_O2_,xg_CO_,xg_S2_,xg_X_,xg_H2_,xg_H2O_,xg_CO2_,xg_OCS_,xg_SO2_,xg_H2S_,xg_CH4_,Xg_t,wm_CO2_,wm_CH4_,wm_CO_,wm_H2O_,wm_H2_,wm_SO3_,wm_S_,wm_H2S_,wm_ST_,wm_X_,xm_H2O_,xm_CO2_,Xm_t,Fe32,Fe3T,S62,S6T = system(guesses) 
        C = ((wt_C/M_C) - (wm_CO2_/M_CO2) - (wm_CH4_/M_CH4) - (wm_CO_/M_CO))/(((xg_CO2_+xg_CO_+xg_CH4_+xg_OCS_)/Xg_t) - (wm_CO2_/M_CO2) - (wm_CH4_/M_CH4) - (wm_CO_/M_CO))
        O = ((wt_O/M_O) - (wm_H2O_/M_H2O) - ((2.0*wm_CO2_)/M_CO2) - (3.0*wm_SO3_/M_SO3) - (wm_CO_/M_CO) - (wt_Fe/M_Fe)*((1.5*Fe32+1.0)/(Fe32+1.0)))/(((2.0*xg_CO2_ + xg_CO_ + 2.0*xg_O2_ + xg_H2O_ + 2.0*xg_SO2_ + xg_OCS_ )/Xg_t) - (wm_H2O_/M_H2O) - ((2.0*wm_CO2_)/M_CO2) - (3.0*wm_SO3_/M_SO3) - (wm_CO_/M_CO))
        H = ((wt_H/(2.0*M_H)) - (wm_H2O_/M_H2O) - (wm_H2_/M_H2) - (2.*wm_CH4_/M_CH4) - (wm_H2S_/M_H2S))/(((xg_H2O_+xg_H2_+2.0*xg_CH4_+xg_H2S_)/Xg_t) - (wm_H2O_/M_H2O) - (wm_H2_/M_H2) - (2.*wm_CH4_/M_CH4) - (wm_H2S_/M_H2S))
        S = ((wt_S/M_S) - (wm_S_/M_S) - (wm_SO3_/M_SO3) - (wm_H2S_/M_H2S))/(((xg_SO2_+2.0*xg_S2_+xg_H2S_+xg_OCS_)/Xg_t) - (wm_S_/M_S) - (wm_SO3_/M_SO3) - (wm_H2S_/M_H2S))
        X = ((wt_X/M_X) - (wm_X_/M_X))/(((xg_X_)/Xg_t) - (wm_X_/M_X))
        return C-O,C-H,C-S,C-X,C

    
    def solve_system(guesses):        
        CO,CH,CS,CX,C = solvers(guesses)
        return [float(CO),float(CH),float(CS),float(CX)]
    
    newguess = optimize.fsolve(solve_system, [float(guessx),float(guessy),float(guessz),float(guessw)])
    xg_O2_,xg_CO_,xg_S2_,xg_X_,xg_H2_,xg_H2O_,xg_CO2_,xg_OCS_,xg_SO2_,xg_H2S_,xg_CH4_,Xg_t,wm_CO2_,wm_CH4_,wm_CO_,wm_H2O_,wm_H2_,wm_SO3_,wm_S_,wm_H2S_,wm_ST_,wm_X_,xm_H2O_,xm_CO2_,Xm_t,Fe32,Fe3T,S62,S6T = system(newguess)
    CO,CH,CS,CX,wt_g = solvers(newguess)
    return xg_O2_,xg_CO_,xg_S2_,xg_X_,xg_H2_,xg_H2O_,xg_CO2_,xg_OCS_,xg_SO2_,xg_H2S_,xg_CH4_,Xg_t,wm_CO2_,wm_CH4_,wm_CO_,wm_H2O_,wm_H2_,wm_SO3_,wm_S_,wm_H2S_,wm_ST_,wm_X_,xm_H2O_,xm_CO2_,Xm_t,Fe32,Fe3T,S62,S6T,wt_g