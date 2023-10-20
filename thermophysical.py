import math as mat


def Pw(tpwc): # Returns the saturated vapor pressure of pure water (Pa), tpwc in °C (http://dx.doi.org/10.1016/j.desal.2016.02.024)
     tpw = tpwc+273.15
     return mat.exp((-5800 / tpw) + 1.3915 - 0.04864 * tpw + 4.1765 * 10 ** -5 * tpw ** 2 - 1.4452 * 10 ** -8 * tpw ** 3 + 6.546 * mat.log(tpw))

def Psw(tPswc, psal): #Calculating the saturated vapor pressure of pure water (Pa) according to salinity psal (g/kg)  for temperature tPswc in °C   (http://dx.doi.org/10.1016/j.desal.2016.02.024)
     pw_t = Pw(tPswc)
     return (mat.exp(-4.5818 * 10 ** -4 * psal - 2.0443 * 10 ** -6 * psal ** 2)) * pw_t

ttd_C=70
saltd=35
RH_ftd=0.9

ttd = ttd_C #+ 273.15
if saltd == 0:
    psw_tss = Pw(ttd)

else:
    psw_tss = Psw(ttd, saltd)

psw_ts = RH_ftd * psw_tss


if ttd >= 0:
    Tdb_asr = 6.54 + 14.526 * mat.log(psw_ts * (10 ** (-3))) + 0.7389 * ((mat.log(psw_ts * (10 ** (-3)))) ** 2) + 0.09486 * ((mat.log(psw_ts * 10 ** (-3))) ** 3) + 0.4569 * ((psw_ts * (10 ** (-3))) ** 0.1984)
else:
    Tdb_asr = 6.09 + 12.608 * mat.log(psw_ts * (10 ** (-3))) + 0.4959 * ((mat.log(psw_ts * (10 ** (-3)))) ** 2)


def Twb_stull(ttw, RH_ftw):
            # Calculating the temperature in °C of the wet bulb according to Stull's equation (https://doi.org/10.1175/JAMC-D-11-0143.1),
            # given the temperature in °C and the humidity is entered as a fraction but in the formula it is transformed to a percentage
       arg1=0.151977 * (RH_ftw * 100 + 8.313659) ** 0.5
       arg2=ttw + RH_ftw * 100
       arg3=RH_ftw * 100 - 1.676331
       arg4=0.023101 * RH_ftw * 100
       return ttw * mat.atan(arg1) + mat.atan(arg2) - mat.atan(arg3) + 0.00391838 * (RH_ftw * 100) ** 1.5 * mat.atan(arg4) - 4.686035

twbul=Twb_stull(70, 0.9)


# Calculating specific heat of water vapor in air (J/kg K) with
# the air temperature tvap in °C
def CpVapor(tvap):
    Cpvap = (1.3605 * (10 ** (3))) + (2.31334 * (tvap + 273.15)) - ((2.46784 * (10 ** (-10))) * ((tvap + 273.15) ** 5))+((5.91332 * (10 ** (-13))) * ((tvap + 273.15) ** 6))
    return Cpvap

# Calculating specific heat of dry air (J/kg K) with
# the air temperature tcpac in °C
def Cp_air(tcpac):
    tcpa = tcpac + 273.15
    Cp_a = ((1.045356 * 10 ** 3)) - ((3.161783 * 10 ** (-1)) * tcpa) + ((7.083814 * 10 ** (-4)) * tcpa ** 2) - ((2.705209 * 10 ** (-7)) * tcpa ** 3)
    return Cp_a

def Cons_air(tgkc, prefg):
               # Calculating the g1-g8 constants, according to Perez-Galindo et al., equations (https://www.techstreet.com/standards/lb-07-047-thermodynamic-properties-for-saturated-air-an-engineering-correlation?product_id=1712682)

        tgk = tgkc + 273.15


        g1 = 2.2770286 - 0.02406584 * tgk + 1.8213945 * 10 ** (-4) * tgk ** 2 - 6.8894708 * 10 ** (
             -7) * tgk ** 3 + 1.297668 * 10 ** (-9) * tgk ** 4 - 9.7078508 * 10 ** (-13) * tgk ** 5 + 3.9945654 * 10 ** (
             -8) * prefg
        g2 = 1.9208388 - 0.018226313 * tgk + 1.4278754 * 10 ** (-4) * tgk ** 2 - 5.5526317 * 10 ** (
             -7) * tgk ** 3 + 1.073933 * 10 ** (-9) * tgk ** 4 - 8.2740746 * 10 ** (-13) * tgk ** 5 + 3.9755302 * 10 ** (
             -9) * prefg
        g3 = 1.0041923 * (tgk - 273.15)
        g4 = 1.8642569 * (tgk - 273.15)
        g5 = 1.0041923 * mat.log(tgk)
        g6 = 1.8642569 * mat.log(tgk)
        g7 = -4.1281106 * 10 ** (-18) * prefg ** 3 + 1.8191667 * 10 ** (-12) * prefg ** 2 - 3.0990568 * 10 ** (
             -6) * prefg + 0.28844128
        g8 = 1.3323852 * 10 ** (-23) * prefg ** 4 - 6.0489891 * 10 ** (-18) * prefg ** 3 + 1.0524313 * 10 ** (
            -12) * prefg ** 2 - 8.897083 * 10 ** (-8) * prefg + 3.2465229 * 10 ** (-3)

        val=[g1,g2,g3,g4,g5,g6,g7,g8]
        return val


def Xsw_sns(txmc,salxm,paxm,RH_xmf):  # Calculating water mole fraction of pure water vapor (salxm=0),  or for seawater (salxm > 0): (ASHRAE Fundamentals Handbook)
    pxm_s = 0
    fxm = 0
    try:
# With respect to salinity salxm g/kg given at temperature txmc in °C
# a pressure paxm (Pa) and relative humidity RH_xmf as a fractional value

# Calculating the value of f. If the relative humidity is 1 f = g1,
# if the relative humidity is between 0 and less than 1, f = 1

      if RH_xmf == 1:
         fxmi = Cons_air(txmc, paxm)
         fxm = fxmi[0]
      elif RH_xmf > 0 and RH_xmf < 1:
         fxm = 1

    # ************************************************
    # Calculating the saturated vapor pressure
      if salxm == 0:
          pxm_s = Pw(txmc)
      elif salxm > 0:
          pxm_s = Psw(txmc, salxm)

      pxm = RH_xmf * pxm_s
      Xsw_sns = fxm * pxm / paxm
      return Xsw_sns

    except ZeroDivisionError:
          print("The values entered generate a division by 0, which is not admissible.")
    except:
        print("Wrong parameters have been entered")




def Xa_sns(taxm,salaxm,paaxm,RH_axmf):  # Calculating air mole fraction of pure water vapor (salxm=0),  or for seawater (salxm > 0): (ASHRAE Fundamentals Handbook)
                                        # aith respect to salinity salaxm g/kg given at temperature taxm in °C
                                        # a pressure paaxm (Pa) and relative humidity RH_axmf as a fractional value

                                        # Calculating the value of f. If the relative humidity is 1 f = g1,
                                        # if the relative humidity is between 0 and less than 1, f = 1.
                                        #  f is the diffusion factor considered in saturated humid air.
      pxma_s=0
      fxma=0
      try:
          if RH_axmf == 1:
              fxmai = Cons_air(taxm, paaxm)
              fxma = fxmai[0]
          elif RH_axmf > 0 and RH_axmf < 1:
              fxma = 1

       # ************************************************
       # Calculating the saturated vapor pressure

          if salaxm == 0:
                pxma_s = Pw(taxm)
          elif salaxm > 0:
                pxma_s = Psw(taxm, salaxm)
          elif salaxm < 0:
                print("Salinity must be a number greater than or equal to zero, given in g/kg in this case.")

          pxma = RH_axmf * pxma_s

          Xa_sns = mat.fabs(paaxm - fxma * pxma) / paaxm
          return Xa_sns
      except ZeroDivisionError:
             print("The values entered generate a division by 0, which is not admissible.")
      except:
             print("Wrong parameters have been entered")


#ddd
def Water_air(tws, salws, paws,RH_fws):  # Calculating absolute humidity or humidity ratio (kg water /kg air) according to salinity salws(g/kg)
                                         # temperature tws °C
                                         # Pressure paws (Pa), Relative humidity RH_fws in fraction (-)
       xsw = Xsw_sns(tws, salws, paws, RH_fws)
       xasw = Xa_sns(tws, salws, paws, RH_fws)
       print("The intermediate values for calculating absolute humidity are: ",xsw," - ",xasw)
       res_air = xsw * 0.621945 / xasw
       return res_air

                                        #Calculating the specific heat capacity of moist air (J/kg K), whether saturated or unsaturated,
                                         # or moisture from pure water or seawater, for a given pressure paws (Pa), salinity salws in (g/kg),
                                         # temperature tws in °C, and relative humidity RH_fws in fraction (-)

def Cp_ma(tws, salws, paws, RHfws):
          Cpa = Cp_air(tws)
          Cpv = (CpVapor(tws))
          Wabs = Water_air(tws, salws, paws, RHfws)
          CapOK = Cpa + Wabs * Cpv
          return CapOK

# Calculating enthalpy of moist air (J/kg) for a given pressure paha (Pa), salinity salha in (g/kg),
# temperature tha in °C, and relative humidity RH_fha in fraction (-)
#Pérez-Galindo et al.,  https://www.techstreet.com/standards/lb-07-047-thermodynamic-properties-for-saturated-air-an-engineering-correlation?product_id=1712682

def H_air(tha, salha, paha, RH_fha):
         tha_c = tha
         Fh_ref = 1551.57881 - 103454.917 * 10 ** (-6) * tha_c - 30029.857 * 10 ** (-6) * tha_c ** 2 + 501.77652 * 10 ** (
             -6) * tha_c ** 3 - 2.83396 * 10 ** (-6) * tha_c ** 4
         gh = Cons_air(tha, paha)
         gh3 = gh[2]
         gh4 = gh[3]
         gh7 = gh[6]
         xha = Xsw_sns(tha, salha, paha, RH_fha)
         xaha = Xa_sns(tha, salha, paha, RH_fha)
         H_a = (gh3 + (xha / xaha) * (gh4 + Fh_ref) + gh7) * 1000
         return H_a

         # Calculating entropy of moist air (J/kg K), either saturated or unsaturated, or moisture from pure water or seawater,
         # for a given pressure pasa (Pa), salinity salsa in (g/kg), temperature tsa in °C and relative humidity RH_fsa in fraction
def S_air(tsa, salsa, pasa, RH_fsa):
         tsa_c = tsa
         Fs_wref = -633481.3313 * 10 ** (-5) + 923.94626 * 10 ** (-5) * tsa_c - 50448.581 * 10 ** (
             -8) * tsa_c ** 2 + 835.5436 * 10 ** (-8) * tsa_c ** 3 - 5.03383 * 10 ** (-8) * tsa_c ** 4
         gs = Cons_air(tsa, pasa)
         gs2 = gs[1]
         gs5 = gs[4]  #Cons_air(tsa, pasa, 5)
         gs6 = gs[5]  #Cons_air(tsa, pasa, 6)
         gs8 = gs[7]  #Cons_air(tsa, pasa, 8)
         xsa = Xsw_sns(tsa, salsa, pasa, RH_fsa)
         xasa = Xa_sns(tsa, salsa, pasa, RH_fsa)
         S_a = (gs5 - 5.63354 + (xsa / xasa) * (gs6 + Fs_wref) - (0.287042 / xasa) * mat.log(
             pasa / 101325) + (0.287042 / xasa) * mat.log(
             gs2 / xasa) + (xsa / xasa) * 0.287042 * mat.log(gs2 / xsa) + gs8)*1000
         return S_a


def Tdb_asrhae(ttd_C, saltd,RH_ftd):  # Calculating the dew temperature °C according to ASRHAE 2017 equation,
                                      # temperature ttd_C in °C and rel. humidity RH_ftd in fraction, according to water salinity (g/kg).
    psw_tss = 0
    Tdb_asr = 0
    if saltd == 0:
        psw_tss = Pw(ttd_C)

    elif saltd > 0:
        psw_tss = Psw(ttd_C, saltd)

    psw_ts = RH_ftd * psw_tss

    if ttd >= 0:
        Tdb_asr = 6.54 + 14.526 * mat.log(psw_ts * 10 ** (-3)) + 0.7389 * (mat.log(psw_ts * 10 ** (-3))) ** 2 + 0.09486 * (
            mat.log(psw_ts * 10 ** (-3))) ** 3 + 0.4569 * (psw_ts * 10 ** (-3)) ** 0.1984

    elif ttd < 0:
        Tdb_asr = 6.09 + 12.608 * mat.log(psw_ts * 10 ** (-3)) + 0.4959 * (mat.log(psw_ts * 10 ** (-3))) ** 2

    return Tdb_asr

    # Functions for calculating the thermophysical properties of water:
    # http://dx.doi.org/10.1016/j.desal.2016.02.024  and http://web.mit.edu/seawater/2017_MIT_Seawater_Property_Tables_r2b.pdf
    # https://doi.org/10.5004/dwt.2010.1079

def hw(tcel):  # Calculation of the enthalpy of pure water (J/kg), temperature tcel in °C is entered.
      return 141.355 + 4202.07 * tcel - 0.535 * tcel ** 2 + 0.004 * tcel ** 3

def hw_s(tcels, sa):  # Calculation of the enthalpy of salt water (J/kg), enter temperature tcels in °C and salinity sa in g/kg
      hw_t = hw(tcels)
      sal = sa / 1000
      return hw_t - sal * (-2.34825 * 10 ** 4 + 3.15183 * 10 ** 5 * sal + 2.80269 * 10 ** 6 * sal ** 2 - 1.44606 * 10 ** 7 * sal ** 3 + 7.82607 * 10 ** 3 * tcels - 4.41733 * 10 ** 1 * tcels ** 2 + 2.1394 * 10 ** (
                -1) * tcels ** 3 - 1.99108 * 10 ** 4 * sal * tcels + 2.77846 * 10 ** 4 * sal ** 2 * tcels + 9.72801 * 10 ** 1 * sal * tcels ** 2)


def Sw(stcel):  # Calculation of entropy of pure water (J/(kg*K)), stcel temperature in °C
      return 0.1543 + 15.383 * stcel - 2.996 * 10 ** -2 * stcel ** 2 + 8.193 * 10 ** -5 * stcel ** 3 - 1.37 * 10 ** -7 * stcel ** 4

def Sw_s(stcels, sa):  # Calculation of the entropy of salt water (J/(kg*K)), enter stcels temperature in °C and salinity sa in g/kg
       sw_t = Sw(stcels)
       ssal = sa / 1000
       return sw_t - ssal * (-4.231 * 10 ** 2 + 1.463 * 10 ** 4 * ssal - 9.88 * 10 ** 4 * ssal ** 2 + 3.095 * ssal ** 3 + 2.562 * 10 * stcels - 1.443 * 10 ** -1 * stcels ** 2 + 5.879 * 10 ** -4 * stcels ** 3 - 6.111 * 10 * ssal * stcels + 8.041 * 10 * ssal ** 2 * stcels + 3.035 * 10 ** -1 * ssal * stcels ** 2)

def Lw(tLw):  # By calculating the latent heat of vaporization of pure water (J/kg) given a temperature in °C
      return ((2.501 * 10 ** 6)) - ((2.369 * 10 ** 3) * tLw) + ((2.678 * 10 ** -1) * tLw ** 2) - (
              (8.103 * 10 ** -3) * tLw ** 3) - ((2.079 * 10 ** -5) * tLw ** 4)

def Lsw(tLsw,  salLswc):  # Calculating the latent heat of vaporization of water (J/kg) with salLswc salinity (g/kg) given a temperature in °C
       lw_st = Lw(tLsw)
       if salLswc == 0:
           hswt=lw_st
       elif salLswc < 0:
           print("Salinity cannot be a value less than zero.")
       else:
           salLsw = salLswc / 1000
           hswt = lw_st * (1 - salLsw)

       return hswt

def Cw_s(tcswc,ssw):  # Calculating the specific heat of water (J/kg K) according to salinity ssw (g/kg) given a temperature tcswc in °C
      tcsw = tcswc + 273.15
      Ac = 5.328 - 9.76 * (10 ** -2) * ssw + 4.04 * (10 ** -4) * ssw ** 2
      Bc = -6.913 * (10 ** -3) + 7.351 * (10 ** -4) * ssw - 3.15 * (10 ** -6) * ssw ** 2
      Cc = 9.6 * (10 ** -6) - 1.927 * (10 ** -6) * ssw + 8.23 * (10 ** -9) * ssw ** 2
      Dc = 2.5 * (10 ** -9) + 1.666 * (10 ** -9) * ssw - 7.125 * (10 ** -12) * ssw ** 2
      Cw_sr = (Ac + Bc * tcsw + Cc * tcsw ** 2 + Dc * tcsw ** 3) * 1000
      return Cw_sr

def tempbulbohumediter(ta,sala,presa,rha):   # Calculate the wet bulb temperature °C by iteration
      ha = H_air(ta,sala,presa,rha)
      wa = Water_air(ta,sala,presa,rha)

      taa = int(((ta + 1)//0.01)+10)
      err = 0
      errc = 1000
      twb = ta

      if rha >= 1:
          twb = ta
      else:
          for it in range(0,taa,1):
              itt=it*0.01
              has = H_air(itt, sala, presa, 1)
              was = Water_air(itt, sala, presa, 1)
              hwatv = hw_s(itt,sala)
              valiz = ha+(was-wa)*hwatv
              valdr = has
              err=(mat.fabs(valdr-valiz))*100/valdr
              if err <= errc:
                  twb = itt
                  errc = err


      return twb


def Iter_H_air(Orig_H_air, IHa_salha, IHa_paha, IHa_RH_fha): # Determine temperature in °C  given a value of specific enthalpy of air Orig_H_air J/kg,
                                                            # salinity of origin IHa_salha g/kg, pressure IHa_paha (Pa) and relative humidity IHa_RH_fha in fraction (-)

     tmax = int((100 // 0.01) + 1)
     terr = 0
     tsal = 0
     terrc = 1000

     for ith in range(0, tmax, 1):
         itht = ith * 0.01
         Cheq_H_air = H_air(itht, IHa_salha, IHa_paha, IHa_RH_fha)
         terr = mat.fabs((Orig_H_air - Cheq_H_air)) * 100 / Orig_H_air
         if terr < terrc:
             tsal = itht
             terrc = terr

     return tsal

def Iter_H_and_w_air(Orig_H_air, IHa_salha, IHa_paha,wref):  # Determine temperature in °C  given a value of specific enthalpy of air Orig_H_air J/kg,
                                                            # salinity of origin IHa_salha g/kg, pressure IHa_paha (Pa) and relative humidity IHa_RH_fha in fraction (-) and absolute humidity wref (kg water / kg air)
        lres=[1,1]
        tmax = int(1 // 0.01)
        terrc = 1000

        for itw in range(0,tmax,1):
            IHa_RH_fha=1-itw*0.01
            itht =Iter_H_air(Orig_H_air,IHa_salha,IHa_paha,IHa_RH_fha)
            Cheq_W_air = Water_air(itht, IHa_salha, IHa_paha, IHa_RH_fha)
            terr = mat.fabs((wref - Cheq_W_air)) * 100 / wref
            if terr < terrc:
                tsal = itht
                rhsal = IHa_RH_fha
                lres.insert(0, tsal)
                lres.insert(1, rhsal)
                lres.pop(2)
                lres.pop(2)
                terrc = terr

        return lres

#@@@@@@@@@@@@@@@@@******************************************************@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Start of calculation with Runge_Kutta For cooling towers (humidifier) according to Pope's method
# https://doi.org/10.1016/j.desal.2015.04.021
# https://doi.org/10.1016/j.ijheatmasstransfer.2004.09.004
# http://hdl.handle.net/10019.1/1476

def Runge_Kutta_merk_pop_basic(twat,ha,wa,sal,pres,tfmw,tfma): # Calculating the basic equations of the cooling tower or humidifier process
        #Runge_Kutta_merk_pop_interm(ltw[git-1],thpas[git-1],twpas[git-1],lsal[git-1],pres,ltmw[git-1],tfma)
        # (twat,tair,rhair,sal,pres,tfmw,tfma)      (twat,rhair,rwair,sal,pres,tfmw,tfma)
        # twat: water temperature °C
        # tair: air temperature °C
        # rhair: relative humidity of air (-)
        # hair:  air enthalpy (J/kg)
        #rwair: absolute humidity (kg water / kg air)
        # sal: salinity of the water that generates the humidity in the analized air (g/kg)
        # pres: air pressure (Pa)
        # tfmw,tfma mass flow rates of water and air respectively (kg/s)

        wsw = Water_air(twat,sal,pres,1)      # This is to calculate the Lewis factor for both saturated and unsaturated air (lef),
                                              # and that of supersaturated air (lefs).

        tPsw = Psw(twat,sal)  # This is to recalculate relative humidity
        rhair = (wa * pres) / (tPsw * (0.62198 + wa))  # This is to recalculate relative humidity too

        # It is now possible to iterate to find the temperature of the air.
        tair = Iter_H_air(ha, sal,pres, rhair)

        wss = Water_air(tair, sal, pres, 1)
        le = (wsw + 0.622) / (wa + 0.622)
        les = (wsw + 0.622) / (wss + 0.622)


        lef = (0.865 ** 0.667) * (le - 1) * ((mat.log(le)) ** -1)        # Factores de Lewis antes referidos
        lefs = (0.865 ** 0.667) * (les - 1) * ((mat.log(les)) ** -1)

        cpw = Cw_s(twat,sal)
        coe1 = tfmw * cpw / tfma
        dwt = wsw - wa
        dwts = wsw - wss
        hasw = H_air(twat, sal, pres, 1)
        has = H_air(tair, sal, pres, 1)
        hv = Lsw(tair, sal) + CpVapor(tair) * tair
        dg = hasw - ha + (lef - 1) * (hasw - ha - dwt * hv) - dwt * cpw * twat
        dgs = hasw - has + (lefs - 1) * (hasw - has - dwts * hv) - dwts * cpw * twat



        vw = coe1 * (dwt / dg)                # Calculating absolute humidity of saturated or unsaturated air
        vws = coe1 * (dwts / dgs)             # Calculating absolute humidity of supersaturated air


        vha = coe1 * (1 + (dwt * cpw * twat / dg))  # Calculating enthalpy of saturated or unsaturated air
        vhas = coe1 * (1 + (dwts * cpw * twat / dgs))  # Calculating enthalpy of supersaturated air

        vme = cpw / dg                  # Calculating Merkel number of saturated or unsaturated air
        vmes = cpw / dgs                # Calculating Merkel number of supersaturated air


        # It is tested whether the air is supersaturated or not.
        tessatur = 0
        twb = tempbulbohumediter(tair,sal,pres,rhair)    #  Wet bulb temperature is determined

        result = [1,1,1]

        if twb < tair:
            result.insert(0,vw)
            result.insert(1, vha)
            result.insert(2, vme)

            result.pop(3)
            result.pop(3)
            result.pop(3)

        elif twb >= tair:
            result.insert(0, vws)
            result.insert(1, vhas)
            result.insert(2, vmes)

            result.pop(3)
            result.pop(3)
            result.pop(3)
        print("The current iteration has been completed")
        return result

# *******************************************************************************************************************


def Runge_Kutta_merk_pop_values(twat,tair,rhair,sal,pres,tfmw,tfma,tfmb,twin,iter): # Se calculan los coeficientes y los pasos dando como resultado los valores de ha - wa - Me

       # twat Water temperature °C
       # tair Air temperature °C
       # rhair Relative humidity (-)
       # sal salinity of the water that generates the humidity in the analized air (g/kg)
       # pres Air pressure (Pa)
       # tfmw,tfma mass flow rates of water and air respectively (kg/s): tfma = Constant, tfmw It changes
       # tfmb brine water flow rate at humidifier outlet (kg/s)
       # twin water temperature at the inlet of the humidifier or cooling tower °C
       # iter proposed number of iterations

       dtw = (twin - twat) / (iter-1)        # Twat has been entered as the water temperature at the outlet of the humidifier,
                                         # thus dtw is determined as the variation of the water temperature at each step.
       salou = (tfmw / tfmb) * sal       # Salinity at the inlet
       lsal = [salou]                    # List to store the water salinity values at each step, here the index value 0 is given.
       dtsal = (sal - salou) / (iter-1)      # Salinity variation at each step
       dtmw = (tfmw - tfmb) / (iter-1)       # The magnitude dtmw refers to the variation of the water flow rate at each step.
       ltw = [twat]                      # List for storing the water temp. values at each step, index 0 is given here.
       ltmw = [tfmb]                     # List to store the water flow rate values at each step, index 0 is given here.
       h0 = H_air(tair, sal, pres, rhair)  # The initial values of enthalpy h0 and absolute humidity w0 are calculated according to air conditions.
       w0 = Water_air(tair, sal, pres, rhair)   # at humidifier inlet, "tair, salt, pres, rhair" mainly tair temperature and relat. humidity rhair
       thpas = [h0]               # List to store the air enthalpy values at each step, index 0 is given here.
       twpas = [w0]               # List to store the absolute humidity values of the air in each step, here the value 0 is given.
       tme = [0]                  # List to store the values of the Merkel number in each step, here index 0 is given a value.
       resfin = {}        # What is expected is to group the results of the different iterations in this dictionary,
                          # de modo tal que cada indice se denomine "clavn"  y se relacione con una lista como registro,
                          # list containing the values of w, h and Me of each iteration ":[wn,hn,Men]" where n goes from 0 to iter-1

       clv = 'clav' + str(0)   # This initializes the dictionary of results.
       resfin[clv] = [twpas[0], thpas[0], tme[0]]

           # For the time being, no pressure drops are considered in this context.

       for iin in range(1,iter,1):
           ltw.insert(iin,twat + iin * dtw)            # Values that the water temperature takes at each step
           ltmw.insert(iin,tfmb + iin * dtmw)          # Values that the water flow rate takes at each step
           lsal.insert(iin,salou - iin * dtsal)        # Values that the salinity takes at each step

       for git in range(1,iter,1):
           valcoe1 = Runge_Kutta_merk_pop_basic(ltw[git-1],thpas[git-1],twpas[git-1],lsal[git-1],pres,ltmw[git-1],tfma)
           jnn1 = dtw * valcoe1[0]
           knn1 = dtw * valcoe1[1]
           lnn1 = dtw * valcoe1[2]

           valcoe2 = Runge_Kutta_merk_pop_basic(ltw[git-1]+dtw/2,thpas[git-1]+knn1/2,twpas[git - 1] + jnn1 / 2,lsal[git-1],pres,ltmw[git-1],tfma)
           jnn2 = dtw * valcoe2[0]
           knn2 = dtw * valcoe2[1]
           lnn2 = dtw * valcoe2[2]

           valcoe3 = Runge_Kutta_merk_pop_basic(ltw[git-1] + dtw / 2, thpas[git-1] + knn2 / 2, twpas[git-1] + jnn2 / 2, lsal[git-1], pres, ltmw[git-1], tfma)
           jnn3 = dtw * valcoe3[0]
           knn3 = dtw * valcoe3[1]
           lnn3 = dtw * valcoe3[2]

           valcoe4 = Runge_Kutta_merk_pop_basic(ltw[git-1] + dtw, thpas[git-1] + knn3, twpas[git-1] + jnn3, lsal[git-1], pres, ltmw[git-1], tfma)
           jnn4 = dtw * valcoe4[0]
           knn4 = dtw * valcoe4[1]
           lnn4 = dtw * valcoe4[2]

           twpasn = twpas[git-1] + (jnn1 + 2 * jnn2 + 2 * jnn3 + jnn4) / 6
           thpasn = thpas[git-1] + (knn1 + 2 * knn2 + 2 * knn3 + knn4) / 6
           tmen = tme[git-1] + (lnn1 + 2 * lnn2 + 2 * lnn3 + lnn4) / 6

           twpas.insert(git,twpasn)
           thpas.insert(git,thpasn)
           tme.insert(git,tmen)

           clv = 'clav' + str(git)
           resfin[clv]=[twpasn,thpasn,tmen]
           print("Absolute humidity values calculated in basic Runge-Kutta",twpas)
       return resfin

#******@@@@@@ End Runge-Kutta***********************@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

def entuhumid(twi,two,tai,tmw,tma,sall,press,trhai,ee):  # Calculating Merkel number by the e-ntu method (Pending of a further validation )
       # twi water temperature at humidifier inlet °C
       # two water temperature at humidifier outlet °C
       # tai air temperature at humidifier inlet °C
       # tmw  mass flow rate of water kg/s
       # tma  mass flow rate of air kg/s
       # sall salinity of the water that generates the humidity in the analized air (g/kg)
       # trhai Relative humidity (-)
       # ee   Humidifier effectiveness (-)

       wCp = Cw_s(twi, sall)

       dcpw = tmw * wCp

       dcpa = tma * Cp_ma(tai, sall, press, trhai)

       haswi = H_air(twi,sall, press, 1)

       haswo = H_air(two, sall, press, 1)


       global cmin,cmax,me
       if dcpw < dcpa:
           cmin=dcpw
           cmax=dcpa
       else:
           cmin = dcpa
           cmax = dcpw

       c = cmin / cmax
       ntu = (1/(1 - c)) * mat.log((1 - ee * c) / (1 - ee))

       dhwtw = (haswi - haswo) / (twi - two)
       if tma > (dcpw / dhwtw):
           me = wCp * ntu / dhwtw
       elif tma < (dcpw / dhwtw):
           me = tma * ntu / tmw

       return me

# *******************************************************************************************************************
# Calculating dehumidifier geometry
def Triboix_ntu_de(twid,taid,tmw,tma,sall, press, trhaid, Efec): # Calculating NTU values for the dehumidifier
                                                                 # https://doi.org/10.1016/j.icheatmasstransfer.2008.10.012
         c = 0.44
         d = 0.92
         n = 1.25
         q = 0.15
         ss = 1.15
         pi = mat.pi         #3.14159265358979

         cpw = Cw_s(twid, sall)
         cpa = Cp_ma(taid, sall, press, trhaid)
         dcpw = tmw * cpw
         dcpa = tma * cpa

         if dcpw > dcpa:
             cmin=dcpa
             cmax=dcpw
         else:
             cmin = dcpw
             cmax = dcpa

         rc = cmin / cmax

         if rc > 0.3 and Efec > 0.52:
             ntu1a = 1/(pi * (rc ** q))
             ntu1b = (1 + c * (1 - rc)) / (1 - Efec + c * (1 - rc))
             ntu1c = (ntu1b ** (2 * n)) - d
             ntu1d = ntu1c ** (1 / n)
             ntu1 = ntu1a * ntu1d

             ntu = ntu1
         else:
             ntu2a = 1+(rc**ss)*(mat.log(1-Efec))
             ntu2 = (-1*mat.log(ntu2a))/(rc**ss)

             ntu = ntu2

         print("The NTU value is: ", ntu)
         return ntu

# Calculation of areas of main parts (Humidifier, dehumidifier, solar water heater)

def AltHumid(n1,n2,n3,n4,tmw,tma,me):         # Calculate the humidifier head according to nx coefficients given by the packed bed manufacturer.
       #n1,n2,n3,n4: The nx coefficients given by the packed bed manufacturer

       #               n2            n4       From this equation,
       # me = n1  *  mr  *  n3  *  Hh         Hh can be solved

       mr=tmw/tma

       Hh=(me/((n3**n4)*n1*(mr**n2)))**(1/n4)
       print("The packed bed height is: ",Hh," m")
       return Hh

def AreaRellEvapHum(Hh, fm, tmsw):  # Calculate the perimeter area of the humidifier (its packed bed).
    # Hh packed bed height
    # fm:   It is the mass flow of water at the inlet of the humidifier,
    # recommended by several authors to be always 1.5 kg/(s*m2).      https://doi.org/10.1016/j.desal.2014.06.016
    Ash = tmsw / fm
    Ath = 4 * Hh * (Ash ** 0.5) + 2 * Ash
    print("The surface area of the humidifier's packed bed is: ", Ath," m2")
    return Ath


def AreaCalSo(Qi, Isp, CoefC):  # Calculating the area of the solar water heater
    # Isp: Average solar radiation intensity (W/m2) depends on the place where the equipment is installed (it can be around 1000 W/m2).
    # CoefC: Coefficient of conversion of solar energy into heat by the device or solar heater (-), can reach at most 60%.
    # Qi: Solar radiation intensity (W), calculated according to system operating parameters,
    # water outlet temperature and humidifier inlet temperature, for a solar water heating system.
    Acs = Qi / (Isp * CoefC)
    print("The effective area of the solar water heater is: ", Acs," m2")
    return Acs


def areacondendeshumid(ntu, u, tmw, tma,twid, taid, sall, press, trhaid):  # Calculating the effective area of the condenser in the dehumidifier
                                                                           # https://doi.org/10.1016/j.desal.2014.06.016
       # ntu: Number of transfer units (-)
       # aed: effective area of the condenser, according to the ntu method effectiveness (m2)
       # u: Overall heat transfer coefficient, given by the characteristics of the dehumidifier heat exchanger (W/(m2*K)).
          # It can be in the range of 5 - 100 W/(m2*K).

       cpw = Cw_s(twid, sall)
       cpa = Cp_ma(taid, sall, press, trhaid)
       dcpw = tmw * cpw
       dcpa = tma * cpa

       if dcpw > dcpa:
           cmin = dcpa
       else:
           cmin = dcpw

       aed=ntu * cmin / u

       print("The effective area of the condenser in the dehumidified system is: ",aed," m2")
       return aed

