
import math

from thermophysical import *
#import numpy as np
import matplotlib.pyplot as plt



class intermediate_calculus:

     def __init__(self, tswid=30, tswih=69, msw=0.0899, RH_ih=0.9, RH_oh=0.9, Efd=0.85, Efh=0.85, Sal=35,
                 Pr_a=101325): # Initialization of variables of intermediate_calculus



          self.tswid = tswid  # Seawater temperature at the dehumidifier inlet (0C)
          self.tswih = tswih  # Seawater temperature at the humidifier inlet (0C)
          self.msw = msw  # Mass flow rate of seawater  (kg/s)
          self.RH_ih = RH_ih  # Relative humidity at the humidifier inlet (-)
          self.RH_oh = RH_oh  # Relative humidity at the humidifier outlet (-)
          self.Efd = Efd  # Effectiveness of the dehumidifier (-)
          self.Efh = Efh  # Effectiveness of the humidifier (-)
          self.Sal = Sal  # Seawater salinity (g/kg)

          self.Pr_a = Pr_a  # Ambient pressure  (Pa)

          self.taid = 75     # Air temperature at the dehumidifier inlet (0C)
          self.taod = 40     # Air temperature at the dehumidifier outlet (0C)
          self.taih = 40     # Air temperature at the humidifier inlet (0C)
          self.taoh = 75     # Air temperature at the humidifier outlet (0C)
          self.tswod = 60     # Seawater temperature at the dehumidifier outlet (0C)
          self.tswoh = 45     # Seawater temperature at the humidifier outlet (0C)

          self.hswid=3000   # Specific entalpy of seawater at the dehumidifier inlet (J/kg)
          self.hswih=8000   # Specific entalpy of seawater at the humidifier inlet (J/kg)
          #self.hswod=6000   # Specific entalpy of seawater at the dehumidifier outlet (J/kg)
          #self.hswoh=4500   # Specific entalpy of seawater at the humidifier outlet (J/kg)
          #self.haid=7500     # Specific entalpy of air at the dehumidifier inlet (J/kg)
          #self.haih=4000     # Specific entalpy of air at the humidifier inlet (J/kg)
          #self.haod=4000     # Specific entalpy of air at the dehumidifier outlet (J/kg)

          self.haoh=7500     # Specific entalpy of air at the humidifier outlet (J/kg)
          self.Woh = 0.04  # Absolute humidity of air at the humidifier outlet (kg air/kg water)
          # self.haoh and self.Woh are intermediate calculus, but they have to be saved in order to the Runge-kutta validation


          #self.Wih=0.02      # Absolute humidity of air at the humidifier inlet (kg air/kg water)

          self.Sswid=300   # Specific entropy of seawater at the dehumidifier inlet [J/(kg*K)]
          self.Sswih=800   # Specific entropy of seawater at the humidifier inlet [J/(kg*K)]
          #self.Sswod=600   # Specific entropy of seawater at the dehumidifier outlet [J/(kg*K)]
          #self.Sswoh=450   # Specific entropy of seawater at the humidifier outlet [J/(kg*K)]
          #self.Said=400     # Specific entropy of air at the dehumidifier inlet [J/(kg*K)]
          #self.Saih=750     # Specific entropy of air at the humidifier inlet [J/(kg*K)]
          #self.Saod=750    # Specific entropy of air at the dehumidifier outlet [J/(kg*K)]
          #self.Saoh=400     # Specific entropy of air at the humidifier outlet [J/(kg*K)]

          #self.hpw_ideal=3000    # Ideal specific enthalpy of pure water if it had the temperature tswid [J/(kg*K)]
          self.haod_ideal=7000    # Ideal specific enthalpy of saturated air at the dehumidifier outlet if it had the temperature tswid [J/(kg*K)]
          #self.tpw=49            # Pure water temperature  (0C)
          #self.hpw=5000          # Specific entalpy of pure water (J/kg)
          #self.Spw=500           # Specific entropy of pure water [J/(kg*K)]
          self.hswod_ideal=6500   # Ideal specific enthalpy of seawater at the dehumidifier outlet if it had the temperature taid (J/kg)
          self.haoh_ideal=4000    # Ideal specific enthalpy of saturated air at the humidifier outlet if it had the temperature tswih
          #self.hswoh_ideal=4300   # Ideal specific enthalpy of seawater at the dehumidifier outlet if it had the temperature taih (J/kg)

          self.ma=0.1                    # Mass flow rate of air  (kg/s)
          self.mpw=self.msw*0.03         # Mass flow rate of pure water  (kg/s)
          #self.mpw_idea=self.msw*0.03   # Ideal mass flow rate of pure water  (kg/s)
          self.mb=self.msw-self.mpw      # Mass flow rate of brine  (kg/s)

          self.ArHu=20            # Perimeter area of the humidifier (m2)
          self.ArEd=10            # Heat transfer area of  dehumidifier (m2)
          self.ArCS=20            # Area of the SWH (m2)
          self.nMe=1              # Merkel number (-)
          self.Hh=1               # Height of the packed bed (m)
          self.qi=2000            # Heat input rate  to the system (W)

          # Variables used to form objective functions

          #********************************************************************

          self.GOR=1              # Gain Output Ratio (-)
          self.RR=0.03            # Recovery Ratio (-)
          self.Atot=35            # Total area of the components (m2)

          self.DEfd=1             # Difference of calculated values of dehumidifier effectiveness by solving system of equations (-)
          self.DEfh=1             # Difference of calculated values of humidifier effectiveness by solving system of equations (-)
          self.pDEfd = 1          # Relative Error: Percent relative to the difference of calculated values of dehumidifier effectiveness by solving system of equations (%)
          self.pDEfh = 1          # Relative Error: Percent relative to the calculated values of humidifier effectiveness by solving system of equations (%)
          self.DEniod = 1         # Difference of calculated values of incoming (i) and outgoing (o) energy in dehumidifier by solving system of equations (W)
          self.DEnioh = 1         # Difference of calculated values of incoming (i) and outgoing (o) energy in humidifier by solving system of equations (W)
          self.pDEniod = 1        # Relative Error: Percent relative to the difference of calculated values of incoming (i) and outgoing (o) energy in dehumidifier by solving system of equations (%)
          self.pDEnioh = 1        # Relative Error: Percent relative to the difference of calculated values of incoming (i) and outgoing (o) energy in humidifier by solving system of equations (%)
          self.DEntiod = 1  # Difference of calculated values of incoming (i) and outgoing (o) specific entropy in dehumidifier by solving system of equations [J/(kg*K)]
          self.DEntioh = 1  # Difference of calculated values of incoming (i) and outgoing (o) specific entropy in humidifier by solving system of equations [J/(kg*K)]
          self.pDEntiod = 1  # Relative Error: Percent relative to the difference of calculated values of incoming (i) and outgoing (o) specific entropy in dehumidifier by solving system of equations (%)
          self.pDEntioh = 1  # Relative Error: Percent relative to the difference of calculated values of incoming (i) and outgoing (o) specific entropy in humidifier by solving system of equations (%)

          self.ValTerm = 0  # Checking the validity of the thermal calculation according to configuration:  elf.ValTerm = 0 non-valid  elf.ValTerm = 1 valid
          self.ValRK = 0    # Checking the validity of the thermal Runge-kutta calculation according to configuration:  self.ValRK = 0 non-valid  self.ValRK = 1 valid
          self.verhp = 0    # Verificar que la altura del relleno evaporador no sea muy grande respecto a el agua pura producida de forma razonable
          self.valt = 2     # Checking the validity of all system calculations:  self.valt = 0 non-valid  self.valt = 1 valid

          # ********************************************************************

     def CheqEnt(self, Ento, Enti):  # Function for: Entropy calculation check in compliance with the laws of physics: (Ent,ou-Ent,in)>= o
               if (Ento - Enti) > 0 and ((mat.fabs(Ento - Enti)) * 100 / Ento) >= 2:
                     ve = 0
               elif (Ento - Enti) > 0 and ((mat.fabs(Ento - Enti)) * 100 / Ento) < 2:
                     ve = 5
               elif (Ento - Enti) == 0:
                     ve = 10
               elif (Ento - Enti) < 0 and ((mat.fabs(Ento - Enti)) * 100 / Ento) < 0.1:
                     ve = 15
               elif (Ento - Enti) < 0 and ((mat.fabs(Ento - Enti)) * 100 / Ento) < 0.25:
                     ve = 25
               elif (Ento - Enti) < 0 and ((mat.fabs(Ento - Enti)) * 100 / Ento) < 0.5:
                     ve = 35
               elif (Ento - Enti) < 0 and ((mat.fabs(Ento - Enti)) * 100 / Ento) < 1:
                     ve = 550
               else:
                   ve = 100000
               return ve

     def CheqEne(self, Eneo, Enei):  # Checking energy and effectiveness to ensure proper conditions:  Eneo reference value and Enei value to compare
         if (Eneo - Enei) == 0:
             vf = 0
         elif (mat.fabs(Eneo - Enei)) * 100 / Eneo < 0.001:
             vf = 0.1
         elif (mat.fabs(Eneo - Enei)) * 100 / Eneo < 0.01:
             vf = 2
         elif (mat.fabs(Eneo - Enei)) * 100 / Eneo < 0.1:
             vf = 10
         elif (mat.fabs(Eneo - Enei)) * 100 / Eneo < 0.25:
             vf = 50
         elif (mat.fabs(Eneo - Enei)) * 100 / Eneo < 0.5:
             vf = 100
         elif (mat.fabs(Eneo - Enei)) * 100 / Eneo < 1:
             vf = 1000
         else:
             vf=6000
         return vf

     def CheqEft(self, Ef1, Ef2, Ef):  # Checking the effectiveness EF is in a given span between Ef1 to Ef2
         if Ef >= Ef1 and Ef <= Ef2:
             rc = 1
         else:
             rc = 0
         return rc

     def calcularhdh(self,vt11,vt21,vt12,vt22,pas,trf1,trf2,trf3,trf4,rrp,fr,csma):

               #  Range of terminal temperature diferences (TTD in 0C)
               # vt11 Lower value of taid
               # vt21 Lower value of taod
               # vt12 Higher value of tswod
               # vt22 Higher value of  tswoh
               # pas: Size of temperature variation steps
               # trf1=taid=taoh -- trf2=taod=taih -- trf3=tswoh -- trf4=tswod. They are used for temperature limiting and to refine the calculation with few iterations.
               # rrp=: Refers to the value + - by which trf* varies
               # fr=0 se trata de una iteración sin acotamiento ----- fr=1 Se está realizando una iteración acotada (se acota con una función externa a esta)
               # calcultermico_fino es la función externa con la que se acota y esta sería una función interna de ella
               # csma: If ma is calculated or is a given value. csma=0 calculated ---- csma=1 a given value
            try:
               sumfin = 100000.00
               self.hswid = hw_s(self.tswid, self.Sal)
               self.hswih = hw_s(self.tswih, self.Sal)
               self.Sswid = Sw_s(self.tswid, self.Sal)
               self.Sswih = Sw_s(self.tswih, self.Sal)

               Lah = Lsw(self.tswid,self.Sal)

               self.haod_ideal = H_air(self.tswid, 0, self.Pr_a, 1)
               self.haoh_ideal = H_air(self.tswih, self.Sal, self.Pr_a, 1)
               thpw_ideal = hw(self.tswid)
               tWih_id = Water_air(self.tswid, 0, self.Pr_a, 1)  # Ideal absolute humidity of the air at the dehumidifier output

               #The first "for" ************************
               #pas=0.25
               vtr1 = int(((vt12-vt11)//pas)+1)
               vtr2 = int(((vt22-vt21)//pas)+1)

               ttaid = 0     # Making these variables as global, they can "enter in and out" of the different indentations.
               ttaod = 0
               ttswohu = 0
               ttswode = 0
               caud_a = 0.01

               for i in range(vt11,vtr1,1):
                   if fr == 0:
                       ttaid = self.tswih - (i*pas)      # Iterations for air temperatures
                   elif fr == 1:
                       ttaid = (trf1 + rrp) - (i * pas)
                   ttaoh = ttaid
                   thaid = H_air(ttaid, self.Sal, self.Pr_a, self.RH_oh)
                   thaoh = thaid
                   tSaid = S_air(ttaid, self.Sal, self.Pr_a, self.RH_oh)
                   tSaoh = tSaid
                   tWoh = Water_air(ttaoh, self.Sal, self.Pr_a, self.RH_oh)
                   thswod_ideal = hw_s(ttaid,self.Sal)  # Ideal enthalpy of water leaving the dehumidifier
                   #End of the first "for"*******

                   #Second "for" ########################333************

                   for j in range(vt21,vtr2,1):
                       if fr == 0:
                           ttaod = self.tswid + (j*pas)  # Iterations for air temperatures
                       elif fr == 1:
                           ttaod = (trf2 - rrp) + (j * pas)
                       ttaih = ttaod
                       thaod = H_air(ttaod, 0, self.Pr_a, self.RH_ih)
                       thaih = thaod  # Enthalpy of air at dehumidifier outlet and humidifier inlet.
                       tSaod = S_air(ttaod, 0, self.Pr_a, self.RH_ih)
                       tSaih = tSaod  # Entropy of air at dehumidifier outlet and humidifier inlet.
                       tWih = Water_air(ttaih, 0, self.Pr_a, self.RH_ih)  # Absolute humidity of the air at the humid inlet.
                       thswoh_ideal = hw_s(ttaih,self.Sal)  # Ideal water enthalpy at humidifier outlet
                                                    # a la entrada del humid.
                       tpw_n = (Tdb_asrhae(ttaid, self.Sal,self.RH_oh) + ttaod) / 2  # The temperature of the pure water obtained is calculated as follows
                                                                          # the average between the dehumidifier inlet dew point temperature
                                                                          # and dehumidifier outlet temperature
                       thpw = hw(tpw_n)  # Taking the temperature of the pure water obtained, entropy and enthalpy of the water are calculated.
                       tSpw = Sw(tpw_n)

                       #End of the second "for" ****************************

                       # Third "for"  ****************************

                       for k in range(vt21,vtr2,1):
                           if fr == 0:
                               ttswohu = ttaih + (k*pas)  # Iterations for salt water temperatures at the humidifier output
                           elif fr == 1:
                               ttswohu = (trf3 - rrp) + (k * pas)
                           thswoh = hw_s(ttswohu,self.Sal)  # This is used to calculate the enthalpy and entropy of the seawater leaving the humidifier.
                           tSswoh = Sw_s(ttswohu, self.Sal)

                           # Fin del tercer for ***#######*********

                           # Cuarto for **************************

                           for l in range(vt11,vtr1,1):
                               if fr == 0:
                                   ttswode = ttaid - (l * pas)       # Iterations for saltwater temperatures at the dehumidifier outlet
                               elif fr == 1:
                                   ttswode = (trf4 + rrp) - (l * pas)
                               thswod = hw_s(ttswode,self.Sal)  # This is used to calculate the enthalpy and entropy of the seawater leaving the dehumidifier.
                               tSswod = Sw_s(ttswode, self.Sal)

                               tDeltHwatMax_de = thswod_ideal-self.hswid  # Maximum enthalpy variation of water in dehumidifier
                               tDeltHairMax_de = thaid-self.haod_ideal-thpw_ideal*(tWoh - tWih_id)  # Maximum air enthalpy variation in the dehumidifier
                               if csma == 0:
                                   caud_a = self.msw * tDeltHwatMax_de / tDeltHairMax_de  # Calculating ma according to the criterion of  HCR=1
                               elif csma == 1:
                                   caud_a = self.ma

                               caud_pw = caud_a * (tWoh - tWih)
                               caud_b = self.msw - caud_pw


                               tDeltHairMax_hu = self.haoh_ideal-thaih     # Maximum air enthalpy variation in humidifier
                               tDeltHwatMax_hu = self.hswih-thswoh_ideal   # Maximum enthalpy variation of water in the humidifier. It is not used.


                               tDeltHair_de = thaid-thaod-thpw*(tWoh - tWih)  # Enthalpy variation of air in the dehumidifier (hot fluid)
                               tDeltHwat_de = thswod-self.hswid  # Enthalpy variation of water in the dehumidifier (cold fluid)

                               tDeltHair_hu = thaoh-thaih         # Enthalpy variation of the air in the humidifier (cold fluid)
                               tDeltHwat_hu = self.hswih-thswoh   # Enthalpy variation of water in the humidifier (hot fluid). It is not used.


                               if (self.msw*tDeltHwatMax_de)<(caud_a*tDeltHairMax_de):
                                   Ef_de_n = tDeltHwat_de/tDeltHwatMax_de
                               else:
                                   Ef_de_n = tDeltHair_de / tDeltHairMax_de

                               if (self.msw*self.hswih-caud_b*thswoh_ideal)<(caud_a*tDeltHairMax_hu):
                                   Ef_hu_n = (self.msw*self.hswih-caud_b*thswoh)/(self.msw*self.hswih-caud_b*thswoh_ideal)
                               else:
                                   Ef_hu_n = caud_a*tDeltHair_hu / caud_a*tDeltHairMax_hu


                               # Calculating the input and output energies of the dehumidifier and humidifier.

                               Enei_de = self.msw * self.hswid + caud_a * thaid
                               Eneo_de = self.msw * thswod + caud_a * thaod + caud_pw * thpw
                               Enei_hu = self.msw * self.hswih + caud_a * thaih
                               Eneo_hu = caud_b * thswoh + caud_a * thaoh

                              # Calculating the input and output entropies of the dehumidifier and humidifier.

                               Enti_de = self.msw * self.Sswid + caud_a * tSaid
                               Ento_de = self.msw * tSswod + caud_a * tSaod + caud_pw * tSpw
                               Enti_hu = self.msw * self.Sswih + caud_a * tSaih
                               Ento_hu = caud_b * tSswoh + caud_a * tSaoh

                               CEn_hu = self.CheqEne(Eneo_hu,Enei_hu)  # Compare input and output energies have a difference
                                                                       # as close as possible to zero.


                               CEn_de = self.CheqEne(Eneo_de, Enei_de)



                               CEt_de = self.CheqEnt(Ento_de, Enti_de)  # Compare input and output entropy greater than or equal to zero
                                                                        # (negative values very close to zero are contemplated)

                               CEt_hu = self.CheqEnt(Ento_hu, Enti_hu)


                               CEf_hu = self.CheqEne(self.Efh,Ef_hu_n)  # Comparing the effectiveness of the solving Ef_hu_n y Ef_de_n with
                                                                        # the one we want Eff_hu y Eff_de

                               CEf_de = self.CheqEne(self.Efd, Ef_de_n)


                               sumat = CEt_de + CEt_hu + CEn_hu + CEn_de + CEf_hu + CEf_de

                               if sumat < sumfin:
                                   self.taid = ttaid
                                   self.taod = ttaod
                                   self.taih = ttaih
                                   self.taoh = ttaoh
                                   self.tswod = ttswode
                                   self.tswoh = ttswohu
                                   self.ma = caud_a
                                   self.mpw = caud_pw
                                   self.mb = caud_b
                                   self.Ef_de = Ef_de_n
                                   self.Ef_hu = Ef_hu_n
                                   self.qi = self.msw*(self.hswih - thswod)
                                   self.GOR = (caud_pw * Lah) / (self.qi + 0.000000000000000000000000000000001)
                                   self.RR = caud_pw / self.msw
                                   if self.GOR >= 1.1 and self.GOR <= 3.5 and self.RR >= 0.01 and self.RR <= 0.049 and CEn_hu <= 10 and CEn_de <= 10 and CEf_hu <= 10 and CEf_de <= 10 and CEt_de <= 10 and CEt_hu <= 10:
                                      self.ValTerm = 1                         # The thermal solution is recorded as valid



                                   # Values indicating quality of the solution of the system of thermal equations
                                   self.DEfd = self.Efd - Ef_de_n
                                   self.DEfh = self.Efh - Ef_hu_n
                                   self.pDEfd = mat.fabs(self.DEfd*100/self.Efd)
                                   self.pDEfh = mat.fabs(self.DEfh*100/self.Efh)
                                   self.DEniod = Eneo_de - Enei_de
                                   self.DEnioh = Eneo_hu - Enei_hu
                                   self.pDEniod = self.DEniod*100/Eneo_de
                                   self.pDEnioh = self.DEnioh*100/Eneo_hu
                                   # Quality of solution of the system of thermal equations with respect to entropy (second law of thermodynamics)
                                   self.DEntiod = Ento_de - Enti_de    # Must be greater than or equal to zero
                                   self.DEntioh = Ento_hu - Enti_hu    # Must be greater than or equal to zero
                                   self.pDEntiod = self.DEntiod*100/Ento_de
                                   self.pDEntioh = self.DEntioh*100/Ento_hu
                                   sumfin = sumat
                                   self.haoh = thaoh
                                   self.Woh = tWoh
               print("The thermal values of DS-HDH have already been calculated, with an iteration step of = ",pas)

               # End of the fourth "for" ######################@@@@@@@@@@***********
     #@@@@@@@@@@@@*********End of the function to calculate HDH thermal variables.*******@@@@@@@@@@@@@@@@@
            except ZeroDivisionError:
                print("The values entered generate a division by 0, which is not admissible.")
            except:
                print("Wrong parameters have been entered to calculate the solar desalination system")
          # finally:
           #    print("System calculations have been completed")

     def calcultermico_fino(self,vt11,vt21,vt12,vt22,vcsma):
             #  Sending to calculate shorter and shorter temperature sections to increase accuracy
             # and decrease number of iterations.
             # vcsma: If ma is calculated or given. vcsma=0 calculated ---- vcsma=1 given
             # vt11 y vt21 start of the temperature variations
             # vt12 y vt22 end of the temperature variations
             pasos=[0.5,0.25,0.1,0.01]
             global vtp11, vtp21, vtp12, vtp22, trp1,trp2,trp3,trp4,rp,mfr
             for nn in pasos:
                 if nn==0.5:
                     vtp11=vt11
                     vtp21=vt21
                     vtp12=vt12
                     vtp22=vt22
                     trp1 = 1
                     trp2 = 1
                     trp3 = 1
                     trp4 = 1
                     rp = 0
                     mfr = 0
                 elif nn==0.25:
                     vtp11 = 0
                     vtp21 = 0
                     vtp12 = 5
                     vtp22 = 5
                     trp1 = self.taid
                     trp2 = self.taod
                     trp3 = self.tswoh
                     trp4 = self.tswod
                     rp = 2.5
                     mfr = 1
                 elif nn==0.1:
                     vtp11 = 0
                     vtp21 = 0
                     vtp12 = 2
                     vtp22 = 2
                     trp1 = self.taid
                     trp2 = self.taod
                     trp3 = self.tswoh
                     trp4 = self.tswod
                     rp = 1
                     mfr = 1
                 elif nn==0.01:
                     vtp11 = 0
                     vtp21 = 0
                     vtp12 = 0.2
                     vtp22 = 0.2
                     trp1 = self.taid
                     trp2 = self.taod
                     trp3 = self.tswoh
                     trp4 = self.tswod
                     rp = 0.1
                     mfr = 1
                 self.calcularhdh(vtp11, vtp21, vtp12, vtp22, nn,trp1,trp2,trp3,trp4,rp,mfr,vcsma)

     def MostrarValores(self):
                   print("GOR=", self.GOR, " - RR=", self.RR, " - Mpw=", self.mpw, "\n",
                         "It is a problem with base values: Msw=", self.msw, " - Ma=", self.ma, " - MR=", self.msw/self.ma, " - TwaterMax=",
                         self.tswih, " - TwaterMin=", self.tswid, "\n",
                         "                           RHih=", self.RH_ih, " - RHoh=", self.RH_oh, " - Efh=", self.Efh,
                         " - Efd=", self.Efd, "\n",
                         "Results of the solution of the system of equations:  Dif.Efect. dehum.=", self.DEfd,
                         " - Dif.Efect. hum.=", self.DEfh, " - % Dif.Efect. dehum.=", self.pDEfd,
                         " - % Dif.Efect. hum.=", self.pDEfh, "\n",
                         "                            Dif.Energ. in-ou dehum.=", self.DEniod,
                         "Dif.Energ. in-ou hum.=", self.DEnioh, "% Dif.Energ. in-ou dehum.=",
                         self.pDEniod, "% Dif.Energ. in-ou hum.=", self.pDEnioh, "\n",
                         "                            Dif.Entrop. in-ou dehum.=", self.DEntiod,
                         " - Dif.Entrop. in-ou hum.=", self.DEntioh,
                         " - % Dif.Entrop. in-ou dehum.=", self.pDEntiod,
                         " - % Dif.Entrop. in-ou hum.=",self.pDEntioh,"\n",
                         "            The calculated air and water temperatures at the extremes are as follows:             \n",
                         "taide  =  taohu            taode  =  taihu            tswohu            tswode \n",
                         "    ",self.taid,"                    ",self.taod,"                ",self.tswoh,"           ",self.tswod)

     def fun_y(self,xv):    # This intermediate function is used to calculate normalized enthalpy to simplify the creation
                            # of temperature vs. enthalpy midpoints.

         ent = (hw_s(xv, self.Sal))/self.ma
         return ent

     def crearpuntoslinrec(self,x1,x2):  # Create the straight line points: water temperature vs. enthalpy in the humidifier

         livar = [x1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, x2]
         # livag = [0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,xx,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32]
         lifun = [self.fun_y(x1), 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                  1, self.fun_y(x2)]

         livar[16] = (livar[0] + livar[32]) / 2
         lifun[16] = (self.fun_y(livar[0]) + self.fun_y(livar[32])) / 2

         livar[8] = (livar[0] + livar[16]) / 2
         lifun[8] = (self.fun_y(livar[0]) + self.fun_y(livar[16])) / 2

         livar[24] = (livar[16] + livar[32]) / 2
         lifun[24] = (self.fun_y(livar[16]) + self.fun_y(livar[32])) / 2

         livar[4] = (livar[0] + livar[8]) / 2
         lifun[4] = (self.fun_y(livar[0]) + self.fun_y(livar[8])) / 2

         livar[12] = (livar[8] + livar[16]) / 2
         lifun[12] = (self.fun_y(livar[8]) + self.fun_y(livar[16])) / 2

         livar[20] = (livar[16] + livar[24]) / 2
         lifun[20] = (self.fun_y(livar[16]) + self.fun_y(livar[24])) / 2

         livar[28] = (livar[24] + livar[32]) / 2
         lifun[28] = (self.fun_y(livar[24]) + self.fun_y(livar[32])) / 2

         livar[2] = (livar[0] + livar[4]) / 2
         lifun[2] = (self.fun_y(livar[0]) + self.fun_y(livar[4])) / 2

         livar[6] = (livar[4] + livar[8]) / 2
         lifun[6] = (self.fun_y(livar[4]) + self.fun_y(livar[8])) / 2

         livar[10] = (livar[8] + livar[12]) / 2
         lifun[10] = (self.fun_y(livar[8]) + self.fun_y(livar[12])) / 2

         livar[14] = (livar[12] + livar[16]) / 2
         lifun[14] = (self.fun_y(livar[12]) + self.fun_y(livar[16])) / 2

         livar[18] = (livar[20] + livar[16]) / 2
         lifun[18] = (self.fun_y(livar[20]) + self.fun_y(livar[16])) / 2

         livar[22] = (livar[24] + livar[20]) / 2
         lifun[22] = (self.fun_y(livar[24]) + self.fun_y(livar[20])) / 2

         livar[26] = (livar[24] + livar[28]) / 2
         lifun[26] = (self.fun_y(livar[24]) + self.fun_y(livar[28])) / 2

         livar[30] = (livar[28] + livar[32]) / 2
         lifun[30] = (self.fun_y(livar[28]) + self.fun_y(livar[32])) / 2

         livar[1] = (livar[0] + livar[2]) / 2
         lifun[1] = (self.fun_y(livar[0]) + self.fun_y(livar[2])) / 2

         livar[3] = (livar[2] + livar[4]) / 2
         lifun[3] = (self.fun_y(livar[2]) + self.fun_y(livar[4])) / 2

         livar[5] = (livar[4] + livar[6]) / 2
         lifun[5] = (self.fun_y(livar[4]) + self.fun_y(livar[6])) / 2

         livar[7] = (livar[6] + livar[8]) / 2
         lifun[7] = (self.fun_y(livar[6]) + self.fun_y(livar[8])) / 2

         livar[9] = (livar[8] + livar[10]) / 2
         lifun[9] = (self.fun_y(livar[8]) + self.fun_y(livar[10])) / 2

         livar[11] = (livar[10] + livar[12]) / 2
         lifun[11] = (self.fun_y(livar[10]) + self.fun_y(livar[12])) / 2

         livar[13] = (livar[12] + livar[14]) / 2
         lifun[13] = (self.fun_y(livar[12]) + self.fun_y(livar[14])) / 2

         livar[15] = (livar[14] + livar[16]) / 2
         lifun[15] = (self.fun_y(livar[14]) + self.fun_y(livar[16])) / 2

         livar[17] = (livar[16] + livar[18]) / 2
         lifun[17] = (self.fun_y(livar[16]) + self.fun_y(livar[18])) / 2

         livar[19] = (livar[18] + livar[20]) / 2
         lifun[19] = (self.fun_y(livar[18]) + self.fun_y(livar[20])) / 2

         livar[21] = (livar[20] + livar[22]) / 2
         lifun[21] = (self.fun_y(livar[20]) + self.fun_y(livar[22])) / 2

         livar[23] = (livar[22] + livar[24]) / 2
         lifun[23] = (self.fun_y(livar[22]) + self.fun_y(livar[24])) / 2

         livar[25] = (livar[24] + livar[26]) / 2
         lifun[25] = (self.fun_y(livar[24]) + self.fun_y(livar[26])) / 2

         livar[27] = (livar[26] + livar[28]) / 2
         lifun[27] = (self.fun_y(livar[26]) + self.fun_y(livar[28])) / 2

         livar[29] = (livar[28] + livar[30]) / 2
         lifun[29] = (self.fun_y(livar[28]) + self.fun_y(livar[30])) / 2

         livar[31] = (livar[30] + livar[32]) / 2
         lifun[31] = (self.fun_y(livar[30]) + self.fun_y(livar[32])) / 2

         parej = {'temp': livar, 'ental': lifun}

         listemp = parej['temp']

         listent = parej['ental']



         return parej

     def calcularpinchpoint(self,tx1,tx2):   # With this function, pairs of temperatures and differences of normalized enthalpies (DEN) are delivered,
                                             # to see if there are no negative differences of DEN.

           dicdev = self.crearpuntoslinrec(tx1,tx2)  # Moving on to generate a list of points on the temp. vs. enthalpy of water line

           listemp = dicdev['temp']  # The temperature list is stored here

           listent = dicdev['ental']  # The list of enthalpies is stored here

           listenta = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1]
                                                      # Initializing enthalpy of air list
           cont = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1]
                                     # Initializing list (cont) that will store differences between normalized water and air enthalpies of the humidifier
           global verc, mindex  # verc, If == 0, there are not negative difference in cont, If == 1, yes
                                # mindex, to put there the indexes of the two smaller values of cont
           verc = [0,0]                # Initialize verc
           mindex = list([0, 0])    # Initialize mindex

           for h_a in range(0, 33, 1):
               listenta[h_a] = (H_air(listemp[h_a], self.Sal, self.Pr_a, 1)) / self.ma

               cont[h_a] = listenta[h_a] - listent[h_a]  # The differences between both normalized enthalpies are saved, that is, hw,a/ma

               if cont[h_a] <= 0:  # This check is in case values less than or equal to zero occur, it is known that conditions are not met.
                   verc[0] = 1
                   print("Thermodynamic laws are not fulfilled.")
                   break
           print("List of normalized enthalpy differences without ordering", cont)
           contor = sorted(cont)  # Se ordena la lista cont para la nueva lista contor
           print("List of sorted normalized enthalpy differences",contor)

           conts = [0, 0]  # conts is initialized to only search for the values of interest
           conts[0] = contor[0]
           conts[1] = contor[1]

           for cc in range(0, 2, 1):
               for jj in range(0, 33, 1):
                   if conts[cc] == int(cont[jj]):
                       mindex[cc] = jj  # Here we compared the two values of conts with respect to all the values of cont,
                                        # that is, to see to which index corresponds the closest values of temp and the difference of both enthalpies.

           tpinch1 = listemp[mindex[0]]
           tpinch2 = listemp[mindex[1]]

           hpinch1 = conts[0]
           hpinch2 = conts[1]



           result = [tpinch1,tpinch2,hpinch1,hpinch2,verc[0],listemp,listent,listenta]    #  It is stored in this list, in the following order:
                                                         # The lowest candidate pinch temperature, the highest,
                                                         # the normalized enthalpy differences corresponding to the temperatures given above
                                                         # and whether there were negative values or not (verc=0 -No- ::: verc=1 -Yes-)
           print("The two temperatures and enthalpy differences",result)

           return result

     def validarseley(self):     # The two pinch temperature values are taken and if: their difference gives a relative value
                                 # less than 0.00001 and it is verified there are no negative values the pinch point is given
                                 # as from the average of those temperatures, and it is validated of course.
             global sx1, sx2, tpinc,expl,cumsley
                                    # If cumsley = 0, is not met, If cumsley = 1  it is met
             sx1 = self.tswoh
             sx2 = self.tswih
             tpinc = [1, 1]
             verif = self.calcularpinchpoint(sx1, sx2)

             listemp = verif[5]
             listent = verif[6]
             listenta = verif[7]
             dicdev1 = self.crearpuntoslinrec(self.tswid,self.tswod)  # Moving on to generate a list of points on the temperature vs. enthalpy of water straight line
             listdesh =   [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1, 1, 1]
             listempd = dicdev1['temp']  # The temperature list is stored here
             listdeshg = dicdev1['ental']  # The list of enthalpies is stored here
             for thb in range(0, 33, 1):
                  listdesh[thb] = listdeshg[thb] / self.ma

             m = 0
             if verif[4] == 0:

                 while m == 0:
                     verif = self.calcularpinchpoint(sx1,sx2)



                     if math.fabs((verif[0] - verif[1]) * 100 / verif[0]) <= 0.00001:
                         tpinc[0] = (verif[0] + verif[1]) / 2
                         hap = (H_air(tpinc[0], self.Sal, self.Pr_a, 1)) / self.ma
                         hwp = (hw_s(tpinc[0], self.Sal)) / self.ma
                         tpinc[1] = hap - hwp
                         m = 1


                     else:
                         sx1 = verif[0]
                         sx2 = verif[1]
                     cumsley = 1
                     expl = "The pinch point, in temperature and pinch enthalpy, respectively, is: "
             else:
                 cumsley = 0
                 expl = "There is no pinch point, as the solution is not valid due to negative values in the normalized temp. vs. enthalp. curve: "
             print(expl,tpinc)
             return cumsley   # Here it provides the pinch point, temperature and the difference of normalized enthalpies.

     def validarseleyextdeshumi(self):     # To validate the graph temperature vs. enthalpy dehumidifier ends
                                           # If comp = 0 one of the ends does not comply with the second law If comp = 1 it does comply with it one
           global comp
           comp = 0

           hapb = (H_air(self.taod, self.Sal, self.Pr_a, 1)) / self.ma
           hwpb = (hw_s(self.tswid, self.Sal)) / self.ma

           hapt = (H_air(self.taid, self.Sal, self.Pr_a, 1)) / self.ma
           hwpt = (hw_s(self.tswod, self.Sal)) / self.ma

           if (hapb - hwpb) > 0 and (hapt - hwpt) > 0:
               comp = 1

           return comp



     def hdhcalcul(self,vt11,vt21,vt12,vt22,iter_rk,uu,n1,n2,n3,n4,iradsol, eficcolsol,calmaf):


      # Function to calculate solar humidifier-dehumidifier desalination system, open water circulation and closed air circulation,
      # water  heating . Includes all parameters, thermal and geometry.

      # fmsw mass flow rate of seawater fed to the system
      # ftswid,ftswih Temperature of the water entering the system through the dehumidifier and temperature of the water entering the humidifier
      # frhih,frhoh Relative humidity of the air, at the humidifier inlet and outlet respectively.
      # fefd,fefh   Effectiveness of dehumidifier and humidifier respectively
      # fsal,fpres  Water salinity and air pressure (ambient pressure in this framework) respectively

      # vt11,vt21,vt12,vt22,vcsma
      # calmaf: If ma is calculated or given. calmaf=0 it is calculated ---- calmaf=1 it is given
      #  vt11 y vt21 Start of the temperature variation
      #vt12 y vt22 End of the temperature variation
      # iter_rk Number of iterations for the Runge-Kutta method of calculating the Merkel number of the packed bed
      # uu   Overall heat transfer coefficient, given by the characteristics of the heat exchanger of dehumidifier(W/(m2*K)).
      #      Range from 5 to 100 W/(m2*K)
      # n1,n2,n3,n4: These are coefficients that depend on the values given by the manufacturer of the evaporator filler.
      # Interpretation of the formula:
      #               n2            n4       From this equation, Hh can be solved as follows
      # me = n1  *  mr  *  n3  *  Hh
      # iradsol  Intensity of solar radiation (W/m2)
     # eficcolsol   Solar heater conversion efficiency (-)



           global hdhret,value,valint,chhp
           value = 0
           hdhret = [1000000000, 1000000000, 1000000000, 2,1,1,1,1,1,1,1,1]
           valint = 2
           sg4 = 0
           chhp=1
           self.calcultermico_fino(vt11,vt21,vt12,vt22,calmaf)    # Calculating thermal design variables of the system


           sx1 = self.tswoh
           sx2 = self.tswih
           verif = self.calcularpinchpoint(sx1, sx2)
           if verif[4] == 0:
               value = 1                                         # If verif[4] = 0  value = 0, the non-overlapping verification of pinch enthalpies in the humidifier is fulfilled.
           value2 = self.validarseleyextdeshumi()                   # If value2 = 0, the non-overlapping verification of pinch enthalpies in the dehumidifier is fulfilled.

           if value == 1 and value2 == 1 and self.ValTerm == 1:
               tntu = Triboix_ntu_de(self.tswid,self.taid,self.msw, self.ma, self.Sal, self.Pr_a, self.RH_oh,self.Efd)
               # Effective area of the condenser (dehumidifier) will be calculated.

               valrungkutt = Runge_Kutta_merk_pop_values(self.tswoh,self.taih,self.RH_ih,self.Sal,self.Pr_a,self.msw,self.ma,self.mb,self.tswih,iter_rk)
               valrungkutt2 = Runge_Kutta_merk_pop_values(self.tswoh, self.taih, self.RH_ih, self.Sal, self.Pr_a,self.msw, self.ma, self.mb, self.tswih, iter_rk+2)
               # Merkel number is calculated (it gives also other values, air enthalpy and absolute humidity output humidifier
               # those data are compared with the thermal calculations and validate this other, they should be as close as possible)

               clv = 'clav' + str(iter_rk-1)
               clv2 = 'clav' + str(iter_rk + 2 - 1)
               print("The values of wou, hou and Me in that order for the packed bed are as follows: ",valrungkutt[clv])       # Esto muestra los valores de la última iteración del método, [0]=wou -- [1]=hou -- [2]=me
               sacarme = valrungkutt[clv]

               print("The Merkel number for this case is me: ",sacarme[2])

               sacarme2 = valrungkutt2[clv2]

               me = sacarme[2]   # This assigns to a variable the number of Merkel obtained
               wrk = sacarme[0]  # This assigns to a variable the value of Woh according to Runge - Kutta
               hrk = sacarme[1]  # This assigns to a variable the value of hoh according to Runge - Kutta

               me1 = sacarme2[2]  # This assigns to a variable the number of Merkel obtained
               wrk1 = sacarme2[0]  # This assigns to a variable the value of Woh according to Runge - Kutta
               hrk1 = sacarme2[1]  # This assigns to a variable the value of hoh according to Runge - Kutta

               come = (abs(me1 - me)) * 100 / me1

               self.nMe = me1
               hre = AltHumid(n1,n2,n3,n4,self.msw, self.ma,self.nMe)

               self.Hh = hre
               relflu = 1.5     #   This value of water mass flow at the humidifier inlet is recommended by several authors.

               asupre = AreaRellEvapHum(self.Hh, relflu, self.msw) # ###@@ Calculating the surface area of the packed bed



               acsol = AreaCalSo(self.qi, iradsol, eficcolsol)    # ###@@ Calculating the area of the solar water heater

               aefdesh = areacondendeshumid(tntu, uu, self.msw, self.ma,self.tswid, self.taid, self.Sal, self.Pr_a, self.RH_oh)
      # ###@@ Calculating the effective area of the dehumidifier condenser.



               self.Atot = asupre + aefdesh + acsol     # ###@@ Calculating the total area (related to the cost) of the system

               print("The relative total area of the system (serves as a reference to its cost) is: ",self.Atot," m2")
               print("The GOR of the system is: ",self.GOR," (-)")
               print("The RR of the system is: : ",self.RR," (-)")

               f1 = self.Atot / self.mpw  ###@@ Value of objective function 1, for this run and this desalinator configuration
               f2 = 1 / self.RR             ###@@ Value of objective function 2, for this run and this desalinator configuration
               f3 = self.Atot / self.GOR  ###@@ Value of objective function 3, for this run and this desalinator configuration

               print(" The ratio of total area (related to cost) of the system with respect to the mass flow rate of pure water produced f1 : ",f1,
                     "\n The ratio of mass flow rates of water fed to the system with respect to potable water produced f2: ",f2,
                     "\n The ratio of total area (related to cost) of the system with respect to GOR f3 : ",f3)

               self.MostrarValores()  # Showing thermal values of the system.

               difw = abs(wrk1-self.Woh)
               porw = difw * 100 / wrk1
               difh = abs(hrk1-self.haoh)
               porh = difh * 100 / hrk1


               chhp = ((self.mpw * 3600) / self.Hh)
               if chhp >= 1:
                   self.verhp = 1

               print(
                   "Runge - Kutta validation \n Absolute humidity value at humidifier outlet according to thermal calculation and then according to calculation by Rung-Kut with percent error in that order",
                   self.Woh, " - ", wrk, " - ", porw,
                   "\n Enthalpy value of the air leaving the humidifier, enthalpy value at the humidifier outlet, but according to Rung-Kut and percent of error in that order",
                   self.haoh, " - ", hrk, " - ", porh)


               if porw < 2 and porh < 2 and come < 0.5:
                   self.ValRK = 1
                   valint = self.ValTerm * self.ValRK * value * value2 * self.verhp
               if valint == 1:
                   self.valt = -1
               hdhret = [f1, f2, f3, self.valt,asupre,aefdesh,acsol,self.Hh,self.Atot,self.GOR,self.RR,self.mpw]
          #             [f1 ,f2, f3, valido, AH,Adesh,Aswh,Hp,At,GOR,RR,mpw]


           print("Transferred for processing Pareto Front",hdhret)
           return hdhret

