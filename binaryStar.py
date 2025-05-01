import math
from astropy.constants import M_sun, L_sun
from astropy import units as u
from dataclasses import dataclass

@dataclass
class Hvezda:
   mag: float = 0.0
   abs_mag: float = 0.0
   zarivy_vykon: float = 0.0
   hmotnost: float = 0.0

# T celkova obezna doba v letech
# M_celkove_Ms Hmotnost soustavy v nasobcich M Slunce
# return: Velikost hlavni poloosy v AU
def hlavni_poloosa_krok2(T_roky, M_celkove_Ms):
    return (T_roky**2 * (M_celkove_Ms))**(1/3)

# a_AU Velikost hlavni poloosy v AU
# uhlova_vterina - uhlovy rozmer hlavni poloosy v uhlovych vterinach
# return: Vzdalenost soustavy od pozorovatele
def vzdalenost_krok3(a_AU, uhlova_vterina):
    # 1 st = 3600 vterin
    # 1 st = PI/180
    # 1 vterina = PI/(180*3600) rad
    return a_AU / math.tan(uhlova_vterina *math.pi/(180*3600))

# mag relativni mag
# d vzdalenost v PC
# return: abs. mag
def abs_magnituda_krok4(mag, d_pc):
    return mag + 5 - 5 * math.log10(d_pc)

# absMag absolutni mag
# return: Zarivy vykon v W
def zarivy_vykon_krok5(absMag):
    # absolutni mag Slunce  4.83 - viz Skripta Uvod to Astronomie 
    return L_sun * 10 **(((4.83 * u.mag).value- absMag) / 2.5)

# L zarivy vykon v W
# return: Hmotnost v kg
def hmotnost_z_zariveho_krok6(L):
    return M_sun * ((L / L_sun) ** (1/3.5))

def main():
    #uvodni inicializace ze zadani
    a_uhlove = 4.5
    b_uhlove = 3.4
    T_castecna = 11
    H1_mag = 3.9
    H2_mag = 5.3
    tolerance = 0.01

    # excentricita
    h = math.sqrt(a_uhlove**2 - b_uhlove**2)
    
    # plocha usece
    epsilon = a_uhlove * b_uhlove * (math.acos(h / a_uhlove) - (h / a_uhlove**2) * math.sqrt(a_uhlove**2 - h**2))

    # plocha cele elipsy
    S = math.pi * a_uhlove * b_uhlove

    #perioda obehu
    T_celkove = S / epsilon * T_castecna

    # ******************************
    # methoda dynamicke paralaxe
    
    # pocatecni hodnoty v nasobcich Ms
    M_celkove_Ms = 2.0
    diff = float('inf')

    H1 = Hvezda(H1_mag)
    H2 = Hvezda(H2_mag)

    while diff > tolerance * M_celkove_Ms:
        a_AU = hlavni_poloosa_krok2(T_celkove, M_celkove_Ms)
        d_AU = vzdalenost_krok3(a_AU, a_uhlove)
        # vzdalenost d_AU na parseky
        d_pc = (d_AU * u.AU).to(u.pc).value
        
        # Absolutni magnituda
        H1.abs_mag = abs_magnituda_krok4(H1.mag, d_pc)
        H2.abs_mag = abs_magnituda_krok4(H2.mag, d_pc)
       
        # zarivy vykon
        H1.zarivy_vykon = zarivy_vykon_krok5(H1.abs_mag)
        H2.zarivy_vykon = zarivy_vykon_krok5(H2.abs_mag)

        # Hmotnosti
        H1.hmotnost = hmotnost_z_zariveho_krok6(H1.zarivy_vykon)
        H2.hmotnost = hmotnost_z_zariveho_krok6(H2.zarivy_vykon)

        nove_M_celkove_Ms = (H1.hmotnost/M_sun) + (H2.hmotnost/M_sun)
        diff = abs(nove_M_celkove_Ms - M_celkove_Ms)
        M_celkove_Ms = nove_M_celkove_Ms

    print ("Obezna perioda = {0} roku".format(T_celkove))
    print ("Vzdalenost soustavy = {0} AU = {1}".format(d_AU, (d_AU * u.AU).to(u.lyr)));
    print ("Abs Mag H1 = {0}".format(H1.abs_mag))
    print ("Abs Mag H2 = {0}".format(H2.abs_mag))
    print ("Hmotnost H1 = {0} Ms".format(H1.hmotnost/M_sun))
    print ("Hmotnost H2 = {0} Ms".format(H2.hmotnost/M_sun))
    
if __name__=="__main__":
    main()  
