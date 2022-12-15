from sympy.concrete.tests.test_sums_products import NS

from sympy.core.singleton import S
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.physics.units import convert_to, coulomb_constant, elementary_charge, gravitational_constant, planck
from sympy.physics.units.definitions.unit_definitions import statcoulomb, coulomb, second, gram, centimeter, erg, \
    newton, joule, dyne, speed_of_light, meter, farad, henry, statvolt, volt, ohm
from sympy.physics.units.systems import SI
from sympy.physics.units.systems.cgs import cgs_gauss


def test_conversion_to_from_si():
    assert convert_to(statcoulomb, coulomb, cgs_gauss) == coulomb/2997924580
    assert convert_to(coulomb, statcoulomb, cgs_gauss) == 2997924580*statcoulomb
    assert convert_to(statcoulomb, sqrt(gram*centimeter**3)/second, cgs_gauss) == centimeter**(S(3)/2)*sqrt(gram)/second
    assert convert_to(coulomb, sqrt(gram*centimeter**3)/second, cgs_gauss) == 2997924580*centimeter**(S(3)/2)*sqrt(gram)/second

    # SI units have an additional base unit, no conversion in case of electromagnetism:
    assert convert_to(coulomb, statcoulomb, SI) == coulomb
    assert convert_to(statcoulomb, coulomb, SI) == statcoulomb

    # SI without electromagnetism:
    assert convert_to(erg, joule, SI) == joule/10**7
    assert convert_to(erg, joule, cgs_gauss) == joule/10**7
    assert convert_to(joule, erg, SI) == 10**7*erg
    assert convert_to(joule, erg, cgs_gauss) == 10**7*erg


    assert convert_to(dyne, newton, SI) == newton/10**5
    assert convert_to(dyne, newton, cgs_gauss) == newton/10**5
    assert convert_to(newton, dyne, SI) == 10**5*dyne
    assert convert_to(newton, dyne, cgs_gauss) == 10**5*dyne


def test_ohm_cgs_gauss():

    assert convert_to(ohm,second/centimeter,cgs_gauss) == 25000*second/(22468879468420441*centimeter)
    assert convert_to(1*ohm,second/centimeter,cgs_gauss) == 25000*second/(22468879468420441*centimeter)
    assert convert_to(2*ohm,second/centimeter,cgs_gauss) == 50000*second/(22468879468420441*centimeter)
    assert NS(convert_to(ohm,second/centimeter,cgs_gauss)) == '1.11265005605362e-12*second/centimeter'
    assert NS(convert_to(2*ohm,second/centimeter,cgs_gauss)) == '2.22530011210724e-12*second/centimeter'

def test_henry_cgs_gauss():
    assert convert_to(henry,second**2/centimeter,cgs_gauss) == 25000*second**2/(22468879468420441*centimeter)
    assert convert_to(1*henry,second**2/centimeter,cgs_gauss) == 25000*second**2/(22468879468420441*centimeter)
    assert convert_to(2*henry,second**2/centimeter,cgs_gauss) == 50000*second**2/(22468879468420441*centimeter)
    assert NS(convert_to(henry,second**2/centimeter,cgs_gauss)) == '1.11265005605362e-12*second**2/centimeter'
    assert NS(convert_to(2*henry,second**2/centimeter,cgs_gauss)) == '2.22530011210724e-12*second**2/centimeter'

def test_volt_cgs_gauss():
    assert convert_to(volt,statvolt,cgs_gauss) == 10**6*statvolt/299792458
    assert convert_to(1*volt,statvolt,cgs_gauss) == 10**6*statvolt/299792458
    assert convert_to(2*volt,statvolt,cgs_gauss) == 2*10**6*statvolt/299792458

def test_farad_cgs_gauss():

   assert convert_to(farad,centimeter,cgs_gauss) == 299792458**2*centimeter/10**5
   assert convert_to(1*farad,centimeter,cgs_gauss) == 299792458**2*centimeter/10**5
   assert convert_to(2*farad,centimeter,cgs_gauss) == 2*299792458**2*centimeter/10**5



def test_cgs_gauss_convert_constants():

    assert convert_to(speed_of_light, centimeter/second, cgs_gauss) == 29979245800*centimeter/second

    assert convert_to(coulomb_constant, 1, cgs_gauss) == 1
    assert convert_to(coulomb_constant, newton*meter**2/coulomb**2, cgs_gauss) == 22468879468420441*meter**2*newton/(2500000*coulomb**2)
    assert convert_to(coulomb_constant, newton*meter**2/coulomb**2, SI) == 22468879468420441*meter**2*newton/(2500000*coulomb**2)
    assert convert_to(coulomb_constant, dyne*centimeter**2/statcoulomb**2, cgs_gauss) == centimeter**2*dyne/statcoulomb**2
    assert convert_to(coulomb_constant, 1, SI) == coulomb_constant
    assert NS(convert_to(coulomb_constant, newton*meter**2/coulomb**2, SI)) == '8987551787.36818*meter**2*newton/coulomb**2'

    assert convert_to(elementary_charge, statcoulomb, cgs_gauss)
    assert convert_to(gravitational_constant, dyne*centimeter**2/gram**2, cgs_gauss)
    assert NS(convert_to(planck, erg*second, cgs_gauss)) == '6.62607015e-27*erg*second'
