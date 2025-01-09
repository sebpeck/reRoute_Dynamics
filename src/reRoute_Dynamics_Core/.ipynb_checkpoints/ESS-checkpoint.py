'''
ESS.py
'''

DEFAULT_MOTOR_EFF = .916 # unknown source, not listed in paper -- tourque and speed depemdamt, next level would be based on tourqe
DEFAULT_INVERTER_EFF =  .971 # unknown source, not listed in paper
DEFAULT_AUX_EFF = .89 # unknown source, Auxulliary system efficiency?
DEFAULT_LOAD_SIMPLE = 7000 # watts, Abdelaty & Mohamed []
DEFAULT_REGEN_EFF = .6 # Gallet et al []
DEFAULT_MAX_REGEN = -100000 # Watts

def calc_instance_power(value,
                        motor_eff=DEFAULT_MOTOR_EFF,
                        invert_eff=DEFAULT_INVERTER_EFF,
                        aux_eff=DEFAULT_AUX_EFF,
                        aux_load=DEFAULT_LOAD_SIMPLE,
                        regen_eff=DEFAULT_REGEN_EFF,
                        max_regen=DEFAULT_MAX_REGEN):
    '''
    calc_instance_power takes in a power value,
    and converts it to the corresponding load on the ESS.
    This is a simple stopgap.

    Parameters:
    value: a power value in Watts, as an int or float.

    Returns:
    converted battery power as a float.
    '''

    # set the battery power to zero.
    bat_pow = 0

    # Including Auxilliary load, though not strictly important at the moment. 
    if (value >= 0):
        # Discharging, converting the needed power into power battery must exert
        bat_pow = value/(motor_eff*invert_eff) + (aux_load/aux_eff)
    elif(value*regen_eff*motor_eff > max_regen):

        #charging, the regenerative braking ALL the time, max regen is 100
        bat_pow = (value*regen_eff*motor_eff) + (aux_load/aux_eff)
        
    else:
        bat_pow = max_regen + (aux_load/aux_eff)


    # Return the battery power.
    return bat_pow

    